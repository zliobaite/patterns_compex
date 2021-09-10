import pdb
import datetime
import re
import os
import numpy
import scipy.stats

from TuningPeg_StronaUG17 import CB, UG, UG_mod
# FF null model, Curveball (Strona et al. 2014)
# unbiased proportional null model by Ulrich and Gotelli (2012)

circ_def = 2*numpy.pi*6370997.

GEN_WITH_INDET = True


def common_preff(pieces):
    if len(pieces) == 0:
        return ""
    elif len(pieces) == 1:
        return pieces[0]
    ps = sorted(pieces, key=lambda x: len(x))
    i = 0
    while i < len(ps[0]) and all([ps[0][i] == pj[i] for pj in ps[1:]]):
        i += 1
    return ps[0][:i]


def combined_preff(pieces):
    if len(pieces) == 0:
        return ""
    elif len(pieces) == 1:
        return pieces[0]
    c = common_preff(pieces)
    return c+"+".join([p[len(c):] for p in pieces])


def ratio_str(v_num, v_den):
    # assume: v_num and v_den both >= 0, and v_num <= v_den
    if v_den == 0:
        return "---"
    elif v_num == v_den:
        return "1"
    else:
        return ("%.3f" % (v_num/v_den))[1:]


def cumcounts_str(marg):
    return " ".join(["% 3d" % v for v in numpy.cumsum(numpy.bincount(marg))])


def make_fld(fld):
    if not os.path.exists(fld):
        os.makedirs(fld)


def get_ranks(vs):
    if len(vs) == 0:
        return numpy.array([])
    if len(vs) == 1:
        return numpy.array([1])
    ord_ids = numpy.argsort(vs)
    ord_vals = vs[ord_ids]
    ord_ranks = 1.+numpy.arange(len(ord_vals), dtype=float)
    counts_ties = []
    if numpy.min(ord_vals[1:] - ord_vals[:-1]) == 0:
        b_start = 0
        while b_start < len(ord_vals):
            b_stop = b_start + 1
            while b_stop < len(ord_vals) and ord_vals[b_start] == ord_vals[b_stop]:
                b_stop += 1
            if b_stop > (b_start + 1):
                m_rnk = (ord_ranks[b_start] + ord_ranks[b_stop-1])/2
                ord_ranks[b_start:b_stop] = m_rnk
                counts_ties.append((m_rnk, b_stop-b_start))
            b_start = b_stop
    ranks = -numpy.ones(len(ord_vals))
    ranks[ord_ids] = ord_ranks
    return ranks, counts_ties


def format_field(d, head, k, fields_format):
    if k in fields_format:
        return fields_format[k](d[head[k]])
    return d[head[k]]

# COMPARING SPECIES


def is_indeterminate_str(v):
    return "." in v  # or re.match("incertae", v)


def is_indeterminate_gs(gen, spc):
    return is_indeterminate_str(gen) or is_indeterminate_str(spc)


def is_indeterminate(sid, species, genus_field, species_field):
    return is_indeterminate_gs(species[sid][genus_field], species[sid][species_field])


def nb_distinct_species(subset, species, genus_field, species_field):
    dets = [(species[s][genus_field], species[s][species_field]) for s in subset if not is_indeterminate_gs(species[s][genus_field], species[s][species_field])]
    nb_indet = len(subset) - len(dets)
    nb_det_distinct = len(set(dets))
    # return nb_det_distinct, max(nb_det_distinct, 1*(nb_indet > 0)), nb_det_distinct+nb_indet
    return nb_det_distinct+nb_indet


def is_pair_same(sidA, sidB, species, genus_field, species_field):
    return species[sidA][genus_field] == species[sidB][genus_field]


def is_pair_indeterminate(sidA, sidB, species, genus_field, species_field):
    return is_indeterminate(sidA, species, genus_field, species_field) or is_indeterminate(sidB, species, genus_field, species_field)


def pair_status(sidA, sidB, species, genus_field, species_field):
    if is_pair_indeterminate(sidA, sidB, species, genus_field, species_field):
        return -1
    elif is_pair_same(sidA, sidB, species, genus_field, species_field):
        return 1
    else:
        return 0


def pair_group(sidA, sidB, species, order_field, orders_to_num):
    if species[sidA][order_field] == species[sidB][order_field]:
        return orders_to_num.get(species[sidA][order_field], -1)
    return -1


def pairA_group(sidA, sidB, species, order_field, orders_to_num):
    return orders_to_num.get(species[sidA][order_field], -1)


def pairB_group(sidA, sidB, species, order_field, orders_to_num):
    return orders_to_num.get(species[sidB][order_field], -1)


# TIME SLICES
def get_slice_info(si, slce, slices_fields, pref="", missing=None, to_str=False):
    if to_str:
        tmp = {pref+"MAX_AGE": "%s" % slce[0], pref+"MIN_AGE": "%s" % slce[1], pref+"NAME": slce[-1], pref+"ID": "%s" % si}
    else:
        tmp = {pref+"MAX_AGE": slce[0], pref+"MIN_AGE": slce[1], pref+"NAME": slce[-1], pref+"ID": si}
    return [tmp.get(f, missing) for f in slices_fields]


def is_in_slice(site_max_age, site_min_age, slice_max_age, slice_min_age, max_ratio_out):
    # if any([sitet[0] <= s[0] and sitet[1] >= s[1] for s in slices]):
    if site_max_age > slice_min_age and site_min_age < slice_max_age:
        outside = max(0, site_max_age - slice_max_age) + max(0, slice_min_age - site_min_age)
        return outside <= max_ratio_out*(site_max_age-site_min_age)
    return False


def get_slice(slices, site, map_fields_sites, max_ratio_out=0):
    cands = [(si, slce) for si, slce in enumerate(slices) if is_in_slice(float(site[map_fields_sites["MAX_AGE"]]), float(site[map_fields_sites["MIN_AGE"]]), slce[0], slce[1], max_ratio_out=max_ratio_out)]
    if len(cands) == 1:
        return cands[0]
    else:
        return (-1, cands)


def assign_sites_to_slices(slices, sites, map_fields_sites, max_ratio_out=0):
    return [[sj for sj, site in enumerate(sites) if is_in_slice(site[map_fields_sites["MAX_AGE"]], site[map_fields_sites["MIN_AGE"]], slce[0], slce[1], max_ratio_out=max_ratio_out)] for si, slce in enumerate(slices)]


def compute_slices_overlap(slices, sites, map_fields_sites, file_out):  # DATA_FLD+"overlaps_matrix.csv"
    contained = []
    O = numpy.zeros((len(slices)+2, len(sites)+2))
    for i, s in enumerate(slices):
        O[i, -1] = s[1]
        O[i, -2] = s[0]
        for j, site in enumerate(sites):
            if i == 0:
                O[-1, j] = site[map_fields_sites["MIN_AGE"]]
                O[-2, j] = site[map_fields_sites["MAX_AGE"]]

            if site[map_fields_sites["MAX_AGE"]] >= s[1] and site[map_fields_sites["MIN_AGE"]] <= s[0]:
                maxa = max(site[map_fields_sites["MAX_AGE"]], s[0])
                mina = min(site[map_fields_sites["MIN_AGE"]], s[1])

                O[i, j] = (s[0]-s[1])/(maxa-mina)

        contained.append([sites_dict[site[map_fields_sites["LIDNUM"]]] for site in sites
                          if site[map_fields_sites["MAX_AGE"]] <= s[0] and site[map_fields_sites["MIN_AGE"]] >= s[1]])

    I = numpy.where(numpy.sum(O[:-2, :-2] > 0, axis=1) > 0)[0]
    J = numpy.where(numpy.sum(O[:-2, :-2] == 1, axis=0) == 0)[0]

    for i, j in zip(*numpy.where(O[:-2, J] > 0.95)):
        print(O[i, J[j]], "\t", i, J[j], "\t", O[[-1, -2], J[j]], "\t", O[i, [-1, -2]])

    numpy.savetxt(file_out, O[numpy.append(I, [-1, -2]), :][:, numpy.append(J, [-1, -2])], fmt='%.3f')


# GEO DISTANCES
def greatCircleAngle(lat, lng):
    # compute great circle distance for degree coordinates (lat, long)
    rd0 = numpy.radians(lat)
    rd1 = numpy.radians(lng)
    return numpy.arccos(numpy.sin(rd0[0]) * numpy.sin(rd1[0])
                        + numpy.cos(rd0[0]) * numpy.cos(rd1[0]) * numpy.cos(rd0[1] - rd1[1]))


def greatCircleDist(lat, lng):
    # compute great circle distance for degree coordinates (lat, long)
    return circ_def * greatCircleAngle(lat, lng)

# INPUT/OUTPUT


def read_slices(slices_file):
    sep = ","
    slices = []
    with open(slices_file) as fp:
        for line in fp:
            parts = line.strip().split(sep)
            if True:  # len(parts) == 3 or (len(parts) > 3 and parts[3] != "Age"):
                slices.append((float(parts[1]), float(parts[2]), parts[0]))
    return slices


def read_centroids(centroids_file, terra_file, min_ratio_land, boundaries):
    C = numpy.loadtxt(centroids_file)
    T = numpy.loadtxt(terra_file, skiprows=1)
    X = numpy.zeros(int(numpy.max(T[:, 0])), dtype=bool)
    X[T[numpy.where(T[:, 1] > min_ratio_land)[0], 0].astype('int')-1] = 1
    within = X & (C[:, 0] >= boundaries["LONG_MIN"]) & (C[:, 0] <= boundaries["LONG_MAX"]) & (C[:, 1] >= boundaries["LAT_MIN"]) & (C[:, 1] <= boundaries["LAT_MAX"])
    return C[within, :]


def get_site_key(site, map_fields_sites, with_lidnum=True):
    if with_lidnum:
        return site[map_fields_sites["LIDNUM"]]
    else:
        return (float(site[map_fields_sites["LAT"]]), float(site[map_fields_sites["LONG"]]), float(site[map_fields_sites["MAX_AGE"]]), float(site[map_fields_sites["MIN_AGE"]]))


def read_fossils(list_file, slices_file, fields_format, fields_sites, fields_species, exclude_indet=False):
    head = None
    data = []
    sep = ","
    with open(list_file) as fp:
        for line in fp:
            parts = line.strip().split(sep)
            if head is None:
                head = dict([(v, k) for (k, v) in enumerate(parts)])
            else:
                # if re.match("mn[0-9\-]+", parts[head["BFA_MIN"]]) and re.match("mn[0-9\-]+", parts[head["BFA_MAX"]]) and parts[head["BFA_MAX"]] == parts[head["BFA_MIN"]]:
                # if not re.match("mn[0-9\-]+", parts[head["BFA_MIN"]]) or not re.match("mn[0-9\-]+", parts[head["BFA_MAX"]]):
                # if re.match("c[0-9]", parts[head["BFA_MIN"]]) and re.match("c[0-9]", parts[head["BFA_MAX"]]):
                if True:  # not (re.match("mn[0-9\-]+", parts[head["BFA_MIN"]]) or re.match("c[0-9]", parts[head["BFA_MIN"]])) and not (re.match("mn[0-9\-]+", parts[head["BFA_MAX"]]) or re.match("c[0-9]", parts[head["BFA_MAX"]])):

                    data.append(tuple(parts))

    sites = sorted(set([tuple([format_field(d, head, k, fields_format) for k in fields_sites]) for d in data]))
    map_fields_sites = dict([(v, k) for (k, v) in enumerate(fields_sites)])
    with_lidnum = ("LIDNUM" in map_fields_sites)
    if "SLICE_ID" in map_fields_sites and slices_file is not None:
        slice_index = map_fields_sites["SLICE_ID"]
        slices = read_slices(slices_file)
    else:
        slice_index = None
        slices = []

    sites_dict = {}
    sliced_sites = [[] for s in slices]
    for k, v in enumerate(sites):
        sk = get_site_key(v, map_fields_sites, with_lidnum)
        sites_dict[sk] = k
        if slice_index is not None and v[slice_index] < len(sliced_sites):
            sliced_sites[v[slice_index]].append(k)

    species = sorted(set([tuple([format_field(d, head, k, fields_format) for k in fields_species]) for d in data]))
    if exclude_indet:
        n = len(species)
        gen_field = fields_species.index("GENUS")
        spc_field = fields_species.index("SPECIES")
        species = [s for s in species if not is_indeterminate_gs(s[gen_field], s[spc_field])]
        print("From %d to %d determinate species..." % (n, len(species)))
    species_dict = dict([(v, k) for k, v in enumerate(species)])
    occurrences = numpy.zeros((len(sites), len(species)), dtype=bool)  # store counts? (few multi-occurrences)
    extras = {}
    for di, d in enumerate(data):
        site_key = get_site_key(d, head, with_lidnum)
        # site_key = d[head[fields_sites[0]]]
        species_key = tuple([d[head[k]] for k in fields_species])
        if species_key in species_dict:
            if occurrences[sites_dict[site_key], species_dict[species_key]] > 0:
                if species_key not in extras:
                    extras[species_key] = []
                extras[species_key].append((site_key, di))
            else:
                occurrences[sites_dict[site_key], species_dict[species_key]] = 1

    # check extras occur in only once at each site
    print("occs", occurrences.shape, "extras", len(extras))
    uniq_xtr = {}
    for xk, xs in extras.items():
        uniq_sites = [set()]
        for x in xs:
            usi = 0
            while usi < len(uniq_sites) and x[0] in uniq_sites[usi]:
                usi += 1
            if usi >= len(uniq_sites):
                uniq_sites.append(set())
            uniq_sites[usi].add(x[0])
        for usi, uss in enumerate(uniq_sites):
            spc_k = tuple([x for x in xk[:-1]]+["%s_%s" % (xk[-1], usi)])
            uniq_xtr[spc_k] = uss
    if len(uniq_xtr) > 0:
        xtr_occurrences = numpy.zeros((len(sites), len(uniq_xtr)), dtype=bool)  # store counts? (few multi-occurrences)
        for xi, (xspc, xsites) in enumerate(sorted(uniq_xtr.items())):
            xtr_occurrences[[sites_dict[site_key] for site_key in xsites], xi] = 1
            species_dict[xspc] = len(species)
            species.append(xspc)
        occurrences = numpy.hstack([occurrences, xtr_occurrences])

    return {"sites": sites, "sites_dict": sites_dict,
            "species": species, "species_dict": species_dict,
            "slices": slices, "sliced_sites": sliced_sites,
            "occurrences": occurrences}


def extract_slices(fossils_data, subset_slices):
    map_slices_names = dict([(v[-1], k) for k, v in enumerate(fossils_data["slices"])])
    sub_slices = []
    sub_sliced = []
    subset_sites = set()
    for i in subset_slices:
        sub_slices.append(fossils_data["slices"][map_slices_names.get(i, i)])
        sub_sliced.append(fossils_data["sliced_sites"][map_slices_names.get(i, i)])
        subset_sites.update(sub_sliced[-1])

    ssites = sorted(subset_sites)
    sites_org_to_new = -numpy.ones(len(fossils_data["sites"]), dtype=int)
    for (k, v) in enumerate(ssites):
        sites_org_to_new[v] = k
    sub_sliced_sites = [[sites_org_to_new[i] for i in ss] for ss in sub_sliced]

    sub_sites = [fossils_data["sites"][i] for i in ssites]
    reverse_sites_dict = dict([(v, k) for (k, v) in fossils_data["sites_dict"].items()])
    sub_sites_dict = dict([(reverse_sites_dict[i], k) for k, i in enumerate(ssites)])

    sspecies = numpy.where(numpy.sum(fossils_data["occurrences"][ssites, :], axis=0))[0]
    species_org_to_new = -numpy.ones(len(fossils_data["species"]), dtype=int)
    for (k, v) in enumerate(sspecies):
        species_org_to_new[v] = k

    sub_species = [fossils_data["species"][i] for i in sspecies]
    reverse_species_dict = dict([(v, k) for (k, v) in fossils_data["species_dict"].items()])
    sub_species_dict = dict([(reverse_species_dict[i], k) for k, i in enumerate(sspecies)])

    xs, ys = numpy.where(fossils_data["occurrences"])
    sub_occurrences = numpy.zeros((len(ssites)+1, len(sspecies)+1), dtype=bool)
    sub_occurrences[sites_org_to_new[xs], species_org_to_new[ys]] = fossils_data["occurrences"][xs, ys]
    return {"sites": sub_sites, "sites_dict": sub_sites_dict,
            "species": sub_species, "species_dict": sub_species_dict,
            "slices": sub_slices, "sliced_sites": sub_sliced_sites,
            "occurrences": sub_occurrences[:-1, :-1]}


def write_fossils(fossils_data, list_file_out, fields_sites, fields_species):
    sites, species, occurrences = (fossils_data["sites"], fossils_data["species"], fossils_data["occurrences"])
    lines = [",".join(["%s" % v for v in sites[x]+species[y]]) for (x, y) in zip(*numpy.where(occurrences))]
    lines.sort()
    lines.insert(0, ",".join(fields_sites+fields_species))
    with open(list_file_out, "w") as fo:
        fo.write("\n".join(lines))


LONG_PREAMBLE = """\\documentclass[10pt]{article}

\\usepackage[utf8]{inputenc}
\\usepackage[margin=1.5cm]{geometry}

\\usepackage{booktabs}
\\usepackage{longtable}
%%\\input{../macros}
\\newcommand{\groupNA}{NA}
\\newcommand{\groupEU}{EU}
\\newcommand{\groupL}{L}
\\newcommand{\groupS}{S}
\\newcommand{\groupC}{C}

\\begin{document}

\\begin{longtable}{@{\\hspace{5ex}}rl@{\\hspace{.8em}}r@{--}l@{\\hspace{2.5em}}r@{~;~}r@{\\hspace{3.em}}r@{\\hspace{1.5em}}r@{/}r@{/}r@{\\hspace{5ex}}}
\\caption{%s} \\label{tab:%s} \\\\

\\toprule \\multicolumn{2}{l}{\\textbf{Locality}} & \\multicolumn{2}{l}{\\textbf{Age (Ma)}} & \\multicolumn{2}{l}{\\textbf{Coordinates}} & \\multicolumn{1}{l}{\\textbf{H}} & \\multicolumn{3}{l}{\\textbf{Nbs}} \\\\[.8em] 
\\endfirsthead

\\multicolumn{10}{c}%%
{{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}} \\\\
\\multicolumn{2}{l}{\\textbf{Locality}} & \\multicolumn{2}{l}{\\textbf{Age}} & \\multicolumn{2}{l}{\\textbf{Coordinates}} & \\multicolumn{1}{l}{\\textbf{H}} & \\multicolumn{3}{l}{\\textbf{Nbs}} \\\\[.8em]
\\endhead

%%\\multicolumn{11}{r}{{Continued on next page}} \\\\
%%\\endfoot

\\bottomrule
\\endlastfoot

"""
LONG_POSTAMBLE = "\\end{longtable}\n\\end{document}"  # \n%%% Local Variables:\n%%% mode: latex\n%%% TeX-master: t\n%%% End:"


def save_outline(fossils_data, fields_species, fields_sites, outline_file, tex_localities_file=None):

    whch, grp = (fossils_data["WHICH"], fossils_data["GROUP"])
    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])
    occurrences = fossils_data["occurrences"]

    ord_field = fields_species.index("ORDER")
    gen_field = fields_species.index("GENUS")
    spc_field = fields_species.index("SPECIES")
    ords = sorted(set([s[ord_field] for s in species]))
    ords_map = dict([(v, k) for (k, v) in enumerate(ords)])
    N_ords = len(ords)
    map_to_ord = numpy.array([ords_map[s[ord_field]] for s in species])

    gens = sorted(set(["%s_%s" % (s[ord_field], s[gen_field]) for s in species if not is_indeterminate_str(s[gen_field])]))
    gens_map = dict([(v, k) for (k, v) in enumerate(gens)])
    N_gens = len(gens)
    indet_gen_id = len(gens)
    map_to_gen = numpy.array([gens_map.get("%s_%s" % (s[ord_field], s[gen_field]), indet_gen_id) for s in species])

    gens_to_ords = numpy.array([ords_map[s.split("_")[0]] for s in gens])

    if tex_localities_file is not None:
        sites_details, sites_by_slices = ({}, {})
        site_id_field, site_slice_field = (fields_sites.index('LIDNUM'), fields_sites.index('SLICE_ID'))

    with open(outline_file, "w") as fo:

        fo.write(",".join(fields_sites + ["nb_orders", "nb_genera", "nb_spec", "nb_spec_all"] + ["%s_genera" % o for o in ords] + ords + gens + ["indet."]) + "\n")

        O = numpy.bincount(map_to_ord, minlength=N_ords)
        G = numpy.bincount(map_to_gen, minlength=N_gens+1)
        Og = numpy.bincount(gens_to_ords, minlength=N_ords)
        T = [numpy.sum(O > 0), numpy.sum(G[:-1] > 0), numpy.sum(G[:-1]), numpy.sum(G)]
        fo.write(",".join(["TOTALS"] + ["" for f in fields_sites[1:]] + ["%d" % v for v in T] + ["%d" % v for v in Og] + ["%d" % v for v in O] + ["%d" % v for v in G]) + "\n")

        for si, site in enumerate(sites):
            mask = numpy.where(occurrences[si, :])[0]

            O = numpy.bincount(map_to_ord[mask], minlength=N_ords)
            G = numpy.bincount(map_to_gen[mask], minlength=N_gens+1)
            Og = numpy.bincount(gens_to_ords[G[:-1] > 0], minlength=N_ords)
            T = [numpy.sum(O > 0), numpy.sum(G[:-1] > 0), numpy.sum(G[:-1]), numpy.sum(G)]
            fo.write(",".join(list(map(str, site)) + ["%d" % v for v in T] + ["%d" % v for v in Og] + ["%d" % v for v in O] + ["%d" % v for v in G]) + "\n")

            if tex_localities_file is not None:
                if site[site_slice_field] not in sites_by_slices:
                    sites_by_slices[site[site_slice_field]] = []
                sites_by_slices[site[site_slice_field]].append(site[site_id_field])
                sites_details[site[site_id_field]] = dict(zip(*[fields_sites + ["nb_orders", "nb_genera", "nb_spec", "nb_spec_all"], list(site)+T]))

    if tex_localities_file is not None:
        with open(tex_localities_file, "w") as fo:

            tlgd = "List of localities per time unit for \\group%s{}, \\group%s{}. The columns indicate the identifier and name of the locality, the age bounds, the coordinates (latitude; longitude), the average hypsodonty (H) as well as the number of distinct orders, of genera and of species recorded (Nbs)." % (whch, grp)
            tlbl = "tab:locs-%s-%s" % (whch, grp)
            fo.write(LONG_PREAMBLE % (tlgd, tlbl))

            for si, slck in enumerate(sorted(sites_by_slices.keys())):
                slc = slices[slck]
                nb_locs = len(sites_by_slices[slck])
                if si > 0:
                    fo.write("[1.2em]")
                fo.write(f"\n\\multicolumn{{2}}{{l}}{{\\qquad \\textbf{{ {slc[-1]} }} }} & {slc[0]:.2f} & {slc[1]:.2f} & \\multicolumn{{2}}{{l}}{{ {nb_locs} localities}} \\\\\n\\midrule\n")
                for sk in sorted(sites_by_slices[slck]):
                    site = sites_details[sk]
                    site["NAMEb"] = site['NAME'].replace("&", "\\&")
                    site["HYP"] = site['MEAN_HYPSODONTY'] if site['MEAN_HYPSODONTY'] != "\\N" else "-"
                    fo.write(f"{site['LIDNUM']} & {site['NAMEb']} & {site['MAX_AGE']:.2f} & {site['MIN_AGE']:.2f} & {site['LAT']:.3f} & {site['LONG']:.3f} & {site['HYP']} & {site['nb_orders']} & {site['nb_genera']} & {site['nb_spec']} \\\\\n")

            fo.write(LONG_POSTAMBLE)


def read_teeth_details(teeth_file):
    data = {}
    with open(teeth_file) as fp:
        for line in fp:
            parts = line.strip().split(",")
            data[tuple(parts[:-1])] = parts[-1]
    return data


def compute_context_stats(fossils_data, map_fields_species, map_fields_sites):

    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])
    occurrences_org = fossils_data["occurrences"]

    gen_field = map_fields_species["GENUS"]
    spc_field = map_fields_species["SPECIES"]
    gens = sorted(set([s[gen_field] for s in species if not is_indeterminate_str(s[gen_field])]))
    gens_map = dict([(v, k) for (k, v) in enumerate(gens)])
    gen_counts = numpy.zeros((len(gens)+1, len(slices)), dtype=int)

    c_stats = []
    for xi, slce_sites in enumerate(sliced_sites_all):

        hpys = numpy.array([float(sites[ssi][map_fields_sites["MEAN_HYPSODONTY"]]) for ssi in slce_sites if sites[ssi][map_fields_sites["MEAN_HYPSODONTY"]] != "\\N"])
        spc_counts = occurrences_org[slce_sites, :].sum(0)

        sliced_sids = numpy.where(spc_counts > 0)[0]
        sliced_species = [species[ssi] for ssi in sliced_sids]

        gen_counts[-1, xi] = len(slce_sites)
        for si in range(len(sliced_sids)):
            if sliced_species[si][gen_field] in gens_map:  # else it's indeterminate
                gen_counts[gens_map[sliced_species[si][gen_field]], xi] += 1

        # nb of species per genus
        gc = numpy.bincount(gen_counts[:-1, xi])[1:]
        gc_avg = numpy.dot(gc, numpy.arange(1, gc.shape[0]+1))/numpy.sum(gc)
        ngc = gc/numpy.sum(gc)
        ngc[gc == 0] = 1
        gc_H = numpy.sum(-ngc*numpy.log2(ngc))

        spc_gen_det = numpy.array([2*is_indeterminate_str(s[gen_field])+1*is_indeterminate_str(s[spc_field]) for s in sliced_species])

        c_stats.append((numpy.mean(hpys), gc_avg, gc_H, len(slce_sites), len(sliced_species), numpy.sum(spc_gen_det > 0)/spc_gen_det.shape[0]))
        lbls = ["avg_hyp", "avg_spc_gen", "ent_spc_gen", "nb_sites", "nb_species", "ratio_indet"]
        # lbls = ["mean hypsondonty", "mean nb of sp per genus", "entropy of species assignment to genera", "nb of localities", "tot nb of species", "ratio of indeterminate"]
    return numpy.array(c_stats).T, lbls


# RANDOMIZING
#####################################
RND_METHODS = {}

# shuffle species labels (within each slice)


def rnd_shuffle_species(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]

    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slce_sites in blocks:
        species_ids = numpy.where(numpy.sum(occurrences_org[slce_sites, :], axis=0) > 0)[0]
        permutation_ids_pairs = zip(species_ids, numpy.random.permutation(species_ids))
        for iorg, irnd in permutation_ids_pairs:
            occurrences[slce_sites, irnd] = occurrences_org[slce_sites, iorg]

    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["shuffle-species"] = rnd_shuffle_species

# shuffle species labels within each order (and slice)


def rnd_shuffle_species_order(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]
    ord_field = fields_species.index("ORDER")
    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slce_sites in blocks:
        groups = {}
        for i in numpy.where(numpy.sum(occurrences_org[slce_sites, :], axis=0) > 0)[0]:
            if fossils_data["species"][i][ord_field] not in groups:
                groups[fossils_data["species"][i][ord_field]] = []
            groups[fossils_data["species"][i][ord_field]].append(i)
        permutation_ids_pairs = []
        for g, species_ids in groups.items():
            permutation_ids_pairs.extend(zip(species_ids, numpy.random.permutation(species_ids)))
        for iorg, irnd in permutation_ids_pairs:
            occurrences[slce_sites, irnd] = occurrences_org[slce_sites, iorg]
    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["shuffle-species-order"] = rnd_shuffle_species_order

# shuffle species labels within each order (and slice)


def rnd_shuffle_species_indet(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]
    ord_field = fields_species.index("ORDER")
    gen_field = fields_species.index("GENUS")
    spc_field = fields_species.index("SPECIES")
    mask_indet_spc = numpy.array([is_indeterminate_str(s[spc_field]) for s in fossils_data["species"]], dtype=bool)
    mask_indet_gen = numpy.array([is_indeterminate_str(s[gen_field]) for s in fossils_data["species"]], dtype=bool)

    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slce_sites in blocks:
        present = numpy.sum(occurrences_org[slce_sites, :], axis=0) > 0
        permutation_ids_pairs = []
        for species_ids in [numpy.where(present & mask_indet_gen)[0], numpy.where(present & ~mask_indet_gen & mask_indet_spc)[0], numpy.where(present & ~mask_indet_spc)[0]]:
            permutation_ids_pairs.extend(zip(species_ids, numpy.random.permutation(species_ids)))
        for iorg, irnd in permutation_ids_pairs:
            occurrences[slce_sites, irnd] = occurrences_org[slce_sites, iorg]
    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["shuffle-species-indet"] = rnd_shuffle_species_indet


# re-sample occurrence sites for species (within each slice)
def rnd_sample_sites(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]
    gen_field = fields_species.index("GENUS")
    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slce_sites in blocks:
        cnts = numpy.sum(occurrences_org[slce_sites, :], axis=0)
        species_ids = numpy.where(cnts > 0)[0]

        groups = {}
        if "effect_genus_same" in rnd_parameters or "effect_genus_other" in rnd_parameters:
            for i in species_ids:
                if fossils_data["species"][i][gen_field] not in groups:
                    groups[fossils_data["species"][i][gen_field]] = []
                groups[fossils_data["species"][i][gen_field]].append(i)

        if rnd_parameters.get("sort_species") == "+":
            species_ids = sorted(species_ids, key=lambda x: cnts[x])
        elif rnd_parameters.get("sort_species") == "-":
            species_ids = sorted(species_ids, key=lambda x: -cnts[x])
        else:
            numpy.random.shuffle(species_ids)

        for species_id in species_ids:
            p_elems = []
            cnts_all = numpy.sum(occurrences[slce_sites, :], axis=1)
            if len(groups) > 0 and numpy.sum(cnts_all) > 0:
                cnts_same = numpy.sum(occurrences[slce_sites, :][:, groups[fossils_data["species"][species_id][gen_field]]], axis=1)
                if "effect_genus_same" in rnd_parameters and numpy.sum(cnts_same) > 0:
                    mask_same = cnts_same > 0
                    dval_same = numpy.bincount(mask_same, minlength=2)
                    w_same_0 = 1/(dval_same[0]+dval_same[1]*rnd_parameters["effect_genus_same"])
                    w_same_1 = rnd_parameters["effect_genus_same"]*w_same_0
                    p_elems.append(~mask_same*w_same_0 + mask_same*w_same_1)
                if "effect_genus_other" in rnd_parameters:
                    mask_other = (cnts_all - cnts_same) > 0
                    dval_other = numpy.bincount(mask_other, minlength=2)
                    w_other_0 = 1/(dval_other[0]+dval_other[1]*rnd_parameters["effect_genus_other"])
                    w_other_1 = rnd_parameters["effect_genus_other"]*w_other_0
                    p_elems.append(~mask_other*w_other_0 + mask_other*w_other_1)

            if rnd_parameters.get("effect_marg_sites", 1) != 1:
                potential = (rnd_parameters["marg_sites"][slce_sites] - cnts_all)
                w_potential = float(rnd_parameters["effect_marg_sites"])**potential
                p_elems.append(w_potential/numpy.sum(w_potential))

            if len(p_elems) > 0:
                p = numpy.prod(numpy.array(p_elems), axis=0)
                p /= numpy.sum(p)
                # print(numpy.min(p), numpy.max(p))
            else:
                p = None
            # rnd_sites_ids = numpy.random.choice(numpy.arange(len(slce_sites)), cnts[species_id], replace=False, p=p)
            # rnd_sites = slce_sites[rnd_sites_ids]
            rnd_sites = numpy.random.choice(slce_sites, cnts[species_id], replace=False, p=p)
            occurrences[rnd_sites, species_id] = 1
    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["sample-sites"] = rnd_sample_sites

# re-sample species for site (within each slice)


def rnd_sample_species(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]
    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slce_sites in blocks:
        species_ids = numpy.where(numpy.sum(occurrences_org[slce_sites, :], axis=0) > 0)[0]
        cnts = numpy.sum(occurrences_org[slce_sites, :], axis=1)
        for sid in sorted(numpy.arange(len(slce_sites)), key=lambda x: -cnts[x]):
            rnd_species = numpy.random.choice(species_ids, cnts[sid], replace=False)
            occurrences[slce_sites[sid], rnd_species] = 1
    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["sample-species"] = rnd_sample_species


def rnd_null_CB(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]
    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slce_sites in blocks:
        species_ids = numpy.where(numpy.sum(occurrences_org[slce_sites, :], axis=0) > 0)[0]
        M = numpy.array(occurrences_org[slce_sites, :][:, species_ids], dtype=float)
        Mrnd = CB(M, rnd_parameters.get("nb_iter", -1))
        for si, ss in enumerate(slce_sites):
            occurrences[ss, species_ids] = (Mrnd[si, :] == 1)
    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["null-CB"] = rnd_null_CB


def rnd_null_UG(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    MAX_NB_TRIES = 5
    MAX_BUDGET_LOOP = 1000
    rnd_details = {}
    if "seed" in rnd_parameters:
        numpy.random.seed(rnd_parameters["seed"])
        rnd_details["seed"] = rnd_parameters["seed"]
    occurrences_org = fossils_data["occurrences"]
    occurrences = numpy.zeros((len(fossils_data["sites"]), len(fossils_data["species"])), dtype=bool)

    if rnd_parameters.get("sliced", True):
        blocks = fossils_data["sliced_sites"]
    else:
        blocks = [numpy.arange(len(fossils_data["sites"]))]

    for slci, slce_sites in enumerate(blocks):
        species_ids = numpy.where(numpy.sum(occurrences_org[slce_sites, :], axis=0) > 0)[0]
        M = numpy.array(occurrences_org[slce_sites, :][:, species_ids], dtype=float)
        nb_tries = 0
        Mrnd = UG_mod(M, max_budget_loop=MAX_BUDGET_LOOP)
        while nb_tries < MAX_NB_TRIES and Mrnd is None:
            nb_tries += 1
            print("Try again to randomize UG (%d)" % nb_tries)
            Mrnd = UG_mod(M, max_budget_loop=MAX_BUDGET_LOOP)
            if Mrnd is not None:
                print("Done with success")
        if Mrnd is not None:
            for si, ss in enumerate(slce_sites):
                occurrences[ss, species_ids] = (Mrnd[si, :] == 1)
        else:
            print("Randomization failed !", slci)
            for si, ss in enumerate(slce_sites):
                occurrences[ss, species_ids] = (M[si, :] == 1)

    return {"sites": fossils_data["sites"], "sites_dict": fossils_data["sites_dict"],
            "species": fossils_data["species"], "species_dict": fossils_data["species_dict"],
            "slices": fossils_data["slices"], "sliced_sites": fossils_data["sliced_sites"],
            "occurrences": occurrences}, rnd_details


RND_METHODS["null-UG"] = rnd_null_UG


PARAMS_FORMAT = {"i": ("nb_iter", int),
                 "x": ("sliced", "+"),
                 "z": ("sort_species", None),
                 "S": ("effect_genus_same", float),
                 "O": ("effect_genus_other", float),
                 "M": ("effect_marg_sites", float, "marg_sites")}


def make_rnd_series(key):
    params = {}
    parts = key.split("_")
    if parts[0] in RND_METHODS:
        params["method"] = parts.pop(0)
    p_str = "_".join(parts)
    for m in re.finditer("(?P<k>[a-zA-Z]+)(?P<val>[0-9\.\-\+]+)", p_str):
        if m.group("k") in PARAMS_FORMAT:
            k = m.group("k")
            if PARAMS_FORMAT[k][1] is None:
                params[PARAMS_FORMAT[k][0]] = m.group("val")
            elif type(PARAMS_FORMAT[k][1]) is type:
                params[PARAMS_FORMAT[k][0]] = PARAMS_FORMAT[k][1](m.group("val"))
            else:
                params[PARAMS_FORMAT[k][0]] = (PARAMS_FORMAT[k][1] == m.group("val"))

            if len(PARAMS_FORMAT[k]) > 2:
                params[PARAMS_FORMAT[k][2]] = None
    return params


def randomize_records(fossils_data, fields_sites, fields_species, rnd_parameters={}):
    if rnd_parameters.get("method") in RND_METHODS:
        return RND_METHODS[rnd_parameters.get("method")](fossils_data, fields_sites, fields_species, rnd_parameters)
    return fossils_data, None

# COMPUTING SCORES (+storage)
#####################################

# COMPUTE DIVERSITY COUNTS


def compute_diverse(fossils_data, fields_species, map_fields_sites):

    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])
    occurrences = fossils_data["occurrences"]

    gen_field = fields_species.index("GENUS")
    spc_field = fields_species.index("SPECIES")
    gens = sorted(set([s[gen_field] for s in species if not is_indeterminate_str(s[gen_field])]))
    gens_map = dict([(v, k) for (k, v) in enumerate(gens)])
    gen_counts = numpy.zeros((len(gens)+1, len(slices)))
    sliced_divs = []
    for xi, s in enumerate(slices):
        sliced_divs.append([])
        sliced_sites = sliced_sites_all[xi]

        if len(sliced_sites) > 0:
            spc_counts = occurrences[sliced_sites, :].sum(0)

            sliced_sids = numpy.where(spc_counts > 0)[0]
            sliced_species = [species[ssi] for ssi in sliced_sids]

            gen_counts[-1, xi] = len(sliced_sites)
            for si in range(len(sliced_sids)):
                # if sliced_species[si][gen_field] in gens_map: ### else it's indeterminate
                if sliced_species[si][gen_field] in gens_map and (GEN_WITH_INDET or not is_indeterminate_gs(sliced_species[si][gen_field], sliced_species[si][spc_field])):
                    gen_counts[gens_map[sliced_species[si][gen_field]], xi] += 1

            # go over each site in the time slice
            for site in sliced_sites:
                sids = sorted(numpy.where(occurrences[site, sliced_sids])[0])

                gen_subsets = {}
                for si in sids:
                    # if sliced_species[si][gen_field] in gens_map: ### else it's indeterminate
                    if sliced_species[si][gen_field] in gens_map and (GEN_WITH_INDET or not is_indeterminate_gs(sliced_species[si][gen_field], sliced_species[si][spc_field])):
                        genk = gens_map[sliced_species[si][gen_field]]
                        if genk not in gen_subsets:
                            gen_subsets[genk] = []
                        gen_subsets[genk].append(si)
                for genk, subset in gen_subsets.items():
                    cdet = nb_distinct_species(subset, sliced_species, gen_field, spc_field)
                    sliced_divs[-1].append((genk, site, cdet))

    return gens, gen_counts, sliced_divs


def save_divs(gens, gen_counts, sliced_divs, gen_counts_file, sliced_divs_file):
    with open(gen_counts_file, "w") as fo:
        num_fmt = ",".join(gen_counts.shape[1]*["%d"])
        fo.write("GENUS,"+(num_fmt % tuple(gen_counts[-1, :]))+"\n")
        for si, s in enumerate(gen_counts[:-1, :]):
            fo.write(gens[si]+","+(num_fmt % tuple(s))+"\n")
    with open(sliced_divs_file, "w") as fo:
        for xi, vs in enumerate(sliced_divs):
            for (genk, site, cdet) in vs:
                fo.write("%d,%d,%d,%d\n" % (xi, genk, site, cdet))
            # for (genk, site, cdet, cmin, cmax) in vs:
            #     fo.write("%d,%d,%d,%d,%d,%d\n" % (xi, genk, site, cdet, cmin, cmax))


def load_divs(gen_counts_file, sliced_divs_file):
    gens = []
    gen_counts = []
    with open(gen_counts_file) as fp:
        for line in fp:
            parts = list(line.strip().split(","))
            if len(gen_counts) > 0:
                gens.append(parts[0])
            gen_counts.append(tuple(map(int, parts[1:])))
    gen_counts = numpy.array(gen_counts[1:]+[gen_counts[0]])
    sliced_divs = [[] for i in range(gen_counts.shape[1])]
    with open(sliced_divs_file) as fp:
        for line in fp:
            parts = list(map(int, line.strip().split(",")))
            sliced_divs[parts[0]].append(tuple(parts[1:]))
    return gens, gen_counts, sliced_divs


def get_counts(occurrences, sliced_sites, det_mask, ord_mask=None):
    collect = []
    for det_only in [True, False]:
        collect.append([])
        for axis in [1, 0]:
            if det_only:
                if ord_mask is not None:
                    M = occurrences[sliced_sites, :][:, det_mask & ord_mask]
                else:
                    M = occurrences[sliced_sites, :][:, det_mask]
            else:
                if ord_mask is not None:
                    M = occurrences[sliced_sites, :][:, ord_mask]
                else:
                    M = occurrences[sliced_sites, :]

            all_counts = M.sum(axis)
            nnz_counts = all_counts[all_counts > 0]
            if len(nnz_counts) > 0:
                collect[-1].extend([nnz_counts.shape[0], numpy.sum(nnz_counts == 1), numpy.median(nnz_counts), numpy.max(nnz_counts)])  # , numpy.sum(nnz_counts)
            else:
                collect[-1].extend([0, 0, 0, 0])
            if axis == 0:
                collect[-1].append(numpy.sum(nnz_counts))
    return numpy.array(collect).T.flatten()


def table_counts(fossils_data, map_fields_species, table_counts=None):

    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])
    occurrences = fossils_data["occurrences"]

    gen_field = map_fields_species["GENUS"]
    spc_field = map_fields_species["SPECIES"]
    fam_field = map_fields_species["FAMILY"]
    ord_field = map_fields_species["ORDER"]

    indeterminate_mask = numpy.array([is_indeterminate_gs(s[gen_field], s[spc_field]) for spci, s in enumerate(species)], dtype=bool)
    ocnt = {}
    for xi, s in enumerate(slices):
        sliced_sites = sliced_sites_all[xi]

        if len(sliced_sites) > 0:
            ord_gen = {}  # dict: for each order, set of genera
            gen_spc = {}  # dict: for each order, then genus, set of species
            cnt_gen = {}  # dict: for each order, then genus, count occurrences
            ord_masks = {}
            spc_counts = occurrences[sliced_sites, :].sum(0)
            for i in numpy.where(spc_counts > 0)[0]:
                gkey = tuple(species[i][ord_field+1:gen_field+1])
                if species[i][ord_field] not in ord_gen:
                    ord_gen[species[i][ord_field]] = [set(), set()]
                    ord_masks[species[i][ord_field]] = numpy.zeros(occurrences.shape[1], dtype=bool)
                ord_gen[species[i][ord_field]][indeterminate_mask[i]].add(gkey)
                ord_masks[species[i][ord_field]][i] = 1

                if species[i][ord_field] not in gen_spc:
                    gen_spc[species[i][ord_field]] = {}
                    cnt_gen[species[i][ord_field]] = {}
                if gkey not in gen_spc[species[i][ord_field]]:
                    gen_spc[species[i][ord_field]][gkey] = [set(), set()]
                    cnt_gen[species[i][ord_field]][gkey] = [0, 0]
                gen_spc[species[i][ord_field]][gkey][indeterminate_mask[i]].add(species[i][spc_field])
                cnt_gen[species[i][ord_field]][gkey][indeterminate_mask[i]] += spc_counts[i]

            for k, lst in ord_gen.items():
                tmp_spc = numpy.array([(len(l[0]), len(l[0])+len(l[1])) for kk, l in gen_spc[k].items()])
                tmp_occ = numpy.array([(l[0], l[0]+l[1]) for kk, l in cnt_gen[k].items()])
                tmp = (len(lst[0]), len(lst[1] | lst[0]),
                       numpy.median(tmp_spc[:, 0]), numpy.median(tmp_spc[:, 1]),
                       numpy.max(tmp_spc[:, 0]), numpy.max(tmp_spc[:, 1]),
                       numpy.median(tmp_occ[:, 0]), numpy.median(tmp_occ[:, 1]),
                       numpy.max(tmp_occ[:, 0]), numpy.max(tmp_occ[:, 1]))

                ocnt[(xi, k)] = numpy.concatenate([get_counts(occurrences, sliced_sites, ~indeterminate_mask, ord_masks[k]), tmp])

            ocnt[(xi, "-- %s" % s[-1])] = get_counts(occurrences, sliced_sites, ~indeterminate_mask)

    # grep "\-\-" stats.txt | sed 's/.*-- //'
    str_all_h = "\\begin{table}\n\\small"
    str_all_h += "\\begin{tabular}{@{\\hspace{1ex}}l@{\\hspace{4ex}}"+4*"r@{}r@{/}r@{}r@{\\hspace{3ex}}r@{}r@{/}r@{}r@{\\hspace{5ex}}"+"r@{}r@{/}r@{}r"+"@{\\hspace{1ex}}}\n\\toprule\n"
    str_all_h += " & ".join(["", "\\multicolumn{8}{c}{sites}", "\\multicolumn{8}{c}{occs/site}",
                             "\\multicolumn{8}{c}{taxa}", "\\multicolumn{8}{c}{occs/taxon}", "\\multicolumn{4}{c}{occs}"])+"\\\\\n"
    str_all_h += " & ".join([""]+["\\multicolumn{4}{l}{%s}" % v for v in ["total", "single", "med", "max", "total", "single", "med", "max", "total"]])+"\\\\\n"
    str_all_h += "\\midrule\n"

    # str_all_h = "\\begin{table}\n\\small"
    # str_all_h += "\\begin{tabular}{@{\\hspace{1ex}}l@{\\hspace{4ex}}"+4*"r@{/}r@{\\hspace{3ex}}r@{/}r@{\\hspace{5ex}}"+"r@{/}r"+"@{\\hspace{1ex}}}\n\\toprule\n"
    # str_all_h += " & ".join(["", "\\multicolumn{4}{c}{sites}", "\\multicolumn{4}{c}{occs/site}",
    #                            "\\multicolumn{4}{c}{taxa}", "\\multicolumn{4}{c}{occs/taxon}", "\\multicolumn{2}{c}{occs}"])+"\\\\\n"
    # str_all_h += " & ".join([""]+["\\multicolumn{2}{l}{%s}" % v for v in ["total", "single", "med", "max", "total", "single", "med", "max", "total"]])+"\\\\\n"
    # str_all_h += "\\midrule\n"

    str_ord_h = "\\begin{table}\n\\small"
    str_ord_h += "\\begin{tabular}{@{\\hspace{1ex}}l@{\\hspace{4ex}}"+4*"r@{/}r@{\\hspace{3ex}}r@{/}r@{\\hspace{5ex}}"+"r@{/}r@{\\hspace{5ex}}" + \
        "r@{/}r@{\\hspace{5ex}}"+"r@{/}r@{\\hspace{3ex}}r@{/}r@{\\hspace{5ex}}"+"r@{/}r@{\\hspace{3ex}}r@{/}r"+"@{\\hspace{1ex}}}\n\\toprule\n"
    str_ord_h += " & ".join(["", "\\multicolumn{4}{c}{sites}", "\\multicolumn{4}{c}{occs/site}",
                             "\\multicolumn{4}{c}{taxa}", "\\multicolumn{4}{c}{occs/taxon}", "\\multicolumn{2}{c}{occs}",
                             "\\multicolumn{2}{c}{gen/ord}", "\\multicolumn{4}{c}{taxa/gen}", "\\multicolumn{4}{c}{occs/gen}"])+"\\\\\n"
    str_ord_h += " & ".join([""]+["\\multicolumn{2}{l}{%s}" % v for v in ["total", "single", "med", "max",
                                                                          "total", "single", "med", "max", "total",
                                                                          "", "med", "max", "med", "max"]])+"\\\\\n"
    str_ord_h += "\\midrule\n"
    str_all = ""
    str_ord = ""
    patt_fmt = {0: "& $%d$", 1: "$%d$ &"}
    for k, vs in ocnt.items():
        if re.match("--", k[1]):
            # str_all += " & ".join(["(%d) %s" % (k[0], k[1])]+["$%d$" % v for v in vs])+" \\\\\n"
            str_all += " & ".join([k[1].strip("-- ")]+[patt_fmt[(vi % 2)] % v for (vi, v) in enumerate(vs)])+" \\\\\n"
            str_ord += "[.5em] \n"
        else:
            str_ord += " & ".join(["(%d) %s" % (k[0], k[1])]+["$%d$" % v for v in vs])+" \\\\\n"
    str_all_t = "\\bottomrule\n\\end{tabular}\n\\end{table}"
    str_ord_t = "\\bottomrule\n\\end{tabular}\n\\end{table}"

    if table_counts is None:
        print(str_all)
        print(str_ord)
    else:
        with open(table_counts, "w") as fo:
            # fo.write("%%% Overall counts\n")
            # fo.write(str_all_h)
            fo.write(str_all)
            # fo.write(str_all_t)
            # fo.write("\n\n%%% Counts by orders\n")
            # fo.write(str_ord_h)
            # fo.write(str_ord)
            # fo.write(str_ord_t)

# COMPUTE TRIPLETS COUNTS


def compute_counts(fossils_data, fields_species, map_fields_sites):

    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])
    occurrences = fossils_data["occurrences"]

    sites_counts = numpy.array([len(sliced_sites_all[xi]) for xi, s in enumerate(slices)])
    spc_counts = numpy.zeros((len(species)+1, len(slices)))
    sliced_counts = []
    for xi, s in enumerate(slices):
        sliced_counts.append({})
        # print("Slice", s)
        sliced_sites = sliced_sites_all[xi]

        if len(sliced_sites) > 0:
            spc_counts[:-1, xi] = occurrences[sliced_sites, :].sum(0)
            spc_counts[-1, xi] = len(sliced_sites)
            # print("Nb species", numpy.sum(spc_counts[:-1, xi] > 0))

            # go over each site in the time slice, and increment co-occ count for each species pair
            for site in sliced_sites:
                sids = sorted(numpy.where(occurrences[site, :])[0])

                for i in range(len(sids)):
                    for j in range(i):
                        sliced_counts[-1][(sids[j], sids[i])] = sliced_counts[-1].get((sids[j], sids[i]), 0) + 1

    return spc_counts, sliced_counts, sites_counts


def save_counts(species, spc_fields, spc_counts, sliced_counts, spc_counts_file, sliced_counts_file):
    with open(spc_counts_file, "w") as fo:
        txt_fmt = ",".join(len(spc_fields)*["%s"])
        num_fmt = ",".join(spc_counts.shape[1]*["%d"])
        fo.write((txt_fmt % tuple(spc_fields))+","+(num_fmt % tuple(spc_counts[-1, :]))+"\n")
        for si, s in enumerate(spc_counts[:-1, :]):
            fo.write((txt_fmt % tuple(species[si]))+","+(num_fmt % tuple(s))+"\n")
    with open(sliced_counts_file, "w") as fo:
        for xi, vs in enumerate(sliced_counts):
            for k, v in vs.items():
                fo.write("%d,%d,%d,%d\n" % (xi, k[0], k[1], v))


def load_counts(spc_counts_file, sliced_counts_file):
    species = []
    spc_fields = None
    spc_counts = []
    with open(spc_counts_file) as fp:
        for line in fp:
            parts = list(line.strip().split(","))
            if spc_fields is None:
                i = 0
                while i < len(parts) and spc_fields is None:
                    try:
                        v = int(parts[i])
                    except ValueError:
                        i += 1
                    else:
                        spc_fields = parts[:i]
                spc_counts.append(tuple(map(int, parts[i:])))
            else:
                species.append(tuple(parts[:i]))
                spc_counts.append(tuple(map(int, parts[i:])))
    spc_counts = numpy.array(spc_counts[1:]+[spc_counts[0]])
    sliced_counts = [{} for i in range(spc_counts.shape[1])]
    with open(sliced_counts_file) as fp:
        for line in fp:
            parts = list(map(int, line.strip().split(",")))
            sliced_counts[parts[0]][(parts[1], parts[2])] = parts[3]
    return species, spc_fields, spc_counts, sliced_counts


def fisher_test_values(N, spc_counts, species, genus_field, species_field):
    # N = sites_counts[xi]
    scores_tmp = []
    sids = sorted(numpy.where(spc_counts[:-1, xi] > 0)[0])
    for i in range(len(sids)):
        for j in range(i):
            pair_st = pair_status(sids[j], sids[i], species, genus_field, species_field)
            cij, ci, cj = (pairs_counts.get((sids[j], sids[i]), 0), spc_counts[sids[i], xi], spc_counts[sids[j], xi])
            cmean = scipy.stats.hypergeom.mean(N, ci, cj)
            cstd = scipy.stats.hypergeom.std(N, ci, cj)
            pmf = scipy.stats.hypergeom.pmf(cij, N, ci, cj)
            cdf = scipy.stats.hypergeom.cdf(cij, N, ci, cj)
            obs, smaller, larger = (pmf, cdf-pmf, 1-cdf)
            zval = (cij-cmean)/cstd

            #     stts = "%.4f,%.4f,%.4f %d:%d+%d/%d (%.3f/%.3f=%.3f)" % (smaller, obs, larger, cij, ci, cj, N, cmean, cstd, zval)
            #     print(f"{out}\t{xi}\t{pair_st} {species[sids[j]][genus_field]} {species[sids[j]][species_field]}\t{species[sids[i]][genus_field]} {species[sids[i]][species_field]} {stts}")
            scores_tmp.append((pair_st, smaller, obs, larger, cij, ci, cj, N, cmean, cstd, zval, sids[j], sids[i], xi))
    return scores_tmp


CDF = {}
PMF = {}
CCS = {}
C = 0


def pair_fisher(N, ci, cj, cij):
    # global C
    if N not in CDF:
        CDF[N] = -numpy.ones((N, N, N))
        PMF[N] = -numpy.ones((N, N, N))
        CCS[N] = numpy.zeros((N, N, N))
    cmin, cmax = (int(min(ci, cj)), int(max(ci, cj)))
    cinter = int(cij)
    if CDF[N][cmin, cmax, cinter] == -1:
        # C = C + 1
        CDF[N][cmin, cmax, cinter] = scipy.stats.hypergeom.cdf(cij, N, ci, cj)
        PMF[N][cmin, cmax, cinter] = scipy.stats.hypergeom.pmf(cij, N, ci, cj)
    CCS[N][cmin, cmax, cinter] += 1
    return CDF[N][cmin, cmax, cinter]-0.5*PMF[N][cmin, cmax, cinter]
    # C = C + 1
    # return scipy.stats.hypergeom.cdf(cij, N, ci, cj) - 0.5*scipy.stats.hypergeom.pmf(cij, N, ci, cj)
# C-score  (c_i - c_ij)(c_j - c_ij)


def pair_cscore(N, ci, cj, cij):
    return (ci - cij)*(cj - cij)
# Jaccard |S_i inter S_j|/|S_i union S_j| = c_ij/(c_i+c_j-c_ij)


def pair_jacc(N, ci, cj, cij):
    if ci == 0 or cj == 0:
        return 0
    return cij/(ci + cj - cij)
# min |S_i inter S_j|/min(|S_i|, |S_j|) = c_ij/min(c_i,c_j)


def pair_min(N, ci, cj, cij):
    if ci == 0 or cj == 0:
        return 0
    return cij/min(ci, cj)


MEASURES = {"fisher": pair_fisher, "cscore": pair_cscore, "jacc": pair_jacc, "min": pair_min}


def compute_pairs_scores(species, spc_fields, spc_counts, sliced_counts, sites_counts, measures):
    # tic = datetime.datetime.now()
    # global C
    # C = 0
    genus_field = spc_fields.index("GENUS")
    species_field = spc_fields.index("SPECIES")
    order_field = spc_fields.index("ORDER")
    ord_counts = {}
    for s in species:
        if s[order_field] not in ord_counts:
            ord_counts[s[order_field]] = 0
        ord_counts[s[order_field]] += 1
    orders_to_num = dict([(v, i) for (i, v) in enumerate(sorted(ord_counts.keys(), key=lambda x: -ord_counts[x]))])
    scores = []
    tots = []
    for xi, pairs_counts in enumerate(sliced_counts):
        N = sites_counts[xi]
        scores_tmp = []
        sids = sorted(numpy.where(spc_counts[:-1, xi] > 0)[0])
        for i in range(len(sids)):
            for j in range(i):
                cij, ci, cj = (pairs_counts.get((sids[j], sids[i]), 0), spc_counts[sids[i], xi], spc_counts[sids[j], xi])
                pair_st = pair_status(sids[j], sids[i], species, genus_field, species_field)
                pair_gA = pairA_group(sids[j], sids[i], species, order_field, orders_to_num)
                pair_gB = pairB_group(sids[j], sids[i], species, order_field, orders_to_num)
                scores_tmp.append([pair_st, pair_gA, pair_gB]+[MEASURES[m](N, ci, cj, cij) for m in measures]+[N, ci, cj, cij, sids[j], sids[i]])
        scores.append(numpy.array(scores_tmp))

        # ONLY COOCURRING PAIRS
        # scores.append(numpy.array([[pair_status(k[0], k[1], species, genus_field, species_field)]+[MEASURES[m](N, spc_counts[k[0], xi], spc_counts[k[1], xi], v) for m in measures] for k, v in pairs_counts.items()]))

        # tots.append([0,0,0])
        # spcs = list(numpy.where(spc_counts[:-1,xi] > 0)[0]) ### last is number of sites
        # for i in range(len(spcs)):
        #     for j in range(i):
        #         tots[-1][pair_status(spcs[i], spcs[j], species, genus_field, species_field)] += 1
    # print("Pairs score [%d], elapsed %s" % (C, datetime.datetime.now() - tic))
    return scores, tots


def get_pair_statusvars():
    return ["pair_status", "gA_status", "gB_status"]


def get_pair_extravars():
    return ["N", "ci", "cj", "cij", "sidj", "sidi"]


CMP_FUN = {"gt": numpy.greater, "lt": numpy.less, "ge": numpy.greater_equal, "le": numpy.less_equal, "ne": numpy.not_equal, "eq": numpy.equal}
VAL_FUN = {"avg": numpy.mean, "std": numpy.std, "var": numpy.var, "qunt": numpy.percentile, "len": len}


def eval_fun(c_fun, x, c_arg):
    if len(x) == 0:
        return numpy.nan
    ratio = False
    if re.match("ratio_", c_fun):
        ratio = True
        c_fun = c_fun.lstrip("ratio_")
        tot = len(x)
    if c_fun in VAL_FUN:
        if c_arg is None:
            v = VAL_FUN[c_fun](x)
        elif type(c_arg) is dict:
            v = VAL_FUN[c_fun](x, **c_arg)
        else:
            v = VAL_FUN[c_fun](x, c_arg)
    if c_fun in CMP_FUN:
        if c_arg is None:
            v = numpy.sum(CMP_FUN[c_fun](x))
        elif type(c_arg) is dict:
            v = numpy.sum(CMP_FUN[c_fun](x, **c_arg))
        else:
            v = numpy.sum(CMP_FUN[c_fun](x, c_arg))
    if ratio:
        v = v/tot
    return v


def compute_pairs_stats(slices, spc_counts, scores, tots, stats_details, grp=None):
    ID_STATUS = 0
    ID_GA = 1
    ID_GB = 2
    collect = []
    for si, slce in enumerate(slices):
        collect.append([])

        X = scores[si]

        mask_all = (X[:, ID_STATUS] > -1)
        mask_same_gen = (X[:, ID_STATUS] == 1)
        mask_diff_gen = (X[:, ID_STATUS] == 0)
        # mask_diff = (X[:,ID_STATUS]==0)
        # mask_same_ord = (X[:,ID_GA] == X[:,ID_GB])
        if grp is not None and grp == "L":
            mask_ordAm = (X[:, ID_GA] == 0) | (X[:, ID_GA] == 1)
            mask_ordBm = (X[:, ID_GB] == 0) | (X[:, ID_GB] == 1)
            mask_ordAo = (X[:, ID_GA] > 1)
            mask_ordBo = (X[:, ID_GB] > 1)
        else:
            mask_ordAm = (X[:, ID_GA] == 0)
            mask_ordBm = (X[:, ID_GB] == 0)
            mask_ordAo = (X[:, ID_GA] > 0)
            mask_ordBo = (X[:, ID_GB] > 0)

        Mmm = (mask_ordAm & mask_ordBm)
        Moo = (mask_ordAo & mask_ordBo)
        Mmo = (mask_ordAm & mask_ordBo) | (mask_ordAo & mask_ordBm)

        masks = [mask_same_gen, mask_diff_gen, mask_all,
                 Mmm & mask_same_gen, Mmm & mask_diff_gen, Mmm & mask_all,
                 Moo & mask_same_gen, Moo & mask_diff_gen, Moo & mask_all,
                 Mmo & mask_all]

        if grp is not None and grp == "L":
            mask_ordA0 = (X[:, ID_GA] == 0)
            mask_ordB0 = (X[:, ID_GB] == 0)
            mask_ordA1 = (X[:, ID_GA] == 1)
            mask_ordB1 = (X[:, ID_GB] == 1)

            M00 = (mask_ordA0 & mask_ordB0)
            M11 = (mask_ordA1 & mask_ordB1)
            M01 = (mask_ordA0 & mask_ordB1) | (mask_ordA1 & mask_ordB0)

            masks.extend([
                M00 & mask_same_gen, M00 & mask_diff_gen, M00 & mask_all,
                M11 & mask_same_gen, M11 & mask_diff_gen, M11 & mask_all,
                M01 & mask_all])

        for stat in stats_details:
            id_score, c_fun, c_arg = stat[:3]
            collect[-1].extend([eval_fun(c_fun, X[msk, id_score], c_arg) for msk in masks])
    return collect


def get_agg_series(grp=None):
    if grp is not None and grp == "L":
        return ["same", "diff", "all", "ordMsame", "ordMdiff", "ordMall", "ordOsame", "ordOdiff", "ordOall", "ordMxO",
                "ord0same", "ord0diff", "ord0all", "ord1same", "ord1diff", "ord1all", "ord0x1"]
    else:
        return ["same", "diff", "all", "ordMsame", "ordMdiff", "ordMall", "ordOsame", "ordOdiff", "ordOall", "ordMxO"]


def compute_marg_stat(fossils_data, fossils_data_org=None):
    collect = []
    for si, ss in enumerate(fossils_data["sliced_sites"]):
        if fossils_data_org is not None:
            mask_spc = numpy.sum(fossils_data_org["occurrences"][ss, :], axis=0) > 0
            dmargR = numpy.sum(fossils_data_org["occurrences"][ss, :][:, mask_spc], axis=1) - numpy.sum(fossils_data["occurrences"][ss, :][:, mask_spc], axis=1)
            dmargC = numpy.sum(fossils_data_org["occurrences"][ss, :][:, mask_spc], axis=0) - numpy.sum(fossils_data["occurrences"][ss, :][:, mask_spc], axis=0)

            collect.append((numpy.mean(dmargC**2), numpy.mean(dmargC), numpy.median(dmargC), numpy.std(dmargC),
                            numpy.mean(dmargR**2), numpy.mean(dmargR), numpy.median(dmargR), numpy.std(dmargR)))
        else:
            collect.append((0, 0, 0, 0, 0, 0, 0, 0))
    return collect


def get_marg_details():
    return [(0, "mrgCsqavg", "Avg. sq. col margins discrepancy", "%.3f", True),
            (1, "mrgCavg", "Average col margins discrepancy", "%.3f", False),
            (2, "mrgCmed", "Median col margins discrepancy", "%.3f", False),
            (3, "mrgCstd", "Stdev col margins discrepancy", "%.3f", False),
            (4, "mrgRsqavg", "Avg. sq. row margins discrepancy", "%.3f", False),
            (5, "mrgRavg", "Average row margins discrepancy", "%.3f", False),
            (6, "mrgRmed", "Median row margins discrepancy", "%.3f", False),
            (7, "mrgRstd", "Stdev row margins discrepancy", "%.3f", False)]


def compute_div_stat(which, slices, gens, gen_counts, sliced_divs):
    eps = 10**-8

    obins = [0, 1, 2, 3, 6, 4000]

    bins = [0, 1, 2, 3, 6, 4000]
    fbins = [0, eps, 0.1, 0.25, 0.5, 0.66, 1-eps, 1]
    if which == "NA":
        cbins = [0, 1, 2, 4, 8, 16, 25, 100]
    else:
        cbins = [0, 1, 2, 5, 10, 25, 50, 100]

    ID_GENK, ID_SITE, ID_CDET = (0, 1, 2)
    store = []
    for si, slce in enumerate(slices):

        # nb of species per genus
        gc, _ = numpy.histogram(gen_counts[:-1, si], bins)
        # print("\tNb genera\t", gc)
        nb_gen = numpy.sum(gc[1:])
        nb_gen_mul = numpy.sum(gc[2:])
        nb_gen_div = 0
        if numpy.sum(gen_counts[:-1, si] > 1) > 0:
            # COLLECT FOR EACH OCCURRING GENUS with multiple species:
            # the number of sites where exactly one of its species occurs
            # the number of sites where one or more of its species occurs
            X = numpy.array(sliced_divs[si])
            collect = []
            for gi in numpy.where(gen_counts[:-1, si] > 1)[0]:
                mask_gen = X[:, ID_GENK] == gi
                h_cdet, _ = numpy.histogram(X[mask_gen, ID_CDET], obins)
                collect.append((gen_counts[gi, si], gi, h_cdet[1], numpy.sum(h_cdet[1:])))

            collect = numpy.array(collect)

            # ratio of species diverse sites
            nb_gen_div = numpy.sum(collect[:, 2] < collect[:, 3])
        store.append((nb_gen, nb_gen_mul, nb_gen_div))
    return store


def get_gen_details():
    return [(0, "nbTot", "Nb. of genera total", "%.0f", False),
            (1, "nbMulti", "Nb. of genera with multiple species", "%.0f", False),
            (2, "nbCooc", "Nb. of genera with co-occurring species", "%.0f", True)]

#####################################
# ALL MEASURES AND NULL MODELS DETAILS


def get_statistics_details(group):
    type_measures = ["gen", "pairs", "marg"]
    # type_measures = ["gen"]

    values = {}
    values["gen"] = get_gen_details()
    values["marg"] = get_marg_details()

    pair_measures = ["fisher", "cscore"]  # "jacc", "min"
    # prefix variables indicating status of pair genera in scores produced by compute_pairs_scores used by compute_pairs_stats, stored in DATA_PROBAS_FILE
    nb_status = len(get_pair_statusvars())
    pm_map = dict([(v, k+nb_status) for (k, v) in enumerate(pair_measures)])
    pair_stats = [(0, "len", None, "nbPairs", "nb. of pairs", "%.0f", False),
                  (pm_map["fisher"], "ratio_gt", 0.95, "FisherP95", "Fraction of pairs Fisher P>0.95", "%.4f", False),
                  (pm_map["fisher"], "ratio_lt", 0.05, "FisherP05", "Fraction of pairs Fisher P<0.05", "%.4f", False),
                  (pm_map["fisher"], "qunt", 95, "FisherQ95", "Fisher Q=0.95", "%.6f", True),
                  (pm_map["fisher"], "qunt", 5, "FisherQ05", "Fisher Q=0.05", "%.6f", True),
                  (pm_map["cscore"], "qunt", 95, "CscoreQ95", "C-score Q=0.95", "%.6f", True),
                  (pm_map["cscore"], "qunt", 5, "CscoreQ05", "C-score Q=0.05", "%.6f", False),
                  (pm_map["cscore"], "avg", None, "CscoreAvg", "Avg. C-score", "%.6f", True),
                  ]

    pair_agg_series = get_agg_series(group)
    values["pairs"] = []
    for si, stat in enumerate(pair_stats):
        values["pairs"].append((len(pair_agg_series)*si, stat[-4], stat[-3], stat[-2], stat[-1]))
    return type_measures, pair_measures, pm_map, values, pair_stats, pair_agg_series


def get_nullmodels_details():
    rnd_subseries = {"UG": [], "CB": [], "shuffle": []}  # , "sample": []}
    for slc in ["+"]:  # , "-"]:
        rnd_subseries["UG"].extend(["null-UG_x%s" % slc])
        rnd_subseries["CB"].extend(["null-CB_x%si100" % slc])
        rnd_subseries["shuffle"].extend(["shuffle-species-indet_x%s" % slc])
        # rnd_subseries["shuffle"].extend(["shuffle-species_x%s" % slc, "shuffle-species-indet_x%s" % slc])
        # for srt in ["0"]:  # , "+", "-"]:
        #     for gs in ["0.01", "0.1", "0.5", "1", "2", "10", "100"]:
        #         rnd_subseries["sample"].append("sample-sites_x%sz%sS%s" % (slc, srt, gs))
        #     for m in ["1.05", "1.5"]:  # , "1.1"
        #         rnd_subseries["sample"].append("sample-sites_x%sz%sM%s" % (slc, srt, m))
        #     for m in ["1.05", "1.5"]:  # , "1.1"
        #         for gs in ["0.1", "10"]:
        #             rnd_subseries["sample"].append("sample-sites_x%sz%sM%sS%s" % (slc, srt, m, gs))
    return rnd_subseries

# SAVING AND LOADING STATS
#####################################


def save_stats_boxes(collect, id_val, lgd_val, fmt_val, store_stats_file):
    S = numpy.array(collect)
    if S.shape[0] > 1:
        numpy.savetxt(store_stats_file, S[1:, :, id_val], fmt=fmt_val, delimiter=',', header=lgd_val)
    else:
        numpy.savetxt(store_stats_file, S[:1, :, id_val], fmt=fmt_val, delimiter=',', header=lgd_val)


def load_stats_boxes(fmt_val, store_stats_file):
    X = numpy.loadtxt(store_stats_file, delimiter=',', skiprows=1)
    if not numpy.any(numpy.isnan(X)) and (re.match("%d", fmt_val) or re.match("%.0f", fmt_val)):
        X = numpy.array(X, dtype=int)
    return X


def collect_values(rndk, focus_values, val_keys, pair_agg_series, values, data_statb_patt):
    if focus_values is None:
        focus_values = [k for (k, v) in val_keys.items() if v[-1] >= 0]

    cv = []
    for vname in focus_values:
        typ, vi, o = val_keys[vname]
        if typ != "marg" or rndk != "original":
            if o == -1:
                for off, gname in enumerate(pair_agg_series):
                    cv.append((vname+"_"+gname, values[typ][vi][3]))
            else:
                cv.append((vname, values[typ][vi][3]))

    collected = {}
    for vname, vfmt in cv:
        try:
            collected[vname] = load_stats_boxes(vfmt, data_statb_patt % (vname, rndk))
        except FileNotFoundError:
            print(">> No such file: ", data_statb_patt % (vname, RNDK))
    # print("Collected", rndk, collected.keys())
    return collected
