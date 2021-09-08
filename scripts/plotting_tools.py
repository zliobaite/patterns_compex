import backend_cust_svg
import pdb
import re
import numpy
import matplotlib.pyplot as plt

import common_tools as ct

import matplotlib.cm
import matplotlib as mpl
mpl.rcParams['font.size'] = 10
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Times']
mpl.rcParams['text.usetex'] = True
# mpl.rcParams['xtick.major.pad'] = 20
# mpl.rcParams['xtick.minor.pad'] = 20
# mpl.rcParams['svg.fonttype'] = 'none'
# ax.tick_params(axis='x', which="both", pad=20.0)

# from matplotlib import rc
# #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('font',**{'size': 10, 'family':'serif','serif':['Times']})
# rc('text', usetex=True)

FSZ_SML = "medium"
FSZ_XSML = "small"

# C_NEUTRAL = "#882E72"
C_NEUTRAL = "#666666"
C_BLACK = "#000000"
C_RED = "#D92120"
C_GREY = "#666666"
# C_DIFF = "#CC6677"
# C_SAME = "#117733"
# C_DIFF = "#4477AA"
# C_SAME = "#DDCC77"
C_DIFF = "#D95F0E"
C_SAME = "#FEC44F"
C_FOUR = "#999933"  # "#88CCEE", "#332288", "#DDCC77", "#999933"

C_POSITIVE = "#999933"
C_NEGATIVE = "#88CCEE"
C_SIGN = "#88CCEE"

ALPH_FULL = "FF"
ALPH_MID = "66"

# DEFINE AN ORDERED LIST OF COLORS AND ASSIGN TO GENERA BY DECREASING FREQ
# COLORS_ORD = ["#117733", "#999933", "#332288", "#AA4466", "#6699CC", "#AA4499"]
COLORS_ORD = ["#7DB874", "#E39C37", "#6699CC", "#AA4466", "#44AA99", "#AA4499"]
ORD_COLORS = {}
for i, v in enumerate(["Artiodactyla", "Perissodactyla", "Proboscidea", "Primates", "Hyracoidea"]):  # L
    ORD_COLORS[v] = COLORS_ORD[i]
for i, v in enumerate(["Rodentia", "Eulipotyphla", "Lagomorpha", "Leptictida", "Cimolesta", "Didelphimorphia"]):  # S
    ORD_COLORS[v] = COLORS_ORD[i]
for i, v in enumerate(["Carnivora", "Creodonta", "Tubulidentata"]):  # C
    ORD_COLORS[v] = COLORS_ORD[i]

PLBLS = {"NA": ["O", "W", "A1", "A2", "A3", "A4", "Hf1", "Hf2", "B1", "B2", "C1", "C2", "C3", "Hp1", "Hp2", "Hp3", "Hp4", "BE"],
         "EU": ["MN2", "MN3", "MN4", "MN5", "MN6", "MN7-8", "MN9", "MN10", "MN11", "MN12", "MN13", "MN14", "MN15", "MN16"]}

LABELS_ON = False
FIGW_LRG, FIGH_LRG = (9, 2.4)
FIGW_SML, FIGH_SML = (4.66, 2.)
ADJN_LEFT, ADJW_LEFT = (0.08, 0.1)
ADJN_BOTTOM, ADJW_BOTTOM = (0.15, 0.25)
ADJN_RIGHT, ADJW_RIGHT = (0.98, 0.95)
ADJN_TOP, ADJW_TOP = (0.98, 0.95)

FIGSIZE = {}
SUB_ADJ = {}

FIGSIZE["counts_orders"] = (FIGW_SML, FIGH_LRG)
SUB_ADJ["counts_orders"] = {"left": ADJN_LEFT, "bottom": ADJN_BOTTOM, "right": ADJN_RIGHT, "top": ADJW_TOP}

FIGSIZE["outline_simple"] = (FIGW_SML, 1.)  # (FIGW_LRG, 1.5)
SUB_ADJ["outline_simple"] = {"left": ADJN_LEFT, "bottom": 0.05, "right": ADJN_RIGHT, "top": ADJW_TOP}

FIGSIZE["div_hists"] = (FIGW_SML, FIGH_SML)
SUB_ADJ["div_hists"] = {"left": ADJN_LEFT, "bottom": ADJN_BOTTOM, "right": ADJN_RIGHT, "top": ADJW_TOP}

FIGSIZE["pairs_hists"] = (FIGW_LRG, FIGH_SML)
SUB_ADJ["pairs_hists"] = {"left": 0.04, "bottom": ADJN_BOTTOM, "right": ADJN_RIGHT, "top": ADJN_TOP}

FIGSIZE["boxes_simple"] = (FIGW_SML, FIGH_SML)
SUB_ADJ["boxes_simple"] = {"left": ADJN_LEFT, "bottom": ADJN_BOTTOM, "right": ADJN_RIGHT, "top": ADJN_TOP}

FIGSIZE["multi_boxes_series"] = (FIGW_SML, FIGH_SML)
SUB_ADJ["multi_boxes_series"] = {"left": ADJN_LEFT, "bottom": ADJN_BOTTOM, "right": ADJN_RIGHT, "top": ADJN_TOP}

#######################
FIGSIZE["outline_time"] = (FIGW_LRG, 1.5)
SUB_ADJ["outline_time"] = {"left": ADJN_LEFT, "bottom": 0.22, "right": ADJN_RIGHT, "top": ADJN_TOP}

FIGSIZE["counts_species"] = (FIGW_LRG, 2.6)
SUB_ADJ["counts_species"] = {"left": ADJN_LEFT, "bottom": 0.3, "right": ADJN_RIGHT, "top": ADJN_TOP, "wspace": 0.01}

FIGSIZE["pairs_hists_cut"] = (FIGW_LRG, FIGH_LRG)
SUB_ADJ["pairs_hists_cut"] = {"left": ADJN_LEFT, "bottom": ADJW_BOTTOM, "right": ADJN_RIGHT, "top": ADJN_TOP}

FIGSIZE["div_complex"] = (FIGW_LRG, 3)
SUB_ADJ["div_complex"] = {"left": ADJN_LEFT, "bottom": 0.20, "right": ADJN_RIGHT, "top": ADJN_TOP}

SUB_ADJ["multi_boxes"] = {"left": ADJN_LEFT, "bottom": ADJW_BOTTOM, "right": ADJW_RIGHT, "top": ADJN_TOP}

FIGSIZE["overlap_areas"] = (FIGW_LRG, 3)
SUB_ADJ["overlap_areas"] = {"left": ADJN_LEFT, "bottom": 0.138, "right": ADJN_RIGHT, "top": ADJN_TOP, "wspace": 0.01}


def prep_fig(pstyle, with_axe=True):
    if pstyle in FIGSIZE:
        main_fig = plt.figure(figsize=FIGSIZE[pstyle])
    else:
        main_fig = plt.figure()
    if with_axe:
        ax1 = plt.gca()
    else:
        ax1 = None
    return main_fig, ax1


def finish_fig(fig, axe, pstyle, fig_file=None, params={}):
    sub_adj = dict(SUB_ADJ.get(pstyle, {}))
    # x-axis
    if "slices" in params:
        s_off = 0
        x_ticks, x_lbls = zip(*[(si+0.5+s_off, slce[-1]) for si, slce in enumerate(params["slices"])])
        if x_lbls[0] == "MN2":
            x_lbls = PLBLS["EU"]
        else:
            x_lbls = PLBLS["NA"]
        axe.tick_params(axis='x', which="both", pad=8.0)
        # axe.tick_params(axis='x', which="both", pad=15.0)
        # sub_adj["bottom"] = 0.22
        axe.set_xticks(x_ticks)
        axe.set_xticklabels(x_lbls, rotation=25, va="center", fontsize=FSZ_XSML)
        if LABELS_ON:
            axe.set_xlabel("Time bins")
        axe.set_xlim(0+s_off, len(x_lbls)+s_off)

    else:
        if "xparams" in params:
            axe.tick_params(axis='x', which="both", **params["xparams"])
        if "xlim" in params:
            axe.set_xlim(params["xlim"][0], params["xlim"][1])
        if "xlabel" in params and LABELS_ON:
            axe.set_xlabel(params["xlabel"], fontsize=FSZ_XSML, **params.get("xlabelparams", {}))
        if "xticks" in params and "xticklabels" in params:
            axe.set_xticks(params["xticks"])
            axe.set_xticklabels(params["xticklabels"], fontsize=FSZ_XSML, **params.get("xtickparams", {}))

    # y-axis
    axs = [axe]
    if params.get("ax2") is not None:
        axs.append(params["ax2"])
    for ax in axs:
        if params.get("yticks") == "int":
            y_ticks = ax.get_yticks()
            y_ticks = sorted(set([int(vy) for vy in y_ticks]))
            y_lbls = ["%d" % abs(vy) for vy in y_ticks]
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_lbls, fontsize=FSZ_XSML, **params.get("ytickparams", {}))
        elif params.get("yticks") == "auto":
            y_ticks = ax.get_yticks()
            fmt = params.get("yticksfmt", "%.1f")
            y_ticks = sorted(set([float(fmt % vy) for vy in y_ticks]))
            y_lbls = [fmt % vy for vy in y_ticks]
            ax.set_yticks(y_ticks)
            ax.set_yticklabels(y_lbls, fontsize=FSZ_XSML, **params.get("ytickparams", {}))
        elif "yticks" in params and "yticklabels" in params:
            ax.set_yticks(params["yticks"])
            ax.set_yticklabels(params["yticklabels"], fontsize=FSZ_XSML, **params.get("ytickparams", {}))

        if "yparams" in params:
            ax.tick_params(axis='y', which="both", **params["yparams"])
        if "ylim" in params:
            ax.set_ylim(params["ylim"][0], params["ylim"][1])
        if "ylabel" in params and LABELS_ON:
            ax.set_ylabel(params["ylabel"], fontsize=FSZ_XSML, **params.get("ylabelparams", {}))

    if len(sub_adj) > 0:
        plt.subplots_adjust(**sub_adj)

    if fig_file is None:
        plt.show()
    else:
        plt.savefig(fig_file)
        plt.close()


# BACKGROUND
#####################################
def plot_counts_species(fossils_data, spc_counts, sliced_counts, map_fields_species, map_fields_sites, fig_outline_file):
    pstyle = "counts_species"
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])

    gen_field = map_fields_species["GENUS"]
    spc_field = map_fields_species["SPECIES"]
    indeterminate_sids_all = set([spci for spci, s in enumerate(species) if ct.is_indeterminate_gs(s[gen_field], s[spc_field])])

    gen_tots = {}
    for ssi, s in enumerate(species):
        if s[gen_field] not in gen_tots:
            gen_tots[s[gen_field]] = 0
        gen_tots[s[gen_field]] += 1
    # map_gens = dict([(v,k) for k,v in enumerate(sorted(gen_tots.keys(), key=lambda x: gen_tots[x]))])
    map_gens = dict([(v, k) for k, v in enumerate(sorted(gen_tots.keys()))])

    main_fig, ax1 = prep_fig(pstyle, with_axe=False)
    splts = main_fig.subplots(1, len(sliced_counts), sharey=True)
    for si, scounts in enumerate(sliced_counts):
        sids = numpy.where(spc_counts[:-1, si] > 0)[0]
        full_mask = numpy.ones(len(sids), dtype=bool)
        det_mask = numpy.array([(i not in indeterminate_sids_all) for i in sids], dtype=bool)

        counts = spc_counts[sids, si]
        sc, _ = numpy.histogram(counts, bins=[1, 2, 3, 4, 5, 100])
        scd, _ = numpy.histogram(counts[det_mask], bins=[1, 2, 3, 4, 5, 100])
        splts[si].bar(numpy.arange(sc.shape[0]), sc, color=C_NEUTRAL+ALPH_MID)
        splts[si].bar(numpy.arange(scd.shape[0]), scd, color=C_NEUTRAL+ALPH_FULL)
        splts[si].tick_params(axis='x', which="both", pad=3.0)
        splts[si].set_xticks([0, 2, 4])
        splts[si].set_xticklabels(["1", "3", "5+"])
        if si == 0:
            splts[si].set_ylabel("Nb. species")
        splts[si].set_xlabel(slices[si][-1], rotation=25, va="center", fontsize=FSZ_SML, labelpad=20)
        # #    splts[si].set_yticks([])
    finish_fig(main_fig, ax1, pstyle, fig_outline_file)


def plot_counts_orders(fossils_data, fields_species, map_fields_sites, fig_outline_file=None):
    pstyle = "counts_orders"
    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])
    occurrences = fossils_data["occurrences"]

    gen_field = fields_species.index("GENUS")
    spc_field = fields_species.index("SPECIES")
    g_field = fields_species.index("ORDER")
    # g_field = fields_species.index("FAMILY")
    indeterminate_sids_mask = numpy.array([ct.is_indeterminate_gs(s[gen_field], s[spc_field]) for spci, s in enumerate(species)])
    g_counts = {}
    for si in numpy.where(~indeterminate_sids_mask)[0]:
        if species[si][g_field] not in g_counts:
            g_counts[species[si][g_field]] = 0
        g_counts[species[si][g_field]] += 1
    gs = sorted(g_counts.keys(), key=lambda x: -g_counts[x])

    store_counts = []
    M_counts = []
    for xi, s in enumerate(slices):
        sliced_sites = sliced_sites_all[xi]
        if len(sliced_sites) > 0:
            spc_counts = occurrences[sliced_sites, :].sum(0)
            mask_spc = (~indeterminate_sids_mask) & (spc_counts > 0)
            sids = numpy.where(mask_spc)[0]
            Xgen_set = set([species[sid][gen_field] for sid in numpy.where(spc_counts > 0)[0]])
            gen_set = set([species[sid][gen_field] for sid in sids])
            M = occurrences[sliced_sites, :][:, mask_spc]
            spc_margs = numpy.sum(M, axis=0)
            # spc_list = ["%s %s" % (species[sids[k]][gen_field], species[sids[k]][spc_field]) for k in numpy.where(spc_margs==1)[0]]
            loc_margs = numpy.sum(M, axis=1)
            locs_uniq_spc = sorted(numpy.where(numpy.sum(M[:, spc_margs == 1], axis=1) > 0)[0], key=lambda x: -loc_margs[x])
            counts_gs = dict([(si, [0, 0, 0]) for si in gs])
            counts_gs["tot_sites"] = [numpy.sum(loc_margs > 1), numpy.sum(loc_margs == 1)]
            counts_gs["tot_species"] = [numpy.sum(spc_margs > 1), numpy.sum(loc_margs == 1)]
            counts_gens = {}
            tot_true = 0
            for si, u in enumerate(spc_margs):
                o = species[sids[si]][g_field]
                counts_gs[o][0] += (u > 1)
                counts_gs[o][1] += (u == 1)
                counts_gs[o][2] += (u == 1) & (numpy.sum(occurrences[:, sids[si]]) == 1)
                tot_true += (u == 1) & (numpy.sum(occurrences[:, sids[si]]) == 1)
            true_unique = numpy.where(numpy.sum(occurrences[:, sids[spc_margs == 1]], axis=0) == 1)[0]
            xxx = numpy.where(spc_margs == 1)[0][true_unique]
            spc_list = sorted([" ".join(species[sids[k]]) for k in xxx])
            #     o = species[sids[si]][gen_field]
            #     if o not in counts_gens:
            #         counts_gens[o] = [0,0]
            #     counts_gens[o][0] += (u>1)
            #     counts_gens[o][1] += (u==1)

            # gen_A = [gen for gen, c in counts_gens.items() if c[0] > 1]
            # gen_B = [gen for gen, c in counts_gens.items() if c[0]+c[1] > 1]
            # #print("%d\t%d\t%d\t%s" % (len(gen_A), len(gen_B), len(counts_gens), gen_A))
            # print("%d\t%d\t%d" % (len(gen_A), len(gen_B), len(counts_gens)))
            # store_counts.append((counts_gs, counts_gens))
            store_counts.append(counts_gs)
            tmp = [0]
            for si in gs:
                tmp.extend(counts_gs[si][:2])
            M_counts.append(tmp)
            # print("%s\tloc:%d/%d=%.3f\tspc:[%d] %d/%d=%.3f\t%s\t%s" % (s[-1], numpy.sum(loc_margs==1), len(loc_margs), numpy.sum(loc_margs==1)/ len(loc_margs), tot_true, numpy.sum(spc_margs==1), len(spc_margs), numpy.sum(spc_margs==1)/ len(spc_margs), "; ".join(["%s [%d] %d/%d=%.3f" % (si, counts_gs[si][-1], counts_gs[si][0], counts_gs[si][0]+counts_gs[si][1], counts_gs[si][0]/(counts_gs[si][0]+counts_gs[si][1]) if (counts_gs[si][0]+counts_gs[si][1]) > 0 else -1) for si in gs]), spc_list))
            # for loc in locs_uniq_spc:
            #     u_sids = numpy.where((M[loc, :] > 0) & (spc_margs ==1))[0]
            #     # spc_list = sorted(["%s %s" % (species[sids[k]][gen_field], species[sids[k]][spc_field]) for k in u_sids])
            #     spc_list = sorted([" ".join(species[sids[k]]) for k in u_sids])
            #     print("At %s\t%d/%d\t%s" % (sites[sliced_sites[loc]][0], len(spc_list), numpy.sum(M[loc, :] > 0), "; ".join(spc_list)))

    M_counts = numpy.array(M_counts)
    CC = numpy.cumsum(M_counts, axis=1)
    # hcmap = matplotlib.cm.get_cmap("tab20")
    # cmap_nbc = 20
    # colors = []
    # while len(colors) < M_counts.shape[1]:
    #     colors.extend([hcmap(i/cmap_nbc) for i in range(cmap_nbc)])
    colors = []
    for g in gs:
        colors.append(ORD_COLORS[g]+ALPH_FULL)
        colors.append(ORD_COLORS[g]+ALPH_MID)
    colors = colors[:M_counts.shape[1]-1]

    main_fig, ax1 = prep_fig(pstyle)
    for i in range(M_counts.shape[0]):
        plt.bar(i*numpy.ones(M_counts.shape[1]-1)+.5, M_counts[i, 1:], 0.85, bottom=CC[i, :-1], color=colors)
    top = numpy.max(CC[:, -1])

    # xlgd_0 = 0.1
    # xx_off = 0.12
    xlgd_0 = 0.1
    xx_off = 0.2*len(slices)/len(PLBLS["EU"])
    if gs[0] == "Rodentia":
        if slices[0][-1] == "MN2":
            xlgd_0 = 7
        else:
            xlgd_0 = 12
    plt.text(xlgd_0-xx_off+0.5, top, "1", verticalalignment='center', horizontalalignment='center', zorder=10)
    plt.text(xlgd_0+xx_off+0.5, top, "+", verticalalignment='center', horizontalalignment='center', zorder=10)
    for gi, g in enumerate(gs):
        ypos = top*(1-0.075*(len(gs)-gi))
        plt.scatter([xlgd_0+xx_off+0.5, xlgd_0-xx_off+0.5], [ypos, ypos], c=colors[2*gi:2*(gi+1)], marker="o", zorder=10)
        # plt.text(xlgd_0+0.28, ypos, "%s" % g, verticalalignment='center', fontsize=FSZ_SML, zorder=10)
        plt.text(xlgd_0+2.2*xx_off+0.5, ypos, "$\\textit{%s}$" % g, verticalalignment='center', fontsize=FSZ_XSML, zorder=10)

    params = {"slices": slices, "ylabel": "Nb. species", "ylim": [0, 1.05*top], "yticks": "int"}
    finish_fig(main_fig, ax1, pstyle, fig_outline_file, params)
    # return store_counts, gs


def plot_bubbles_species(fossils_data, spc_counts, sliced_counts, map_fields_species, map_fields_sites, fig_outline_file, teeth_file=None):
    pstyle = "bubble_species"
    out_svg = fig_outline_file is not None and re.search("svg$", fig_outline_file) is not None
    if out_svg:
        # pstyle = "bubble_species_svg"
        store_ff = mpl.rcParams['font.family']
        mpl.rcParams['font.family'] = 'sans-serif'

        def dot_size(vs):
            # return 2*numpy.log2(1+vs)
            return 2.5*(2+vs)
    else:
        def dot_size(vs):
            # return 2*numpy.log2(1+vs)
            return 2+vs

    if teeth_file is not None and fig_outline_file is not None and re.search("genD", fig_outline_file):
        teeth_data = read_teeth_details(teeth_file)
    else:
        teeth_data = {}

    # MARKERS = {None: "o", "bra": ">", "hyp": "<", "hys": "^", "mes": "v"}
    MARKERS = {None: "o", "bra": "^", "hyp": "s", "hys": "s", "mes": "s"}

    cmapC = "tab20"
    # cmapC = "rainbow"
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])

    gen_field = map_fields_species["GENUS"]
    spc_field = map_fields_species["SPECIES"]
    fam_field = map_fields_species["FAMILY"]
    ord_field = map_fields_species["ORDER"]
    if fig_outline_file is not None and re.search("genF", fig_outline_file):
        color_field = map_fields_species["FAMILY"]
    else:
        color_field = map_fields_species["ORDER"]
    indeterminate_sids_all = set([spci for spci, s in enumerate(species) if ct.is_indeterminate_gs(s[gen_field], s[spc_field])])
    filtr = None
    if re.search("genR", fig_outline_file):
        filtr = (ord_field, "Rodentia")
        color_field = map_fields_species["FAMILY"]

    gen_to_color = {}
    gen_to_ord = {}
    collect = []
    gen_tots = {}

    # slices = slices[:2]
    # sliced_counts = sliced_counts[:2]
    for si, scounts in enumerate(sliced_counts):
        sids = numpy.where(spc_counts[:-1, si] > 0)[0]

        gen_counts = {}
        for s in sids:
            if s not in indeterminate_sids_all and (filtr is None or species[s][filtr[0]] == filtr[1]):
                if species[s][gen_field] not in gen_counts:
                    gen_counts[species[s][gen_field]] = []
                kt = tuple([species[s][xx] for xx in [fam_field, gen_field, spc_field]])
                xtr = MARKERS[teeth_data.get(kt, None)]
                gen_counts[species[s][gen_field]].append((spc_counts[s, si], xtr, "LBL:%s %s (%d)" % (species[s][gen_field], species[s][spc_field], spc_counts[s, si])))
                if species[s][gen_field] not in gen_to_color:
                    if re.match("incert", species[s][color_field]):
                        print(species[s])
                    gen_to_color[species[s][gen_field]] = species[s][color_field]
                    gen_to_ord[species[s][gen_field]] = species[s][ord_field]

        for k in gen_counts.keys():
            gen_counts[k].sort()
            if k not in gen_tots:
                gen_tots[k] = []
            gen_tots[k].append((si, len(gen_counts[k])))
        collect.append(gen_counts)

    for k in gen_tots.keys():
        gen_tots[k] = list(zip(*gen_tots[k]))
    colors = sorted(set(gen_to_color.values()))
    colors_num = dict([(k, i/len(colors)) for (i, k) in enumerate(colors)])
    gen_max = dict([(k, max(cs[1])) for (k, cs) in gen_tots.items()])
    gen_multi = [k for (k, m) in gen_max.items() if m > 1]
    gen_pos = dict([(k, i) for (i, k) in enumerate(sorted(gen_multi, key=lambda x: gen_tots[x]))])
    max_nb = max(gen_max.values()) + 1
    top = max(gen_pos.values())+1

    if out_svg:
        fmt_params = {"fig_h": 0.12*top, "fig_w": 13.5, "lgd_fontsize": "large",
                      "adj_l": 0.12, "adj_b": 9/top, "adj_r": 0.88, "adj_t": 0.98}
    else:
        fmt_params = {"fig_h": 0.075*top, "fig_w": 9, "lgd_fontsize": FSZ_SML,
                      "adj_l": 0.165, "adj_b": 9/top, "adj_r": 0.835, "adj_t": 0.98}
    main_fig = plt.figure(figsize=(fmt_params["fig_w"], fmt_params["fig_h"]))

    y_lbls, y_ticks = zip(*gen_pos.items())
    for yy in y_ticks[1::2]:
        plt.plot([0, len(slices)], [yy, yy], "--", c="#DDDDDD", linewidth=0.1, zorder=0)
    for yy in y_ticks[::2]:
        plt.plot([-1, len(slices)-0.048], [yy, yy], "--", c="#AAAAAA", linewidth=0.1, zorder=0)

    bottom = 0
    # splts = main_fig.subplots(1, len(sliced_counts), sharey=True)
    for i, gs in enumerate(collect):
        cts_singles = []
        gen_singles = []
        lbl_singles = []
        mrk_singles = []
        for j, (g, xxx) in enumerate(gs.items()):
            cts, mrks, lbls = zip(*xxx)
            if g in gen_pos:
                xs = numpy.arange(len(cts))/max_nb + i
                ys = gen_pos[g]*numpy.ones(len(cts))
                for mrk in set(mrks):
                    cis = [mi for mi, m in enumerate(mrks) if m == mrk]
                    urls = [lbls[mi] for mi, m in enumerate(mrks) if m == mrk]
                    # plt.scatter(xs[cis], ys[cis], s=dot_size(numpy.array(cts)[cis]), c=numpy.ones(len(cis))*colors_num[gen_to_color[g]], vmin=0, vmax=1, cmap=cmapC, marker=mrk, zorder=5)
                    plt.scatter(xs[cis], ys[cis], s=dot_size(numpy.array(cts)[cis]), c=ORD_COLORS[gen_to_ord[g]], marker=mrk, zorder=5)
                    if out_svg:
                        plt.scatter(xs[cis], ys[cis], s=dot_size(numpy.array(cts)[cis]), c="#ffffff01", marker='o', zorder=10, urls=urls)
            else:
                cts_singles.append(cts[0])
                mrk_singles.append(mrks[0])
                lbl_singles.append(lbls[0])
                gen_singles.append(g)

        ys = numpy.floor_divide(numpy.arange(len(cts_singles)), max_nb)
        xs = (numpy.arange(len(cts_singles))/(max_nb) - ys) + i
        cs = numpy.array([colors_num[gen_to_color[g]] for g in gen_singles])
        clrs = [ORD_COLORS[gen_to_ord[g]] for g in gen_singles]
        for mrk in set(mrk_singles):
            cis = [mi for mi, m in enumerate(mrk_singles) if m == mrk]
            urls = [lbl_singles[mi] for mi, m in enumerate(mrk_singles) if m == mrk]
            # plt.scatter(xs[cis], -(2+ys[cis]), s=dot_size(numpy.array(cts_singles)[cis]), c=cs[cis], vmin=0, vmax=1, cmap=cmapC, marker=mrk, zorder=5)
            plt.scatter(xs[cis], -(2+ys[cis]), s=dot_size(numpy.array(cts_singles)[cis]), c=[clrs[ii] for ii in cis], marker=mrk, zorder=5)  # , urls=urls)
            if out_svg:
                plt.scatter(xs[cis], -(2+ys[cis]), s=dot_size(numpy.array(cts_singles)[cis]), c="#ffffff01", marker='o', zorder=10, urls=urls)
        bottom = max(bottom, numpy.max(ys))
    bottom = -(2+bottom)

    # legend
    for oi, (order, num) in enumerate(colors_num.items()):
        # plt.scatter([0.2], [top-2.7*(oi)-1], c=[num], vmin=0, vmax=1, cmap=cmapC, marker="o", zorder=10)
        plt.scatter([0.35], [top-2.5*(oi)-1], c=ORD_COLORS[order], marker="o", zorder=10)
        # plt.text(0.6, top-2.5*(oi)-1, "%s" % order, verticalalignment='center', zorder=10, fontsize=fmt_params["lgd_fontsize"])
        plt.text(0.6, top-2.5*(oi)-1, "$\\textit{%s}$" % order, verticalalignment='center', zorder=10, fontsize=fmt_params["lgd_fontsize"])
    for oi, num in enumerate([1, 2, 4, 8, 16]):
        plt.scatter([4.2], [top-2.5*(oi)-1], s=dot_size(num), c="#222222", marker="o", zorder=10)
        plt.text(4.4, top-2.5*(oi)-1, "%d" % num, verticalalignment='center', zorder=10, fontsize=fmt_params["lgd_fontsize"])
    # for oi, (th, mrk) in enumerate(MARKERS.items()):
    if len(teeth_data) > 0:
        for oi, (th, mrk) in enumerate([(None, MARKERS[None]), ("brachydont", MARKERS["bra"]), ("other", MARKERS["hys"])]):
            if th is None:
                th = "?"
            plt.scatter([5.2], [top-2.5*(oi)-1], c="#222222", marker=mrk, zorder=10)
            plt.text(5.4, top-2.5*(oi)-1, th, verticalalignment='center', zorder=10, fontsize=fmt_params["lgd_fontsize"])

    ax1 = plt.gca()
    ax2 = ax1.twinx()

    ax1.set_ylim([bottom-.75, top+0.75])
    ax2.set_ylim([bottom-.75, top+0.75])
    ax1.set_yticks(y_ticks[::2])
    # ax1.set_yticklabels(["%s" % l for l in y_lbls[::2]], fontsize=fmt_params["lgd_fontsize"])
    ax1.set_yticklabels(["$\\textit{%s}$" % l for l in y_lbls[::2]], fontsize=fmt_params["lgd_fontsize"])

    ax2.set_yticks(y_ticks[1::2])
    # ax2.set_yticklabels(["%s" % l for l in y_lbls[1::2]], fontsize=fmt_params["lgd_fontsize"])
    ax2.set_yticklabels(["$\\textit{%s}$" % l for l in y_lbls[1::2]], fontsize=fmt_params["lgd_fontsize"])

    for si in range(len(slices)):
        plt.plot([si-.05, si-.05], [bottom-1, top+1], 'k:', linewidth=0.75, color="#333333")
    x_ticks, x_lbls = zip(*[(si+.495, slce[-1]) for si, slce in enumerate(slices)])
    ax1.tick_params(axis='x', which="both", pad=20.0)
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_lbls, rotation=25, va="center", fontsize=fmt_params["lgd_fontsize"])
    # plt.text(si+0.5, 1.05, slce[-1], horizontalalignment="center", verticalalignment="bottom", rotation=60)
    plt.xlim([-0.048, len(slices)])

    plt.subplots_adjust(left=fmt_params["adj_l"], bottom=fmt_params["adj_b"], right=fmt_params["adj_r"], top=fmt_params["adj_t"])
    # ax1.set_xlabel("Time slices")
    finish_fig(main_fig, ax1, pstyle, fig_outline_file)
    if out_svg:
        mpl.rcParams['font.family'] = store_ff


def plot_outline_simple(fossils_data, map_fields_species, map_fields_sites, fig_outline_file):
    pstyle = "outline_simple"
    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])

    D = numpy.array([(int(ss[map_fields_sites["SLICE_ID"]]),
                      -float(ss[map_fields_sites["MIN_AGE"]]), -float(ss[map_fields_sites["MAX_AGE"]]),
                      float(ss[map_fields_sites["MEAN_HYPSODONTY"]]) if ss[map_fields_sites["MEAN_HYPSODONTY"]] != "\\N" else 0,
                      int(ss[map_fields_sites["LIDNUM"]]), ssi) for ssi, ss in enumerate(sites)])
    ID_SLICE, ID_MINAGE, ID_MAXAGE, ID_HYP, ID_NUM, ID_ORD = range(6)

    sort_ids = sorted(range(D.shape[0]), key=lambda x: (D[x, ID_SLICE], D[x, ID_HYP], D[x, ID_MAXAGE], D[x, ID_MINAGE], D[x, ID_NUM]))
    # sort_ids = sorted(range(D.shape[0]), key=lambda x: (D[x, ID_SLICE], D[x, ID_MAXAGE], D[x, ID_MINAGE], D[x, ID_HYP], D[x, ID_NUM]))

    hcmap = matplotlib.cm.get_cmap("copper")
    colors = [hcmap(D[i, ID_HYP]/3) if D[i, ID_HYP] > 0 else (.4, .4, .4, .4) for i in sort_ids]

    offZ = 0.5
    bc = numpy.bincount([int(ss[map_fields_sites["SLICE_ID"]]) for ssi, ss in enumerate(sites)])
    max_in_slice = numpy.max(bc)
    top_line = numpy.max(bc)+offZ

    ys = [0]
    step = -1
    for ii in range(len(sort_ids)):
        if ii > 0 and D[sort_ids[ii], ID_SLICE] != D[sort_ids[ii-1], ID_SLICE]:
            # change slice
            ys.append(0)
        else:
            ys.append(ys[-1]+step)
    ys = numpy.array(ys[1:])

    min_x_org = -slices[0][0]
    max_x_org = -slices[-1][1]
    min_x_map = 0
    max_x_map = len(slices)

    main_fig, ax1 = prep_fig(pstyle)
    # ws = (D[sort_ids, ID_MAXAGE]-D[sort_ids, ID_MINAGE])* (max_x_map-min_x_map)/(max_x_org-min_x_org)
    # ls = (D[sort_ids, ID_MINAGE]-min_x_org)*(max_x_map-min_x_map)/(max_x_org-min_x_org)

    ls = D[sort_ids, ID_SLICE]+0.5

    # XY = numpy.dstack([D[:,[1,2]], numpy.vstack([numpy.arange(D.shape[0]),numpy.arange(D.shape[0])]).T])
    # plt.plot(XY[:,:,0].T, XY[:,:,1].T)#, color=colors)
    # plt.barh(ys, ws, height=1., left=ls, color=colors)
    plt.bar(ls, 1., 0.85, bottom=ys, color=colors, align='center', linewidth=0)
    ticks = []
    hypnz_ids = (D[:, ID_HYP] > 0)
    for si, slce in enumerate(slices):
        sids = (D[:, ID_SLICE] == si)
        min_v = numpy.min(D[hypnz_ids & sids, ID_HYP])
        max_v = numpy.max(D[hypnz_ids & sids, ID_HYP])
        avg_v = numpy.mean(D[hypnz_ids & sids, ID_HYP])

        # plt.plot([si, si], [0, -1.15*top_line], 'k:', linewidth=0.75, color="#333333")

        ticks.append((si+0.5, ""))
        # ticks.append((si, slce[0]))
        # plt.text(p, top_line+2, slce[0], color="k", verticalalignment="bottom", horizontalalignment="center", rotation=25)
        # plt.text(si+0.5, -1.05*top_line, "%.2f" % min_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
        plt.text(si+0.5, -1.05*top_line, "%.2f" % avg_v, color="k", rotation=25, verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_XSML)
        # plt.text(si+0.5, -1.15*top_line, "%.2f" % max_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)

    # ticks.append((len(slices), slices[-1][1]))
    # plt.text(p, top_line+2, slices[-1][1], color="k", verticalalignment="bottom", horizontalalignment="center", rotation=25)
    x_ticks, x_lbls = zip(*ticks)
    params = {"xlim": [0, len(slices)], "xticks": x_ticks, "xticklabels": x_lbls,
              "xparams": {"pad": 8.0, "top": True, "bottom": False, "labeltop": True, "labelbottom": False},
              "ylim": [-1.06*top_line, 0], "ylabel": "Nb. localities", "yticks": "int"}
    finish_fig(main_fig, ax1, pstyle, fig_outline_file, params)


def plot_outline_time(fossils_data, map_fields_species, map_fields_sites, fig_outline_file):
    pstyle = "outline_time"
    sites, sites_dict = (fossils_data["sites"], fossils_data["sites_dict"])
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])
    slices, sliced_sites_all = (fossils_data["slices"], fossils_data["sliced_sites"])

    D = numpy.array([(int(ss[map_fields_sites["SLICE_ID"]]),
                      -float(ss[map_fields_sites["MIN_AGE"]]), -float(ss[map_fields_sites["MAX_AGE"]]),
                      float(ss[map_fields_sites["MEAN_HYPSODONTY"]]) if ss[map_fields_sites["MEAN_HYPSODONTY"]] != "\\N" else 0,
                      int(ss[map_fields_sites["LIDNUM"]]), ssi) for ssi, ss in enumerate(sites)])
    ID_SLICE, ID_MINAGE, ID_MAXAGE, ID_HYP, ID_NUM, ID_ORD = range(6)

    # sort_ids = sorted(range(D.shape[0]), key=lambda x: (D[x, ID_SLICE], D[x, ID_HYP], D[x, ID_MAXAGE], D[x, ID_MINAGE], D[x, ID_NUM]))
    sort_ids = sorted(range(D.shape[0]), key=lambda x: (D[x, ID_SLICE], D[x, ID_MAXAGE], D[x, ID_MINAGE], D[x, ID_HYP], D[x, ID_NUM]))

    hcmap = matplotlib.cm.get_cmap("copper")
    colors = [hcmap(D[i, ID_HYP]/3) for i in sort_ids]

    offZ = 0.5
    bc = numpy.bincount([int(ss[map_fields_sites["SLICE_ID"]]) for ssi, ss in enumerate(sites)])
    max_in_slice = numpy.max(bc)
    top_line = numpy.max(bc*(numpy.mod(numpy.arange(bc.shape[0]), 2) == 0))+offZ
    bottom_line = -numpy.max(bc*(numpy.mod(numpy.arange(bc.shape[0]), 2) == 1))-offZ

    ys = [offZ]
    step = 1
    for ii in range(len(sort_ids)):
        if ii > 0 and D[sort_ids[ii], ID_SLICE] != D[sort_ids[ii-1], ID_SLICE]:
            # change slice
            step = -1*step
            ys.append(numpy.sign(step)*offZ)
        else:
            ys.append(ys[-1]+step)
    ys = numpy.array(ys[1:])

    min_x_org = -slices[0][0]
    max_x_org = -slices[-1][1]
    min_x_map = 0
    max_x_map = len(slices)

    top_plt = top_line+.12*(top_line-bottom_line+2)
    top_pltx = top_line+.08*(top_line-bottom_line+2)

    main_fig, ax1 = prep_fig(pstyle)
    ws = (D[sort_ids, ID_MAXAGE]-D[sort_ids, ID_MINAGE]) * (max_x_map-min_x_map)/(max_x_org-min_x_org)
    ls = (D[sort_ids, ID_MINAGE]-min_x_org)*(max_x_map-min_x_map)/(max_x_org-min_x_org)
    # ws = numpy.ones(D.shape[0])
    # ls = D[sort_ids, ID_SLICE]

    # XY = numpy.dstack([D[:,[1,2]], numpy.vstack([numpy.arange(D.shape[0]),numpy.arange(D.shape[0])]).T])
    # plt.plot(XY[:,:,0].T, XY[:,:,1].T)#, color=colors)
    plt.barh(ys, ws, height=0.8, left=ls, color=colors)
    ticks = []
    hypnz_ids = (D[:, ID_HYP] > 0)
    for si, slce in enumerate(slices):
        sids = (D[:, ID_SLICE] == si)
        min_v = numpy.min(D[hypnz_ids & sids, ID_HYP])
        max_v = numpy.max(D[hypnz_ids & sids, ID_HYP])
        avg_v = numpy.mean(D[hypnz_ids & sids, ID_HYP])

        stat_line = numpy.min([-bottom_line, top_line])

        # p = si
        # p_up = si+1
        # plt.plot([si, si, si, si], [bottom_line-1, top_line+1, top_pltx, top_plt], 'k:', linewidth=0.75, color="#333333")

        p = (-slce[0]-min_x_org)*(max_x_map-min_x_map)/(max_x_org-min_x_org)
        p_up = (-slce[1]-min_x_org)*(max_x_map-min_x_map)/(max_x_org-min_x_org)
        plt.plot([p, p, si, si], [bottom_line-1, top_line+1, top_pltx, top_plt], 'k:', linewidth=0.75, color="#333333")

        ticks.append((p, slce[0]))
        # plt.text(p, top_line+2, slce[0], color="k", verticalalignment="bottom", horizontalalignment="center", rotation=25)
        if si % 2 == 1:
            # plt.text((p+p_up)/2, 0.66*stat_line, "%.2f" % min_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
            plt.text((p+p_up)/2, 0.5*stat_line, "%.2f" % avg_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
            # plt.text((p+p_up)/2, 0.33*stat_line, "%.2f" % max_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
        else:
            # plt.text((p+p_up)/2, -0.33*stat_line, "%.2f" % min_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
            plt.text((p+p_up)/2, -0.5*stat_line, "%.2f" % avg_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
            # plt.text((p+p_up)/2, -0.66*stat_line, "%.2f" % max_v, color="k", verticalalignment="bottom", horizontalalignment="center", fontsize=FSZ_SML)
    # plt.text(p, top_line+2, slices[-1][1], color="k", verticalalignment="bottom", horizontalalignment="center", rotation=25)
    p = (-slices[-1][1]-min_x_org)*(max_x_map-min_x_map)/(max_x_org-min_x_org)
    plt.plot([p, p, len(slices)], [bottom_line-1, top_line+1, top_plt], 'k:', linewidth=0.75, color="#333333")

    ticks.append((p, slices[-1][1]))
    x_ticks, x_lbls = zip(*ticks)
    params = {"xlim": [0, len(slices)], "xticks": x_ticks, "xticklabels": x_lbls,
              "xtickparams": {"verticalalignment": "top", "horizontalalignment": "center", "rotation": 40},
              "ylim": [bottom_line, top_plt], "yticks": [], "yticklabels": []}
    finish_fig(main_fig, ax1, pstyle, fig_outline_file, params)


def plot_slices_times(slices, sites, which, map_fields_sites, fig_tslices_file=None, max_ratio_out=0):
    pstyle = "slices_times"
    times_dict = {}
    for site in sites:
        sk = (site[map_fields_sites["MAX_AGE"]], site[map_fields_sites["MIN_AGE"]])
        if sk not in times_dict:
            times_dict[sk] = []
        times_dict[sk].append(site[map_fields_sites["LIDNUM"]])
    site_times = sorted(set([(site[map_fields_sites["MAX_AGE"]], site[map_fields_sites["MIN_AGE"]]) for site in sites]), reverse=True)
    mxt = max([s[0] for s in site_times])
    mnt = min([s[1] for s in site_times])

    main_fig, ax1 = prep_fig(pstyle)
    for si, slce in enumerate(slices):
        if slce[1] > mxt or slce[0] < mnt:
            continue
        if "-" in slce[-1] and which == "NA":
            plt.plot([-slce[0], -slce[0]], [-1, len(site_times)+1], 'k:', linewidth=0.75, color="#333333")
        else:
            plt.plot([-slce[0], -slce[0]], [-1, len(site_times)+1], 'k-', linewidth=0.75, color="#333333")
        plt.text(-(slce[0]+slce[1])/2, -5, slce[-1], color="k", verticalalignment="bottom", horizontalalignment="center", rotation=60)

    ml = max([len(tt) for tt in times_dict.values()])
    for si, sitet in enumerate(site_times):
        l = len(times_dict[sitet])
        color = "r"
        if any([ct.is_in_slice(sitet[0], sitet[1], s[0], s[1], max_ratio_out) for s in slices]):
            # if any([sitet[0] <= s[0] and sitet[1] >= s[1] for s in slices]):
            color = "g"
        plt.fill([-sitet[0], -sitet[1], -sitet[1], -sitet[0]], [si, si, si+0.33+0.6*l/ml, si+0.33+0.6*l/ml], color)
    finish_fig(main_fig, ax1, pstyle, fig_tslices_file)


# CO_OCCURRENCE EXPERIMENTS
#####################################

# PAIRS HISTOGRAMS


def plot_pairs_hists(slices, spc_counts, scores, tots, fig_hist_file=None, id_score=1, ylbl="Occurrences overlap"):
    pstyle = "pairs_hists"
    ID_STATUS, ID_SCORE = (0, id_score)
    all_min = numpy.min([numpy.min(s[:, ID_SCORE]) for s in scores])
    all_max = numpy.max([numpy.max(s[:, ID_SCORE]) for s in scores])

    unit_range = (all_min >= 0 and all_max <= 1)

    # bin_step = 2**-6
    bin_step = 2**-4
    eps = 10**-8
    # eps = 10**-4
    bins = numpy.arange(-bin_step, 1+2*bin_step, bin_step)
    if unit_range:
        bins[0] = 0
        bins[1] = eps
        # bins = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1, 1.05]
        hbins = bins
    else:
        # bb = []
        # for bv in bins[1:-2]:
        #     bh = numpy.percentile(numpy.concatenate([s[:, ID_SCORE] for s in scores])/all_max, 100*bv)
        #     if len(bb) == 0 or (bb[-1]- bh)**2 > 10**-5:
        #         bb.append(bh)
        # bin_step = 1/len(bb)
        # bins = numpy.concatenate([bins[:1], numpy.arange(0, 1, bin_step), bins[-2:]])
        # hbins = numpy.concatenate([bins[:1], bb, bins[-2:]])
        nb_bins = 28
        bin_step = 1/nb_bins
        bins = numpy.array([i*bin_step for i in range(nb_bins+2)])

    main_fig, ax1 = prep_fig(pstyle)
    # plt.scatter(X[:, ID_TSLC]+0.1+0.45*(1-X[:,ID_STATUS])+0.3*numpy.random.random(X.shape[0]), X[:,ID_SCORE], c=X[:,ID_STATUS], vmin=0, vmax=1, s=3, cmap="PiYG")
    for si, slce in enumerate(slices):

        nb_sites = spc_counts[-1, si]
        X = scores[si]

        plt.plot([si, si], [-0.1, 1.1], "k:", linewidth=0.75, color="#333333")
        plt.plot([si+0.5, si+0.5], [0., 1.05], "-", linewidth=0.1, color="#666666")

        mask_same = (X[:, ID_STATUS] == 1)
        mask_diff = (X[:, ID_STATUS] == 0)

        mm = -1
        if not unit_range:
            mm = int(numpy.max(X[mask_same | mask_diff, ID_SCORE]))
            stp = int(numpy.ceil(mm/nb_bins))
            hbins = numpy.array([i*stp for i in range(nb_bins+2)])
            # V = X[:, ID_SCORE]/all_max
            if stp > 1:
                plt.text(si+0.5, 1.115, "$\\times %d$" % stp, horizontalalignment="left", verticalalignment="top", fontsize=FSZ_XSML)

        V = X[:, ID_SCORE]
        c_same, _ = numpy.histogram(V[mask_same], hbins)
        c_diff, _ = numpy.histogram(V[mask_diff], hbins)
        # print("------", mm, "\t", hbins)
        # print(hbins)
        # print(c_diff)
        # print(c_same)

        non_zeros = (numpy.sum(mask_diff)-c_diff[0], numpy.sum(mask_same)-c_same[0])

        if len(tots) > 0:
            tot = tots[si]
            zeros = [tot[0]-non_zeros[0], tot[1]-non_zeros[1]]
        else:
            tot = [numpy.sum(mask_diff), numpy.sum(mask_same)]
            zeros = [c_diff[0], c_same[0]]

        if unit_range and (c_diff[0]/numpy.sum(c_diff) > 0.5):
            low_id = 1
        else:
            low_id = 0

        h_same = c_same[low_id:]/((2.1*numpy.max(c_same[low_id:])) if numpy.max(c_same[low_id:]) > 0 else 1)
        h_diff = c_diff[low_id:]/((2.1*numpy.max(c_diff[low_id:])) if numpy.max(c_diff[low_id:]) > 0 else 1)

        plt.barh(bins[low_id:-1], -h_diff, 0.95*bin_step, si+0.5, align="edge", color=C_DIFF)
        plt.barh(bins[low_id:-1], h_same, 0.95*bin_step, si+0.5, align="edge", color=C_SAME)

        if True:  # id_score == 3:
            # for sgn, c_which in [(-1, c_diff), (1, c_same)]:
            #     mmx = 0.01*numpy.max(c_which)
            #     xs = 0.5*numpy.arange(0,1,1/(mmx+1))[1:]
            #     for few_id in numpy.where((c_which > 0) & (c_which < mmx))[0]:
            #         plt.plot(0.5+si+sgn*xs[:c_which[few_id]], 0.5*bin_step+bins[few_id]*numpy.ones(c_which[few_id]), ".r", markersize=3)

            mmx = 10
            xs = 0.5*numpy.arange(0, 1, 1/(mmx+1))[1:]
            for few_id in numpy.where((c_diff < mmx) & (c_same < mmx))[0]:
                for sgn, c_which, ccl in [(-1, c_diff, C_DIFF), (1, c_same, C_SAME)]:
                    if c_which[few_id] > 0:
                        plt.plot(0.5+si+sgn*xs[:c_which[few_id]], 0.5*bin_step+bins[few_id]*numpy.ones(c_which[few_id]), ".", markersize=2, color=ccl)

        # plt.text(si+0.28, -.05, ratio_str(zeros[0], tot[0]), horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)
        # plt.text(si+0.72, -.05, ratio_str(zeros[1], tot[1]), horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)

        # plt.text(si+0.28, 1.11, ratio_str(c_diff[-1], tot[0]), horizontalalignment="center", verticalalignment="top", fontsize=FSZ_SML)
        # plt.text(si+0.72, 1.11, ratio_str(c_same[-1], tot[1]), horizontalalignment="center", verticalalignment="top", fontsize=FSZ_SML)

        # plt.text(si+0.5, -.158, "%d" % nb_sites, horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_XSML)

        plt.text(si+0.05, -.11, "%d" % tot[0], horizontalalignment="left", verticalalignment="bottom", fontsize=FSZ_XSML)
        # plt.text(si+0.95, -.08, "%d" % tot[1], horizontalalignment="right", verticalalignment="bottom", fontsize=FSZ_XSML)
        plt.text(si+0.95, -.09, "%d" % tot[1], horizontalalignment="right", verticalalignment="top", fontsize=FSZ_XSML)

    if unit_range:
        # y_ticks, y_lbls = [-.1, 0., 0.25, 0.5, 0.75, 1], ["Nb. pairs", "0.00", "0.25", "0.50", "0.75", "1.00"]
        y_ticks, y_lbls = [0., 0.25, 0.5, 0.75, 1], ["0.00", "0.25", "0.50", "0.75", "1.00"]
    else:
        # y_ticks = []
        # y_lbls = []
        # print(all_max, hbins)
        # for hi, hb in enumerate(hbins):
        #     y_ticks.append(bins[hi])
        #     y_lbls.append(hb*all_max)
        ####
        # vticks = bins[1:-1]
        # vlbls = ["%s" % (all_max*h) for h in hbins[1:-1]]
        # if len(vticks) > 12:
        #     y_ticks = [-.12, -.06]+[vticks[0]]+list(vticks[2:-1:2])+[vticks[-1]]
        #     y_lbls = ["sites", "total"]+list([vlbls[0]]+vlbls[2:-1:2])+[vlbls[-1]]
        # else:
        #     y_ticks = [-.12, -.06]+list(vticks)
        #     y_lbls = ["sites", "total"]+list(vlbls)
        sstp = int(numpy.ceil(nb_bins/12))
        y_ticks = [v/nb_bins for v in range(0, nb_bins+1, sstp)]
        y_lbls = ["%d" % v for v in range(0, nb_bins, sstp)]

    params = {"slices": slices,
              "yticks": y_ticks, "yticklabels": y_lbls, "ylim": [-.17, 1.125],
              "ylabel": ylbl, "ylabelparams": {"labelpad": -10}}
    finish_fig(main_fig, ax1, pstyle, fig_hist_file, params)


def plot_pairs_hists_cut(slices, spc_counts, scores, tots, fig_hist_file=None, id_score=1, ylbl="Occurrences overlap"):
    pstyle = "pairs_hists_cut"
    # bin_step = 2**-6
    bin_step = 2**-4
    eps = 10**-8
    # eps = 10**-4
    bins = numpy.arange(-bin_step, 1+2*bin_step, bin_step)
    bins[0] = 0
    bins[1] = eps
    # bins = [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1, 1.05]
    ID_STATUS, ID_SCORE = (0, id_score)

    INDET = False
    main_fig, ax1 = prep_fig(pstyle)
    # plt.scatter(X[:, ID_TSLC]+0.1+0.45*(1-X[:,ID_STATUS])+0.3*numpy.random.random(X.shape[0]), X[:,ID_SCORE], c=X[:,ID_STATUS], vmin=0, vmax=1, s=3, cmap="PiYG")
    for si, slce in enumerate(slices):

        nb_sites = spc_counts[-1, si]
        X = scores[si]
        tot = tots[si]
        mask_same = (X[:, ID_STATUS] == 1)
        mask_diff = (X[:, ID_STATUS] == 0)
        mask_indet = (X[:, ID_STATUS] == -1)

        plt.plot([si, si], [-0.1, 1.1], "k:", linewidth=0.75, color="#333333")
        plt.plot([si+0.5, si+0.5], [0., 1.05], "-", linewidth=0.1, color="#666666")

        c_same, _ = numpy.histogram(X[mask_same, ID_SCORE], bins)
        c_diff, _ = numpy.histogram(X[mask_diff, ID_SCORE], bins)
        c_indet, _ = numpy.histogram(X[mask_indet, ID_SCORE], bins)
        # h_same = c_same[1:]/(2*numpy.sum(c_same[1:]))
        # h_diff = c_diff[1:]/(2*numpy.sum(c_diff[1:]))
        h_same = c_same[1:]/((2.1*numpy.max(c_same[1:])) if numpy.max(c_same[1:]) > 0 else 1)
        h_diff = c_diff[1:]/((2.1*numpy.max(c_diff[1:])) if numpy.max(c_diff[1:]) > 0 else 1)
        if INDET and tot[-1] > 0:
            h_indet = c_indet[1:]/((2.1*numpy.max(c_diff[1:])) if numpy.max(c_indet[1:]) > 0 else 1)

        plt.barh(bins[1:-1], -h_diff, 0.95*bin_step, si+0.5, align="edge", color=C_DIFF)
        plt.barh(bins[1:-1], h_same, 0.95*bin_step, si+0.5, align="edge", color=C_SAME)
        if INDET and tot[-1] > 0:
            plt.barh(bins[1:-1], -h_indet, 0.95*bin_step, si+0.5, align="edge", color="#AAAAAA", alpha=0.6)

        non_zeros = (numpy.sum(mask_diff)-c_diff[0], numpy.sum(mask_same)-c_same[0], numpy.sum(mask_indet)-c_indet[0])
        zeros = [tot[0]-non_zeros[0], tot[1]-non_zeros[1], tot[-1]-non_zeros[-1]]

        if INDET and tot[-1] > 0:
            plt.text(si+0.33, 1.14, ratio_str(c_indet[-1], tot[-1]), horizontalalignment="center", verticalalignment="top", fontsize=FSZ_SML)
        plt.text(si+0.28, 1.11, ratio_str(c_diff[-1], tot[0]), horizontalalignment="center", verticalalignment="top", fontsize=FSZ_SML)
        plt.text(si+0.72, 1.11, ratio_str(c_same[-1], tot[1]), horizontalalignment="center", verticalalignment="top", fontsize=FSZ_SML)

        plt.text(si+0.5, -.15, "%d" % nb_sites, horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)

        if INDET and tot[-1] > 0:
            plt.text(si+0.33, -.06, ratio_str(zeros[-1], tot[-1]), horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)
        plt.text(si+0.28, -.05, ratio_str(zeros[0], tot[0]), horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)
        plt.text(si+0.72, -.05, ratio_str(zeros[1], tot[1]), horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)

        if INDET and tot[-1] > 0:
            plt.text(si+0.33, -.12, "%d" % tot[-1], horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)
        plt.text(si+0.28, -.10, "%d" % tot[0], horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)
        plt.text(si+0.72, -.10, "%d" % tot[1], horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)

    params = {"slices": slices,
              "yticks": [-.12, -.06, 0., 0.25, 0.5, 0.75, 1],
              "yticklabels": ["sites", "total", "0.00", "0.25", "0.50", "0.75", "1.00"],
              "ylim": [-.16, 1.125],
              "ylabel": ylbl, "xlabel": "Time slices"}
    finish_fig(main_fig, ax1, pstyle, fig_hist_file, params)


def plot_div_complex(which, slices, gens, gen_counts, sliced_divs, fig_div_file=None):
    pstyle = "div_complex"
    eps = 10**-8

    obins = [0, 1, 2, 3, 6, 4000]

    bins = [0, 1, 2, 3, 6, 4000]
    # bcmap = matplotlib.cm.get_cmap("viridis")
    bcmap = matplotlib.cm.get_cmap("Blues")
    bcolors = [bcmap(1-x/bins[-2]) for x in bins[:-1]]

    fbins = [0, eps, 0.1, 0.25, 0.5, 0.66, 1-eps, 1]
    fcmap = matplotlib.cm.get_cmap("inferno")
    fcolors = [fcmap(xi/(len(fbins)-1)) for xi, x in enumerate(fbins)]

    if which == "NA":
        cbins = [0, 1, 2, 4, 8, 16, 25, 100]
    else:
        cbins = [0, 1, 2, 5, 10, 25, 50, 100]
    # cbins = [0, 1, 5, 25, 100]
    ccmap = matplotlib.cm.get_cmap("Greens")
    # ccmap = matplotlib.cm.get_cmap("copper")
    ccolors = [ccmap(1-xi/(len(cbins)-1)) for xi, x in enumerate(cbins)]
    # ccolors = [ccmap(x/cbins[-1]) for x in cbins]

    ID_GENK, ID_SITE, ID_CDET = (0, 1, 2)

    main_fig, ax1 = prep_fig(pstyle)
    max_non_z = numpy.max(numpy.sum(gen_counts[:-1, :] > 0, axis=0))
    top = 1.1*max_non_z
    for si, slce in enumerate(slices):
        # print(">>> SLICE", slce[-1])
        plt.plot([si, si], [-0.1, top], "k:", linewidth=0.75, color="#333333")

        # nb of species per genus
        gc, _ = numpy.histogram(gen_counts[:-1, si], bins)
        # print("\tNb genera\t", gc)
        nb_gen = numpy.sum(gc[1:])
        nb_gen_mul = numpy.sum(gc[2:])
        nb_gen_div = 0
        gc = gc[::-1]
        # gc =  gc/max_non_z
        xs = (si+0.15)*numpy.ones(len(bins)-2)
        bs = numpy.concatenate([[0], numpy.cumsum(gc)[:-2]])
        plt.bar(xs, gc[:-1], width=0.2, bottom=bs, align="center", color=bcolors, alpha=0.6)

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
            h_cc, _ = numpy.histogram(1 - collect[:, 2]/collect[:, 3], fbins)
            # print("\tRatio diverse\t", h_cc)
            h_cc = h_cc[::-1]
            xs = (si+0.5)*numpy.ones(len(fbins)-1)
            bs = numpy.concatenate([[0], numpy.cumsum(h_cc)[:-1]])
            plt.bar(xs, h_cc, width=0.4, bottom=bs, align="center", color=fcolors[1:])

            # nb of occurrence sites
            h_cc, _ = numpy.histogram(collect[:, 3], cbins)
            # print("\tNb sites\t", h_cc)
            h_cc = h_cc[::-1]
            xs = (si+0.85)*numpy.ones(len(cbins)-1)
            bs = numpy.concatenate([[0], numpy.cumsum(h_cc)[:-1]])
            plt.bar(xs, h_cc, width=0.2, bottom=bs, align="center", color=ccolors[1:], alpha=0.6)

        plt.text(si+0.05, nb_gen+.5, "%d" % nb_gen, horizontalalignment="left", verticalalignment="bottom", fontsize=FSZ_SML)
        plt.text(si+0.95, nb_gen_mul+.5, "[%d]" % gen_counts[-1, si], horizontalalignment="right", verticalalignment="bottom", fontsize=FSZ_SML)
        plt.text(si+0.5, nb_gen_div+.5, "%d" % nb_gen_div, horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_SML)

    # LEGEND
    hb = 0.4*max_non_z/len(slices)
    nb_lgd = 3
    marg = 0.05
    marg_w = 0.05*len(slices)
    width = (len(slices) - (nb_lgd+1)*marg_w)/nb_lgd

    # nb species
    # bws = numpy.diff(bins)
    bws = numpy.ones(len(bins)-1)
    btot = numpy.sum(bws[:-1])
    ws = width*bws[:-1]/btot
    bs = marg_w+numpy.concatenate([[0], numpy.cumsum(ws)])
    plt.barh(1.2*max_non_z, ws, height=hb, left=bs[:-1], align="center", color=bcolors[:-1][::-1], alpha=0.6)
    for vi, v in enumerate(bins[1:-1]):
        plt.text(bs[vi], 1.2*max_non_z-hb, "%s" % v, verticalalignment="top", horizontalalignment="center")
    plt.text(.5*width+1*marg_w, 1.2*max_non_z+hb, "Nb. species", verticalalignment="bottom", horizontalalignment="center")

    # diverse
    bws = numpy.ones(len(fbins)-1)
    bws[[0, -1]] = 0.1
    btot = numpy.sum(bws)
    ws = width*bws/btot
    bs = width+2*marg_w+numpy.concatenate([[0], numpy.cumsum(ws)])
    plt.barh(1.2*max_non_z, ws, height=hb, left=bs[:-1], align="center", color=fcolors[::-1])
    for vi, v in enumerate(fbins):
        if vi != 1 and vi != (len(fbins)-2):
            plt.text(bs[vi], 1.2*max_non_z-hb, "%.2f" % v, verticalalignment="top", horizontalalignment="center")
    plt.text(1.5*width+2*marg_w, 1.2*max_non_z+hb, "Ratio of occurrence sites with multiple species", verticalalignment="bottom", horizontalalignment="center")

    # nb sites
    # bws = numpy.diff(cbins)
    bws = numpy.ones(len(cbins)-1)
    btot = numpy.sum(bws[:-1])
    ws = width*bws[:-1]/btot
    bs = 2*width+3*marg_w+numpy.concatenate([[0], numpy.cumsum(ws)])
    plt.barh(1.2*max_non_z, ws, height=hb, left=bs[:-1], align="center", color=ccolors[:-1][::-1], alpha=0.6)
    for vi, v in enumerate(cbins[1:-1]):
        plt.text(bs[vi], 1.2*max_non_z-hb, "%s" % v, verticalalignment="top", horizontalalignment="center")
    plt.text(2.5*width+3*marg_w, 1.2*max_non_z+hb, "Nb. occurrence sites", verticalalignment="bottom", horizontalalignment="center")

    params = {"slices": slices,
              "yticks": "int",
              "ylim": [-.125, 1.3*max_non_z],
              # "xlabel": "Time slices",
              "ylabel": "Nb. genera"}
    finish_fig(main_fig, ax1, pstyle, fig_div_file, params)


def plot_div_hists(which, slices, gens, gen_counts, sliced_divs, fig_div_file=None):
    pstyle = "div_hists"
    eps = 10**-8

    obins = [0, 1, 2, 3, 6, 4000]
    bins = [0, 1, 2, 3, 6, 4000]
    fbins = [0, eps, 0.1, 0.25, 0.5, 0.66, 1-eps, 1]
    fcmap = matplotlib.cm.get_cmap("inferno")
    fcolors = [fcmap(xi/(len(fbins)-1)) for xi, x in enumerate(fbins)]

    ID_GENK, ID_SITE, ID_CDET = (0, 1, 2)

    main_fig, ax1 = prep_fig(pstyle)
    max_non_z = numpy.max(numpy.sum(gen_counts[:-1, :] > 0, axis=0))
    top = 1.1*max_non_z
    for si, slce in enumerate(slices):

        # nb of species per genus
        gc, _ = numpy.histogram(gen_counts[:-1, si], bins)
        # print("\tNb genera\t", gc)
        nb_gen = numpy.sum(gc[1:])
        nb_gen_mul = numpy.sum(gc[2:])
        nb_gen_div = 0
        gc = gc[::-1]
        # gc =  gc/max_non_z
        xs = (si+0.15)*numpy.ones(len(bins)-2)
        bs = numpy.concatenate([[0], numpy.cumsum(gc)[:-2]])

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
            nb_gen_tot = sum(gen_counts[:-1, si] > 0)
            nb_gen_div = numpy.sum(collect[:, 2] < collect[:, 3])
            h_cc, _ = numpy.histogram(1 - collect[:, 2]/collect[:, 3], fbins)
            # print("\tRatio diverse\t", h_cc)
            h_cc = h_cc[::-1]
            xs = (si+0.5)*numpy.ones(len(fbins)-1)
            bs = numpy.concatenate([[0], numpy.cumsum(h_cc)[:-1]])

            plt.bar([si+0.5], [nb_gen_tot], width=0.6, bottom=0, align="center", color=C_GREY+ALPH_MID)
            plt.bar(xs, h_cc, width=0.6, bottom=bs, align="center", color=fcolors[1:])

        plt.text(si+0.5, nb_gen_tot+.5, "%d" % nb_gen_tot, horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_XSML)
        plt.text(si+0.5, nb_gen_div+.5, "%d" % nb_gen_div, horizontalalignment="center", verticalalignment="bottom", fontsize=FSZ_XSML)

    # LEGEND
    hb = 0.4*max_non_z/len(slices)
    nb_lgd = 3
    marg = 0.05
    marg_w = 0.05*len(slices)
    width = len(slices)/2

    # diverse
    bws = numpy.ones(len(fbins)-1)
    bws[[0, -1]] = 0.1
    btot = numpy.sum(bws)
    ws = width*bws/btot
    bs = 0.5*width+numpy.concatenate([[0], numpy.cumsum(ws)])
    plt.barh(1.25*max_non_z, ws, height=hb, left=bs[:-1], align="center", color=fcolors[::-1])
    for vi, v in enumerate(fbins):
        if vi != 1 and vi != (len(fbins)-2):
            plt.text(bs[vi], 1.25*max_non_z-hb, "%.2f" % v, verticalalignment="top", horizontalalignment="center", fontsize=FSZ_XSML)
    # plt.text(width, max_non_z+hb, "Ratio of occurrence sites with multiple species", verticalalignment="bottom", horizontalalignment="center")

    params = {"slices": slices,
              "yticks": "int",
              "ylim": [0, 1.3*max_non_z],
              # "xlabel": "Time slices",
              "ylabel": "Nb. genera"}
    finish_fig(main_fig, ax1, pstyle, fig_div_file, params)


def plot_boxes_simple(collect, id_val, lgd_val, slices, fig_boxes_file=None):
    pstyle = "boxes_simple"
    nb_series = 1
    S = numpy.array(collect)

    ftsize = 9
    bwidth = 0.45
    dxs = [0]
    dxls = [0]
    colors = [C_GREY]
    xcolor = C_RED

    main_fig, ax1 = prep_fig(pstyle)
    for off in range(nb_series):
        v_org = S[0, :, id_val+off]
        v_rnds = S[1:, :, id_val+off]
        ax1.boxplot(v_rnds, positions=numpy.arange(v_rnds.shape[1])+.5+dxs[off], widths=bwidth,  # showfliers=False,
                    boxprops={"color": colors[off]}, medianprops={"linewidth": 2., "color": colors[off]},
                    flierprops={"marker": ".", "markeredgecolor": colors[off]},
                    whiskerprops={"color": colors[off]}, capprops={"color": colors[off]})
        ax1.plot(numpy.arange(v_rnds.shape[1])+.5+dxs[off], v_org, "x", markeredgecolor=xcolor, markerfacecolor=xcolor)

    bottom, top = ax1.get_ylim()
    new_bottom = bottom-0.05*(top-bottom)
    new_top = top+0.05*(top-bottom)
    line_bottom = bottom-0.045*(top-bottom)
    line_top = top+0.045*(top-bottom)

    params = {"slices": slices,
              "yticks": "int",
              "ylim": [new_bottom, new_top],
              "ylabel": lgd_val}
    finish_fig(main_fig, ax1, pstyle, fig_boxes_file, params)


def plot_multi_boxes_series(collect, id_val, lgd_val, slices, fig_multib_file=None, nb_series=1, lbl_series=None):
    pstyle = "multi_boxes_series"
    if lbl_series is not None:
        nb_series = len(lbl_series)

    S = numpy.array(collect)
    nb_rnd = len(collect)-1
    nb_rnds = [None for i in range(nb_series)]
    org_rnks = [None for i in range(nb_series)]
    org_zs = [None for i in range(nb_series)]
    if nb_series == 1:
        ftsize = 9
        bwidth = 0.45
        dxs = [0]
        dxls = [0]
        colors = [C_GREY]
        xcolor = C_RED
    elif nb_series == 2:
        ftsize = 9
        bwidth = 0.2
        dxs = [-0.15, 0.15]
        dxls = [-0.25, 0.25]
        colors = [C_DIFF, C_SAME]
        xcolor = C_BLACK
    elif nb_series == 3:
        ftsize = 8
        bwidth = 0.15
        dxs = [-0.2, 0, 0.2]
        dxls = [-0.25, 0, 0.25]
        colors = [C_GREY, C_DIFF, C_SAME]
        xcolor = C_BLACK
    elif nb_series == 4:
        ftsize = 8
        bwidth = 0.12
        dxs = [-0.3, -0.1, 0.1, 0.3]
        dxls = [-0.3, -0.1, 0.1, 0.3]
        # colors = ["#117733", "#498CC2", "#404096", "#666666"]
        colors = [C_GREY, C_DIFF, C_SAME, C_FOUR]
        xcolor = C_BLACK

    main_fig, ax1 = prep_fig(pstyle)
    if slices[-1][-1] == "Blancan-Early" and lbl_series is not None and lbl_series[0] == "CscoreAvg_all":  # cut off the top with outlier
        ax2 = ax1.twinx()
    else:
        ax2 = None

    for off in range(nb_series):
        v_org = S[0, :, id_val+off]
        org_rnks[off] = [numpy.where(numpy.argsort(S[:, xi, id_val+off]) == 0)[0][0] for xi in range(S.shape[1])]
        org_zs[off] = ["%+.2f" % ((v_org[xi]-numpy.mean(S[1:, xi, id_val+off]))/numpy.std(S[1:, xi, id_val+off])) if numpy.unique(S[1:, xi, id_val+off]).shape[0] > 1 else "---" for xi in range(S.shape[1])]
        if numpy.isnan(S[:, :, id_val+off]).any():
            nb_rnds[off] = [numpy.sum(~numpy.isnan(S[:, j, id_val+off]))-1 for j in range(S.shape[1])]
            v_rnds = [S[~numpy.isnan(S[:, j, id_val+off]), j, id_val+off][1:] for j in range(S.shape[1])]
        else:
            v_rnds = S[1:, :, id_val+off]

        if ax2 is not None:
            ax2.boxplot(v_rnds[:, [-1]], positions=[v_rnds.shape[1]-.5+dxs[off]], widths=bwidth,  # showfliers=False,
                        boxprops={"color": colors[off]}, medianprops={"linewidth": 2., "color": colors[off]},
                        flierprops={"marker": ".", "markeredgecolor": colors[off]},
                        whiskerprops={"color": colors[off]}, capprops={"color": colors[off]})
            ax2.plot(v_rnds.shape[1]-.5+dxs[off], v_org[-1], "x", markeredgecolor=xcolor, markerfacecolor=xcolor)
            v_rnds = v_rnds[:, :-1]
            v_org = v_org[:-1]

        ax1.boxplot(v_rnds, positions=numpy.arange(v_rnds.shape[1])+.5+dxs[off], widths=bwidth,  # showfliers=False,
                    boxprops={"color": colors[off]}, medianprops={"linewidth": 2., "color": colors[off]},
                    flierprops={"marker": ".", "markeredgecolor": colors[off]},
                    whiskerprops={"color": colors[off]}, capprops={"color": colors[off]})
        ax1.plot(numpy.arange(v_rnds.shape[1])+.5+dxs[off], v_org, "x", markeredgecolor=xcolor, markerfacecolor=xcolor)

    bottom, top = ax1.get_ylim()
    # if slices[-1][-1] == "Blancan-Early" and lbl_series is not None and lbl_series[0] == "CscoreAvg_all": ## cut off the top with outlier
    #     top = 7.2
    if lbl_series is not None and lbl_series[0] == "FisherQ95_all":  # cut off the top with outlier
        bottom = 0.75
    new_bottom = bottom-0.05*(top-bottom)
    new_top = top+0.05*(top-bottom)
    line_bottom = bottom-0.045*(top-bottom)
    line_top = top+0.045*(top-bottom)

    if ax2 is not None:
        ax1.plot([len(slices)-1-.05, len(slices)-1+.05]*20, numpy.linspace(new_bottom, new_top, 40), ":", color="#888888", linewidth=.5)

    # ax1.text(.5, top, "/%d" % nb_rnd, verticalalignment="center", horizontalalignment="center")
    # for i, ors in enumerate(org_rnks):
    #     for xi, v in enumerate(ors):
    #         plt.text(xi+1+dxls[i], line_bottom, "%d" % v, verticalalignment="bottom", horizontalalignment="center", fontsize=ftsize, color=colors[i])
    for i, ozs in enumerate(org_zs):
        # for xi, v in enumerate(ozs):
        #     plt.text(xi+1+dxls[i], line_bottom, v, verticalalignment="bottom", horizontalalignment="center", fontsize=ftsize, color=colors[i], rotation=90)
        if nb_rnds[i] is not None:
            for j, n in enumerate(nb_rnds[i]):
                ax1.text(j+.5+dxls[i], top, "/%d" % n, verticalalignment="center", horizontalalignment="center", fontsize=ftsize, color=colors[i])

    # if lbl_series is not None:
    #     for i, lbl in enumerate(lbl_series):
    #         plt.text(((.5+i)/len(lbl_series))*len(slices)+0.5, new_top, lbl, verticalalignment="bottom", horizontalalignment="center", fontsize=9, color=colors[i])
    #     new_top = new_top+0.06*(top-bottom)

    params = {"slices": slices,
              "yticks": "auto",
              "ax2": ax2,
              "ylim": [new_bottom, new_top],
              "ylabel": lgd_val}
    finish_fig(main_fig, ax1, pstyle, fig_multib_file, params)


# EXTRAS
#####################################
def plot_trend_values(ax, xvals, yvals, xlims=None, ylims=None, xlbl=None, ylbl=None, plbl=None):

    ax.plot(xlims, [0, 0], ":", color=ct.C_GREY, linewidth=.5)
    for si in range(len(xvals)-1):
        ax.arrow(xvals[si], yvals[si], xvals[si+1]-xvals[si], yvals[si+1]-yvals[si], color=ct.C_RED,
                 length_includes_head=True, width=.0004, head_width=0.035, head_length=0.08, capstyle="round")
    ax.plot(xvals, yvals, "o", color=ct.C_RED, markersize=1)
    ax.tick_params(axis='x', which="both", pad=3.0)

    # for si in range(len(xvals)):
    #     if plbl is None:
    #         ax.text(xvals[si], yvals[si], "%d" % si, verticalalignment="bottom", horizontalalignment="center", fontsize=ct.FSZ_XSML)
    #     else:
    #         ax.text(xvals[si], yvals[si], plbl[si], verticalalignment="bottom", horizontalalignment="center", fontsize=ct.FSZ_XSML)

    # ax.text(xl+.98*(xu-xl), yl+.01*(yu-yl), "%.3f" % numpy.corrcoef(xvals, yvals)[0, 1], color=ct.C_GREY, verticalalignment="bottom", horizontalalignment="right", fontsize=ct.FSZ_SML)  # (yu+yl)/2
    # if plbl is not None and "O" in plbl:
    #     ax.text(xl+.02*(xu-xl), yl+.98*(yu-yl), "%.3f" % numpy.corrcoef(xvals[:10], yvals[:10])[0, 1], color=ct.C_GREY, verticalalignment="top", horizontalalignment="left", fontsize=ct.FSZ_XSML)  # (yu+yl)/2
    #     ax.text(xl+.98*(xu-xl), yl+.98*(yu-yl), "%.3f" % numpy.corrcoef(xvals[8:], yvals[8:])[0, 1], color=ct.C_GREY, verticalalignment="top", horizontalalignment="right", fontsize=ct.FSZ_XSML)  # (yu+yl)/2

    if xlims is not None:
        ax.set_xlim(xlims)
    if ylims is not None:
        nylims = (ylims[0]-0.07*(ylims[1]-ylims[0]), ylims[1]+0.07*(ylims[1]-ylims[0]))
        ax.set_ylim(nylims)
    xl, xu = ax.get_xlim()
    yl, yu = ax.get_ylim()

    if xlbl is not None:
        # ax.set_title(xlbl)
        ax.set_xlabel(xlbl)
    if ylbl is not None:
        ax.set_ylabel(ylbl, multialignment='center')


def plot_tbars(ax, xvals, yvals, ylims=None, ylbl=None, xlbl=None):

    def get_bar_colors(v):
        if v < -3:
            return ct.C_SIGN
        elif v > 3:
            return ct.C_SIGN
        return ct.C_NEUTRAL
    xs = numpy.arange(len(xvals))
    cs = [get_bar_colors(v) for v in yvals]

    ax.bar(xs, yvals, color=cs)  # , ".", color=ct.C_GREY, markersize=1)
    if ylims is not None:
        for si, v in enumerate(yvals):
            if v < ylims[0]:
                ax.text(si, ylims[0], "%.1f" % v, rotation=25, verticalalignment="bottom", horizontalalignment="center", fontsize=ct.FSZ_XSML)
            elif v > ylims[1]:
                ax.text(si, ylims[1], "%.1f" % v, rotation=25, verticalalignment="top", horizontalalignment="center", fontsize=ct.FSZ_XSML)

    ax.set_xticks(xs)
    ax.set_xticklabels(xvals, rotation=25, va="center", fontsize=ct.FSZ_SML)
    ax.tick_params(axis='x', which="both", pad=8.0)
    ax.set_xlim([-0.8, len(xvals)-.2])
    if ylims is not None:
        nylims = (ylims[0]-0.07*(ylims[1]-ylims[0]), ylims[1]+0.07*(ylims[1]-ylims[0]))
        ax.set_ylim(nylims)
    if xlbl is not None:
        # ax.set_title(xlbl)
        ax.set_xlabel(xlbl)
    if ylbl is not None:
        ylbl_folded = ylbl.replace("genera with", "genera\nwith").replace("vs. null", "vs.\nnull")
        ax.set_ylabel(ylbl_folded, multialignment='center')


def plot_overlap_areas(fossils_data, spc_counts, sliced_counts, map_fields_species, map_fields_sites):
    pstyle = "overlap_areas"
    species, species_dict = (fossils_data["species"], fossils_data["species_dict"])

    gen_field = map_fields_species["GENUS"]
    spc_field = map_fields_species["SPECIES"]
    indeterminate_sids_all = set([spci for spci, s in enumerate(species) if ct.is_indeterminate_gs(s[gen_field], s[spc_field])])

    gen_tots = {}
    for ssi, s in enumerate(species):
        if s[gen_field] not in gen_tots:
            gen_tots[s[gen_field]] = 0
        gen_tots[s[gen_field]] += 1

    masks_styles = [("#B178A6", 1.), ("#882E72", 1.), ("#90C987", 1.), ("#4EB265", 1.)]

    # map_gens = dict([(v,k) for k,v in enumerate(sorted(gen_tots.keys(), key=lambda x: gen_tots[x]))])
    map_gens = dict([(v, k) for k, v in enumerate(sorted(gen_tots.keys()))])
    main_fig, ax1 = prep_fig(pstyle)
    splts = main_fig.subplots(1, len(sliced_counts), sharey=True)
    for si, scounts in enumerate(sliced_counts):
        sids = numpy.where(spc_counts[:-1, si] > 0)[0]
        full_mask = numpy.ones(len(sids), dtype=bool)
        det_mask = numpy.array([(i not in indeterminate_sids_all) for i in sids], dtype=bool)

        counts = spc_counts[sids, si]
        sc, _ = numpy.histogram(counts, bins=[1, 2, 3, 4, 5, 100])
        scd, _ = numpy.histogram(counts[det_mask], bins=[1, 2, 3, 4, 5, 100])
        splts[si].bar(numpy.arange(sc.shape[0]), sc, alpha=0.7, color="#882E72")
        splts[si].bar(numpy.arange(scd.shape[0]), scd, color="#882E72")
        splts[si].set_xticks([0, 2, 4])
        splts[si].set_xticklabels(["1", "3", "5+"])
        # gen_counts = {}
        # for ssi, sid in enumerate(sids):
        #     if species[sid][gen_field] not in gen_counts:
        #         gen_counts[species[sid][gen_field]] = numpy.zeros(len(sids), dtype=bool)
        #     gen_counts[species[sid][gen_field]][ssi] = 1
        # all_nb = len(sids)
        # map_sids = dict([(v,k) for (k,v) in enumerate(sids)])
        # C = numpy.zeros((len(sids), len(sids)), dtype=int)
        # for (i,j), c in scounts.items():
        #     C[map_sids[i], map_sids[j]] = c
        #     C[map_sids[j], map_sids[i]] = c

        # areas = numpy.zeros((C.shape[0], 10))
        # for ci in range(C.shape[0]):
        #     ccts_store = []
        #     sums_store = []
        #     areas[ci, -1] = det_mask[ci]
        #     areas[ci, -2] = map_gens[species[sids[ci]][gen_field]]
        #     sg_mask = gen_counts[species[sids[ci]][gen_field]]
        #     masks = [full_mask, det_mask, sg_mask, det_mask & sg_mask]

        #     for mi, mask in enumerate(masks):
        #         ccts = numpy.cumsum(numpy.bincount(C[ci, mask])[:0:-1])[::-1]
        #         sums_store.append(numpy.sum(mask))
        #         ccts_store.append(ccts)
        #         if ccts.shape[0] > 0:
        #             dts = ccts/(spc_counts[sids[ci], si]*(numpy.sum(mask)-mask[ci]))
        #             dtsR = ccts/(spc_counts[sids[ci], si]*all_nb)
        #             areas[ci, mi] = numpy.sum(dts)
        #             areas[ci, 4+mi] = numpy.sum(dtsR)
        #         elif (numpy.sum(mask)-mask[ci]) == 0:
        #             areas[ci, mi] = -.1
        #             areas[ci, 4+mi] = -.1

        #     if areas[ci,2] > 0.75 and spc_counts[sids[ci], si] == 1:
        #         area_str = " / ".join(["%.3f" % v for v in areas[ci, :4]])
        #         areaR_str = " / ".join(["%.3f" % v for v in areas[ci, 4:-2]])
        #         sums_str = " / ".join(["%d" % v for v in sums_store])
        #         print("%s %s (%d):\t%s\t%s\t%s" % (species[sids[ci]][gen_field], species[sids[ci]][spc_field], si, area_str, areaR_str, sums_str))
        #     elif areas[ci,2] > 2: #0.75:
        #         plt.figure(figsize=(3.8, 3.8))
        #         for mi, ccts in enumerate(ccts_store):
        #             plt.bar(numpy.arange(ccts.shape[0]), ccts, width=1, color=masks_styles[mi][0], align="edge", alpha=masks_styles[mi][1])
        #         plt.xlim([0, spc_counts[sids[ci], si]])
        #         plt.ylim([0, all_nb])
        #         plt.xticks(numpy.arange(spc_counts[sids[ci], si]+1))
        #         area_str = " / ".join(["%.3f" % v for v in areas[ci, :4]])
        #         areaR_str = " / ".join(["%.3f" % v for v in areas[ci, 4:-2]])
        #         sums_str = " / ".join(["%d" % v for v in sums_store])
        #         plt.text(spc_counts[sids[ci], si]/2, 0.952*all_nb, area_str, verticalalignment="bottom", horizontalalignment="center")
        #         plt.text(spc_counts[sids[ci], si]/2, 0.948*all_nb, areaR_str, verticalalignment="top", horizontalalignment="center")
        #         plt.text(spc_counts[sids[ci], si]/2, 0.9*all_nb, sums_str, verticalalignment="top", horizontalalignment="center")
        #         plt.xlabel("Nb. of localities")
        #         plt.ylabel("Nb. of species")
        #         plt.title("%s %s (%d)" % (species[sids[ci]][gen_field], species[sids[ci]][spc_field], si))
        #         plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.94)

        # splts[si].scatter(areas[:,0], areas[:,2], s=10+5*areas[:,-1], c=areas[:,-2], alpha=0.6, cmap="rainbow", vmin=0, vmax=len(map_gens)-1)
        # splts[si].set_ylim([-.2, 1.1])
        if si == 0:
            splts[si].set_ylabel("Nb. species")
        # #    splts[si].set_yticks([])
    finish_fig(main_fig, ax1, pstyle)

# ### GEO POLYGONS
# from shapely.geometry import Point, MultiPoint
# def plot_species_polygon(sid, occurrences, sites_coords, C , C_points):
#     spc_coords = sites_coords[numpy.where(occurrences[:,sid])[0],:]
#     P = MultiPoint(spc_coords).convex_hull
#     contained = [ci for ci, c in enumerate(C_points) if P.contains(c)]

#     plt.plot(C[:,0], C[:,1], 'k.')
#     plt.plot(spc_coords[:,0], spc_coords[:,1], 'bo')
#     plt.plot(C[contained,0], C[contained,1], 'r.')
