import argparse
import numpy

import common_tools as ct

import pdb

# Example usage: run randomization experiments and store the stats for North America, large mammals, 1000 data copies randomized with the Curveball algorithm, computing the number of genera with co-occurring species
# python run_rnd.py  --continent NA --group L --rnd_nb 1000 --null_model CB --measure gen

SEED_MAX = 2**32

params = {"continent": ["NA", "EU"],
          "group": ["L", "S", "C"],
          "null_model": ["CB", "UG", "shuffle"],
          "measure": ["gen", "pairs"]}
hlp = {"continent": "continent data subset",
       "group": "ecological group data subset",
       "null_model": "null-model randomization algorithm",
       "measure": "co-occurrence measure type"}

parser = argparse.ArgumentParser(description='Perform co-occurrence randomization experiments.')
parser.add_argument("-r", "--rnd_nb", type=int, help="number of randomized data copies to generate", default=10)
parser.add_argument("-k", "--keep_nb", type=int, help="number of randomized data copies to generate to store", default=2)
parser.add_argument("-d", "--data_folder", type=str, help="folder containing the data subsets", default="../prepared_data/")
parser.add_argument("-t", "--time_folder", type=str, help="folder containing the time bins specifications", default="../times/")
parser.add_argument("-x", "--xps_folder", type=str, help="folder to store the raw results", default="../xps_rnd/")
parser.add_argument("-z", "--xps_suffix", type=str, help="series suffix for the raw results", default="")
for p in params.keys():
    parser.add_argument("-"+p[0], "--"+p, type=str, choices=list(params[p]), action='append', default=argparse.SUPPRESS, help=hlp[p])

pargs = vars(parser.parse_args())
for k, v in pargs.items():
    params[k] = v

RND_FLD = params["xps_folder"]
ct.make_fld(RND_FLD)

for WHICH in params["continent"]:
    for GROUP in params["group"]:
        print("### %s-%s" % (WHICH, GROUP))

        #######################################
        # OUTPUT FILES
        SSUFF = "%s-%s%s" % (WHICH, GROUP, params["xps_suffix"])

        DATA_STATB_FILE = RND_FLD+"statB-%s_" + SSUFF+"_%s.csv"
        DATA_RND_FILE = RND_FLD+"fossils_"+SSUFF+"_%s_%d.csv"
        DATA_PROBAS_FILE = None  # RND_FLD+"pairs_"+SSUFF+"_%s.csv"
        SEED_STORE_FILE = RND_FLD+"seeds_"+SSUFF+".csv"

        #######################################
        # LOADING THE DATA
        DATA_FLD = params["data_folder"]
        TIMES_FLD = params["time_folder"]

        if WHICH == "NA":
            LIST_FILE = DATA_FLD+"fossils_NorthAmerica-%s.csv" % GROUP
            SLICES_FILE = TIMES_FLD+"time-slices_america_filter.csv"
            BOUNDARIES = {"LAT_MIN": 19, "LAT_MAX": 84,
                          "LONG_MIN": -140, "LONG_MAX": -70,
                          "MIN_AGE_MIN": 2.5, "MAX_AGE_MAX": 34}
            MAX_RATIO_OUT = 0.1
        else:
            WHICH = "EU"
            LIST_FILE = DATA_FLD+"fossils_Europe-wide-%s.csv" % GROUP
            SLICES_FILE = TIMES_FLD+"time-slices_europe_filter.csv"
            BOUNDARIES = {"LAT_MIN": 14, "LAT_MAX": 82,
                          "LONG_MIN": -24, "LONG_MAX": 75,
                          "MIN_AGE_MIN": 2.5, "MAX_AGE_MAX": 23}
            MAX_RATIO_OUT = 0.1

        FIELDS_SITES = ["LIDNUM", "NAME", "LAT", "LONG", "MAX_AGE", "MIN_AGE", "SLICE_ID", "SLICE_NAME", "MEAN_HYPSODONTY"]
        # FIELDS_SITES = ["LAT","LONG","MAX_AGE","MIN_AGE","SLICE_ID","SLICE_NAME"]
        MAP_FIELDS_SITES = dict([(v, k) for (k, v) in enumerate(FIELDS_SITES)])
        # FIELDS_SPECIES = ["ORDER","FAMILY","GENUS","SPECIES"] #, "UNIQUE"]
        FIELDS_SPECIES = ["ORDER", "SUBORDERORSUPERFAMILY", "FAMILY", "SUBFAMILY", "GENUS", "SPECIES"]
        MAP_FIELDS_SPECIES = dict([(v, k) for (k, v) in enumerate(FIELDS_SPECIES)])
        FIELDS_FORMAT = {"LAT": float, "LONG": float, "MAX_AGE": float, "MIN_AGE": float, "SLICE_ID": int}

        fossils_data = ct.read_fossils(LIST_FILE, SLICES_FILE, FIELDS_FORMAT, FIELDS_SITES, FIELDS_SPECIES)
        marg_species = numpy.sum(fossils_data["occurrences"], axis=0)
        sort_species = numpy.argsort(marg_species)
        marg_sites = numpy.sum(fossils_data["occurrences"], axis=1)
        sort_sites = numpy.argsort(marg_sites)

        #######################################
        # CHOICE OF CO-OCCURRENCE STATISTICS
        TYPE_MEASURES, PAIR_MEASURES, PAIR_MEASURES_MAP, VALUES, PAIR_STATS, PAIR_AGG_SERIES = ct.get_statistics_details(GROUP)
        type_measure = params["measure"]

        #######################################
        # CHOICE NULL OF MODELS
        nb_rnd = params["rnd_nb"]
        rnd_subseries = ct.get_nullmodels_details()
        rnd_series = []
        for r in params["null_model"]:
            rnd_series.extend(rnd_subseries[r])

        org_stats = {}
        spc_counts, sliced_counts = None, None
        for ti, typ in enumerate(type_measure):

            if typ == "gen":
                # ORIGINAL DATA, GENERA-BASED
                gens, gen_counts, sliced_divs = ct.compute_diverse(fossils_data, FIELDS_SPECIES, MAP_FIELDS_SITES)
                if nb_rnd > 0:
                    org_stats[typ] = ct.compute_div_stat(WHICH, fossils_data["slices"], gens, gen_counts, sliced_divs)
                    for vi, vname, vlgd, vfmt, plot_box in VALUES[typ]:
                        ct.save_stats_boxes([org_stats[typ]], vi, vlgd, vfmt, DATA_STATB_FILE % (vname, "original"))

            elif typ == "pairs":
                # ORIGINAL DATA, PAIR-BASED
                if spc_counts is None:
                    spc_counts, sliced_counts, sites_counts = ct.compute_counts(fossils_data, FIELDS_SPECIES, MAP_FIELDS_SITES)
                scores, tots = ct.compute_pairs_scores(fossils_data["species"], FIELDS_SPECIES, spc_counts, sliced_counts, sites_counts, PAIR_MEASURES)

                if DATA_PROBAS_FILE is not None:
                    for xi in range(len(scores)):
                        scs = scores[xi]
                        rids = numpy.argsort(scs[:, 3])
                        with open(DATA_PROBAS_FILE % xi, "w") as fo:
                            fo.write(",".join(ct.get_pair_statusvars()+PAIR_MEASURES+ct.get_pair_extravars()[:-2]+["genj", "spcj", "geni", "spci"])+"\n")
                            for rid in rids:
                                fo.write(",".join(["%s" % s for s in scs[rid, :-2]] +
                                                  [fossils_data["species"][int(scs[rid, -2])][MAP_FIELDS_SPECIES["GENUS"]],
                                                   fossils_data["species"][int(scs[rid, -2])][MAP_FIELDS_SPECIES["SPECIES"]],
                                                   fossils_data["species"][int(scs[rid, -1])][MAP_FIELDS_SPECIES["GENUS"]],
                                                   fossils_data["species"][int(scs[rid, -1])][MAP_FIELDS_SPECIES["SPECIES"]]])+"\n")

                if nb_rnd > 0:
                    org_stats[typ] = ct.compute_pairs_stats(fossils_data["slices"], spc_counts, scores, tots, stats_details=PAIR_STATS, grp=GROUP)
                    for vi, vname, vlgd, vfmt, plot_box in VALUES[typ]:
                        for off, gname in enumerate(PAIR_AGG_SERIES):
                            ct.save_stats_boxes([org_stats[typ]], vi+off, vlgd, vfmt, DATA_STATB_FILE % (vname+"_"+gname, "original"))

        # randomization part
        if nb_rnd > 0 and len(type_measure) > 0:

            # generate and store seeds
            seeds_series = numpy.random.randint(SEED_MAX, size=len(rnd_series))
            with open(SEED_STORE_FILE, "w") as fo:
                fo.write("\n".join(["%s\t%d" % (rndk, seeds_series[rndi]) for rndi, rndk in enumerate(rnd_series)]))

            for rndi, rndk in enumerate(rnd_series):
                numpy.random.seed(seeds_series[rndi])

                rnd_parameters = ct.make_rnd_series(rndk)
                print("--- RNDS", "/".join([x for x in type_measure]), SSUFF, rndk)

                if rnd_parameters.get("marg_sites", 0) is None:
                    rnd_parameters["marg_sites"] = marg_sites

                collect = {}
                for typ in type_measure:
                    if typ == "marg":
                        collect[typ] = [ct.compute_marg_stat(fossils_data)]
                    else:
                        collect[typ] = [org_stats[typ]]

                for i in range(nb_rnd):
                    if i % 10 == 0:
                        print("-- Random %d" % i)

                    fossils_data_rnd, rnd_details = ct.randomize_records(fossils_data, FIELDS_SITES, FIELDS_SPECIES, rnd_parameters)
                    if i < params["keep_nb"]:
                        ct.write_fossils(fossils_data_rnd, DATA_RND_FILE % (rndk, i), FIELDS_SITES, FIELDS_SPECIES)

                    for ti, typ in enumerate(type_measure):

                        if typ == "marg":
                            collect[typ].append(ct.compute_marg_stat(fossils_data_rnd, fossils_data))

                        ### RND, GENERA-BASED
                        if typ == "gen":
                            gens, gen_counts, sliced_divs = ct.compute_diverse(fossils_data_rnd, FIELDS_SPECIES, MAP_FIELDS_SITES)
                            collect[typ].append(ct.compute_div_stat(WHICH, fossils_data["slices"], gens, gen_counts, sliced_divs))

                        ### RND, PAIR-BASED
                        elif typ == "pairs":
                            spc_counts, sliced_counts, sites_counts = ct.compute_counts(fossils_data_rnd, FIELDS_SPECIES, MAP_FIELDS_SITES)
                            scores, tots = ct.compute_pairs_scores(fossils_data["species"], FIELDS_SPECIES, spc_counts, sliced_counts, sites_counts, PAIR_MEASURES)
                            collect[typ].append(ct.compute_pairs_stats(fossils_data["slices"], spc_counts, scores, tots, stats_details=PAIR_STATS, grp=GROUP))

                # SAVE AND PLOT RND STATS
                for typ in type_measure:
                    if typ in ["gen", "marg"]:
                        for vi, vname, vlgd, vfmt, plot_box in VALUES[typ]:
                            ct.save_stats_boxes(collect[typ], vi, vlgd, vfmt, DATA_STATB_FILE % (vname, rndk))

                    elif typ == "pairs":
                        for vi, vname, vlgd, vfmt, plot_box in VALUES[typ]:
                            for off, gname in enumerate(PAIR_AGG_SERIES):
                                ct.save_stats_boxes(collect[typ], vi+off, vlgd, vfmt, DATA_STATB_FILE % (vname+"_"+gname, rndk))
