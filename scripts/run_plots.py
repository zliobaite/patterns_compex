import argparse
import numpy
import matplotlib.pyplot as plt

import common_tools as ct
import plotting_tools as pt

import pdb

# Example usage: prepare figures for North America, large mammals, randomization with Curveball, computing the number of genera with co-occurring species, experiments steps 1 and 2 (histograms and boxplots)
# python run_plots.py  --continent NA --group L --null_model CB --measure gen --step 1 --step 2
# other example
# python run_plots.py --data_folder ../prepared_1125/ --xps_folder ../rnd_1000-LSC/ --figs_folder ../figs_1000-LSC/

nb_steps = 3
params = {"continent": ["NA", "EU"],
          "group": ["L", "S", "C"],
          "null_model": ["CB", "UG", "shuffle"],
          "measure": ["gen", "pairs"],
          "vars_context": ["avg_hyp", "nb_species"]}  # , "avg_spc_gen", "ent_spc_gen", "nb_sites", "ratio_indet"]

hlp = {"continent": "continent data subset",
       "group": "ecological group data subset",
       "null_model": "null-model randomization algorithm",
       "measure": "co-occurrence measure type",
       "vars_context": "environmental context variables"}

parser = argparse.ArgumentParser(description='Perform co-occurrence randomization experiments.')
parser.add_argument("-d", "--data_folder", type=str, help="folder containing the data subsets", default="../prepared_data/")
parser.add_argument("-t", "--time_folder", type=str, help="folder containing the time bins specifications", default="../times/")
parser.add_argument("-x", "--xps_folder", type=str, help="folder containing the raw results", default="../xps_rnd/")
parser.add_argument("-z", "--xps_suffix", type=str, help="series suffix for the raw results", default="")
parser.add_argument("-f", "--figs_folder", type=str, help="folder to store the figures", default="../figs/")
parser.add_argument("-s", "--step", type=int, choices=list(range(nb_steps+1)), action='append', default=argparse.SUPPRESS, help="experiments step(s) figures should be plotted for")
group = parser.add_mutually_exclusive_group()
group.add_argument('--bbl', action='store_true', dest="plot_bbl", help="make bbl plots (step 0)", default=argparse.SUPPRESS)
group.add_argument('--no-bbl', action='store_false', dest="plot_bbl", default=argparse.SUPPRESS)
group = parser.add_mutually_exclusive_group()
group.add_argument('--svg', action='store_true', dest="svg_bbl", help="make svg bbl plots (step 0)", default=argparse.SUPPRESS)
group.add_argument('--no-svg', action='store_false', dest="svg_bbl", default=argparse.SUPPRESS)

for p in params.keys():
    parser.add_argument("-"+p[0], "--"+p, type=str, choices=list(params[p]), action='append', default=argparse.SUPPRESS, help=hlp[p])

pargs = vars(parser.parse_args())
for k, v in pargs.items():
    params[k] = v

DO_STEP = {}
for i in range(nb_steps+1):
    DO_STEP[i] = i in params.get("step", [i])

RND_FLD = params["xps_folder"]
FIGS_FLD = params["figs_folder"]
ct.make_fld(FIGS_FLD)

rnd_subseries = ct.get_nullmodels_details()
rnd_series = []
for r in params["null_model"]:
    rnd_series.extend(rnd_subseries[r])
context_vars = params["vars_context"]

#######################################
# CHOICE OF CO-OCCURRENCE STATISTICS
TYPE_MEASURES, PAIR_MEASURES, PAIR_MEASURES_MAP, VALUES, PAIR_STATS, PAIR_AGG_SERIES = ct.get_statistics_details("L")
VAL_KEYS = {}
for typ in TYPE_MEASURES:
    for vi, v in enumerate(VALUES[typ]):
        VAL_KEYS[v[1]] = (typ, vi, -1*(typ == "pairs"))
for off, gname in enumerate(PAIR_AGG_SERIES):
    for vi, v in enumerate(VALUES["pairs"]):
        VAL_KEYS[v[1]+"_"+gname] = ("pairs", vi, off)

VALUES_LBLS = {"min": "Occurrences overlap [min]", "jacc": "Occurrences overlap [Jacc]",
               "fisher": "Fisher mid-P val", "cscore": "C-score",
               "CscoreAvg": "C-score [average]", "FisherQ95": "Fisher mid-P val [95 percentile]",
               "avg_hyp": "Mean hyp", "nb_species": "Nb. species",
               "avg_spc_gen": "Avg. number of species per genus",
               "ent_spc_gen": "Entropy of the distribution of species over genera",
               "nb_sites": "Number of localities",
               "ratio_indet": "Ratio of undeterminate species"}

params_series = {}
params_series["box_gen"] = {}
params_series["box_gen"]["rnd_series"] = rnd_series
params_series["box_gen"]["boxes_values"] = ["nbCooc"]

params_series["trend_gen"] = {}
params_series["trend_gen"]["rnd_series"] = rnd_series
params_series["trend_gen"]["trends_rows"] = [["nbCooc"]]
params_series["trend_gen"]["trends_cols"] = context_vars

params_series["box_pairs"] = {}
params_series["box_pairs"]["rnd_series"] = rnd_series
params_series["box_pairs"]["boxes_values"] = [
    ["CscoreAvg", "CscoreAvg_diff", "CscoreAvg_same"],
    ["FisherQ95", "FisherQ95_diff", "FisherQ95_same"]]
# params_series["box_pairs"]["boxes_values"] = [
#     ["CscoreAvg", "CscoreAvg_all", "CscoreAvg_diff", "CscoreAvg_same"],
#     ["FisherQ95", "FisherQ95_all", "FisherQ95_diff", "FisherQ95_same"]]

params_series["trend_pairs"] = {}
params_series["trend_pairs"]["rnd_series"] = rnd_series
params_series["trend_pairs"]["trends_rows"] = [["CscoreAvg_diff", "CscoreAvg_same"], ["FisherQ95_diff", "FisherQ95_same"]]
params_series["trend_pairs"]["trends_cols"] = context_vars

params_series_keys = []
for pk in params_series.keys():
    ptyp, mtyp, *_ = pk.split("_")
    if (DO_STEP[2] and ptyp == "box") or (DO_STEP[3] and ptyp == "trend"):
        if mtyp in params["measure"]:
            params_series_keys.append(pk)

c_org = {}
for WHICH in params["continent"]:
    for GROUP in params["group"]:
        print("### %s-%s" % (WHICH, GROUP))

        #######################################
        # OUTPUT FILES
        SSUFF = "%s-%s%s" % (WHICH, GROUP, params["xps_suffix"])

        CSV_OUTLINE_FILE = FIGS_FLD+"outline_"+SSUFF+".csv"
        TEX_LOCALITIES_FILE = FIGS_FLD+"localities_"+SSUFF+".tex"
        PLOTLY_VARS_FILE = FIGS_FLD+"extra_vars.js"
        FIG_OSIMPLE_FILE = FIGS_FLD+"outline_simple_"+SSUFF+".pdf"
        # TAB_CNTO_FILE = FIGS_FLD+"cnt_"+SSUFF+".tex"
        TAB_CNTO_FILE = FIGS_FLD+"counts_"+SSUFF+".tex"
        FIG_CNTS_FILE = FIGS_FLD+"cnt-species_"+SSUFF+".pdf"
        FIG_CNTO_FILE = FIGS_FLD+"cnt-order_"+SSUFF+".pdf"
        FIG_BBL_FILE = FIGS_FLD+"bbl-genO_"+SSUFF+".pdf"
        FIG_BBLD_FILE = FIGS_FLD+"bbl-genD_"+SSUFF+".pdf"
        FIG_BBL_SVG = FIGS_FLD+"bbl-genO_"+SSUFF+".svg"

        FIG_HIST_FILE = FIGS_FLD+"hist-%s_"+SSUFF+".pdf"

        #######################################
        # LOADING THE DATA
        DATA_FLD = params["data_folder"]
        TIMES_FLD = params["time_folder"]
        TEETH_FILE = DATA_FLD+"teeth_details.csv"  # for bubble plots indicating the hypsodonty type (bbl-genD), ignored otherwise

        if WHICH == "NA":
            LIST_FILE = DATA_FLD+"fossils_NorthAmerica-%s.csv" % GROUP
            SLICES_FILE = TIMES_FLD+"time-slices_america_filter.csv"
            BOUNDARIES = {"LAT_MIN": 19, "LAT_MAX": 84,
                          "LONG_MIN": -140, "LONG_MAX": -70,
                          "MIN_AGE_MIN": 2.5, "MAX_AGE_MAX": 34}
            MAX_RATIO_OUT = 0.1
        elif WHICH == "EU":
            LIST_FILE = DATA_FLD+"fossils_Europe-wide-%s.csv" % GROUP
            SLICES_FILE = TIMES_FLD+"time-slices_europe_filter.csv"
            BOUNDARIES = {"LAT_MIN": 14, "LAT_MAX": 82,
                          "LONG_MIN": -24, "LONG_MAX": 75,
                          "MIN_AGE_MIN": 2.5, "MAX_AGE_MAX": 23}
            MAX_RATIO_OUT = 0.1
        else:
            exit()

        FIELDS_SITES = ["LIDNUM", "NAME", "LAT", "LONG", "MAX_AGE", "MIN_AGE", "SLICE_ID", "SLICE_NAME", "MEAN_HYPSODONTY"]
        # FIELDS_SITES = ["LAT","LONG","MAX_AGE","MIN_AGE","SLICE_ID","SLICE_NAME"]
        MAP_FIELDS_SITES = dict([(v, k) for (k, v) in enumerate(FIELDS_SITES)])
        # FIELDS_SPECIES = ["ORDER","FAMILY","GENUS","SPECIES"] #, "UNIQUE"]
        FIELDS_SPECIES = ["ORDER", "SUBORDERORSUPERFAMILY", "FAMILY", "SUBFAMILY", "GENUS", "SPECIES"]
        MAP_FIELDS_SPECIES = dict([(v, k) for (k, v) in enumerate(FIELDS_SPECIES)])
        FIELDS_FORMAT = {"LAT": float, "LONG": float, "MAX_AGE": float, "MIN_AGE": float, "SLICE_ID": int}

        fossils_data = ct.read_fossils(LIST_FILE, SLICES_FILE, FIELDS_FORMAT, FIELDS_SITES, FIELDS_SPECIES)
        fossils_data.update({"WHICH": WHICH, "GROUP": GROUP})
        marg_species = numpy.sum(fossils_data["occurrences"], axis=0)
        sort_species = numpy.argsort(marg_species)
        marg_sites = numpy.sum(fossils_data["occurrences"], axis=1)
        sort_sites = numpy.argsort(marg_sites)

        if DO_STEP[0]:
            print("--- bckg")
            #########################################################
            # FIGURES AND TABLES OF DATA CHARACTERISTICS
            ct.save_outline(fossils_data, FIELDS_SPECIES, FIELDS_SITES, BOUNDARIES, CSV_OUTLINE_FILE, tex_localities_file=TEX_LOCALITIES_FILE, plotly_vars_file=PLOTLY_VARS_FILE)
            # OUTLINES OF SITES SPANS IN TIME WITH HYP INFO, NOT SCALED TO TIME
            pt.plot_outline_simple(fossils_data, MAP_FIELDS_SPECIES, MAP_FIELDS_SITES, FIG_OSIMPLE_FILE)
            # VARIOUS COUNTS
            pt.plot_counts_orders(fossils_data, FIELDS_SPECIES, MAP_FIELDS_SITES, FIG_CNTO_FILE)
            if params.get("plot_bbl", True):
                spc_counts, sliced_counts, sites_counts = ct.compute_counts(fossils_data, FIELDS_SPECIES, MAP_FIELDS_SITES)
                pt.plot_bubbles_species(fossils_data, spc_counts, sliced_counts, MAP_FIELDS_SPECIES, MAP_FIELDS_SITES, FIG_BBL_FILE, TEETH_FILE)
                if params.get("svg_bbl", True):
                    pt.plot_bubbles_species(fossils_data, spc_counts, sliced_counts, MAP_FIELDS_SPECIES, MAP_FIELDS_SITES, FIG_BBL_SVG, TEETH_FILE)
            # pt.plot_bubbles_species(fossils_data, spc_counts, sliced_counts, MAP_FIELDS_SPECIES, MAP_FIELDS_SITES, FIG_BBLD_FILE, TEETH_FILE)
            # pt.plot_counts_species(fossils_data, spc_counts, sliced_counts, MAP_FIELDS_SPECIES, MAP_FIELDS_SITES, FIG_CNTS_FILE)
            ct.table_counts(fossils_data, MAP_FIELDS_SPECIES, TAB_CNTO_FILE)

        if DO_STEP[1]:
            #########################################################
            # FIGURES OF CO-OCCURRENCE STATISTICS IN ORIGINAL DATA
            spc_counts, sliced_counts = None, None
            for ti, typ in enumerate(params["measure"]):
                print("--- hist_%s" % typ)
                if typ == "gen":
                    # ORIGINAL, GENERA-BASED
                    gens, gen_counts, sliced_divs = ct.compute_diverse(fossils_data, FIELDS_SPECIES, MAP_FIELDS_SITES)
                    pt.plot_div_hists(WHICH, fossils_data["slices"], gens, gen_counts, sliced_divs, FIG_HIST_FILE % typ)

                elif typ == "pairs":
                    # ORIGINAL, PAIR-BASED
                    if spc_counts is None:
                        spc_counts, sliced_counts, sites_counts = ct.compute_counts(fossils_data, FIELDS_SPECIES, MAP_FIELDS_SITES)
                    scores, tots = ct.compute_pairs_scores(fossils_data["species"], FIELDS_SPECIES, spc_counts, sliced_counts, sites_counts, PAIR_MEASURES)
                    for measure in PAIR_MEASURES:
                        pt.plot_pairs_hists(fossils_data["slices"], spc_counts, scores, tots, FIG_HIST_FILE % measure, id_score=PAIR_MEASURES_MAP[measure], ylbl=VALUES_LBLS[measure])

        if DO_STEP[2] or DO_STEP[3]:
            #########################################################
            # COLLECTING VALUES FROM RANDOMIZED EXPERIMENTS FOR PLOTS
            c_org[(WHICH, GROUP)] = {}
            c_org[(WHICH, GROUP)]["context"], c_org[(WHICH, GROUP)]["context_lbls"] = ct.compute_context_stats(fossils_data, MAP_FIELDS_SPECIES, MAP_FIELDS_SITES)
            c_org[(WHICH, GROUP)]["slices"] = fossils_data["slices"]

#########################################################
# PLOTTING RESULTS FROM RANDOMIZED EXPERIMENTS (normalized across experiments)
for k in params_series_keys:
    print("---", k)
    # Check which stats need to be loaded
    BOXES_VALUES = params_series[k].get("boxes_values", [])
    TRENDS_ROWS = params_series[k].get("trends_rows", [])
    TRENDS_COLS = params_series[k].get("trends_cols", [])
    RND_SERIES = params_series[k].get("rnd_series", [])

    c_rnd_series = set(RND_SERIES)
    f_values = set().union(*TRENDS_ROWS)

    for v in BOXES_VALUES:
        if type(v) is list:
            f_values.update(v[1:])
        else:
            f_values.add(v)

    collected = {}
    collected_org = {}
    for w in params["continent"]:
        for g in params["group"]:

            SSUFF = "%s-%s%s" % (w, g, params["xps_suffix"])

            # to load randomized results
            DATA_STATB_PATT = RND_FLD+"statB-%s_"+SSUFF+"_%s.csv"
            # to save the resulting plots
            FIG_STATB_FILE = FIGS_FLD+"box-%s_"+SSUFF+"_%s.pdf"

            collected_org[(w, g)] = ct.collect_values("original", f_values, VAL_KEYS, PAIR_AGG_SERIES, VALUES, DATA_STATB_PATT)
            collected_org[(w, g)].update(c_org[(w, g)])
            for r in c_rnd_series:
                vals = ct.collect_values(r, f_values, VAL_KEYS, PAIR_AGG_SERIES, VALUES, DATA_STATB_PATT)
                if vals is not None:
                    collected[(w, g, r)] = vals

                    ##########################################
                    # BOX PLOTS, DIFFERENT MEASURES/AGGREGATIONS
                    for v in BOXES_VALUES:
                        vname = v
                        if type(v) is list:
                            series = v[1:]
                            vname = v[0]
                            leg = VALUES_LBLS.get(vname, "")
                        elif VAL_KEYS[vname][-1] == -1:
                            series = [vname+"_"+vg[1] for vg in PAIR_AGG_SERIES]
                            leg = VALUES[VAL_KEYS[vname][0]][VAL_KEYS[vname][1]][2]
                        else:
                            series = [vname]
                            leg = VALUES[VAL_KEYS[vname][0]][VAL_KEYS[vname][1]][2]

                        values = []
                        lbl_series = []
                        for vgname in series:
                            if vgname in collected[(w, g, r)]:
                                lbl_series.append(vgname)
                                if vgname in collected_org[(w, g)]:
                                    values.append(numpy.vstack([[collected_org[(w, g)][vgname]], collected[(w, g, r)][vgname]]))
                                else:
                                    values.append(numpy.vstack([[numpy.zeros(collected[(w, g, r)][vgname].shape[1])], collected[(w, g, r)][vgname]]))

                        if len(values) > 0:
                            values = numpy.transpose(numpy.array(values), (1, 2, 0))
                            if len(lbl_series) == 1:
                                pt.plot_boxes_simple(values, 0, leg, collected_org[(w, g)]["slices"], FIG_STATB_FILE % (vname, r))
                            else:
                                # pdb.set_trace()
                                pt.plot_multi_boxes_series(values, 0, leg, collected_org[(w, g)]["slices"], FIG_STATB_FILE % (vname, r), nb_series=values.shape[2], lbl_series=lbl_series)

    #########################################################
    # PLOTTING RESULTS FROM RANDOMIZED EXPERIMENTS (normalized across experiments)
    # TREND SCATTER PLOTS
    if len(TRENDS_ROWS) > 0 and len(TRENDS_COLS) > 0:
        zvals = {}
        min_max = {}
        kvs = set(TRENDS_COLS).union(*TRENDS_ROWS)
        for kv in kvs:
            zvals[kv] = {}
            for (w, g, r) in collected.keys():
                if r not in zvals[kv]:
                    zvals[kv][r] = {}
                # context variable
                if kv in collected_org[w, g]["context_lbls"]:
                    ii = collected_org[w, g]["context_lbls"].index(kv)
                    zvals[kv][r][(w, g)] = collected_org[w, g]["context"][ii, :]
                elif kv in collected_org[(w, g)]:
                    ms = numpy.mean(collected[(w, g, r)][kv], axis=0)
                    if numpy.sum(ms == 0) > 0:
                        pdb.set_trace()
                    zvals[kv][r][(w, g)] = (collected_org[(w, g)][kv]-ms)/ms
                else:
                    # Value is not available for original data
                    pdb.set_trace()

                # # constant value across
                # elif numpy.sum(numpy.std(collected[(w, g, r)][kv], axis=0)) == 0:
                #     zvals[kv][r][(w, g)] = collected[(w, g, r)][kv][0, :]
                # elif kv in collected_org[(w, g)]:
                #     stds = numpy.std(collected[(w, g, r)][kv], axis=0)
                #     stds[stds < 10**-10] = 1
                #     zvals[kv][r][(w, g)] = (collected_org[(w, g)][kv]-numpy.mean(collected[(w, g, r)][kv], axis=0))/stds
                # else:
                #     stds = numpy.std(collected[(w, g, r)][kv], axis=0)
                #     stds[stds == 0] = 1
                #     zvals[kv][r][(w, g)] = numpy.mean(collected[(w, g, r)][kv], axis=0)/stds

            for r, rzs in zvals[kv].items():
                zvs_mms = numpy.array([[numpy.min(zvs), numpy.max(zvs)] for zvs in rzs.values()])
                zmin, zmax = (numpy.min(zvs_mms[:, 0]), numpy.max(zvs_mms[:, 1]))
                ymin, ymax = (zmin-0.05*(zmax-zmin), zmax+0.1*(zmax-zmin))
                min_max[(kv, r)] = (ymin, ymax)
                # if kv == "FisherQ95_all":
                #     min_max[(kv, r)] = None  # (numpy.maximum(-20, ymin), ymax)
                if kv == "avg_hyp":
                    min_max[(kv, r)] = (0.8, 3.2)
                elif kv == "nb_species":
                    min_max[(kv, r)] = (0, min_max[(kv, r)][1])
                # else:
                #     min_max[(kv, r)] = (numpy.maximum(-6.5, ymin), numpy.minimum(6.5, ymax))

        for (w, g, r) in collected.keys():

            SSUFF = "%s-%s%s" % (w, g, params["xps_suffix"])
            TIKZ_TRENDS_FILE = FIGS_FLD+"trend-%s_"+SSUFF+"_%s.tikz"

            for kj, kvj in enumerate(TRENDS_COLS):
                for ti, trow in enumerate(TRENDS_ROWS):

                    rcomb = ct.combined_preff(trow)
                    ymin = numpy.min([min_max[(kvi, r)][0] for ki, kvi in enumerate(trow)])
                    ymax = numpy.max([min_max[(kvi, r)][1] for ki, kvi in enumerate(trow)])

                    fig_str = "%%%%%% trend-%s_%s-%s%s_%s\n" % (rcomb+"-"+kvj, w, g, params["xps_suffix"], r)
                    fig_str += "\\pgfmathtruncatemacro{\\nbs}{%d}%%\n" % len(trow)  # number of series
                    fig_str += "\\pgfmathtruncatemacro{\\nbt}{%d}%%\n" % len(zvals[kvj][r][(w, g)])  # number of time bins, points per series
                    labels = ""  # "xlabel={%s}, ylabel={%s},%%\n" % (VALUES_LBLS[kvj], "z-score" if g == "L" else "")
                    fig_str += "\\begin{axis}[fbyf,%%\n%sxmin=%f, xmax=%f, ymin=%f, ymax=%f]%%\n" % \
                        (labels, min_max[(kvj, r)][0], min_max[(kvj, r)][1], ymin, ymax)
                    fig_str += "\\addplot[hrz] coordinates {(%f,0) (%f,0)};" % (min_max[(kvj, r)][0], min_max[(kvj, r)][1])

                    for ki, kvi in enumerate(trow):
                        fig_str += "\n%%%%%% series %s\n" % kvi
                        fig_str += "\n".join(["\\node[tsc node] (s%d_n%d) at (axis cs:%f,%f) {};" % (ki+1, i+1, x, y)
                                              for i, (x, y) in enumerate(zip(zvals[kvj][r][(w, g)], zvals[kvi][r][(w, g)]))])
                    fig_str += "\n\\end{axis}"

                    if False:
                        print(fig_str+"\n")
                    else:
                        with open(TIKZ_TRENDS_FILE % (rcomb+"-"+kvj, r), "w") as fo:
                            fo.write(fig_str)
