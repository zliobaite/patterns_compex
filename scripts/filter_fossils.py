import argparse

import common_tools as ct

# awk -F'\t' '{ if ( $59 != "\\N") {print $35 "," $36 "," $37 "," $59} }' NOW_20200930_public.csv | sort | uniq > ../prepared_0930/teeth_details.csv

# cut -d ',' -f 3-8,10 fossils_NorthAmerica-L.csv | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sitesOrders_NorthAmerica-L.csv
# cut -d ',' -f 3-8,10 fossils_Europe-wide-C.csv | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sitesOrders_Europe-wide-C.csv

# cut -d ',' -f 3-8,10 fossils_NorthAmerica.csv | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sitesOrders_NorthAmerica.csv
# cut -d ',' -f 3-8 fossils_NorthAmerica.csv | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sites_NorthAmerica.csv

# awk '{ print $0 "," ARGIND }' fossils_NorthAmerica-L.csv fossils_NorthAmerica-S.csv fossils_NorthAmerica-C.csv | cut -d ',' -f 3-8,17 | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sites_NorthAmerica-LSC.csv
# sed 's/ *\"\([A-Z][a-z][a-z]\)\(.*\)\": \"#\([A-F0-9]*\)\".*/\\textcolor{c\1}{\1\2} \\\\/' ../colors_LSC.txt

# awk '{ print $0 "," ARGIND }' fossils_Europe-wide-L.csv fossils_Europe-wide-S.csv fossils_Europe-wide-C.csv | cut -d ',' -f 3-8,17 | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sites_Europe-wide-LSC.csv

# awk '{ print $0 "," ARGIND }' fossils_Europe-drop-L.csv fossils_Europe-drop-S.csv fossils_Europe-drop-C.csv | cut -d ',' -f 3-8,17 | sort -n | uniq -c | sed 's/^ *\([0-9]*\) \(.*\)$/\2,\1/' | awk '{print FNR-1","$0}' > sites_Europe-drop-LSC.csv


parser = argparse.ArgumentParser(description='Prepare data subsets for co-occurrence experiments from a NOW database csv dump.')
# parser.add_argument("-i", "--db_dump", type=argparse.FileType('r'), help="NOW database csv dump", default="../db_dumps/NOW_latest_public.csv")
parser.add_argument("-i", "--db_dump", type=str, help="NOW database csv dump", default="../db_dumps/NOW_latest_public.csv")
parser.add_argument("-d", "--data_folder", type=str, help="folder to store the data subsets", default="../prepared_data/")
parser.add_argument("-t", "--time_folder", type=str, help="folder containing the time bins specifications", default="../times/")
parser.add_argument("--in_sep", type=str, help="input field separator", default=",")
parser.add_argument("--out_sep", type=str, help="output field separator", default=",")
parser.add_argument("--missing", type=str, help="missing value token", default="\\N")
args = parser.parse_args()

fmt_args = {"in_sep": args.in_sep,
            "out_sep": args.out_sep,
            "missing": args.missing}

FOSSILS_FILE = args.db_dump
TIMES_FLD = args.time_folder

DATA_FLD = args.data_folder
ct.make_fld(DATA_FLD)


OUT_PATT = DATA_FLD+"fossils_%s-%s.csv"
KEEP_FIELDS = ["LIDNUM", "NAME", "LAT", "LONG",
               "MAX_AGE", "MIN_AGE", "BFA_MAX", "BFA_MIN", "FRAC_MAX", "FRAC_MIN", "COUNTRY", "STATE",
               "MEAN_HYPSODONTY", "TCRWNHT", "HORIZODONTY",
               "SIDNUM", "ORDER", "FAMILY", "GENUS", "SPECIES",
               "SUBCLASSORSUPERORDER", "SUBORDERORSUPERFAMILY", "SUBFAMILY",
               "UNIQUE", "TAXON_STATUS", "ID_STATUS", "ADD_INFO", "SOURCE_NAME", "SPCOMMENT", "SYNONYMS"]
SLICES_FIELDS = ["SLICE_ID", "SLICE_NAME"]

ORD_SMALL_HERB = set(["Cimolesta", "Dermoptera", "Didelphimorphia", "Eulipotyphla", "Lagomorpha", "Leptictida", "Macroscelidea", "Multituberculata", "Proteutheria", "Rodentia"])
ORD_LARGE_HERB = set(["Arctostylopida", "Artiodactyla", "Dinocerata", "Embrithopoda", "Hyracoidea", "Pantodonta", "Perissodactyla", "Primates", "Proboscidea", "Taeniodonta", "Tillodontia"])
ORD_LARGE_CARN = set(["Carnivora", "Condylarthra", "Creodonta", "Mesonychia", "Ptolemaiida", "Tubulidentata"])
ORD_OTHER = set(["Cetacea", "Chiroptera", "Cingulata", "indet.", "Marsupialia", "Palaeanodonta", "Pholidota", "Pilosa", "Pinnipedia"])
ORD_UNKNOWN = set(["incertae sedis"])

EU_COUNTRIES = set(["Armenia", "Austria", "Azerbaijan", "Belgium", "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Czech Republic", "France", "Georgia", "Germany", "Greece", "Hungary", "Iran", "Iraq", "Italy", "Kazakhstan", "Moldova", "North Macedonia", "Poland", "Portugal", "Romania", "Russia", "Serbia", "Serbia and Montenegro", "Slovakia", "Spain", "Switzerland", "Syria", "Turkey", "Ukraine"])


GROUPS = {"S": ORD_SMALL_HERB, "L": ORD_LARGE_HERB, "C": ORD_LARGE_CARN,
          "O": ORD_OTHER, "X": ORD_UNKNOWN}
MAP_GROUPS = {}
for k, vs in GROUPS.items():
    for v in vs:
        MAP_GROUPS[v] = k

subsets_defs = [{"name": "Europe-wide",
                 "slices_file": TIMES_FLD+"time-slices_europe_filter.csv",
                 "slice_max_ratio_out": 0.1,
                 "tests": [("LAT", "MIN", 14),
                           ("LAT", "MAX", 82),
                           ("LONG", "MIN", -24),
                           ("LONG", "MAX", 75),
                           ("COUNTRY", "IN", EU_COUNTRIES),
                           ("MIN_AGE", "MIN", 2.5),
                           ("MAX_AGE", "MAX", 23)]
                 },
                {"name": "NorthAmerica",
                 "slices_file": TIMES_FLD+"time-slices_america_filter.csv",
                 "slice_max_ratio_out": 0.1,
                 "tests": [("LAT", "MIN", 19),
                           ("LAT", "MAX", 84),
                           ("LONG", "MIN", -140),
                           ("LONG", "MAX", -70),
                           ("COUNTRY", "IN", None),
                           ("MIN_AGE", "MIN", 2.5),
                           ("MAX_AGE", "MAX", 34)]
                 }
                ]


def filter_fossils(traits_file, keep_fields, subsets_defs, slices_fields=[], in_sep=",", out_sep=",", missing="\\N"):
    head_traits = None
    with open(traits_file) as fp:
        for line in fp:
            parts = line.strip().split(in_sep)
            if head_traits is None:
                head_traits = dict([(v, k) for (k, v) in enumerate(parts)])
                row = out_sep.join([k for k in keep_fields if k in head_traits]+slices_fields) + "\n"
                for subset_def in subsets_defs:
                    for kg, vg in GROUPS.items():
                        subset_def["fo_%s" % kg].write(row)
            else:
                row_pieces = [parts[head_traits[k]] for k in keep_fields if k in head_traits]
                for subset_def in subsets_defs:
                    include = parts[head_traits["ORDER"]] in MAP_GROUPS
                    for test in subset_def["tests"]:
                        if test[2] is not None:
                            if test[1] == "MAX":
                                include &= float(parts[head_traits[test[0]]]) <= test[2]
                            if test[1] == "MIN":
                                include &= float(parts[head_traits[test[0]]]) >= test[2]
                            if test[1] == "IN":
                                include &= parts[head_traits[test[0]]] in test[2]
                            if test[1] == "NOT IN":
                                include &= parts[head_traits[test[0]]] not in test[2]

                    if include:
                        if subset_def.get("slices") is not None:
                            si, slce = ct.get_slice(subset_def["slices"], parts, head_traits, max_ratio_out=subset_def.get("slice_max_ratio_out", 0))
                            if si != -1:
                                subset_def["fo_%s" % MAP_GROUPS[parts[head_traits["ORDER"]]]].write(out_sep.join(row_pieces+ct.get_slice_info(si, slce, slices_fields, pref="SLICE_", missing=missing, to_str=True)) + "\n")
                        else:
                            subset_def["fo_%s" % MAP_GROUPS[parts[head_traits["ORDER"]]]].write(out_sep.join(row_pieces+[missing for f in slices_fields]) + "\n")


for s in subsets_defs:
    if s.get("slices_file") is not None:
        s["slices"] = ct.read_slices(s["slices_file"])
    for kg, vg in GROUPS.items():
        s["fo_%s" % kg] = open(OUT_PATT % (s["name"], kg), "w")

filter_fossils(FOSSILS_FILE, KEEP_FIELDS, subsets_defs, SLICES_FIELDS, **fmt_args)

for s in subsets_defs:
    for kg, vg in GROUPS.items():
        s["fo_%s" % kg].close()
