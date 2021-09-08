import sys, re, os
import numpy
import matplotlib.pyplot as plt
import matplotlib.cm
from matplotlib.backends.backend_pdf import PdfPages


import pdb

header = ["smaller", "obs", "larger", "cij", "ci", "cj", "N", "cmean", "cstd", "zval", "pair_status", "genj", "spcj", "geni", "spci", "continent", "group", "slice"]
map_header = dict([(v,k) for (k,v) in enumerate(header)])

suff = "EU-L"
if len(sys.argv) > 2:
    suff = "%s-%s" % (sys.argv[1], sys.argv[2])

slc_lbls = ["MN2", "MN3", "MN4", "MN5", "MN6", "MN7-8", "MN9", "MN10", "MN11", "MN12", "MN13", "MN14", "MN15", "MN16", "all"]
if suff[:2] == "NA":
    slc_lbls = ["Orellan", "Whitneyan", "Arikareean-1", "Arikareean-2", "Arikareean-3", "Arikareean-4", "Hemingfordian-1", "Hemingfordian-2", "Barstovian-1", "Barstovian-2", "Clarendonian-1", "Clarendonian-2", "Clarendonian-3", "Hemphillian-1", "Hemphillian-2", "Hemphillian-3", "Hemphillian-4", "Blancan-Early", "all"]



    
FILE_PAIRS = "../oprob/all_pairs-10000.csv"
# FILE_PAIRS = "../oprob/all_pairs.csv"
FILE_PAIRS = "../oprob/all_pairs_det_"+suff+".csv"

FIG_OUT = "../oprob/dist_pairsO_det_"+suff+".pdf"

CONTINENTS = {b"EU": 0, b"NA": 1}
GROUPS = {b"L": 0, b"S": 1, b"C": 2}
def convert_continent(s):
    return CONTINENTS.get(s, -1)
def convert_group(s):
    return GROUPS.get(s, -1)
def convert_dum(s):
    return 0

converters = {map_header["continent"]: convert_continent, map_header["group"]: convert_group}
for v in ["genj", "spcj", "geni", "spci"]:
    converters[map_header[v]] = convert_dum

bin_step = 2**-4
eps = 10**-8
bins = numpy.arange(-bin_step,1+2*bin_step,bin_step)    
bins[0] = 0
bins[1] = eps
bbs = numpy.concatenate([[-bin_step, 0], bins[2:-1]])

X_all = numpy.loadtxt(FILE_PAIRS, delimiter=",", converters=converters)
# ctop = numpy.max(X_all[:, [map_header["cj"], map_header["ci"]]])
xvals = X_all[:, map_header["cij"]]+0.33*numpy.random.random(X_all.shape[0])
xtop = numpy.max(xvals)

yvals = X_all[:, map_header["cmean"]]
ytop = numpy.max(yvals)

cvals = numpy.abs(X_all[:, map_header["cmean"]] - X_all[:, map_header["cij"]])/X_all[:, map_header["cstd"]]
ctop = numpy.max(cvals)

pmids = X_all[:, map_header["larger"]]+0.5*X_all[:, map_header["obs"]]

mask_same = X_all[:, map_header["pair_status"]] == 1

with PdfPages(FIG_OUT) as pdf:
# if True:
    for slc in range(-1, int(numpy.max(X_all[:, map_header["slice"]]))):
        if slc == -1:
            mask_slice = X_all[:, map_header["slice"]] > -1
        else:
            mask_slice = X_all[:, map_header["slice"]]==slc
            
        main_fig = plt.figure(figsize=(9, 6))
        splts = main_fig.subplots(1, 2)

        h, _ = numpy.histogram(pmids[mask_slice], bins)
        splts[0].barh(bbs, h, 0.95*bin_step, 0, align="edge", color="#888888")
        splts[0].set_xscale('log')
        splts[0].set_ylim(-1.5*bin_step,1+1.5*bin_step)
        splts[0].set_ylabel("P")
        splts[0].set_xlabel("nb pairs")
            
        splts[1].scatter(xvals[mask_slice & ~mask_same], yvals[mask_slice & ~mask_same], s=5, c=cvals[mask_slice & ~mask_same], vmin=-1, vmax=ctop+1, cmap="gray") #"#CCCCCC")
        splts[1].scatter(xvals[mask_slice & mask_same], yvals[mask_slice & mask_same], s=3, c="red") #c="#117733")

        splts[1].plot([0, xtop+1], [0, xtop+1], ":", c="#AAAAAA")
        splts[1].set_ylim(0,ytop+1)
        splts[1].set_xlim(0,xtop+1)
        splts[1].set_xlabel("Observed overlap")
        splts[1].set_ylabel("Expected overlap")
        
        # splts[si, sj+1].set_xlabel(y)
        # splts[si, sj+1].set_ylim(-1.5*bin_step,1+1.5*bin_step)
        # if ytop is not None:
        #     splts[si, sj+1].set_xlim(0,ytop+1)
        plt.subplots_adjust(left=0.08, bottom=0.08, right=0.95, top=0.95)        
        plt.suptitle("%s %s" % (suff, slc_lbls[slc]))        
        pdf.savefig()
        # # plt.savefig(FIG_OUT % slc)
        # plt.show()
        # plt.close()

