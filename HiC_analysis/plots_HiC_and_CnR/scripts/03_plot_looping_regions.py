import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os
import glob
import cooler
import cooltools.lib.plotting

from matplotlib.ticker import EngFormatter
from matplotlib.colors import LogNorm
import matplotlib.patches as patches

import os.path
import sys

### plot the raw and corrected data in logscale ###
from mpl_toolkits.axes_grid1 import make_axes_locatable

bp_formatter = EngFormatter('b')

# create a functions that would return a series of rectangles around called dots
# in a specific region, and exposing importnat plotting parameters
def rectangles_around_dots(dots_df, region, loc="upper", lw=1, ec="cyan", fc="none"):
    """
    yield a series of rectangles around called dots in a given region
    """
    # select dots from the region:
    df_reg = bioframe.select(
        bioframe.select(dots_df, region, cols=("chrom1","start1","end1")),
        region,
        cols=("chrom2","start2","end2"),
    )
    rectangle_kwargs = dict(lw=lw, ec=ec, fc=fc)
    # draw rectangular "boxes" around pixels called as dots in the "region":
    for s1, s2, e1, e2 in df_reg[["start1", "start2", "end1", "end2"]].itertuples(index=False):
        width1 = e1 - s1
        width2 = e2 - s2
        if loc == "upper":
            yield patches.Rectangle((s2, s1), width2, width1, **rectangle_kwargs)
        elif loc == "lower":
            yield patches.Rectangle((s1, s2), width1, width2, **rectangle_kwargs)
        else:
            raise ValueError("loc has to be uppper or lower")

        
def format_ticks(ax, x=True, y=True, rotate=True):
    if y:
        ax.yaxis.set_major_formatter(bp_formatter)
    if x:
        ax.xaxis.set_major_formatter(bp_formatter)
        ax.xaxis.tick_bottom()
    if rotate:
        ax.tick_params(axis='x',rotation=45)

# download test data
# this file is 145 Mb, and may take a few seconds to download
import cooltools

mainDir = "./"

exten = 4

targetResolution = 5000
conditions       = ["WD_NOnub","WD_SUMOnub"]

outDir = str(sys.argv[1])
cond   = str(sys.argv[2])
conditions = [cond]
print(outDir)

inFile = outDir + ".bedpe"

fpInput = open(inFile, "r")

coolFiles = []

chromSizes = {"chr3L" : 28110227,
              "chr3R"	: 32079331,
              "chr2L"	: 23513712,
              "chr2R"	: 25286936,
              "chrX"	: 23542271,
              "chrY"	: 3667352,
              "chrM"	: 19524,
              "chr4"	: 1348131}


norm = LogNorm(vmax=0.1)
norm_raw = LogNorm(vmin=1, vmax=10_000)

for condition in conditions:
    coolFile = glob.glob("../01_cool_files/hic*contacts*%s_*merge*NA*5000*mcool" % (condition))[0]
    coolFiles.append(coolFile)
    
print(coolFiles)
#exit(1)
resolutions = [400,800,1000,2000,4000,8000,10000,20000,40000]
n = 0

fpInput.seek(0)
for loop in fpInput.readlines():
    loop = loop.strip().split()
    if len(loop) == 0:
        continue
    n += 1
    print(loop[0],loop[1],loop[2],loop[3],loop[4],loop[5])
    for resolution in [5000]:

        #res1 = int(float(loop[2]))-int(float(loop[1]))
        #res2 = int(float(loop[5]))-int(float(loop[4]))
        #if res1 in resolutions and res2 == res1:
        #    resolution = res1
        #else:
        #    resolution = 4000
        print("Resolution = ",resolution)

        chrom1 = loop[0]
        start1 = int(int(float(loop[1])/resolution)*resolution)
        end1   = int(int(float(loop[2])/resolution)*resolution)
        chrom2 = loop[3]
        start2 = int(int(float(loop[4])/resolution)*resolution)
        end2   = int(int(float(loop[5])/resolution)*resolution)
        #name   = loop[6]
        name   = str(n) + "_res" + str(resolution)
        nameOld = "res" + str(resolution)
        
        
        if (chromSizes[chrom1] < (end1+exten*resolution)):
            end1   = int(int(float(chromSizes[chrom1])/resolution)*resolution) - exten*resolution
            start1 = end1 - (exten+1)*resolution
        if (chromSizes[chrom2] < (end2+exten*resolution)):
            end2   = int(int(float(chromSizes[chrom2])/resolution)*resolution) - exten*resolution
            start2 = end2 - (exten+1)*resolution 
        
        print(chrom1,start1,end1,chrom2,start2,end2,chromSizes[chrom1])
        
        #if abs(start1-start2) > 150000  and resolution <= 2000:
        #    continue
        #if abs(start1-start2) > 1000000 and resolution <= 5000:
        #    continue
        #if abs(start1-start2) > 5000000 and resolution <= 10000:
        #    continue
        #if abs(start1-start2) > 10000000 and resolution <= 20000:
        #    continue                        

        outFileOld = outDir + "/" + "Loop_%s_%s_%s_%s_%s_%s_loop_%s_flanking_%s_at_%sbp.pdf" % (chrom1, start1, end1, chrom2, start2, end2, nameOld, exten, resolution)        
        outFile = outDir + "/" + "Loop_%s_%s_%s_%s_%s_%s_loop_%s_%s_flanking_%s_at_%sbp.pdf" % (chrom1, start1, end1, chrom2, start2, end2, name,condition, exten, resolution)
        if (os.path.isfile(outFile)):
            continue

        if (os.path.isfile(outFileOld)):
            os.rename(outFileOld, outFile)
            continue

        if start1-exten*resolution < 0:
            print("Redefining start1")
            start1 = (exten+1)*resolution
            extension = (start1-exten*resolution, end1+exten*resolution, end2+exten*resolution, start2-exten*resolution)

        print(loop)
        region = (chrom1, start1-exten*resolution, end2+exten*resolution)
        print(region)
        extension = (start1-exten*resolution, end1+exten*resolution, end2+exten*resolution, start2-exten*resolution)
        #print("Extension of the plot %s" % region)
        

        plt_width=4
        f, axs = plt.subplots(
            figsize=( (plt_width+1)*len(coolFiles), plt_width),
            ncols=int(2*len(coolFiles)),
            nrows=1,
            gridspec_kw={'height_ratios':[4],"wspace":0.01,'width_ratios':[1,.05]*len(coolFiles)},
            constrained_layout=True
        )

        #ax.set_title(f'{chrom1}:{start1:,}-{end1:,} - {chrom2}:{start2:,}-{end2:,}')
        nPlot = 0
        for icoolFile in range(len(coolFiles)):

            coolFile  = coolFiles[icoolFile]
            condition = conditions[icoolFile]
            print(coolFile,condition)


            clr = cooler.Cooler(mainDir + coolFile + '::resolutions/%d' % resolution)            
            #print(f'chromosomes: {clr.chromnames}, binsize: {clr.binsize}')

            sizeX, sizeY=clr.matrix(balance=True).fetch(region).shape
            maxX=int(2*exten)+1
            #print(sizeX,sizeY)
            #print(clr.matrix().fetch(region))

            #print(clr.matrix().fetch(region).shape)            
            #print(clr.matrix().fetch(region)[0:maxX,0:maxX])

            #print(clr.matrix().fetch(region).shape)            
            #print(clr.matrix().fetch(region)[sizeY-maxX:sizeY,sizeY-maxX:sizeY])            
            #print(clr.matrix().fetch(region)[0:maxX,sizeY-maxX:sizeY])
            #exit(1)
            
            ax = axs[nPlot]
            im = ax.matshow(
                clr.matrix().fetch(region),
                norm=norm,
                cmap='fall',
                extent=extension
            );
            #ax.xaxis.set_visible(False)
            
            #ax.add_patch(
            #    patches.Rectangle(
            #        (int(int(start1/resolution)*resolution-0.5*float(resolution)), int(int(start2/resolution)*resolution-0.5*float(resolution))),
            #        2*resolution,
            #        2*resolution,
            #        edgecolor='black',
            #        fill=False,
            #        linestyle='dashed',
            #        lw=2
            #    ) )
            
            nPlot += 1
            
            #ax.xaxis.set_visible(False)
            ax.set_title(condition)
            
            cax = axs[nPlot]
            plt.colorbar(im, cax=cax, label='ICE-balanced counts', fraction=0.046)
            nPlot += 1
            
        plt.savefig(outFile)
        plt.close()

