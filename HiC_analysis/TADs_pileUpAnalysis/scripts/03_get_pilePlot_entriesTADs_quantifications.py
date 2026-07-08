import coolpuppy.lib.io
import glob

import os
import sys

mypath = "./"
inFiles = glob.glob("hic*rescaled*.clpy")
print(inFiles)

for inFile in inFiles:

    print(inFile)
    inFileName = os.path.basename(inFile)
    inFileName = os.path.splitext(inFileName)[0]
    print(inFileName)

    pileup_df = coolpuppy.lib.io.load_pileup_df(inFile)

    # I use the 'resolution' field to get the resolution of the map
    resolution = pileup_df['resolution'][0]
    print("Coolpuppy pile-up matrix at",resolution,"bp resolution")

    # I use the 'data' field, which contains a NxN, where N is the size of the resulting pileup squared matrix
    # In this case the matrix is 21x21, where 21 comes from a region of (100kb (left-flanking) +
    # 10kb (region of interest) + 100kb (right-flanking)) at 10kb resolution
    # n is the counter of the Coolpuppy pile-up matrices contained in inFile
    n = 1
    for matrix in pileup_df['data']:
        print(type(matrix))
        print(matrix)
        
        size = len(matrix)
        print("Size of the Coolpuppy pile-up matrix",size)
        print("Matrix indexes go from 0 to",size-1)

        windowTAD = size / 3.0
        print("Window of the central part is 1/3 of the size: the entire rescaled TAD:", windowTAD)

        print("Central part")
        outFile   = inFileName + "_central_%d.txt" % windowTAD 
        fpoutFile = open(outFile,"w") 
        iStartTAD = int((size-1)/2)-int(windowTAD/2)
        iEndTAD   = int(iStartTAD + windowTAD)-1
        print("The",windowTAD,"x",windowTAD,"sub-matrix of interest goes from index",iStartTAD,"to index",iEndTAD,"included")

        for i in range(iStartTAD,iEndTAD+1):
            for j in range(iStartTAD,iEndTAD+1):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))
                print(i,j,matrix[i][j])
        fpoutFile.close()
                
        print("Top-left part on diagonal")
        window = 66
        outFile = inFileName + "_topLeftTAD_%d.txt" % window
        fpoutFile = open(outFile,"w")         
        yiStart = int(iStartTAD-window/2)                
        yiEnd   = int(iStartTAD+window/2)        
        xiStart = int(iStartTAD-window/2)                
        xiEnd   = int(iStartTAD+window/2)        
        print("The",window,"x",window,"sub-matrix of interest goes from index",xiStart,"to index",xiEnd,"included in x")
        print("The",window,"x",window,"sub-matrix of interest goes from index",yiStart,"to index",yiEnd,"included in y")        

        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                print(i,j,matrix[i][j])

        
        for i in range(xiStart,xiEnd):
            for j in range(yiStart,yiEnd):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))
                print(i,j,matrix[i][j])                
        fpoutFile.close()
              
        print("Bottom-right part on diagonal")
        window = 66
        outFile = inFileName + "_bottomRightTAD_%d.txt" % window
        fpoutFile = open(outFile,"w")                 
        yiStart = int(iEndTAD-window/2)                
        yiEnd   = int(iEndTAD+window/2)        
        xiStart = int(iEndTAD-window/2)                
        xiEnd   = int(iEndTAD+window/2)        
        print("The",window,"x",window,"sub-matrix of interest goes from index",xiStart,"to index",xiEnd,"included in x")
        print("The",window,"x",window,"sub-matrix of interest goes from index",yiStart,"to index",yiEnd,"included in y")        

        for i in range(xiStart,xiEnd):
            for j in range(yiStart,yiEnd):
                fpoutFile.write("%s %d %d %f\n" % (inFile,i,j,matrix[i][j]))
                print(i,j,matrix[i][j])
        fpoutFile.close()
