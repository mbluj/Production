#!/usr/bin/env python                                                                                                                                                                       

import os, re
import commands
import math

inputDir="/home/akalinow/scratch/CMS/HiggsCP/Data/NTUPLES_07_10_2016/"

fileList = os.listdir(inputDir)
        
initString = "inputFile = "

for aFile in fileList:
    if aFile.find("SUSY")!=-1: 
        continue
    initString+=inputDir+"/"+aFile+", "


print initString
