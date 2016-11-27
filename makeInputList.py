#!/usr/bin/env python                                                                                                                                                                       

import os, re
import commands
import math

inputDir="/scratch_local/akalinow/CMS/HiggsCP/Data/NTUPLES_26_11_2016/MT/"

fileList = os.listdir(inputDir)
        
initString = "inputFile = "

for aFile in fileList:
    if aFile.find("SUSY")!=-1: 
        continue
    initString+=inputDir+"/"+aFile+", "


print initString
