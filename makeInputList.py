#!/usr/bin/env python                                                                                                                                                                       

import os, re
import commands
import math

inputDir="/home/akalinow/scratch/CMS/HiggsCP/Data/NTUPLES_07_09_2016/"

fileList = os.listdir(inputDir)
        
initString = "inputFile = "

for aFile in fileList:
    if aFile.find("SUSY")!=-1 or aFile.find("VBF")!=-1: 
        continue
    initString+=inputDir+"/"+aFile+", "


print initString
