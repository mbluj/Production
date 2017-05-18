#!/usr/bin/env python

import os, re
import commands
import math

inputDir="/scratch_local/akalinow/CMS/HiggsCP/Data/NTUPLES_17_05_2017/TT/"

fileList = os.listdir(inputDir)

initString = "inputFiles = "

for aFile in fileList:
    if aFile.find("SUSY")!=-1:
        continue
    initString+=aFile+", "


print initString
