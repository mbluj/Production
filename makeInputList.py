#!/usr/bin/env python

import os, re
import commands
import math

inputDir="/scratch_local/akalinow/CMS/HiggsCP/Data/WAWNTuples/2017/NTUPLES_26_10_2017/MT"

fileList = os.listdir(inputDir)

initString = "inputFiles = "

for aFile in fileList:
    if aFile.find("SUSY")!=-1:
        continue
    initString+=aFile+", "


print initString
