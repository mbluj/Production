#!/usr/bin/env python

import os, re
import commands
import math

inputDir="/home/akalinow/scratch/CMS/HiggsCP/Data/WAWNTuples/Summer17_SVFit_v11/MT/"

fileList = os.listdir(inputDir)

initString = "inputFiles = "

for aFile in fileList:
    if aFile.find("SUSY")!=-1:
        continue
    initString+=aFile+", "


print initString
