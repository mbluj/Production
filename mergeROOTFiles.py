#!/usr/bin/env python

import os, re
import commands

#########################################
#########################################
def mergeDataset(dataset, publish_data_suffix, outputDir):

    workdir = publish_data_suffix
    shortName = dataset.split("/")[1]
    pathPart1 = shortName
    if dataset.split("/")[2].find("Run201")!=-1:
        shortName += "_"+dataset.split("/")[2]
        pathPart1 = "SingleMuon"
    
    dataDirectory =  "Data/"+publish_data_suffix+"/"+pathPart1

    shortName = shortName.replace("-","_")
    shortName = shortName.split("_")[0]+shortName.split("_")[1]+shortName.split("_")[2]
    shortName+="_"+publish_data_suffix
    if dataset.find("ext")!=-1:
        shortName+= "_"+dataset[dataset.find("ext"):dataset.find("ext")+4]

    if dataset.find("part")!=-1:
        shortName+= "_"+dataset[dataset.find("part"):dataset.find("part")+6]

    outputFileName = outputDir+"/"+shortName+".root"

    dataDirectory += "/"+shortName
    dataDirectory+="/"+os.listdir(dataDirectory)[0]

    for directory in os.listdir(dataDirectory):
        tmpDir = dataDirectory+"/"+os.listdir(dataDirectory)[0]
        command = "hadd -f "+outputFileName+" "+tmpDir+"/*.root"
        os.system(command)
#########################################
#########################################


