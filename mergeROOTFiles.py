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
        pathPart1 = dataset.split("/")[1]
    
    dataDirectory =  "Data/"+publish_data_suffix+"/"+pathPart1

    shortName = dataset.split("/")[1]
    if dataset.split("/")[2].find("Run201")!=-1:
        shortName += "_"+dataset.split("/")[2]

    shortName = shortName.replace("-","_")
    shortName = shortName.split("_")[0]+shortName.split("_")[1]+shortName.split("_")[2]

    if dataset.find("PromptReco-v")!=-1:
        shortName+= "_v"+dataset[dataset.find("PromptReco-v")+12:dataset.find("PromptReco-v")+13]
        
    if dataset.find("ext")!=-1:
        shortName+= "_"+dataset[dataset.find("ext"):dataset.find("ext")+4]

    if dataset.find("part")!=-1:
        shortName+= "_"+dataset[dataset.find("part"):dataset.find("part")+6]

    if dataset.find("t-channel")!=-1:
        shortName+= "_"+dataset[dataset.find("channel")+7:dataset.find("channel")+15]
        
    shortName = shortName.rstrip("-")
    shortName+="_"+publish_data_suffix

    outputFileNameMT = outputDir+"/"+shortName+".root"
    outputFileNameTT = outputDir+"/TT_"+shortName+".root"
    outputFileNameMM = outputDir+"/MM_"+shortName+".root"
    

    dataDirectory += "/"+shortName
    dataDirectory+="/"+os.listdir(dataDirectory)[0]

    command = "mkdir -p "+outputDir
    os.system(command)

    command = "hadd -f "+outputFileNameMT+" "+dataDirectory+"/*/WAWMT*.root"
    os.system(command)
    command = "hadd -f "+outputFileNameTT+" "+dataDirectory+"/*/WAWTT*.root"
    os.system(command)
    command = "hadd -f "+outputFileNameMM+" "+dataDirectory+"/*/WAWMM*.root"
    os.system(command)
#########################################
#########################################


