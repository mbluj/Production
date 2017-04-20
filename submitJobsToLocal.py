#!/usr/bin/env python

import os, re
import commands
import math, time

from analyzerMC import *
########################################################
########################################################
def submitJob(aFile, dataPath, back):
    #Update the CMSSW configuration
    process.source.fileNames =  cms.untracked.vstring()
    process.source.fileNames.append('file:'+ dataPath + aFile)
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))
    ##Prepare working directory
    workdir = dataPath.split("Data/")[1] + aFile.split(".root")[0]
    command = "mkdir -p "+workdir

    os.system(command)
    command = "cp job.sh "+workdir
    os.system(command)
    out = open(workdir+'/'+'tmpConfig.py','w')
    out.write(process.dumpPython())
    out.close()
    command = "cd "+workdir+"; ./job.sh >& a.out"
    if back==True:
        command+=" &"
    os.system(command)
########################################################
########################################################
dataPaths = ["/scratch/cms/akalinow/CMS/HiggsCP/Data/SUSYGluGluToBBHToTauTau_M-1000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM/"]

########################################################
for path in dataPaths:
    command = "ls "+path
    fileList = commands.getoutput(command).split("\n")
    index = 0
    for aFile in fileList:
        command = "ps aux | grep cmsRun"
        while commands.getoutput(command).count("cmsRun tmpConfig.py")>=1:
            time.sleep(150)
        submitJob(aFile,dataPath=path, back=True)
        index+=1
########################################################
