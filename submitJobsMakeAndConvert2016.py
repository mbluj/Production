#!/usr/bin/env python

import os, re
import commands
import math

from crab3 import *
#########################################
#########################################
def prepareCrabCfg(dataset,
                   crabCfgName,
                   eventsPerJob,
                   jsonFile,
                   storage_element, 
                   publish_data_suffix):

    workdir = publish_data_suffix
    shortName = dataset.split("/")[1]
    if dataset.split("/")[2].find("Run201")!=-1:
        shortName += "_"+dataset.split("/")[2]

    shortName = shortName.replace("-","_")
    shortName = shortName.split("_")[0]+shortName.split("_")[1]+shortName.split("_")[2]
    shortName+="_"+publish_data_suffix
    if dataset.find("ext")!=-1:
        shortName+= "_"+dataset[dataset.find("ext"):dataset.find("ext")+4]

    if dataset.find("part")!=-1:
        shortName+= "_"+dataset[dataset.find("part"):dataset.find("part")+6]
        
    shortName = shortName.rstrip("-")

    ##Modify CRAB3 configuration
    config.JobType.psetName = 'analyzerMC.py'

    config.JobType.disableAutomaticOutputCollection = True
    config.JobType.scriptExe = 'makeAndConvert.py'
    config.JobType.outputFiles = ['WAW_HTauTauAnalysis.root']
    config.JobType.inputFiles = ['HTauTauTree_v2.C', 'HTauTauTree_v2.h', 'HTTEvent.cxx', 'HTTEvent.h', 'PropertyEnum.h', 'TriggerEnum.h']
    
    config.Site.storageSite = storage_element
    config.General.requestName = shortName

    config.Data.inputDataset = dataset
    config.Data.outLFNDirBase = '/store/user/akalinow/'+publish_data_suffix+"/"
    config.Data.outputDatasetTag = shortName
    config.Data.inputDBS = 'global'
    config.Data.splitting = 'EventAwareLumiBased'
    config.Data.unitsPerJob = eventsPerJob
    config.Data.totalUnits =  -1
    config.Data.lumiMask=""
    if dataset.split("/")[2].find("Run201")!=-1:
        config.Data.lumiMask=jsonFile
        config.JobType.psetName = 'analyzerData.py'
    out = open('crabTmp.py','w')
    out.write(config.pythonise_())
    out.close()        
    os.system("crab submit -c crabTmp.py")
#########################################
#########################################
eventsPerJob = 200000

datasets = [
    #Data
    "/SingleMuon/Run2016B-PromptReco-v2/MINIAOD",
    "/SingleMuon/Run2016C-PromptReco-v2/MINIAOD",
    "/SingleMuon/Run2016D-PromptReco-v2/MINIAOD",
    "/SingleMuon/Run2016E-PromptReco-v2/MINIAOD",
    #Signal SM
    "/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM",
    "/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM",
    #Signal MSSM
    "/SUSYGluGluToHToTauTau_M-90_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM",
    "/SUSYGluGluToHToTauTau_M-120_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM",
    "/SUSYGluGluToHToTauTau_M-130_TuneCUETP8M1_13TeV-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14-v1/MINIAODSIM",
    #DY
    "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/MINIAODSIM", 
    #Wjets
    "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext1-v1/MINIAODSIM", 
    #TT
    "/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16MiniAODv2-PUSpring16RAWAODSIM_reHLT_80X_mcRun2_asymptotic_v14_ext3-v1/MINIAODSIM", 
    ##   
]
##TEST
#datasets = [ "/SingleMuon/Run2016C-PromptReco-v2/MINIAOD",
#             "/SingleMuon/Run2016D-PromptReco-v2/MINIAOD"]
###############

jsonFile2016 = "/home/akalinow/scratch/CMS/HiggsCP/Prod/Crab/JSON/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt"
########################################################
for dataset in datasets:
    prepareCrabCfg(crabCfgName="crab3.py",
                   dataset=dataset,
                   eventsPerJob=eventsPerJob,
                   jsonFile=jsonFile2016,
                   storage_element="T2_PL_Swierk",
                   publish_data_suffix = "v11")
########################################################


