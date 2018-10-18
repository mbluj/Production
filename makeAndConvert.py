#!/usr/bin/env python

import os

from ROOT import gSystem, TChain, TSystem, TFile

doSvFit = True
if doSvFit :
    print "Run with SVFit computation"

#Checkout DNN training files
command = "cd $CMSSW_BASE/src; "
command += "git clone https://github.com/cms-tau-pog/RecoTauTag-TrainingFiles -b master RecoTauTag/TrainingFiles/data; "
command += "cd -;"
os.system(command)

#Some system have problem runnig compilation (missing glibc-static library?).
#First we try to compile, and only then we start time consuming cmssw
status = gSystem.CompileMacro('HTTEvent.cxx')
gSystem.Load('$CMSSW_BASE/lib/$SCRAM_ARCH/libTauAnalysisClassicSVfit.so')
status *= gSystem.CompileMacro('HTauTauTreeBase.C')
status *= gSystem.CompileMacro('HTauMuTauHTree.C')
status *= gSystem.CompileMacro('HTauhTauhTree.C')
status *= gSystem.CompileMacro('HMuMuTree.C')

print "Compilation status: ",status
if status==0:
    exit(-1)

#Produce framework report required by CRAB
command = "cmsRun -j FrameworkJobReport.xml -p PSet.py"
os.system(command)

from ROOT import HTauMuTauHTree
from ROOT import HTauhTauhTree
from ROOT import HMuMuTree

fileNames = []
#aFile = "file://./HTauTauAnalysis_GluGluHToTauTauM125.root"
#aFile = "file://./HTauTauAnalysis_DYJetsToLLM50.root"
aFile = "file://./HTauTauAnalysis.root"
fileNames.append(aFile)
print "Adding file: ",aFile
print "Making the mu*tau tree"
aROOTFile = TFile.Open(aFile)
aTree = aROOTFile.Get("HTauTauTree/HTauTauTree")
print "TTree entries: ",aTree.GetEntries()
HTauMuTauHTree(aTree,doSvFit).Loop()
print "Making the tau*tau tree"
aROOTFile = TFile.Open(aFile)
aTree = aROOTFile.Get("HTauTauTree/HTauTauTree")
HTauhTauhTree(aTree,doSvFit).Loop()
print "Making the mu*mu tree with"
aROOTFile = TFile.Open(aFile)
aTree = aROOTFile.Get("HTauTauTree/HTauTauTree")
HMuMuTree(aTree, doSvFit).Loop()
