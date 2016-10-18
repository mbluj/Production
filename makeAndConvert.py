#!/usr/bin/env python

import os

from ROOT import gSystem, TChain, TSystem, TFile
from PSet import process

#Produce framework report required by CRAB
command = "cmsRun -j FrameworkJobReport.xml -p PSet.py"
os.system(command)

gSystem.CompileMacro('HTauTauTree.C')
from ROOT import HTauTauTree

gSystem.CompileMacro('HTauhTauhTree.C')
from ROOT import HTauhTauhTree

fileNames = []
aFile = "file://./HTauTauAnalysis.root"
fileNames.append(aFile)

print "Adding file: ",aFile
print "Making the mu*tau tree"
aROOTFile = TFile.Open(aFile)
aTree = aROOTFile.Get("HTauTauTree/HTauTauTree")
print "TTree entries: ",aTree.GetEntries()
HTauTauTree(aTree).Loop()

print "Making the tau*tau tree"
aROOTFile = TFile.Open(aFile)
aTree = aROOTFile.Get("HTauTauTree/HTauTauTree")
HTauhTauhTree(aTree).Loop()


