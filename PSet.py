import FWCore.ParameterSet.Config as cms

process = cms.Process("TESTPROD")

process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(1)
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50000

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                'file:///home/akalinow/scratch/CMS/HiggsCP/Prod/Crab/ntuples5.0/Enriched_miniAOD1.root'
                                #'/store/user/akalinow/EnrichedMiniAOD_v5.1/GluGluHToTauTau_M125_13TeV_powheg_pythia8_v5/GluGluHToTauTau_M125_13TeV_powheg_pythia8/GluGluHToTauTau_M125_13TeV_powheg_pythia8_v5/160524_140947/0000/Enriched_miniAOD_1.root',
                                #'/store/user/davignon/EnrichedMiniAOD/KLUB_v5/SingleMuon_Run2015D_KLUB_v5_09_06_16/SingleMuon/SingleMuon_Run2015D_KLUB_v5_09_06_16/160608_180848/0000/Enriched_miniAOD_299.root',
                                #'/store/user/tstreble/EnrichedMiniAOD/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/WJets_LO_08_06_16/160608_175058/0000/Enriched_miniAOD_1.root'
                                #'/store/user/davignon/EnrichedMiniAOD/KLUB_v5/SingleMuon_Run2015D_KLUB_v5_09_06_16/SingleMuon/SingleMuon_Run2015D_KLUB_v5_09_06_16/160608_180848/0000/Enriched_miniAOD_2.root'
                            )
)

process.p = cms.Path()

process.schedule = cms.Schedule(process.p)
