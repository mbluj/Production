#-----------------------------------------
#
#Producer controller
#
#-----------------------------------------
import os, re
PyFilePath = os.environ['CMSSW_BASE']+"/src/LLRHiggsTauTau/NtupleProducer/"

#samples list (it could be moved to a cfg file for better reading
#samples = [
#]
#apply corrections?
APPLYMUCORR=False
APPLYELECORR=True

USEMVAMET=False
APPLYMETCORR=False # flag to enable (True) and disable (False) Z-recoil corrections
USE_NOHFMET = False # True to exclude HF and run on silver json

BUILDONLYOS=False #If true don't create the collection of SS candidates (and thus don't run SV fit on them)
APPLYTESCORRECTION=False # shift the central value of the tau energy scale before computing up/down variations

IsMC=True
Is25ns=True
HLTProcessName='HLT' #Different names possible, check e.g. at https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD.
if not IsMC:
    HLTProcessName='HLT' #It always 'HLT' for real data
print "HLTProcessName: ",HLTProcessName

#relaxed sets for testing purposes
TAUDISCRIMINATOR="byIsolationMVA3oldDMwoLTraw"
#PVERTEXCUT="!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2" #cut on good primary vertexes
PVERTEXCUT=""#No vertex selection in baseline selection HiggsToTauTauWorking2016#Vertices
MUCUT="isLooseMuon && pt>5"
ELECUT="pt>7"#"gsfTrack.hitPattern().numberOfHits(HitPattern::MISSING_INNER_HITS)<=1 && pt>10"
TAUCUT="(tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 || tauID('byIsolationMVArun2v1DBoldDMwLTraw')>-0.999) && pt>18" #miniAOD tau from hpsPFTauProducer have pt>18 and decaymodefinding ID
JETCUT="pt>10"
LLCUT="mass>0"
BCUT="pt>5"

# ------------------------
DO_ENRICHED=False # do True by default, both ntuples and enriched outputs are saved!
STORE_ENRICHEMENT_ONLY=True # When True and DO_ENRICHED=True only collection additional to MiniAOD standard are stored. They can be used to reproduce ntuples when used together with orygnal MiniAOD with two-file-solution
# ------------------------

is80X = True if 'CMSSW_8' in os.environ['CMSSW_VERSION'] else False# True to run in 80X (2016), False to run in 76X (2015)
print "is80X: " , is80X
is92X = True if 'CMSSW_9_2' in os.environ['CMSSW_VERSION'] else False# True to run in 9XY (2017), False to run in 76X (2015) or 80X (2016)
print "is92X: " , is92X
is94X = True if 'CMSSW_9_4' in os.environ['CMSSW_VERSION'] else False# True to run in 9XY (2017), False to run in 76X (2015) or 80X (2016)
print "is94X: " , is94X
##
## Standard sequence
##

if is80X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_80X.py")
elif is92X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_92X.py")
elif is94X:
    execfile(PyFilePath+"python/HiggsTauTauProducer_94X.py")    
else :
    execfile(PyFilePath+"python/HiggsTauTauProducer.py")

### ----------------------------------------------------------------------
### Source, better to use sample to run on batch
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/home/akalinow/scratch/CMS/HiggsCP/Data/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM/0498CD6A-CC42-E811-95D3-008CFA1CB8A8.root'
    )
)

#Limited nEv for testing purposes. -1 to run all events
process.maxEvents.input = 1000
#process.source.eventsToProcess = cms.untracked.VEventRange('1:2795204') #MET test

# JSON mask for data --> defined in the lumiMask file
# from JSON file
#if not IsMC:
#  execfile(PyFilePath+"python/lumiMask.py")
#  process.source.lumisToProcess = LUMIMASK

##
## Output file
##

process.TFileService=cms.Service('TFileService',fileName=cms.string('HTauTauAnalysis.root'))

if DO_ENRICHED:
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string('Enriched_miniAOD.root'),
        outputCommands = cms.untracked.vstring('keep *'),
        fastCloning     = cms.untracked.bool(False),
        #Compression settings from MiniAOD allowing to save about 10% of disc space compared to defults ->
        compressionAlgorithm = cms.untracked.string('LZMA'),
        compressionLevel = cms.untracked.int32(4),
        dropMetaData = cms.untracked.string('ALL'),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        overrideInputFileSplitLevels = cms.untracked.bool(True)
        # <-
    )
    if STORE_ENRICHEMENT_ONLY:
        # Store only additional collections compared to MiniAOD necessary to reproduce ntuples (basically MVAMET, lepton pairs with SVFit and corrected jets)
        # Size of about 10% of full EnrichedMiniAOD
        process.out.outputCommands.append('drop *')
        process.out.outputCommands.append('keep *_SVllCand_*_*')
        process.out.outputCommands.append('keep *_SVbypass_*_*')
        process.out.outputCommands.append('keep *_barellCand_*_*')
        process.out.outputCommands.append('keep *_corrMVAMET_*_*')
        process.out.outputCommands.append('keep *_MVAMET_*_*')
        process.out.outputCommands.append('keep *_jets_*_*')
        process.out.outputCommands.append('keep *_patJetsReapplyJEC_*_*')
        process.out.outputCommands.append('keep *_softLeptons_*_*')
        process.out.outputCommands.append('keep *_genInfo_*_*')
        #process.out.fileName = 'EnrichementForMiniAOD.root' #FIXME: change name of output file?
    process.end = cms.EndPath(process.out)

#process.options = cms.PSet(skipEvent =  cms.untracked.vstring('ProductNotFound')),
#process.p = cms.EndPath(process.HTauTauTree)
process.p = cms.Path(process.Candidates)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.categories.append('onlyError')
#process.MessageLogger.cerr.onlyError=cms.untracked.PSet(threshold  = cms.untracked.string('ERROR'))
#process.MessageLogger.cerr.threshold='ERROR'
#process.MessageLogger = cms.Service("MessageLogger",
#	destinations = cms.untracked.vstring('log.txt')
#)
#process.MessageLogger.threshold = cms.untracked.string('ERROR')

#processDumpFile = open('process.dump' , 'w')
#print >> processDumpFile, process.dumpPython()

if process.source.fileNames[0].find("92X_upgrade2017_realistic_v10")>-1:
    process.GlobalTag.globaltag = '92X_upgrade2017_realistic_v10'
    print "Switching GlobalTag to", process.GlobalTag.globaltag
