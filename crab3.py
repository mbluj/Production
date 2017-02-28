from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'WplusHToTauTau_M120'
config.General.workArea = 'v1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.sendExternalFolder = True

config.section_("Data")
config.Data.inputDataset = '/WplusHToTauTau_M120_13TeV_powheg_pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2000 #number of files per jobs
config.Data.totalUnits =  -1 #number of event
config.Data.outLFNDirBase = '/store/user/akalinow/EnrichedMiniAOD/WplusHToTauTau_M120_v1/'
config.Data.publication = False
config.Data.outputDatasetTag = 'WplusHToTauTau_M120_v1'



config.section_("Site")
config.Site.storageSite = 'T2_PL_Swierk'
config.Site.blacklist = ['T2_KR_*','T2_CN_*','T2_BR_*','T2_US_Florida','T2_US_UCSD']
