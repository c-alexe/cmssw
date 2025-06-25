import CRABClient
from CRABClient.UserUtilities import config 

config = config()

config.General.requestName = 'refit_Data_F_notFinishedLumis'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'refitFromMINIAOD.py'
config.JobType.outputFiles = ['globalcor_0.root']

config.Data.inputDataset = '/Muon/Run2022F-22Sep2023-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'crab_projects/crab_refit_Data_F/results/notFinishedLumis.json'
# 'https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions22/Cert_Collisions2022_eraF_360390_362167_Muon.json'
config.Data.publication = False
config.Data.outputDatasetTag = 'Run2022F-22Sep2023-v2-with-CVH'

config.Site.storageSite = 'T2_IT_Pisa'
