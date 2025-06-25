import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('RECO2',Run3)

'''
options = VarParsing.VarParsing('analysis')
options.register('globalTag',
                 "auto:run2_mc", # default value
                 VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                 VarParsing.VarParsing.varType.string, # string, int, or float
                 "globalTag")
options.parseArguments()
'''
# fix options if running on CRAB
globalTag = '130X_dataRun3_PromptAnalysis_v1'
inputFiles = 'root://xrootd-cms.infn.it//store/data/Run2022F/Muon/MINIAOD/22Sep2023-v2/2550000/00f90c0a-f66e-4bf8-9985-377caf239002.root'
maxEvents = 100

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
 
    
process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring(inputFiles),
                            secondaryFileNames = cms.untracked.vstring()
                            )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(maxEvents)
)

#process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step1 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('step1_RECO.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.RECOSIMoutput.outputCommands = cms.untracked.vstring("keep *_myRefittedTracks_*_*")

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globalTag, '')

# load geometry needed by g4e propagator (for data it must be the one from the GT of the MC)
process.GlobalTag.toGet = cms.VPSet(
                            cms.PSet(
                                record = cms.string("GeometryFileRcd"),tag = cms.string("GeometryFileExtended2021_132DD4hepV3"),label = cms.untracked.string("Extended"),
                                ),
                            )

####################################################################
# load 3d field map and use it for g4e propagator, geant4 internals via geometry producer and a few other places related to the track refit
####################################################################
from MagneticField.ParametrizedEngine.parametrizedMagneticField_PolyFit3D_cfi import ParametrizedMagneticFieldProducer as PolyFit3DMagneticFieldProducer
process.PolyFit3DMagneticFieldProducer = PolyFit3DMagneticFieldProducer
fieldlabel = "PolyFit3DMf"
process.PolyFit3DMagneticFieldProducer.label = fieldlabel
process.stripCPEESProducer.MagneticFieldLabel = fieldlabel
process.StripCPEfromTrackAngleESProducer.MagneticFieldLabel = fieldlabel
process.TransientTrackBuilderESProducer.MagneticFieldLabel = fieldlabel 
process.siPixelTemplateDBObjectESProducer.MagneticFieldLabel = fieldlabel
process.templates.MagneticFieldLabel = fieldlabel

# track refit stuff
from TrackPropagation.Geant4e.geantRefit_cff import geopro
process.load("TrackPropagation.Geant4e.geantRefit_cff")
#process.Geant4eTrackRefitter.src = cms.InputTag("ALCARECOTkAlZMuMu") # for ALCARECO
#process.Geant4eTrackRefitter.usePropagatorForPCA = cms.bool(True) # for ALCARECO

process.geopro.MagneticFieldLabel = fieldlabel
process.Geant4ePropagator.MagneticFieldLabel = fieldlabel

process.Geant4ePropagator.ForCVH = True
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *

import RecoTracker.TrackProducer.trackProducerFromPatMuons_cfi
process.tracksFromMuons  = RecoTracker.TrackProducer.trackProducerFromPatMuons_cfi.trackProducerFromPatMuons.clone(
    src = "slimmedMuons",
    innerTrackOnly = True,
    ptMin = 10. 
)

process.trackrefit = cms.EDProducer('ResidualGlobalCorrectionMakerG4e',
                                    #src = cms.InputTag("ALCARECOTkAlZMuMu"),
                                    src = cms.InputTag("tracksFromMuons"),
                                    muons = cms.InputTag("slimmedMuons"),
                                    #genParticles = cms.InputTag("packedGenParticles"), # uncomment line for MC
                                    #pileupInfo = cms.InputTag("slimmedAddPileupInfo"), # uncomment line for MC
                                    fitFromGenParms = cms.bool(False),
                                    fitFromSimParms = cms.bool(False),
                                    fillTrackTree = cms.bool(True),
                                    fillGrads = cms.bool(False),
                                    fillJac = cms.bool(False),
                                    fillRunTree = cms.bool(False),
                                    doGen = cms.bool(False), # change for MC
                                    doSim = cms.bool(False),
                                    requireGen = cms.bool(False),
                                    doMuons = cms.bool(True),
                                    doMuonAssoc = cms.bool(True),
                                    doTrigger = cms.bool(True),
                                    triggers = cms.vstring("HLT_IsoMu24"),
                                    doRes = cms.bool(False),
                                    bsConstraint = cms.bool(False),
                                    applyHitQuality = cms.bool(True),
                                    useIdealGeometry = cms.bool(False),
                                    corFiles = cms.vstring(),
                                    MagneticFieldLabel = cms.string("PolyFit3DMf"),
                                    )

process.trackrefitideal = cms.EDProducer('ResidualGlobalCorrectionMakerG4e',
                                         #src = cms.InputTag("ALCARECOTkAlZMuMu"),
                                         src = cms.InputTag("tracksFromMuons"),
                                         muons = cms.InputTag("slimmedMuons"),
                                         #genParticles = cms.InputTag("packedGenParticles"), # uncomment line for MC
                                         #pileupInfo = cms.InputTag("slimmedAddPileupInfo"), # uncomment line for MC
                                         fitFromGenParms = cms.bool(False),
                                         fitFromSimParms = cms.bool(False),
                                         fillTrackTree = cms.bool(True),
                                         fillGrads = cms.bool(False),
                                         fillJac = cms.bool(False),
                                         fillRunTree = cms.bool(False),
                                         doGen = cms.bool(False), # change for MC
                                         doSim = cms.bool(False),
                                         requireGen = cms.bool(False),
                                         doMuons = cms.bool(True),
                                         doMuonAssoc = cms.bool(True),
                                         doTrigger = cms.bool(True),
                                         triggers = cms.vstring("HLT_IsoMu24"),
                                         doRes = cms.bool(False),
                                         bsConstraint = cms.bool(False),
                                         applyHitQuality = cms.bool(True),
                                         useIdealGeometry = cms.bool(True),
                                         corFiles = cms.vstring(),
                                         MagneticFieldLabel = cms.string(""),
                                         )

process.mergedGlobalIdxs = cms.EDProducer("GlobalIdxProducer",
                                          src0 = cms.InputTag("trackrefit", "globalIdxs"),
                                          src1 = cms.InputTag("trackrefitideal", "globalIdxs")
                                          )

# import RecoTracker.TrackProducer.TrackRefitter_cfi
# process.myRefittedTracks = RecoTracker.TrackProducer.TrackRefitter_cfi.TrackRefitter.clone(
#     src = 'tracksFromMuons',
#     NavigationSchool = '',
#     Fitter = 'FlexibleKFFittingSmoother'
# )
                                         
# Path and EndPath definitions
process.reconstruction_step = cms.Path(process.geopro+process.tracksFromMuons*process.trackrefit, #process.myRefittedTracks, 
                                       cms.Task(process.TTRHBuilderAngleAndTemplate, process.templates, process.SiPixelTemplateStoreESProducer))

process.endjob_step = cms.EndPath(process.endOfProcess)
#process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
process.schedule = cms.Schedule(process.reconstruction_step,process.endjob_step)
