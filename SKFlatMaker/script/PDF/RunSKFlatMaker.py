import FWCore.ParameterSet.Config as cms
import sys

#####################
# -- set by hand -- #
#####################
isMC = True
isPrivateSample = False

GT_MC = '94X_mc2017_realistic_v12' # -- 2017 Nov MC
GT_DATA = '94X_dataRun2_ReReco_EOY17_v2' # -- 2017 prompt reco v1

TESTFILE_MC = '/store/user/jskim/MiniAOD/TestMiniAOD_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/005DC030-D3F4-E711-889A-02163E01A62D.root'

if len(sys.argv)>2:
  TESTFILE_MC = sys.argv[2]
print TESTFILE_MC
#TESTFILE_MC = '/store/user/jskim/MiniAOD/TestMiniAOD_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/005DC030-D3F4-E711-889A-02163E01A62D.root'
#TESTFILE_MC = '/store/user/jskim/MiniAOD/TestMiniAOD_DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/006F49D4-E5FD-E711-A911-008CFAE45240.root'

TESTFILE_DATA = ''

####################################################################################################################

process = cms.Process("NTUPLE")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( 
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
  wantSummary = cms.untracked.bool(True) 
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1

## Source
FileName = TESTFILE_DATA
if isMC: FileName = TESTFILE_MC

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( FileName ),
  #skipEvents=cms.untracked.uint32(5),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

# -- Geometry and Detector Conditions (needed for a few patTuple production steps) -- #
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# -- Global Tags -- #
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

if isMC == True:
  process.GlobalTag.globaltag = cms.string(GT_MC)
else:
  process.GlobalTag.globaltag = cms.string(GT_DATA) #prompt-reco global tag


process.TFileService = cms.Service("TFileService",
  fileName = cms.string('ntuple_skim_corrected.root')
)

# -- FastFilters -- //
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
   # src = cms.InputTag("offlinePrimaryVertices"),
   src = cms.InputTag("offlineSlimmedPrimaryVertices"), # -- miniAOD -- #
   cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.FastFilters = cms.Sequence( process.goodOfflinePrimaryVertices )

#################
# -- DY Tree -- #
#################
from SKFlatMaker.SKFlatMaker.SKFlatMaker_cfi import *

process.recoTree = SKFlatMaker.clone()
process.recoTree.DebugLevel = cms.untracked.int32(1)
process.recoTree.StoreHLTObjectFlag = False

# -- Objects without Corrections -- # 
process.recoTree.Muon = cms.untracked.InputTag("slimmedMuons") # -- miniAOD -- #
process.recoTree.Electron = cms.untracked.InputTag("slimmedElectrons") # -- miniAOD -- # before smearing
process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
process.recoTree.FatJet = cms.untracked.InputTag("slimmedJetsAK8")
process.recoTree.MET = cms.InputTag("slimmedMETs")
process.recoTree.GenParticle = cms.untracked.InputTag("prunedGenParticles") # -- miniAOD -- #

if "MLM" in TESTFILE_MC:
  process.recoTree.PDFOrder = cms.string("LO")

if len(sys.argv)>3:
  str_this_shift = sys.argv[3]
  str_this_shift = str_this_shift.replace("M","-")
  str_this_shift = str_this_shift.replace("P","+")
  process.recoTree.PDFIDShift = cms.int32(int(str_this_shift))

if isPrivateSample:
  process.recoTree.LHEEventProduct = cms.untracked.InputTag("source")
  process.recoTree.LHERunInfoProduct = cms.untracked.InputTag("source")

process.recoTree.rho = cms.untracked.InputTag("fixedGridRhoFastjetAll")
process.recoTree.conversionsInputTag = cms.untracked.InputTag("reducedEgamma:reducedConversions") # -- miniAOD -- #

# -- for Track & Vertex -- #
process.recoTree.PrimaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices") # -- miniAOD -- #

# -- Else -- #
process.recoTree.PileUpInfo = cms.untracked.InputTag("slimmedAddPileupInfo")

# -- Filters -- #
process.recoTree.ApplyFilter = False

# -- Store Flags -- #
process.recoTree.StoreMuonFlag = True
process.recoTree.StoreElectronFlag = True
process.recoTree.StorePhotonFlag = True # -- photon part should be updated! later when it is necessary -- #
process.recoTree.StoreJetFlag = True
process.recoTree.StoreMETFlag = True
process.recoTree.StoreGENFlag = isMC
process.recoTree.KeepAllGen = isMC
process.recoTree.StoreLHEFlag = isMC

####################
# -- Let it run -- #
####################
process.p = cms.Path(
  process.FastFilters *
  process.recoTree
)

