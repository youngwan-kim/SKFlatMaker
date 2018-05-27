import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('sampletype', "DATA", VarParsing.multiplicity.singleton, VarParsing.varType.string, "sampletype: DATA/MC/PrivateMC")
options.register('PDFIDShift', "0", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDFIDShift: 0/M1/P1/..")
options.register('PDFOrder', "NLO", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDFOrder: LO/NLO/..")
options.register('PDFType', "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDFType: powheg/madgraph0/madgraph1000")
options.parseArguments()



isMC = True
if ("DATA" in options.sampletype) or ("data" in options.sampletype) or ("Data" in options.sampletype):
  isMC = False
if ("MC" in options.sampletype) or ("mc" in options.sampletype):
  isMC = True
isPrivateSample = False
if ("Private" in options.sampletype) or ("private" in options.sampletype):
  isPrivateSample = True

PDFIDShift = options. PDFIDShift
PDFOrder = options. PDFOrder
PDFType = options. PDFType

print 'isMC = '+str(isMC)
print 'isPrivateSample = '+str(isPrivateSample)
print 'PDFIDShift = '+PDFIDShift
print 'PDFOrder = '+PDFOrder
print 'PDFType = '+PDFType

#####################
# -- set by hand -- #
#####################

GT_MC = '94X_mc2017_realistic_v14'
GT_DATA = '94X_dataRun2_v6'

TESTFILE_MC = '/store/user/jskim/MiniAOD/TestMiniAOD_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/0085B26F-3642-E811-8286-008CFAE45400.root'
TESTFILE_DATA = '/store/user/jskim/MiniAOD/TestMiniAOD_SingleMuon_periodB/31Mar2018-v1/001642F1-6638-E811-B4FA-0025905B857A.root'

####################################################################################################################

process = cms.Process("NTUPLE")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( 
  SkipEvent = cms.untracked.vstring('ProductNotFound'),
  wantSummary = cms.untracked.bool(True) 
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## Source
FileName = TESTFILE_DATA
if isMC: FileName = TESTFILE_MC

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( FileName ),
  #skipEvents=cms.untracked.uint32(5),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

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
  fileName = cms.string('SKFlatNtuple.root')
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
process.recoTree.DebugLevel = cms.untracked.int32(0)
process.recoTree.StoreHLTObjectFlag = False ##FIXME

# -- Objects without Corrections -- # 
process.recoTree.Muon = cms.untracked.InputTag("slimmedMuons") # -- miniAOD -- #
process.recoTree.Electron = cms.untracked.InputTag("slimmedElectrons") # -- miniAOD -- #
process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
process.recoTree.FatJet = cms.untracked.InputTag("slimmedJetsAK8")
process.recoTree.MET = cms.InputTag("slimmedMETs")
process.recoTree.GenParticle = cms.untracked.InputTag("prunedGenParticles") # -- miniAOD -- #

process.recoTree.PDFOrder = cms.string(PDFOrder)
process.recoTree.PDFType = cms.string(PDFType)

str_this_shift = PDFIDShift
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
