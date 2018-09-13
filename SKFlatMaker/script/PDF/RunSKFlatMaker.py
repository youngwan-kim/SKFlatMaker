import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('sampletype', "DATA", VarParsing.multiplicity.singleton, VarParsing.varType.string, "sampletype: DATA/MC/PrivateMC")
options.register('PDFIDShift', "0", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDFIDShift: 0/M1/P1/..")
options.register('PDFOrder', "NLO", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDFOrder: LO/NLO/..")
options.register('PDFType', "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDFType: powheg/madgraph0/madgraph1000")
options.register('year',-1, VarParsing.multiplicity.singleton, VarParsing.varType.int, "year: Which year")
options.parseArguments()

import sys

Is2016 = False
Is2017 = False
if options.year==2016:
  Is2016 = True
elif options.year==2017:
  Is2017 = True
else:
  ErrorMgs = "year is not correct; "+str(options.year)
  sys.exit(ErrorMgs)

isMC = True
if ("DATA" in options.sampletype) or ("data" in options.sampletype) or ("Data" in options.sampletype):
  isMC = False
if ("MC" in options.sampletype) or ("mc" in options.sampletype):
  isMC = True
isPrivateSample = False
if ("Private" in options.sampletype) or ("private" in options.sampletype):
  isPrivateSample = True

options.outputFile = "SKFlatNtuple.root"
if len(options.inputFiles)==0:
  if Is2016:
    if isMC:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/TprimeBToTH_M-1700_RH_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/80000/F0074325-18A2-E811-8F6A-0CC47A57CBCC.root')
      options.outputFile = "SKFlatNtuple_2016_MC.root"
    else:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleMuon/MINIAOD/17Jul2018_ver1-v1/80000/306DAB6C-068C-E811-9E30-0242AC1C0501.root')
      options.outputFile = "SKFlatNtuple_2016_DATA.root"
  if Is2017:
    if isMC:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIIFall17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/40000/D87C6B2A-5C42-E811-8FD7-001E677926A8.root')
      options.outputFile = "SKFlatNtuple_2017_MC.root"
    else:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleMuon/MINIAOD/31Mar2018-v1/80000/54F30BE9-423C-E811-A315-0CC47A7C3410.root')
      options.outputFile = "SKFlatNtuple_2017_DATA.root"

PDFIDShift = options.PDFIDShift
PDFOrder = options.PDFOrder
PDFType = options.PDFType

print 'isMC = '+str(isMC)
print 'isPrivateSample = '+str(isPrivateSample)
print 'PDFIDShift = '+PDFIDShift
print 'PDFOrder = '+PDFOrder
print 'PDFType = '+PDFType


#### Global Tag
#### https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#2017_and_2016_re_miniAOD_94X_ver

GT_MC = '94X_mc2017_realistic_v14'
GT_DATA = '94X_dataRun2_v6'
if Is2016:
  GT_DATA = '94X_dataRun2_v10'
  GT_MC = '94X_mcRun2_asymptotic_v3'


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

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( options.inputFiles ),
  #skipEvents=cms.untracked.uint32(5),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

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
  fileName = cms.string( options.outputFile )
)

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

if Is2017:
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process,
                         runVID=False, #saves CPU time by not needlessly re-running VID
                         era='2017-Nov17ReReco')
  #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

####################
# -- Let it run -- #
####################
if Is2016:
  process.p = cms.Path(
    process.recoTree
  )
if Is2017:
  process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.recoTree
  )
