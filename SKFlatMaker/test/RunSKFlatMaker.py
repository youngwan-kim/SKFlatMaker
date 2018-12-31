import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('sampletype', "DATA", VarParsing.multiplicity.singleton, VarParsing.varType.string, "sampletype: DATA/MC/PrivateMC")
options.register('ScaleIDRange', "-999,-999", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDF Scale ID range: 1,9")
options.register('PDFErrorIDRange', "-999,-999", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDF Error ID range: 1001,1100")
options.register('PDFAlphaSIDRange', "-999,-999", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDF AlphaS ID range: 1101,1102")
options.register('PDFAlphaSScaleValue', "-999,-999", VarParsing.multiplicity.singleton, VarParsing.varType.string, "PDF AlphaS Scale values: 1.5,1.5")
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
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1/120000/80DBA5F3-16BE-E811-854E-A0369FC5E530.root')
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


ScaleIDRange = [int(options.ScaleIDRange.split(',')[0]), int(options.ScaleIDRange.split(',')[1])]
PDFErrorIDRange = [int(options.PDFErrorIDRange.split(',')[0]), int(options.PDFErrorIDRange.split(',')[1])]
PDFAlphaSIDRange = [int(options.PDFAlphaSIDRange.split(',')[0]), int(options.PDFAlphaSIDRange.split(',')[1])]
PDFAlphaSScaleValue = [float(options.PDFAlphaSScaleValue.split(',')[0]), float(options.PDFAlphaSScaleValue.split(',')[1])]

print 'isMC = '+str(isMC)
print 'isPrivateSample = '+str(isPrivateSample)
print 'ScaleIDRange = ',
print ScaleIDRange
print 'PDFErrorIDRange = ',
print PDFErrorIDRange
print 'PDFAlphaSIDRange = ',
print PDFAlphaSIDRange
print 'PDFAlphaSScaleValue = ',
print PDFAlphaSScaleValue
print 'year = '+str(options.year)


#### Global Tag
#### https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#2017_and_2016_re_miniAOD_94X_ver

GT_MC = '94X_mc2017_realistic_v14_fixJECJER'
GT_DATA = '94X_dataRun2_v6_fixJECJER'

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
process.recoTree.DataYear = cms.untracked.int32(options.year)
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

process.recoTree.ScaleIDRange = cms.untracked.vint32(ScaleIDRange)
process.recoTree.PDFErrorIDRange = cms.untracked.vint32(PDFErrorIDRange)
process.recoTree.PDFAlphaSIDRange = cms.untracked.vint32(PDFAlphaSIDRange)
process.recoTree.PDFAlphaSScaleValue = cms.untracked.vdouble(PDFAlphaSScaleValue)

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

#### EGamma ####

myEleID =  [
'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff',
'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
]

if Is2016:

  process.recoTree.electron_IDtoSave = cms.untracked.vstring(
'cutBasedElectronID-Summer16-80X-V1-veto',
'cutBasedElectronID-Summer16-80X-V1-loose',
'cutBasedElectronID-Summer16-80X-V1-medium',
'cutBasedElectronID-Summer16-80X-V1-tight',
'mvaEleID-Spring16-GeneralPurpose-V1-wp80',
'mvaEleID-Spring16-GeneralPurpose-V1-wp90',
'mvaEleID-Spring16-HZZ-V1-wpLoose',
'mvaEleID-Spring16-HZZ-V1-wpLoose', ### DUMMY
'mvaEleID-Spring16-GeneralPurpose-V1-wp80',  ### DUMMY
'mvaEleID-Spring16-GeneralPurpose-V1-wp90', ### DUMMY
'mvaEleID-Spring16-HZZ-V1-wpLoose', ### DUMMY
"heepElectronID-HEEPV70",
  )

  process.p = cms.Path(
    process.recoTree
  )

if Is2017:
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process,
                         runVID=True, #saves CPU time by not needlessly re-running VID
                         era='2017-Nov17ReReco',
                         eleIDModules=myEleID)
  #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

  runMetCorAndUncFromMiniAOD (
          process,
          isData = (not isMC), # false for MC
          fixEE2017 = True,
          fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
          postfix = "ModifiedMET"
  )
  process.recoTree.MET = cms.InputTag("slimmedMETsModifiedMET")

  #################
  ### Reapply JEC
  #################

  from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

  #### AK4

  updateJetCollection(
     process,
     jetSource = cms.InputTag('slimmedJets'),
     labelName = 'UpdatedJECslimmedJets',
     jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
  )
  process.recoTree.Jet = cms.untracked.InputTag("updatedPatJetsUpdatedJECslimmedJets")

  #### AK8

  updateJetCollection(
     process,
     jetSource = cms.InputTag('slimmedJetsAK8'),
     labelName = 'UpdatedJECslimmedJetsAK8',
     jetCorrections = ('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual corrections (always set to 1)
  )
  process.recoTree.FatJet = cms.untracked.InputTag("updatedPatJetsUpdatedJECslimmedJetsAK8")

  #### JEC Sequence

  process.jecSequence = cms.Sequence(
    process.patJetCorrFactorsUpdatedJECslimmedJets *
    process.updatedPatJetsUpdatedJECslimmedJets *
    process.patJetCorrFactorsUpdatedJECslimmedJetsAK8 *
    process.updatedPatJetsUpdatedJECslimmedJetsAK8
  )

  ##########
  #### JER 
  ##########

  process.recoTree.AK4Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK4PFchs.txt')
  process.recoTree.AK4Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Fall17_V3_MC/Fall17_V3_MC_SF_AK4PFchs.txt')
  process.recoTree.AK8Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Fall17_V3_MC/Fall17_V3_MC_PtResolution_AK8PFPuppi.txt')
  process.recoTree.AK8Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Fall17_V3_MC/Fall17_V3_MC_SF_AK8PFPuppi.txt')

  #######################################
  #### ecalBadCalibReducedMINIAODFilter
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_ecal_BadCalibReducedM
  ########################################

  process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

  baddetEcallist = cms.vuint32(
    [872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313]
  )

  process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
    "EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
    ecalMinEt        = cms.double(50.),
    baddetEcal    = baddetEcallist, 
    taggingMode = cms.bool(True),
    debug = cms.bool(False)
  )

  ###########
  #### Path
  ###########

  process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.jecSequence *
    process.fullPatMetSequenceModifiedMET *
    process.ecalBadCalibReducedMINIAODFilter *
    process.recoTree
  )
