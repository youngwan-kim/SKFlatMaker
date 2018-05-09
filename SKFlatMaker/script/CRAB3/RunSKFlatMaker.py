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

GT_MC = '94X_mc2017_realistic_v14' # -- 2017 Nov MC
GT_DATA = '94X_dataRun2_v6' # -- 2017 prompt reco v1

TESTFILE_MC = '/store/user/jskim/SKFlat/TestMiniAOD_DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/0085B26F-3642-E811-8286-008CFAE45400.root'
TESTFILE_DATA = '/stroe/user/jskim/SKFlat/TestMiniAOD_SingleMuon_periodB/31Mar2018-v1/001642F1-6638-E811-B4FA-0025905B857A.root'

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

##################################
# For E/gamma correction and VID #
##################################
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
#setupEgammaPostRecoSeq(process,applyEnergyCorrections=False,
#                       applyVIDOnCorrectedEgamma=False,
#                       isMiniAOD=True)


######################
# MET Phi Correction #
######################
#from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#runMetCorAndUncFromMiniAOD(process, isData=True )  #For MC isData=False

#################
# -- DY Tree -- #
#################
from SKFlatMaker.SKFlatMaker.SKFlatMaker_cfi import *
from SKFlatMaker.SKFlatMaker.PUreweight2012_cff import *

process.recoTree = SKFlatMaker.clone()
process.recoTree.DebugLevel = cms.untracked.int32(0)
process.recoTree.StoreHLTObjectFlag = False ##FIXME

# -- Objects without Corrections -- # 
process.recoTree.Muon = cms.untracked.InputTag("slimmedMuons") # -- miniAOD -- #
process.recoTree.Electron = cms.untracked.InputTag("slimmedElectrons") # -- miniAOD -- # before smearing
process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
#process.recoTree.SmearedElectron = "calibratedPatElectrons" # -- Smeared Electron
#process.recoTree.SmearedPhoton = "calibratedPatPhotons" # -- Smeared Photon
process.recoTree.SmearedElectron = cms.untracked.InputTag("slimmedElectrons")
process.recoTree.SmearedPhoton = cms.untracked.InputTag("slimmedPhotons")
#process.recoTree.SmearedElectron = cms.untracked.InputTag("calibratedPatElectrons") # -- Smeared Electron
#process.recoTree.SmearedPhoton = cms.untracked.InputTag("calibratedPatPhotons") # -- Smeared Photon
#process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
#process.recoTree.MET = cms.untracked.InputTag("slimmedMETs") # -- miniAOD -- #
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

# -- for electrons -- # chaged to 2017 ID Map
process.recoTree.rho = cms.untracked.InputTag("fixedGridRhoFastjetAll")
process.recoTree.conversionsInputTag = cms.untracked.InputTag("reducedEgamma:reducedConversions") # -- miniAOD -- #

# -- for photons -- #
process.recoTree.effAreaChHadFile = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt")
process.recoTree.effAreaNeuHadFile= cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt")
process.recoTree.effAreaPhoFile   = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt")




# -- Corrections -- #

# -- JEC
JEC_files = ('Fall17_17Nov2017BCDEF_V6_DATA', 'Fall17_17Nov2017_V6_MC')
if isMC:
  jecFile = JEC_files[1]
else:
  jecFile = JEC_files[0]

from CondCore.CondDB.CondDB_cfi import CondDB
if hasattr(CondDB, 'connect'): delattr(CondDB, 'connect')
process.jec = cms.ESSource("PoolDBESSource",CondDB,
    connect = cms.string('sqlite_fip:SKFlatMaker/SKFlatMaker/data/JEC/db/%s.db'%jecFile),            
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_%s_AK4PF"%jecFile),
            label= cms.untracked.string("AK4PF")),
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFchs"%jecFile),
            label= cms.untracked.string("AK4PFchs")),
        cms.PSet(
            record = cms.string("JetCorrectionsRecord"),
            tag = cms.string("JetCorrectorParametersCollection_%s_AK4PFPuppi"%jecFile),
            label= cms.untracked.string("AK4PFPuppi")),
    )
)
process.es_prefer_jec = cms.ESPrefer("PoolDBESSource","jec")
print "JEC based on", process.jec.connect
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
if isMC:
  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    )
else :
  updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
    )

process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
process.recoTree.FatJet = cms.untracked.InputTag("slimmedJetsAK8")

# -- MET Correction
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, isData= not isMC, electronColl=cms.InputTag('slimmedElectrons'), jetCollUnskimmed=cms.InputTag('slimmedJets'))
process.recoTree.MET = cms.InputTag("slimmedMETs")

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
  #process.egammaPostRecoSeq *
  #process.calibratedPatElectrons *
  #process.calibratedPatPhotons *  
  #process.selectedElectrons *
  #process.selectedPhotons *
  #process.egmPhotonIDSequence *
  #process.egmGsfElectronIDSequence *
  #process.electronIDValueMapProducer *
  #process.fullPatMetSequence *  #This is the phi corrections part
  process.recoTree
)
