import os
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register('sampletype', "DATA", VarParsing.multiplicity.singleton, VarParsing.varType.string, "sampletype: DATA/MC")
options.register('weightmap', "", VarParsing.multiplicity.singleton, VarParsing.varType.string, "weightmap file path or string: Scale[1,9],PDF[1001,1100],AlphaS[1101,1102],AlphaSScale[1.5,1.5],PSSyst[1,45]")
options.register('era',-1, VarParsing.multiplicity.singleton, VarParsing.varType.string, "era: Which era? 2016preVFP, 2016postVFP, 2017, 2018")
options.register('year',-1, VarParsing.multiplicity.singleton, VarParsing.varType.string, "Deprecated. Use 'era'")
options.register('debug',0, VarParsing.multiplicity.singleton, VarParsing.varType.int, "print level")
options.setDefault('outputFile','SKFlatNtuple.root')
options.parseArguments()

import sys

if options.year!=-1:
  sys.exit("Option 'year' is deprecated. Use 'era'")  

Is2016preVFP = False
Is2016postVFP = False
Is2017 = False
Is2018 = False
if options.era=='2016preVFP':
  Is2016preVFP = True
elif options.era=='2016postVFP':
  Is2016postVFP = True
elif options.era=='2017':
  Is2017 = True
elif options.era=='2018':
  Is2018 = True
else:
  ErrorMgs = "era is not correct; "+str(options.era)
  sys.exit(ErrorMgs)

isMC = True
if "data" in options.sampletype.lower():
  isMC = False
if "mc" in options.sampletype.lower():
  isMC = True

if len(options.inputFiles)==0:
  if Is2016preVFP:
    if isMC:
      #options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL16MiniAODAPV/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v8-v1/270000/19589E57-61A2-7049-9343-788B6D00A07D.root') #MiniAODv1
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/008E4139-6019-CE4C-B83C-A849F56F57B3.root') #MiniAODv2
      options.outputFile = options.outputFile.replace(".root","_2016preVFP_MC.root")
    else:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleMuon/MINIAOD/21Feb2020_ver2_UL2016_HIPM-v1/70000/6805D8D5-79EB-0346-9E07-B267BE5BD848.root')
      options.outputFile = options.outputFile.replace(".root","_2016preVFP_DATA.root")
  if Is2016postVFP:
    if isMC:
      #options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL16MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v13-v2/110000/F0A9DC9B-B0D6-804A-BD80-57F7FA06FACB.root') #MiniAODv1
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/2520000/04A698D5-2AF9-B548-9A6D-DB5AFE92F0A6.root') #MiniAODv2
      options.outputFile = options.outputFile.replace(".root","_2016postVFP_MC.root")
    else:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/data/Run2016H/SingleMuon/MINIAOD/21Feb2020_UL2016-v1/250000/D6CA5522-DB00-CE45-BF90-06176D4DBC86.root')
      options.outputFile = options.outputFile.replace(".root","_2016postVFP_DATA.root")
  elif Is2017:
    if isMC:
      #options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL17MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v6-v2/260000/04B074CC-43FC-6045-ACD2-D49AC762E5A6.root') #MiniAODv1
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL17MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mc2017_realistic_v9-v2/100000/00827E5C-253F-F942-9751-3F3277340A21.root') #MiniAODv2
      options.outputFile = options.outputFile.replace(".root","_2017_MC.root")
    else:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleMuon/MINIAOD/09Aug2019_UL2017-v1/210000/49E82EBF-E7AF-8645-B5C7-6F2B208F3F5F.root')
      options.outputFile = options.outputFile.replace(".root","_2017_DATA.root")
  elif Is2018:
    if isMC:
      #options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer19UL18MiniAOD/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/70000/EA2F219D-8534-7B4D-AF83-5D91AF448EC6.root') #MiniAODv1
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/100000/4E295BA9-D9F7-6643-B993-57789E70C0CB.root') #MiniAODv2
      options.outputFile = options.outputFile.replace(".root","_2018_MC.root")
    else:
      options.inputFiles.append('root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/MINIAOD/12Nov2019_UL2018_rsb-v1/240000/FE143ADE-E9E4-3143-8754-C9ECA89F3541.root')
      options.outputFile = options.outputFile.replace(".root","_2018_DATA.root")

weightmap={}
if options.weightmap!="":
  if os.path.exists(options.weightmap):
    with open(options.weightmap) as f:
      for line in f.readlines():
        line=line.split("#",1)[0]
        words=line.split()
        if len(words)==2:
          weightmap[words[0]]=words[1].split(",")
  else:
    for line in options.weightmap.replace("],","@").rstrip("]").split("@"):
      words=line.replace("[","@").split("@")
      weightmap[words[0]]=words[1].split(",")
  for key in weightmap:
    if key=="AlphaSScale":
      weightmap[key]=[float(i) for i in weightmap[key]]
    else:
      weightmap[key]=[int(i) for i in weightmap[key]]
    if key in ["PDF","Scale", "PSSyst"] and len(weightmap[key])==2:
      weightmap[key]=range(weightmap[key][0],weightmap[key][1]+1)
      
#### Global Tag
#### https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis

GT_MC = ''
GT_DATA = ''

if Is2016preVFP:
  GT_MC = '106X_mcRun2_asymptotic_preVFP_v11'
  GT_DATA = '106X_dataRun2_v35'
if Is2016postVFP:
  GT_MC = '106X_mcRun2_asymptotic_v17'
  GT_DATA = '106X_dataRun2_v35'
elif Is2017:
  GT_MC = '106X_mc2017_realistic_v9'
  GT_DATA = '106X_dataRun2_v35'
elif Is2018:
  GT_MC = '106X_upgrade2018_realistic_v16_L1v1'
  GT_DATA = '106X_dataRun2_v35'
## if MiniAODv1: v9->v8 (17,MC), v16->v15 (18,MC)


print 'GT_MC = '+GT_MC
print 'GT_DATA = '+GT_DATA
print 'isMC = '+str(isMC)
print 'era = '+str(options.era)
for key in weightmap:
  if len(weightmap[key])>5:
    print "weight_{} = [{},{},...,{},{}]".format(key,weightmap[key][0],weightmap[key][1],weightmap[key][-2],weightmap[key][-1])
  else:
    print "weight_{} =".format(key),weightmap[key]

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

# -- Geometry and Detector Conditions (needed for a few patTuple production steps) -- #
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# -- Global Tags -- #
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JecGlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

if isMC == True:
  process.GlobalTag.globaltag = cms.string(GT_MC)
else:
  process.GlobalTag.globaltag = cms.string(GT_DATA)


process.TFileService = cms.Service("TFileService",
  fileName = cms.string( options.outputFile )
)

#################
# -- DY Tree -- #
#################
from SKFlatMaker.SKFlatMaker.SKFlatMaker_cfi import *

process.recoTree = SKFlatMaker.clone()
process.recoTree.DataYear = cms.untracked.int32(int(options.era[0:4]))
process.recoTree.DebugLevel = cms.untracked.int32(options.debug)
process.recoTree.StoreHLTObjectFlag = False ##FIXME

# -- Objects without Corrections -- # 
process.recoTree.Muon = cms.untracked.InputTag("slimmedMuons") # -- miniAOD -- #
process.recoTree.Electron = cms.untracked.InputTag("slimmedElectrons") # -- miniAOD -- #
process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
process.recoTree.FatJet = cms.untracked.InputTag("slimmedJetsAK8")
process.recoTree.MET = cms.InputTag("slimmedMETs")
process.recoTree.PuppiMET = cms.InputTag("slimmedMETsPuppi")
process.recoTree.GenParticle = cms.untracked.InputTag("prunedGenParticles") # -- miniAOD -- #

for key in weightmap:
  if key=="AlphaSScale":
    setattr(process.recoTree,"weight_"+key,cms.untracked.vdouble(weightmap[key]))
  else:
    setattr(process.recoTree,"weight_"+key,cms.untracked.vint32(weightmap[key]))

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
process.recoTree.StorePuppiMETFlag = True
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
myPhoID =  [
'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff',
]


#### GenHFHadronMatcher ####
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
  genParticles = cms.InputTag("prunedGenParticles"),
  jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
)

from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenCHadron
process.matchGenCHadron = matchGenCHadron.clone(
  genParticles = cms.InputTag("prunedGenParticles"),
  jetFlavourInfos = cms.InputTag("slimmedGenJetsFlavourInfos"),
)

process.matchGenSeq = cms.Sequence(process.matchGenBHadron + process.matchGenCHadron)

#### BJetEnergyCorrection ####
from SKFlatMaker.SKFlatMaker.BJetEnergyCorr_cff import bJetNN
from SKFlatMaker.SKFlatMaker.BJetEnergyCorr_cff import cJetNN
process.load('SKFlatMaker.SKFlatMaker.BJetEnergyCorr_cff')

if Is2016preVFP:

  print "################"
  print "Running 2016preVFP"
  print "################"

  ###########################
  #### Rerun EGammaPostReco
  #### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#Recipe_for_running_Scales_and_sm
  ###########################
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process,
                         runVID=True,
                         era='2016preVFP-UL',
                         eleIDModules=myEleID,
                         phoIDModules=myPhoID,
                         runEnergyCorrections=True
  )
  #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

  #################
  ### Reapply JEC
  ### https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
  ### https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  ### 106X_dataRun2_v35 -> Summer19UL16_V7_DATA
  ### 106X_mcRun2_asymptotic_preVFP_v11 -> Summer19UL16APV_V7_MC
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
  #### https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#2016_data
  ##########
  process.recoTree.AK4Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs.txt')
  process.recoTree.AK4Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.txt')
  process.recoTree.AK8Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_PtResolution_AK8PFPuppi.txt')
  process.recoTree.AK8Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16APV_JRV3_MC/Summer20UL16APV_JRV3_MC_SF_AK8PFPuppi.txt')

  #################
  #### BJetEnergyRegression
  #################
  bJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2016.pb")
  bJetNN.outputFormulas = cms.vstring(["at(0)*0.31976690888404846+1.047176718711853","0.5*(at(2)-at(1))*0.31976690888404846"])
  cJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2016.pb")
  cJetNN.outputFormulas = cms.vstring(["at(0)*0.28862622380256653+0.9908722639083862","0.5*(at(2)-at(1))*0.28862622380256653"])

  #################
  #### Update MET
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription
  #### This is independent of jecSequence, but it rather reapply JEC/JER using GT withing this MET corrector module
  #################
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  runMetCorAndUncFromMiniAOD(process,
                           isData=(not isMC),
                           )

  from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
  makePuppiesFromMiniAOD( process, True );
  runMetCorAndUncFromMiniAOD(process,
                             isData=is_data,
                             metType="Puppi",
                             postfix="Puppi",
                             jetFlavor="AK4PFPuppi",
                             )

  #########################
  #### L1Prefire reweight
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  #########################
  from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
  process.prefiringweight = l1PrefiringWeightProducer.clone(
    TheJets = cms.InputTag("updatedPatJetsUpdatedJECslimmedJets"), #this should be the slimmedJets collection with up to date JECs !
    DataEraECAL = cms.string("UL2016preVFP"),
    DataEraMuon = cms.string("2016preVFP"),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
  )

  ###########################
  #### Rochester correction
  ###########################
  process.recoTree.roccorPath = cms.string('SKFlatMaker/SKFlatMaker/data/roccor.Run2.v5/RoccoR2016aUL.txt')

  ###########################
  #### BadPFMuonDz filter
  ###########################
  from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
  process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
  )

  ###########
  #### Path
  ###########

  process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.jecSequence *
    process.fullPatMetSequence *
    process.egmPhotonIDSequence *
    process.puppiMETSequence *
    process.fullPatMetSequencePuppi
  )
  process.p *= process.BadPFMuonFilterUpdateDz
  if isMC:
    process.recoTree.StoreL1PrefireFlag = cms.untracked.bool(True)
    process.p *= process.prefiringweight

    ## GenHFHadron
    process.p *= process.matchGenSeq

  ## BJetEnergyRegression
  process.p *= process.bJetCorrSeq

  process.p *= process.recoTree

if Is2016postVFP:

  print "################"
  print "Running 2016postVFP"
  print "################"

  ###########################
  #### Rerun EGammaPostReco
  #### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#Recipe_for_running_Scales_and_sm
  ###########################
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process,
                         runVID=True,
                         era='2016postVFP-UL',
                         eleIDModules=myEleID,
                         phoIDModules=myPhoID,
                         runEnergyCorrections=True
  )
  #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

  #################
  ### Reapply JEC
  ### https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
  ### https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  ### 106X_dataRun2_v35 -> Summer19UL16_V7_DATA
  ### 106X_mcRun2_asymptotic_v17 -> Summer19UL16_V7_MC
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
  #### https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#2016_data
  ##########
  process.recoTree.AK4Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.txt')
  process.recoTree.AK4Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_SF_AK4PFchs.txt')
  process.recoTree.AK8Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_PtResolution_AK8PFPuppi.txt')
  process.recoTree.AK8Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer20UL16_JRV3_MC/Summer20UL16_JRV3_MC_SF_AK8PFPuppi.txt')

  #################
  #### BJetEnergyRegression
  #################
  bJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2016.pb")
  bJetNN.outputFormulas = cms.vstring(["at(0)*0.31976690888404846+1.047176718711853","0.5*(at(2)-at(1))*0.31976690888404846"])
  cJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2016.pb")
  cJetNN.outputFormulas = cms.vstring(["at(0)*0.28862622380256653+0.9908722639083862","0.5*(at(2)-at(1))*0.28862622380256653"])

  #################
  #### Update MET
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription
  #### This is independent of jecSequence, but it rather reapply JEC/JER using GT withing this MET corrector module
  #################
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  runMetCorAndUncFromMiniAOD(process,
                           isData=(not isMC),
                           )

  from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
  makePuppiesFromMiniAOD( process, True );
  runMetCorAndUncFromMiniAOD(process,
                             isData=is_data,
                             metType="Puppi",
                             postfix="Puppi",
                             jetFlavor="AK4PFPuppi",
                             )

  #########################
  #### L1Prefire reweight
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  #########################
  from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
  process.prefiringweight = l1PrefiringWeightProducer.clone(
    TheJets = cms.InputTag("updatedPatJetsUpdatedJECslimmedJets"), #this should be the slimmedJets collection with up to date JECs !
    DataEraECAL = cms.string("UL2016postVFP"),
    DataEraMuon = cms.string("2016postVFP"),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
  )

  ###########################
  #### Rochester correction
  ###########################
  process.recoTree.roccorPath = cms.string('SKFlatMaker/SKFlatMaker/data/roccor.Run2.v5/RoccoR2016bUL.txt')

  ###########################
  #### BadPFMuonDz filter
  ###########################
  from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
  process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
  )

  ###########
  #### Path
  ###########

  process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.jecSequence *
    process.fullPatMetSequence *
    process.egmPhotonIDSequence *
    process.puppiMETSequence *
    process.fullPatMetSequencePuppi
  )
  process.p *= process.BadPFMuonFilterUpdateDz
  if isMC:
    process.recoTree.StoreL1PrefireFlag = cms.untracked.bool(True)
    process.p *= process.prefiringweight

    ## GenHFHadron
    process.p *= process.matchGenSeq

  ## BJetEnergyRegression
  process.p *= process.bJetCorrSeq

  process.p *= process.recoTree

elif Is2017:

  print "################"
  print "Running 2017"
  print "################"

  ###########################
  #### Rerun EGammaPostReco
  #### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#Recipe_for_running_Scales_and_sm
  ###########################
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process,
                         runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                         era='2017-UL',
                         eleIDModules=myEleID,
                         phoIDModules=myPhoID,
                         runEnergyCorrections=True
  )
  #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

  #################
  ### Reapply JEC
  ### https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
  ### https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  ### 106X_dataRun2_v33 -> Summer19UL17_V5_DATA
  ### 106X_mc2017_realistic_v8 -> Summer19UL17_V5_MC
  ### AK8PFPuppi corrections are copied from prelegacy
  #################

  #### Update jets
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
  #### https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#2017_data
  ##########
  process.recoTree.AK4Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt')
  process.recoTree.AK4Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt')
  # No JER for AK8Jet. Use old one temporary.
  process.recoTree.AK8Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_PtResolution_AK8PFPuppi.txt')
  process.recoTree.AK8Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_SF_AK8PFPuppi.txt')

  #################
  #### BJetEnergyRegression
  #################
  bJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2017.pb")
  bJetNN.outputFormulas = cms.vstring(["at(0)*0.28225210309028625+1.055067777633667","0.5*(at(2)-at(1))*0.28225210309028625"])
  cJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2017.pb")
  cJetNN.outputFormulas = cms.vstring(["at(0)*0.24718524515628815+0.9927206635475159","0.5*(at(2)-at(1))*0.24718524515628815"])

  #################
  #### Update MET
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription
  #### This is independent of jecSequence, but it rather reapply JEC/JER using GT withing this MET corrector module
  #################
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  runMetCorAndUncFromMiniAOD(process,
                           isData=(not isMC),
                           )

  from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
  makePuppiesFromMiniAOD( process, True );
  runMetCorAndUncFromMiniAOD(process,
                             isData=is_data,
                             metType="Puppi",
                             postfix="Puppi",
                             jetFlavor="AK4PFPuppi",
                             )

  #########################
  #### L1Prefire reweight
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  #########################
  from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
  process.prefiringweight = l1PrefiringWeightProducer.clone(
    TheJets = cms.InputTag("updatedPatJetsUpdatedJECslimmedJets"), #this should be the slimmedJets collection with up to date JECs !
    DataEraECAL = cms.string("UL2017BtoF"),
    DataEraMuon = cms.string("20172018"),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
  )

  ###########################
  #### Rochester correction
  ###########################
  process.recoTree.roccorPath = cms.string('SKFlatMaker/SKFlatMaker/data/roccor.Run2.v5/RoccoR2017UL.txt')

  ###########################
  #### BadPFMuonDz filter
  ###########################
  from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
  process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
  )

  ###########
  #### Path
  ###########

  process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.jecSequence *
    process.fullPatMetSequence *
    process.egmPhotonIDSequence *
    process.puppiMETSequence *
    process.fullPatMetSequencePuppi
  )
  process.p *= process.BadPFMuonFilterUpdateDz
  if isMC:
    process.recoTree.StoreL1PrefireFlag = cms.untracked.bool(True)
    process.p *= process.prefiringweight

    ## GenHFHadron
    process.p *= process.matchGenSeq

  ## BJetEnergyRegression
  process.p *= process.bJetCorrSeq

  process.p *= process.recoTree

elif Is2018:

  print "################"
  print "Running 2018"
  print "################"

  ###########################
  #### Rerun EGammaPostReco
  #### https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018#Recipe_for_running_Scales_and_sm
  ###########################
  from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
  setupEgammaPostRecoSeq(process,
                         runVID=True, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                         era='2018-UL',
                         eleIDModules=myEleID,
                         phoIDModules=myPhoID,
                         runEnergyCorrections=True
  )  
  #a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

  #################
  ### Reapply JEC
  ### https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#CorrPatJets
  ### https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  ### 106X_dataRun2_v33 -> Summer19UL18_V5_DATA
  ### 106X_upgrade2018_realistic_v15_L1v1 -> Summer19UL18_V5_MC
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
  #### https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#2018_data
  ##########
  process.recoTree.AK4Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.txt')
  process.recoTree.AK4Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK4PFchs.txt')
  process.recoTree.AK8Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_PtResolution_AK8PFPuppi.txt')
  process.recoTree.AK8Jet_JER_SF_filepath    = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL18_JRV2_MC/Summer19UL18_JRV2_MC_SF_AK8PFPuppi.txt')

  #################
  #### BJetEnergyRegression
  #################
  bJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/breg_training_2018.pb")
  bJetNN.outputFormulas = cms.vstring(["at(0)*0.27912887930870056+1.0545977354049683","0.5*(at(2)-at(1))*0.27912887930870056"])
  cJetNN.weightFile =  cms.FileInPath("PhysicsTools/NanoAOD/data/creg_training_2018.pb")
  cJetNN.outputFormulas = cms.vstring(["at(0)*0.24325256049633026+0.993854820728302","0.5*(at(2)-at(1))*0.24325256049633026"])

  #################
  #### Update MET
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription
  #### This is independent of jecSequence, but it rather reapply JEC/JER using GT withing this MET corrector module
  #################
  from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
  runMetCorAndUncFromMiniAOD(process,
                           isData=(not isMC),
                           )

  from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
  makePuppiesFromMiniAOD( process, True );
  runMetCorAndUncFromMiniAOD(process,
                             isData=is_data,
                             metType="Puppi",
                             postfix="Puppi",
                             jetFlavor="AK4PFPuppi",
                             )

  #########################
  #### L1Prefire reweight
  #### https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe
  #########################
  from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
  process.prefiringweight = l1PrefiringWeightProducer.clone(
    TheJets = cms.InputTag("updatedPatJetsUpdatedJECslimmedJets"), #this should be the slimmedJets collection with up to date JECs !
    DataEraECAL = cms.string("None"),
    DataEraMuon = cms.string("20172018"),
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
  )

  ###########################
  #### Rochester correction
  ###########################
  process.recoTree.roccorPath = cms.string('SKFlatMaker/SKFlatMaker/data/roccor.Run2.v5/RoccoR2018UL.txt')

  ###########################
  #### BadPFMuonDz filter
  ###########################
  from RecoMET.METFilters.BadPFMuonDzFilter_cfi import BadPFMuonDzFilter
  process.BadPFMuonFilterUpdateDz=BadPFMuonDzFilter.clone(
    muons = cms.InputTag("slimmedMuons"),
    vtx   = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    minDzBestTrack = cms.double(0.5),
    taggingMode    = cms.bool(True)
  )

  ###########
  #### Path
  ###########

  process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.jecSequence *
    process.fullPatMetSequence *
    process.egmPhotonIDSequence *
    process.puppiMETSequence *
    process.fullPatMetSequencePuppi
  )
  process.p *= process.BadPFMuonFilterUpdateDz
  
  if isMC:
    process.recoTree.StoreL1PrefireFlag = cms.untracked.bool(True)
    process.p *= process.prefiringweight
    
    ## GenHFHadron 
    process.p *= process.matchGenSeq
 
  ## BJetEnergyRegression
  process.p *= process.bJetCorrSeq
  
  process.p *= process.recoTree
