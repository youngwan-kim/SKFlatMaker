import FWCore.ParameterSet.Config as cms

SKFlatMaker = cms.EDAnalyzer("SKFlatMaker",

  DataYear = cms.untracked.int32(-1),
  processName = cms.untracked.string("HLT"),
  DebugLevel = cms.untracked.int32(0),

  # -- Object Tags -- #
  Muon = cms.untracked.InputTag("slimmedMuons"),
  Tau = cms.untracked.InputTag("slimmedTaus"),
  Electron = cms.untracked.InputTag("slimmedElectrons"),
  Photon = cms.untracked.InputTag("slimmedPhotons"),
  Jet = cms.untracked.InputTag("slimmedJets"),
  GenJet = cms.untracked.InputTag("slimmedGenJets"),
  FatJet = cms.untracked.InputTag("slimmedJetsAK8"),
  GenFatJet = cms.untracked.InputTag("slimmedGenJetsAK8"),
  MET = cms.InputTag("slimmedMETs"),
  PuppiMET = cms.InputTag("slimmedMETsPuppi"),
  LHEEventProduct = cms.untracked.InputTag("externalLHEProducer"),
  LHERunInfoProduct = cms.untracked.InputTag("externalLHEProducer"),
  GenParticle = cms.untracked.InputTag("genParticles"),

  # -- electron information -- #
  rho = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
  rhoNC = cms.untracked.InputTag("fixedGridRhoFastjetCentralNeutral"),
  conversionsInputTag = cms.untracked.InputTag("allConversions"),
  GsfTrack = cms.untracked.InputTag("electronGsfTracks"),
  electron_EA_NHandPh_file = cms.untracked.FileInPath("RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt"),
  electron_IDtoSave = cms.untracked.vstring(
"cutBasedElectronID-Fall17-94X-V2-veto",
"cutBasedElectronID-Fall17-94X-V2-loose",
"cutBasedElectronID-Fall17-94X-V2-medium",
"cutBasedElectronID-Fall17-94X-V2-tight",
'mvaEleID-Fall17-iso-V2-wp80' ,
'mvaEleID-Fall17-iso-V2-wp90' ,
'mvaEleID-Fall17-iso-V2-wpHZZ' ,
'mvaEleID-Fall17-iso-V2-wpLoose',
'mvaEleID-Fall17-noIso-V2-wp80' ,
'mvaEleID-Fall17-noIso-V2-wp90' ,
'mvaEleID-Fall17-noIso-V2-wpLoose',
"heepElectronID-HEEPV70",
  ),

  #### Rochestor
  roccorPath = cms.string('SKFlatMaker/SKFlatMaker/data/roccor.Run2.v5/RoccoR2017UL.txt'),

  # -- photon information -- #
  photon_EA_CH_file = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt"),
  photon_EA_HN_file = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_V2.txt"),
  photon_EA_Ph_file = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased_V2.txt"),

  # -- Jet information -- #
  AK4Jet_payloadName = cms.string('AK4PFchs'),
  AK8Jet_payloadName = cms.string('AK8PFPuppi'),
  AK4Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.txt'),
  AK4Jet_JER_SF_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_SF_AK4PFchs.txt'),
  AK8Jet_JER_PtRes_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_PtResolution_AK8PFPuppi.txt'),
  AK8Jet_JER_SF_filepath = cms.string('SKFlatMaker/SKFlatMaker/data/JRDatabase/textFiles/Summer19UL17_JRV2_MC/Summer19UL17_JRV2_MC_SF_AK8PFPuppi.txt'),

  # -- MET information -- #
  METFilterResults_PAT = cms.InputTag("TriggerResults", "", "PAT"),
  METFilterResults_RECO = cms.InputTag("TriggerResults", "", "RECO"),
     
  # -- Trigger -- #
  TriggerResults = cms.untracked.InputTag("TriggerResults", "", "HLT"),
  TriggerResultsPAT = cms.untracked.InputTag("TriggerResults", "", "PAT"),
  ##TriggerObject = cms.untracked.InputTag("selectedPatTrigger"),
  TriggerObject = cms.untracked.InputTag("slimmedPatTrigger"),     
 
  # -- Else -- #
  GenEventInfo = cms.untracked.InputTag("generator"),
  BeamSpot = cms.untracked.InputTag("offlineBeamSpot"),
  PrimaryVertex = cms.untracked.InputTag("offlinePrimaryVerticesWithBS"),
  Track = cms.untracked.InputTag("generalTracks"),
  PileUpInfo = cms.untracked.InputTag("addPileupInfo"),

  # -- Store Flags -- #
  StoreMuonFlag = cms.untracked.bool(True),
  StoreTauFlag = cms.untracked.bool(True),
  StoreElectronFlag = cms.untracked.bool(True),
  StoreCalibElectronFlag = cms.untracked.bool(True),
  StorePhotonFlag = cms.untracked.bool(True),
  StoreJetFlag = cms.untracked.bool(True),
  StoreFatJetFlag = cms.untracked.bool(True),
  StoreMETFlag = cms.untracked.bool(True),
  StorePuppiMETFlag = cms.untracked.bool(True),
  StoreLHEFlag = cms.untracked.bool(True),
  StoreGENFlag = cms.untracked.bool(True),
  KeepAllGen = cms.untracked.bool(True),   
  StorePriVtxFlag = cms.untracked.bool(True),
  StoreHLTReportFlag = cms.untracked.bool(True),
  StoreHLTObjectFlag = cms.untracked.bool(False),
  StoreL1PrefireFlag = cms.untracked.bool(False),

  # -- Filters -- #
  ApplyFilter = cms.untracked.bool(False),
  FilterType = cms.untracked.int32(0),

  # -- GenHFHadronMatcher  -- #
  genJets = cms.InputTag("slimmedGenJets"),
  genBHadJetIndex = cms.InputTag("matchGenBHadron", "genBHadJetIndex"),
  genBHadFlavour = cms.InputTag("matchGenBHadron", "genBHadFlavour"),
  genBHadFromTopWeakDecay = cms.InputTag("matchGenBHadron", "genBHadFromTopWeakDecay"),
  genBHadPlusMothers = cms.InputTag("matchGenBHadron", "genBHadPlusMothers"),
  genBHadPlusMothersIndices = cms.InputTag("matchGenBHadron", "genBHadPlusMothersIndices"),
  genBHadIndex = cms.InputTag("matchGenBHadron", "genBHadIndex"),
  genBHadLeptonHadronIndex = cms.InputTag("matchGenBHadron", "genBHadLeptonHadronIndex"),
  genBHadLeptonViaTau = cms.InputTag("matchGenBHadron", "genBHadLeptonViaTau"),
  genCHadJetIndex = cms.InputTag("matchGenCHadron", "genCHadJetIndex"),
  genCHadFlavour = cms.InputTag("matchGenCHadron", "genCHadFlavour"),
  genCHadFromTopWeakDecay = cms.InputTag("matchGenCHadron", "genCHadFromTopWeakDecay"),
  genCHadBHadronId = cms.InputTag("matchGenCHadron", "genCHadBHadronId"),
  genCHadIndex = cms.InputTag("matchGenCHadron", "genCHadIndex"),

  # -- bJetEnergyCorrectionNN -- #
  bJetNNCorr = cms.InputTag("bJetNN", "bJetNNCorr"),
  bJetNNRes = cms.InputTag("bJetNN", "bJetNNRes"),
  cJetNNCorr = cms.InputTag("cJetNN", "cJetNNCorr"),
  cJetNNRes = cms.InputTag("cJetNN", "cJetNNRes"),
)
