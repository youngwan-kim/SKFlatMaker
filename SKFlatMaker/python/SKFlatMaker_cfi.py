import FWCore.ParameterSet.Config as cms

SKFlatMaker = cms.EDAnalyzer("SKFlatMaker",
  isMC = cms.untracked.bool(True),
  processName = cms.untracked.string("HLT"),
  DebugLevel = cms.untracked.int32(0),

  # -- Object Tags -- #
  Muon = cms.untracked.InputTag("slimmedMuons"),
  Electron = cms.untracked.InputTag("slimmedElectrons"),
  SmearedElectron = cms.untracked.InputTag("slimmedElectrons"),
  Photon = cms.untracked.InputTag("slimmedPhotons"),
  SmearedPhoton = cms.untracked.InputTag("slimmedPhotons"),
  Jet = cms.untracked.InputTag("slimmedJets"),
  FatJet = cms.untracked.InputTag("slimmedJetsAK8"),
  GenFatJet = cms.untracked.InputTag("slimmedGenJetsAK8"),
  MET = cms.InputTag("slimmedMETs"),
  LHEEventProduct = cms.untracked.InputTag("externalLHEProducer"),
  LHERunInfoProduct = cms.untracked.InputTag("externalLHEProducer"),
  GenParticle = cms.untracked.InputTag("genParticles"),

	# -- electron information -- #
	rho = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
        mvaIsoValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
        mvaNoIsoValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
	conversionsInputTag = cms.untracked.InputTag("allConversions"),
	GsfTrack = cms.untracked.InputTag("electronGsfTracks"),

	# -- photon information -- #
	effAreaChHadFile = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_TrueVtx.txt"),
	effAreaNeuHadFile= cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased_TrueVtx.txt"),
	effAreaPhoFile   = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_90percentBased_TrueVtx.txt"),

  # -- Jet information -- #
  BDiscriminant_tcheff = cms.untracked.double(0.7),
  BDiscriminant_tchpur = cms.untracked.double(0.7),
  BDiscriminant_ssv = cms.untracked.double(2.05),

  # -- MET information -- #
  METFilterResults_PAT = cms.InputTag("TriggerResults", "", "PAT"),
  METFilterResults_RECO = cms.InputTag("TriggerResults", "", "RECO"),
  pfMET = cms.untracked.InputTag("pfMet"),
     
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
  StoreElectronFlag = cms.untracked.bool(True),
  StoreCalibElectronFlag = cms.untracked.bool(True),
  StorePhotonFlag = cms.untracked.bool(True),
  StoreJetFlag = cms.untracked.bool(True),
  StoreMETFlag = cms.untracked.bool(True),
  StoreLHEFlag = cms.untracked.bool(True),
  StoreGENFlag = cms.untracked.bool(True),
  KeepAllGen = cms.untracked.bool(True),   
  StorePriVtxFlag = cms.untracked.bool(True),
  StoreTTFlag = cms.untracked.bool(False),
  StoreHLTReportFlag = cms.untracked.bool(True),

  # -- Filters -- #
  ApplyFilter = cms.untracked.bool(False),
  FilterType = cms.untracked.int32(0),
)
