//-------------------------------------------------
//
//   Class: SKFlatMaker
//
//   Description: Flat Ntuple maker for SNUCMS group data analysis framework
//
//
//   Authors of DYNtupleMaker (https://github.com/KyeongPil-Lee/NtupleMaker)
//   H.D.Yoo           Purdue University
//   K.P.Lee           Seoul National University
//   D.M.Pai           Seoul National University (EGM corrections, emu variables, and Gen-level neutrinos)
//
//   Revised By
//   S.B. Oh           Seoul National University
//   J.S. Kim          Seoul National University
//
//--------------------------------------------------

#include "SKFlatMaker/SKFlatMaker/interface/SKFlatMaker.h"

//
// class decleration
//

using namespace std;
using namespace reco;
using namespace edm;
using namespace pat;
using namespace isodeposit;


// -- Constructor -- //
SKFlatMaker::SKFlatMaker(const edm::ParameterSet& iConfig):
// -- object tokens -- //
MuonToken                           ( consumes< std::vector<pat::Muon> >                    (iConfig.getUntrackedParameter<edm::InputTag>("Muon")) ),
ElectronToken                       ( consumes< edm::View<pat::Electron> >                  (iConfig.getUntrackedParameter<edm::InputTag>("Electron")) ),
PhotonToken                         ( consumes< edm::View<pat::Photon> >                    (iConfig.getUntrackedParameter<edm::InputTag>("Photon")) ),
JetToken                            ( consumes< std::vector<pat::Jet> >                     (iConfig.getUntrackedParameter<edm::InputTag>("Jet")) ),
genJetToken                         ( consumes< reco::GenJetCollection >                    (iConfig.getUntrackedParameter<edm::InputTag>("GenJet")) ),
FatJetToken                         ( consumes< std::vector<pat::Jet> >                     (iConfig.getUntrackedParameter<edm::InputTag>("FatJet")) ),
genFatJetToken                      ( consumes< reco::GenJetCollection >                    (iConfig.getUntrackedParameter<edm::InputTag>("GenFatJet")) ),
MetToken                            ( consumes< std::vector<pat::MET> >                     (iConfig.getParameter<edm::InputTag>("MET")) ),

LHEEventProductToken                ( consumes< LHEEventProduct >                           (iConfig.getUntrackedParameter<edm::InputTag>("LHEEventProduct")) ),
LHEEventProductSourceToken          ( consumes< LHEEventProduct >                           (edm::InputTag("source")) ),
LHERunInfoProductToken              ( consumes< LHERunInfoProduct,edm::InRun >              (iConfig.getUntrackedParameter<edm::InputTag>("LHERunInfoProduct")) ),
LHERunInfoProductSourceToken        ( consumes< LHERunInfoProduct,edm::InRun >              (edm::InputTag("source")) ),
mcLabel_                            ( consumes< reco::GenParticleCollection>                (iConfig.getUntrackedParameter<edm::InputTag>("GenParticle"))  ),

// -- MET Filter tokens -- //
METFilterResultsToken_PAT           ( consumes<edm::TriggerResults>                         (iConfig.getParameter<edm::InputTag>("METFilterResults_PAT"))),
METFilterResultsToken_RECO          ( consumes<edm::TriggerResults>                         (iConfig.getParameter<edm::InputTag>("METFilterResults_RECO"))),

// -- Electron tokens -- //
RhoToken                            ( consumes< double >                                    (iConfig.getUntrackedParameter<edm::InputTag>("rho")) ),
ConversionsToken                    ( consumes< std::vector<reco::Conversion> >             (iConfig.getUntrackedParameter<edm::InputTag>("conversionsInputTag")) ),
GsfTrackToken                       ( consumes< std::vector< reco::GsfTrack > >             (iConfig.getUntrackedParameter<edm::InputTag>("GsfTrack")) ),

// -- Trigger Token -- //
TriggerToken                        ( consumes< edm::TriggerResults >                       (iConfig.getUntrackedParameter<edm::InputTag>("TriggerResults")) ),
TriggerTokenPAT                     ( consumes< edm::TriggerResults >                       (iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsPAT")) ),
TriggerObjectToken                  ( consumes< std::vector<pat::TriggerObjectStandAlone> > (iConfig.getUntrackedParameter<edm::InputTag>("TriggerObject")) ),
// -- Else -- //
GenEventInfoToken                   ( consumes< GenEventInfoProduct >                       (iConfig.getUntrackedParameter<edm::InputTag>("GenEventInfo")) ),
BeamSpotToken                       ( consumes< reco::BeamSpot >                            (iConfig.getUntrackedParameter<edm::InputTag>("BeamSpot")) ),
PrimaryVertexToken                  ( consumes< reco::VertexCollection >                    (iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex")) ),
TrackToken                          ( consumes< edm::View<reco::Track> >                    (iConfig.getUntrackedParameter<edm::InputTag>("Track")) ),
PileUpInfoToken                     ( consumes< std::vector< PileupSummaryInfo > >          (iConfig.getUntrackedParameter<edm::InputTag>("PileUpInfo")) ),
GenLumiInfoHeaderToken              ( consumes< GenLumiInfoHeader,edm::InLumi >             (edm::InputTag("generator")) )
{

  DataYear                          = iConfig.getUntrackedParameter<int>("DataYear");
  if(DataYear<0){
    cout << "DataYear is not set : DataYear = " << DataYear << endl;
    exit(EXIT_FAILURE);
  }
  else{
    cout << "[SKFlatMaker::SKFlatMaker] DataYear = " << DataYear << endl;
  }

  //==== Flag_BadPFMuonDzFilter
  BadPFMuonDzFilter_token = consumes< bool >(edm::InputTag("BadPFMuonFilterUpdateDz"));

  //==== L1 Prefireing
  prefweight_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProb"));
  prefweightup_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbUp"));
  prefweightdown_token = consumes< double >(edm::InputTag("prefiringweight:nonPrefiringProbDown"));

  theDebugLevel                     = iConfig.getUntrackedParameter<int>("DebugLevel", 0);

  if(theDebugLevel) cout << "[SKFlatMaker::SKFlatMaker] Constructor called" << endl;

  nEvt = 0;
  
  processName                       = iConfig.getUntrackedParameter<string>("processName", "HLT");

  electron_IDtoSave                 = iConfig.getUntrackedParameter< std::vector<std::string> >("electron_IDtoSave");
  cout << "[SKFlatMaker::SKFlatMaker]" << endl;
  cout << "############# electron_IDtoSave #############" << endl;
  for(unsigned int i=0; i<electron_IDtoSave.size(); i++){
    cout << electron_IDtoSave.at(i) << endl;
  }
  cout << "#############################################" << endl;
  electron_EA_NHandPh_file          = iConfig.getUntrackedParameter<edm::FileInPath>( "electron_EA_NHandPh_file" );
  photon_EA_CH_file                 = iConfig.getUntrackedParameter<edm::FileInPath>( "photon_EA_CH_file" );
  photon_EA_HN_file                 = iConfig.getUntrackedParameter<edm::FileInPath>( "photon_EA_HN_file" );
  photon_EA_Ph_file                 = iConfig.getUntrackedParameter<edm::FileInPath>( "photon_EA_Ph_file" );

  jet_payloadName_                  = iConfig.getParameter<std::string>("AK4Jet_payloadName");
  fatjet_payloadName_               = iConfig.getParameter<std::string>("AK8Jet_payloadName");

  jet_jecUnc = NULL;
  fatjet_jecUnc = NULL;
  
  // -- Store Flags -- //
  theStorePriVtxFlag                = iConfig.getUntrackedParameter<bool>("StorePriVtxFlag", true);
  theStoreJetFlag                   = iConfig.getUntrackedParameter<bool>("StoreJetFlag", true);
  theStoreFatJetFlag                = iConfig.getUntrackedParameter<bool>("StoreFatJetFlag", true);
  theStoreMETFlag                   = iConfig.getUntrackedParameter<bool>("StoreMETFlag", true);
  theStoreHLTReportFlag             = iConfig.getUntrackedParameter<bool>("StoreHLTReportFlag", true);
  theStoreHLTObjectFlag             = iConfig.getUntrackedParameter<bool>("StoreHLTObjectFlag", true);
  theStoreMuonFlag                  = iConfig.getUntrackedParameter<bool>("StoreMuonFlag", true);
  theStoreElectronFlag              = iConfig.getUntrackedParameter<bool>("StoreElectronFlag", true);
  theStoreLHEFlag                   = iConfig.getUntrackedParameter<bool>("StoreLHEFlag", false);
  theStoreGENFlag                   = iConfig.getUntrackedParameter<bool>("StoreGENFlag", true);
  theKeepAllGen                     = iConfig.getUntrackedParameter<bool>("KeepAllGen", true);
  theStorePhotonFlag                = iConfig.getUntrackedParameter<bool>("StorePhotonFlag", true);
  theStoreL1PrefireFlag             = iConfig.getUntrackedParameter<bool>("StoreL1PrefireFlag",true);

  cout << "[SKFlatMaker::SKFlatMaker] Rochester correction file from : " << edm::FileInPath( iConfig.getParameter<std::string>("roccorPath") ).fullPath() << endl;
  rc.init(edm::FileInPath( iConfig.getParameter<std::string>("roccorPath") ).fullPath());

  //==================
  //==== Prepare JER
  //==================

  //==== this method uses Spring16_25nsV6 (80X, 2016, ICHEP dataset) DATA/MC SFs ...
  //jet_resolution = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
  //jet_resolution_sf = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");

  jet_resolution = JME::JetResolution( edm::FileInPath( iConfig.getParameter<std::string>("AK4Jet_JER_PtRes_filepath") ).fullPath()  );
  jet_resolution_sf = JME::JetResolutionScaleFactor( edm::FileInPath( iConfig.getParameter<std::string>("AK4Jet_JER_SF_filepath") ).fullPath() );
  fatjet_resolution = JME::JetResolution( edm::FileInPath( iConfig.getParameter<std::string>("AK8Jet_JER_PtRes_filepath") ).fullPath()  );
  fatjet_resolution_sf = JME::JetResolutionScaleFactor( edm::FileInPath( iConfig.getParameter<std::string>("AK8Jet_JER_SF_filepath") ).fullPath() );

  //==== Weights
  WeightMap["Scale"]={};
  WeightMap["PDF"]={};
  WeightMap["AlphaS"]={};
  WeightMap["PSSyst"]={};
  PDFAlphaSScaleValue_={1.};
  for(auto pname:iConfig.getParameterNames()){
    TString wname=pname;
    if(wname.Contains("weight_")){
      wname.ReplaceAll("weight_","");
      if(wname=="AlphaSScale"){
	PDFAlphaSScaleValue_ = iConfig.getUntrackedParameter< std::vector<double> >(pname);
      }else{
	WeightMap[wname]=iConfig.getUntrackedParameter<std::vector<int>>(pname);
	cout << "[SKFlatMaker::SKFlatMaker] WeightMap["+wname+"] = " << WeightMap[wname].at(0)<<" - "<<WeightMap[wname].back() << endl;
	weight_[wname]={};
      }
    }
  }
  if(WeightMap["AlphaS"].size()){
    if(PDFAlphaSScaleValue_.size()<WeightMap["AlphaS"].size()){
      PDFAlphaSScaleValue_.resize(WeightMap["AlphaS"].size(),PDFAlphaSScaleValue_.back());
    }
    cout << "[SKFlatMaker::SKFlatMaker] PDFAlphaSScaleValue = " << PDFAlphaSScaleValue_.at(0)<<" - "<<PDFAlphaSScaleValue_.back() << endl;
  }

  if(theDebugLevel) cout << "[SKFlatMaker::SKFlatMaker] Constructor finished" << endl;

}

SKFlatMaker::~SKFlatMaker() { }

//
// member functions
//

// ------------ method called to for each event  ------------ //
void SKFlatMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] called" << endl;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
  
  ///////////////////////////////////////////
  // -- initialize for ntuple variables -- //
  ///////////////////////////////////////////
  
  Flag_goodVertices = false;
  Flag_globalSuperTightHalo2016Filter = false;
  Flag_HBHENoiseFilter = false;
  Flag_HBHENoiseIsoFilter = false;
  Flag_EcalDeadCellTriggerPrimitiveFilter = false;
  Flag_BadPFMuonFilter = false;
  Flag_BadPFMuonDzFilter = false;
  //Flag_hfNoisyHitsFilter = false; //in MiniAODv2
  Flag_BadChargedCandidateFilter = false;
  Flag_eeBadScFilter = false;
  Flag_ecalBadCalibFilter = false;

  nPV = -1;
  PVtrackSize = -1;
  PVchi2 = -1;
  PVndof = -1;
  PVnormalizedChi2 = -1;
  PVx = -1000;
  PVy = -1000;
  PVz = -1000;
  PVprob = -1;
  
  // -- PF iso deposits -- // 
  sumEt = 0;
  photonEt = 0;
  chargedHadronEt = 0;
  neutralHadronEt = 0;
  Rho = 0;

  //==== MET

  pfMET_pt=-999;
  pfMET_phi=-999;
  pfMET_SumEt=-999;
  pfMET_Type1_pt=-999;
  pfMET_Type1_phi=-999;
  pfMET_Type1_SumEt=-999;
  pfMET_Type1_PhiCor_pt=-999;
  pfMET_Type1_PhiCor_phi=-999;
  pfMET_Type1_PhiCor_SumEt=-999;
  pfMET_pt_shifts.clear();
  pfMET_phi_shifts.clear();
  pfMET_SumEt_shifts.clear();
  pfMET_Type1_pt_shifts.clear();
  pfMET_Type1_phi_shifts.clear();
  pfMET_Type1_SumEt_shifts.clear();
  pfMET_Type1_PhiCor_pt_shifts.clear();
  pfMET_Type1_PhiCor_phi_shifts.clear();
  pfMET_Type1_PhiCor_SumEt_shifts.clear();

  //==== Trigger (object)

  HLT_TriggerName.clear();
  HLT_TriggerFilterName.clear();
  HLTObject_pt.clear();
  HLTObject_eta.clear();
  HLTObject_phi.clear();
  HLTObject_FiredFilters.clear();
  HLTObject_FiredPaths.clear();

  //==== LHE

  LHE_Px.clear();
  LHE_Py.clear();
  LHE_Pz.clear();
  LHE_E.clear();
  LHE_Status.clear();
  LHE_ID.clear();
  for(auto& it:weight_){
    weight_[it.first].clear();
  }

  //==== GEN

  gen_phi.clear();
  gen_eta.clear();
  gen_pt.clear();
  gen_mass.clear();
  gen_charge.clear();
  gen_mother_index.clear();
  gen_status.clear();
  gen_PID.clear();
  gen_isPrompt.clear();
  gen_isPromptFinalState.clear();
  gen_isTauDecayProduct.clear();
  gen_isPromptTauDecayProduct.clear();
  gen_isDirectPromptTauDecayProductFinalState.clear();
  gen_isHardProcess.clear();
  gen_isLastCopy.clear();
  gen_isLastCopyBeforeFSR.clear();
  gen_isPromptDecayed.clear();
  gen_isDecayedLeptonHadron.clear();
  gen_fromHardProcessBeforeFSR.clear();
  gen_fromHardProcessDecayed.clear();
  gen_fromHardProcessFinalState.clear();
  gen_isMostlyLikePythia6Status3.clear();
  gen_weight=-999;
  genWeight_Q=-999;
  genWeight_X1=-999;
  genWeight_X2=-999;
  genWeight_id1=-999;
  genWeight_id2=-999;
  genWeight_alphaQCD=-999;
  genWeight_alphaQED=-999;

  //==== Electron
  electron_MVAIso.clear();
  electron_MVANoIso.clear();
  electron_Energy.clear();
  electron_Energy_Scale_Up.clear();
  electron_Energy_Scale_Down.clear();
  electron_Energy_Smear_Up.clear();
  electron_Energy_Smear_Down.clear();
  electron_eta.clear();
  electron_phi.clear();
  electron_charge.clear();
  electron_gsfpt.clear();
  electron_gsfEta.clear();
  electron_gsfPhi.clear();
  electron_gsfCharge.clear();
  electron_scEta.clear();
  electron_scPhi.clear();
  electron_etaWidth.clear();
  electron_phiWidth.clear();
  electron_dEtaIn.clear();
  electron_dEtaInSeed.clear();
  electron_dPhiIn.clear();
  electron_sigmaIEtaIEta.clear();
  electron_Full5x5_SigmaIEtaIEta.clear();
  electron_e2x5OverE5x5.clear();
  electron_e1x5OverE5x5.clear();
  electron_HoverE.clear();
  electron_fbrem.clear();
  electron_eOverP.clear();
  electron_InvEminusInvP.clear();
  electron_dxyVTX.clear();
  electron_dxyerrVTX.clear();
  electron_dzVTX.clear();
  electron_dzerrVTX.clear();
  electron_3DIPVTX.clear();
  electron_3DIPerrVTX.clear();
  electron_dxy.clear();
  electron_sigdxy.clear();
  electron_dz.clear();
  electron_dxyBS.clear();
  electron_dzBS.clear();
  electron_chIso03.clear();
  electron_nhIso03.clear();
  electron_phIso03.clear();
  electron_puChIso03.clear();
  electron_passConversionVeto.clear();
  electron_isGsfCtfScPixChargeConsistent.clear();
  electron_isGsfScPixChargeConsistent.clear();
  electron_isGsfCtfChargeConsistent.clear();
  electron_mHits.clear();
  electron_ecalDriven.clear();
  electron_r9.clear();
  electron_scEnergy.clear();
  electron_scPreEnergy.clear();
  electron_scRawEnergy.clear();
  electron_scEt.clear();
  electron_E15.clear();
  electron_E25.clear();
  electron_E55.clear();
  electron_RelPFIso_dBeta.clear();
  electron_RelPFIso_Rho.clear();
  electron_IDBit.clear();
  electron_IDCutBit.clear();
  electron_EnergyUnCorr.clear();
  electron_chMiniIso.clear();
  electron_nhMiniIso.clear();
  electron_phMiniIso.clear();
  electron_puChMiniIso.clear();
  electron_trackIso.clear();
  electron_dr03EcalRecHitSumEt.clear();
  electron_dr03HcalDepth1TowerSumEt.clear();
  electron_dr03HcalTowerSumEt.clear();
  electron_dr03TkSumPt.clear();
  electron_ecalPFClusterIso.clear();
  electron_hcalPFClusterIso.clear();
  electron_pathbits.clear();
  electron_filterbits.clear();

  //==== Muon
  muon_PfChargedHadronIsoR04.clear();
  muon_PfNeutralHadronIsoR04.clear();
  muon_PfGammaIsoR04.clear();
  muon_PFSumPUIsoR04.clear();
  muon_PfChargedHadronIsoR03.clear();
  muon_PfNeutralHadronIsoR03.clear();
  muon_PfGammaIsoR03.clear();
  muon_PFSumPUIsoR03.clear();
  muon_TypeBit.clear();
  muon_IDBit.clear();
  muon_ishighpt.clear();
  muon_dB.clear();
  muon_phi.clear();
  muon_eta.clear();
  muon_pt.clear();
  muon_mass.clear();
  muon_trkiso.clear();
  muon_hcaliso.clear();
  muon_ecaliso.clear();
  muon_trkisoR05.clear();
  muon_hcalisoR05.clear();
  muon_ecalisoR05.clear();
  muon_charge.clear();
  muon_nChambers.clear();
  muon_matchedstations.clear();
  muon_stationMask.clear();
  muon_nSegments.clear();
  muon_normchi.clear();
  muon_validhits.clear();
  muon_trackerHits.clear();
  muon_pixelHits.clear();
  muon_validmuonhits.clear();
  muon_trackerLayers.clear();
  muon_qoverp.clear();
  muon_theta.clear();
  muon_lambda.clear();
  muon_dxy.clear();
  muon_d0.clear();
  muon_dsz.clear();
  muon_dz.clear();
  muon_dxyBS.clear();
  muon_dzBS.clear();
  muon_dszBS.clear();
  muon_dxyVTX.clear();
  muon_dxyerrVTX.clear();
  muon_dzVTX.clear();
  muon_dzerrVTX.clear();
  muon_3DIPVTX.clear();
  muon_3DIPerrVTX.clear();
  muon_dszVTX.clear();
  muon_vx.clear();
  muon_vy.clear();
  muon_vz.clear();
  muon_Best_pt.clear();
  muon_Best_ptError.clear();
  muon_Best_eta.clear();
  muon_Best_phi.clear();
  muon_Inner_pt.clear();
  muon_Inner_ptError.clear();
  muon_Inner_eta.clear();
  muon_Inner_phi.clear();
  muon_Outer_pt.clear();
  muon_Outer_ptError.clear();
  muon_Outer_eta.clear();
  muon_Outer_phi.clear();
  muon_GLB_pt.clear();
  muon_GLB_ptError.clear();
  muon_GLB_eta.clear();
  muon_GLB_phi.clear();
  muon_TuneP_pt.clear();
  muon_TuneP_ptError.clear();
  muon_TuneP_eta.clear();
  muon_TuneP_phi.clear();
  muon_TuneP_charge.clear();
  muon_roch_sf.clear();
  muon_roch_sf_up.clear();
  muon_PfChargedHadronMiniIso.clear();
  muon_PfNeutralHadronMiniIso.clear();
  muon_PfGammaMiniIso.clear();
  muon_PFSumPUMiniIso.clear();
  muon_MVA.clear();
  muon_lowptMVA.clear();
  muon_softMVA.clear();
  muon_jetPtRatio.clear();
  muon_jetPtRel.clear();
  muon_chi2LocalPosition.clear();
  muon_trkKink.clear();
  muon_segmentCompatibility.clear();
  muon_trackerValidFraction.clear();
  muon_simType.clear();
  muon_simExtType.clear();
  muon_simFlavour.clear();
  muon_simHeaviestMotherFlavour.clear();
  muon_simPdgId.clear();
  muon_simMotherPdgId.clear();
  muon_simMatchQuality.clear();
  muon_pathbits.clear();
  muon_filterbits.clear();

  //==== Jet
  jet_pt.clear();
  jet_eta.clear();
  jet_phi.clear();
  jet_charge.clear();
  jet_area.clear();
  jet_partonFlavour.clear();
  jet_hadronFlavour.clear();
  jet_CSVv2.clear();
  jet_DeepCSV.clear();
  jet_CvsL.clear();
  jet_CvsB.clear();
  jet_DeepFlavour_b.clear();
  jet_DeepFlavour_bb.clear();
  jet_DeepFlavour_lepb.clear();
  jet_DeepFlavour_c.clear();
  jet_DeepFlavour_uds.clear();
  jet_DeepFlavour_g.clear();
  jet_DeepCvsL.clear();
  jet_DeepCvsB.clear();
  jet_chargedHadronEnergyFraction.clear();
  jet_neutralHadronEnergyFraction.clear();
  jet_neutralEmEnergyFraction.clear();
  jet_chargedEmEnergyFraction.clear();
  jet_muonEnergyFraction.clear();
  jet_chargedMultiplicity.clear();
  jet_neutralMultiplicity.clear();
  jet_tightJetID.clear();
  jet_tightLepVetoJetID.clear();
  jet_partonPdgId.clear();
  jet_vtxNtracks.clear();
  jet_m.clear();
  jet_energy.clear();
  jet_PileupJetId.clear();
  jet_shiftedEnUp.clear();
  jet_shiftedEnDown.clear();
  jet_smearedRes.clear();
  jet_smearedResUp.clear();
  jet_smearedResDown.clear();
  jet_JECL1FastJet.clear();
  jet_JECFull.clear();

  //==== FatJet
  fatjet_pt.clear();
  fatjet_eta.clear();
  fatjet_phi.clear();
  fatjet_charge.clear();
  fatjet_area.clear();
  fatjet_partonFlavour.clear();
  fatjet_hadronFlavour.clear();
  fatjet_CSVv2.clear();
  fatjet_DeepCSV.clear();
  fatjet_CvsL.clear();
  fatjet_CvsB.clear();
  fatjet_DeepFlavour_b.clear();
  fatjet_DeepFlavour_bb.clear();
  fatjet_DeepFlavour_lepb.clear();
  fatjet_DeepFlavour_c.clear();
  fatjet_DeepFlavour_uds.clear();
  fatjet_DeepFlavour_g.clear();
  fatjet_DeepCvsL.clear();
  fatjet_DeepCvsB.clear();
  fatjet_tightJetID.clear();
  fatjet_tightLepVetoJetID.clear();
  fatjet_partonPdgId.clear();
  fatjet_vtxNtracks.clear();
  fatjet_m.clear();
  fatjet_energy.clear();
  fatjet_puppi_tau1.clear();
  fatjet_puppi_tau2.clear();
  fatjet_puppi_tau3.clear();
  fatjet_puppi_tau4.clear();
  fatjet_softdropmass.clear();
  fatjet_chargedHadronEnergyFraction.clear();
  fatjet_neutralHadronEnergyFraction.clear();
  fatjet_neutralEmEnergyFraction.clear();
  fatjet_chargedEmEnergyFraction.clear();
  fatjet_muonEnergyFraction.clear();
  fatjet_chargedMultiplicity.clear();
  fatjet_neutralMultiplicity.clear();
  fatjet_shiftedEnUp.clear();
  fatjet_shiftedEnDown.clear();
  fatjet_smearedRes.clear();
  fatjet_smearedResUp.clear();
  fatjet_smearedResDown.clear();
  fatjet_LSF.clear();
  fatjet_LSFlep_PID.clear();
  fatjet_LSFlep_Pt.clear();
  fatjet_LSFlep_Eta.clear();
  fatjet_LSFlep_Phi.clear();

  //==== Photon
  photon_Energy.clear();
  photon_EnergyUnCorr.clear();
  photon_eta.clear();
  photon_phi.clear();
  photon_scEta.clear();
  photon_scPhi.clear();
  photon_HoverE.clear();
  photon_hasPixelSeed.clear();
  photon_Full5x5_SigmaIEtaIEta.clear();
  photon_ChIso.clear();
  photon_NhIso.clear();
  photon_PhIso.clear();
  photon_ChIsoWithEA.clear();
  photon_NhIsoWithEA.clear();
  photon_PhIsoWithEA.clear();
  photon_passMVAID_WP80.clear();
  photon_passMVAID_WP90.clear();
  photon_passLooseID.clear();
  photon_passMediumID.clear();
  photon_passTightID.clear();

  // cout << "[SKFlatMaker::analyze] Varialbe intilization done" << endl;
  
  nEvt++;
  // -- run number & event number -- //
  runNum = iEvent.id().run();
  evtNum = iEvent.id().event();
  lumiBlock = iEvent.id().luminosityBlock();
  
  //get the geometry
  edm::ESHandle<GlobalTrackingGeometry> glbTrackingGeometry;
  iSetup.get<GlobalTrackingGeometryRecord>().get(glbTrackingGeometry);
  
  // -- PileUp Reweighting -- //
  IsData = iEvent.isRealData();
  if( !IsData ){
    edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByToken(PileUpInfoToken, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;

    int npv = -1;
    //int npvin = -1;
      
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI){
      int BX = PVI->getBunchCrossing();
    
      if( BX == 0 ){
        // npvin = PVI->getPU_NumInteractions(); // in time only
        npv = PVI->getTrueNumInteractions(); // in and out of time
        continue;
      }
    
    }
      
    nPileUp = npv;
      
  }

  //====================
  //==== L1 Prefireing
  //====================

  L1PrefireReweight_Central = -999;
  L1PrefireReweight_Up = -999;
  L1PrefireReweight_Down = -999;
  if(theStoreL1PrefireFlag){

    edm::Handle< double > theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight ) ;
    L1PrefireReweight_Central =(*theprefweight);

    edm::Handle< double > theprefweightup;
    iEvent.getByToken(prefweightup_token, theprefweightup ) ;
    L1PrefireReweight_Up =(*theprefweightup);

    edm::Handle< double > theprefweightdown;
    iEvent.getByToken(prefweightdown_token, theprefweightdown ) ;
    L1PrefireReweight_Down =(*theprefweightdown);

  }

  //==================
  //==== Prepare JEC
  //==================

  //==== 1) AK4

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(jet_payloadName_,JetCorParColl);
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  if(jet_jecUnc) delete jet_jecUnc;
  jet_jecUnc = new JetCorrectionUncertainty(JetCorPar);

  //==== For cross check
  //==== https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
  //if(jet_jecUnc_methodB) delete jet_jecUnc_methodB;
  //jet_jecUnc_methodB = new JetCorrectionUncertainty("/u/user/jskim/scratch/CMSSW_9_4_2/src/SKFlatMaker/SKFlatMaker/data/JEC/MC/Fall17_17Nov2017_V6_MC/Fall17_17Nov2017_V6_MC_Uncertainty_AK4PFchs.txt");

  //==== 2) AK8

  edm::ESHandle<JetCorrectorParametersCollection> FatJetCorParColl;
  iSetup.get<JetCorrectionsRecord>().get(fatjet_payloadName_,FatJetCorParColl);
  JetCorrectorParameters const & FatJetCorPar = (*FatJetCorParColl)["Uncertainty"];
  if(fatjet_jecUnc) delete fatjet_jecUnc;
  fatjet_jecUnc = new JetCorrectionUncertainty(FatJetCorPar);



  //==== For JER, get genJets
  iEvent.getByToken(genJetToken, m_genJets);
  iEvent.getByToken(genFatJetToken, m_genFatJets);

  //==== Event varialbes
  edm::Handle< double > rhoHandle;
  iEvent.getByToken(RhoToken,rhoHandle);
  Rho = *rhoHandle;

  // fills
  
  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreHLTReportFlag" << endl;
  if( theStoreHLTReportFlag ) hltReport(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStorePriVtxFlag" << endl;
  if( theStorePriVtxFlag ) fillPrimaryVertex(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreJetFlag" << endl;
  if( theStoreJetFlag ) fillJet(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreFatJetFlag" << endl;
  if( theStoreFatJetFlag ) fillFatJet(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreMETFlag" << endl;
  if( theStoreMETFlag ) fillMET(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreLHEFlag" << endl;
  //if( !IsData && theStoreLHEFlag ) fillLHEInfo(iEvent);
  if(theStoreLHEFlag ) fillLHEInfo(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreGENFlag" << endl;
  //if( !IsData && theStoreGENFlag ) fillGENInfo(iEvent);
  if(theStoreGENFlag ) fillGENInfo(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStorePhotonFlag" << endl;
  if( theStorePhotonFlag ) fillPhotons(iEvent);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreMuonFlag" << endl;
  if( theStoreMuonFlag ) fillMuons(iEvent, iSetup);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] theStoreElectronFlag" << endl;
  if( theStoreElectronFlag ) fillElectrons(iEvent, iSetup);

  if(theDebugLevel) cout << "[SKFlatMaker::analyze] Tree Fill" << endl;
  DYTree->Fill();
  if(theDebugLevel) cout << "[SKFlatMaker::analyze] Tree Fill finished" << endl;

}

// ------------ method called once each job just before starting event loop  ------------ //
void SKFlatMaker::beginJob()
{

  if(theDebugLevel) cout << "[SKFlatMaker::beginJob] called" << endl;

  edm::Service<TFileService> fs;
  fs->file().SetCompressionLevel(9);
  DYTree = fs->make<TTree>("SKFlat","SKFlat");

  // -- global event variables -- //
  DYTree->Branch("IsData",&IsData,"IsData/O");
  DYTree->Branch("nTotal",&nEvt,"nTotal/I");
  DYTree->Branch("run",&runNum,"runNum/I");
  DYTree->Branch("event",&evtNum,"evtNum/l");
  DYTree->Branch("lumi",&lumiBlock,"lumiBlock/I");
  DYTree->Branch("Rho",&Rho,"Rho/F");
  DYTree->Branch("nPV",&nPV,"nPV/I");

  if(theStoreL1PrefireFlag){
    DYTree->Branch("L1PrefireReweight_Central",& L1PrefireReweight_Central,"L1PrefireReweight_Central/F");
    DYTree->Branch("L1PrefireReweight_Up",& L1PrefireReweight_Up,"L1PrefireReweight_Up/F");
    DYTree->Branch("L1PrefireReweight_Down",& L1PrefireReweight_Down,"L1PrefireReweight_Down/F");
  }

  //MET Filters 2017
  DYTree->Branch("Flag_goodVertices",&Flag_goodVertices,"Flag_goodVertices/O");
  DYTree->Branch("Flag_globalSuperTightHalo2016Filter",&Flag_globalSuperTightHalo2016Filter,"Flag_globalSuperTightHalo2016Filter/O");
  DYTree->Branch("Flag_HBHENoiseFilter",&Flag_HBHENoiseFilter,"Flag_HBHENoiseFilter/O");
  DYTree->Branch("Flag_HBHENoiseIsoFilter",&Flag_HBHENoiseIsoFilter,"Flag_HBHENoiseIsoFilter/O");
  DYTree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter",&Flag_EcalDeadCellTriggerPrimitiveFilter,"Flag_EcalDeadCellTriggerPrimitiveFilter/O");
  DYTree->Branch("Flag_BadPFMuonFilter",&Flag_BadPFMuonFilter,"Flag_BadPFMuonFilter/O");
  DYTree->Branch("Flag_BadPFMuonDzFilter",&Flag_BadPFMuonDzFilter,"Flag_BadPFMuonDzFilter/O");
  //DYTree->Branch("Flag_hfNoisyHitsFilter",&Flag_hfNoisyHitsFilter,"Flag_hfNoisyHitsFilter/O");
  DYTree->Branch("Flag_BadChargedCandidateFilter",&Flag_BadChargedCandidateFilter,"Flag_BadChargedCandidateFilter/O");
  DYTree->Branch("Flag_eeBadScFilter",&Flag_eeBadScFilter,"Flag_eeBadScFilter/O");
  DYTree->Branch("Flag_ecalBadCalibFilter",&Flag_ecalBadCalibFilter,"Flag_ecalBadCalibFilter/O");

  
  
  if(theStorePriVtxFlag){
    DYTree->Branch("PVtrackSize", &PVtrackSize,"PVtrackSize/I");
    DYTree->Branch("PVchi2", &PVchi2,"PVchi2/F");
    DYTree->Branch("PVndof", &PVndof,"PVndof/F");
    DYTree->Branch("PVnormalizedChi2", &PVnormalizedChi2,"PVnormalizedChi2/F");
    DYTree->Branch("vertex_X", &PVx,"PVx/F");
    DYTree->Branch("vertex_Y", &PVy,"PVy/F");
    DYTree->Branch("vertex_Z", &PVz,"PVz/F");
  }

  if(theStoreHLTReportFlag){

    DYTree->Branch("HLT_TriggerName", "vector<string>", &HLT_TriggerName);

    if(theStoreHLTObjectFlag){
      DYTree->Branch("HLT_TriggerFilterName", "vector<string>", &HLT_TriggerFilterName);
      DYTree->Branch("HLTObject_pt", "vector<float>", &HLTObject_pt);
      DYTree->Branch("HLTObject_eta", "vector<float>", &HLTObject_eta);
      DYTree->Branch("HLTObject_phi", "vector<float>", &HLTObject_phi);
      DYTree->Branch("HLTObject_FiredFilters", "vector<string>", &HLTObject_FiredFilters);
      DYTree->Branch("HLTObject_FiredPaths", "vector<string>", &HLTObject_FiredPaths);
    }

  }

  if(theStoreJetFlag){

    DYTree->Branch("jet_pt", "vector<float>", &jet_pt);
    DYTree->Branch("jet_eta", "vector<float>", &jet_eta);
    DYTree->Branch("jet_phi", "vector<float>", &jet_phi);
    DYTree->Branch("jet_charge", "vector<float>", &jet_charge);
    DYTree->Branch("jet_area", "vector<float>", &jet_area);
    DYTree->Branch("jet_partonFlavour", "vector<int>", &jet_partonFlavour);
    DYTree->Branch("jet_hadronFlavour", "vector<int>", &jet_hadronFlavour);
    DYTree->Branch("jet_CSVv2", "vector<float>", &jet_CSVv2);
    DYTree->Branch("jet_DeepCSV", "vector<float>", &jet_DeepCSV);
    DYTree->Branch("jet_CvsL", "vector<float>", &jet_CvsL);
    DYTree->Branch("jet_CvsB", "vector<float>", &jet_CvsB);
    DYTree->Branch("jet_DeepFlavour_b", "vector<float>", &jet_DeepFlavour_b);
    DYTree->Branch("jet_DeepFlavour_bb", "vector<float>", &jet_DeepFlavour_bb);
    DYTree->Branch("jet_DeepFlavour_lepb", "vector<float>", &jet_DeepFlavour_lepb);
    DYTree->Branch("jet_DeepFlavour_c", "vector<float>", &jet_DeepFlavour_c);
    DYTree->Branch("jet_DeepFlavour_uds", "vector<float>", &jet_DeepFlavour_uds);
    DYTree->Branch("jet_DeepFlavour_g", "vector<float>", &jet_DeepFlavour_g);
    DYTree->Branch("jet_DeepCvsL", "vector<float>", &jet_DeepCvsL);
    DYTree->Branch("jet_DeepCvsB", "vector<float>", &jet_DeepCvsB);
    DYTree->Branch("jet_chargedHadronEnergyFraction", "vector<float>", &jet_chargedHadronEnergyFraction);
    DYTree->Branch("jet_neutralHadronEnergyFraction", "vector<float>", &jet_neutralHadronEnergyFraction);
    DYTree->Branch("jet_neutralEmEnergyFraction", "vector<float>", &jet_neutralEmEnergyFraction);
    DYTree->Branch("jet_chargedEmEnergyFraction", "vector<float>", &jet_chargedEmEnergyFraction);
    DYTree->Branch("jet_muonEnergyFraction", "vector<float>", &jet_muonEnergyFraction);
    DYTree->Branch("jet_chargedMultiplicity", "vector<int>", &jet_chargedMultiplicity);
    DYTree->Branch("jet_neutralMultiplicity", "vector<int>", &jet_neutralMultiplicity);
    DYTree->Branch("jet_tightJetID", "vector<bool>", &jet_tightJetID);
    DYTree->Branch("jet_tightLepVetoJetID", "vector<bool>", &jet_tightLepVetoJetID);
    DYTree->Branch("jet_partonPdgId", "vector<int>", &jet_partonPdgId);
    DYTree->Branch("jet_vtxNtracks", "vector<int>", &jet_vtxNtracks);
    DYTree->Branch("jet_m", "vector<float>", &jet_m);
    DYTree->Branch("jet_energy", "vector<float>", &jet_energy);
    DYTree->Branch("jet_PileupJetId", "vector<float>", &jet_PileupJetId);
    DYTree->Branch("jet_shiftedEnUp", "vector<float>", &jet_shiftedEnUp);
    DYTree->Branch("jet_shiftedEnDown", "vector<float>", &jet_shiftedEnDown);
    DYTree->Branch("jet_smearedRes", "vector<float>", &jet_smearedRes);
    DYTree->Branch("jet_smearedResUp", "vector<float>", &jet_smearedResUp);
    DYTree->Branch("jet_smearedResDown", "vector<float>", &jet_smearedResDown);
    DYTree->Branch("jet_JECL1FastJet", "vector<float>", &jet_JECL1FastJet);
    DYTree->Branch("jet_JECFull", "vector<float>", &jet_JECFull);

  }
  
  if(theStoreFatJetFlag){
    DYTree->Branch("fatjet_pt", "vector<float>", &fatjet_pt);
    DYTree->Branch("fatjet_eta", "vector<float>", &fatjet_eta);
    DYTree->Branch("fatjet_phi", "vector<float>", &fatjet_phi);
    DYTree->Branch("fatjet_charge", "vector<float>", &fatjet_charge);
    DYTree->Branch("fatjet_area", "vector<float>", &fatjet_area);
    DYTree->Branch("fatjet_partonFlavour", "vector<int>", &fatjet_partonFlavour);
    DYTree->Branch("fatjet_hadronFlavour", "vector<int>", &fatjet_hadronFlavour);
    DYTree->Branch("fatjet_CSVv2", "vector<float>", &fatjet_CSVv2);
    DYTree->Branch("fatjet_DeepCSV", "vector<float>", &fatjet_DeepCSV);
    DYTree->Branch("fatjet_DeepFlavour_b", "vector<float>", &fatjet_DeepFlavour_b);
    DYTree->Branch("fatjet_DeepFlavour_bb", "vector<float>", &fatjet_DeepFlavour_bb);
    DYTree->Branch("fatjet_DeepFlavour_lepb", "vector<float>", &fatjet_DeepFlavour_lepb);
    DYTree->Branch("fatjet_DeepFlavour_c", "vector<float>", &fatjet_DeepFlavour_c);
    DYTree->Branch("fatjet_DeepFlavour_uds", "vector<float>", &fatjet_DeepFlavour_uds);
    DYTree->Branch("fatjet_DeepFlavour_g", "vector<float>", &fatjet_DeepFlavour_g);
    DYTree->Branch("fatjet_CvsL", "vector<float>", &fatjet_CvsL);
    DYTree->Branch("fatjet_CvsB", "vector<float>", &fatjet_CvsB);
    DYTree->Branch("fatjet_DeepCvsL", "vector<float>", &fatjet_DeepCvsL);
    DYTree->Branch("fatjet_DeepCvsB", "vector<float>", &fatjet_DeepCvsB);
    DYTree->Branch("fatjet_tightJetID", "vector<bool>", &fatjet_tightJetID);
    DYTree->Branch("fatjet_tightLepVetoJetID", "vector<bool>", &fatjet_tightLepVetoJetID);
    DYTree->Branch("fatjet_partonPdgId", "vector<int>", &fatjet_partonPdgId);
    DYTree->Branch("fatjet_vtxNtracks", "vector<int>", &fatjet_vtxNtracks);
    DYTree->Branch("fatjet_m", "vector<float>", &fatjet_m);
    DYTree->Branch("fatjet_energy", "vector<float>", &fatjet_energy);
    DYTree->Branch("fatjet_puppi_tau1", "vector<float>", &fatjet_puppi_tau1);
    DYTree->Branch("fatjet_puppi_tau2", "vector<float>", &fatjet_puppi_tau2);
    DYTree->Branch("fatjet_puppi_tau3", "vector<float>", &fatjet_puppi_tau3);
    DYTree->Branch("fatjet_puppi_tau4", "vector<float>", &fatjet_puppi_tau4);
    DYTree->Branch("fatjet_softdropmass", "vector<float>", &fatjet_softdropmass);
    DYTree->Branch("fatjet_chargedHadronEnergyFraction", "vector<float>", &fatjet_chargedHadronEnergyFraction);
    DYTree->Branch("fatjet_neutralHadronEnergyFraction", "vector<float>", &fatjet_neutralHadronEnergyFraction);
    DYTree->Branch("fatjet_neutralEmEnergyFraction", "vector<float>", &fatjet_neutralEmEnergyFraction);
    DYTree->Branch("fatjet_chargedEmEnergyFraction", "vector<float>", &fatjet_chargedEmEnergyFraction);
    DYTree->Branch("fatjet_muonEnergyFraction", "vector<float>", &fatjet_muonEnergyFraction);
    DYTree->Branch("fatjet_chargedMultiplicity", "vector<int>", &fatjet_chargedMultiplicity);
    DYTree->Branch("fatjet_neutralMultiplicity", "vector<int>", &fatjet_neutralMultiplicity);
    DYTree->Branch("fatjet_shiftedEnUp", "vector<float>", &fatjet_shiftedEnUp);
    DYTree->Branch("fatjet_shiftedEnDown", "vector<float>", &fatjet_shiftedEnDown);
    DYTree->Branch("fatjet_smearedRes", "vector<float>", &fatjet_smearedRes);
    DYTree->Branch("fatjet_smearedResUp", "vector<float>", &fatjet_smearedResUp);
    DYTree->Branch("fatjet_smearedResDown", "vector<float>", &fatjet_smearedResDown);
    DYTree->Branch("fatjet_LSF", "vector<float>", &fatjet_LSF);
    DYTree->Branch("fatjet_LSFlep_PID", "vector<float>", &fatjet_LSFlep_PID);
    DYTree->Branch("fatjet_LSFlep_Pt", "vector<float>", &fatjet_LSFlep_Pt);
    DYTree->Branch("fatjet_LSFlep_Eta", "vector<float>", &fatjet_LSFlep_Eta);
    DYTree->Branch("fatjet_LSFlep_Phi", "vector<float>", &fatjet_LSFlep_Phi);
  }

  // Electron
  if( theStoreElectronFlag ){
    DYTree->Branch("electron_MVAIso", "vector<float>", &electron_MVAIso);
    DYTree->Branch("electron_MVANoIso", "vector<float>", &electron_MVANoIso);
    DYTree->Branch("electron_Energy", "vector<float>", &electron_Energy);
    DYTree->Branch("electron_Energy_Scale_Up", "vector<float>", &electron_Energy_Scale_Up);
    DYTree->Branch("electron_Energy_Scale_Down", "vector<float>", &electron_Energy_Scale_Down);
    DYTree->Branch("electron_Energy_Smear_Up", "vector<float>", &electron_Energy_Smear_Up);
    DYTree->Branch("electron_Energy_Smear_Down", "vector<float>", &electron_Energy_Smear_Down);
    DYTree->Branch("electron_eta", "vector<float>", &electron_eta);
    DYTree->Branch("electron_phi", "vector<float>", &electron_phi);
    DYTree->Branch("electron_charge", "vector<int>", &electron_charge);
    DYTree->Branch("electron_gsfpt", "vector<float>", &electron_gsfpt);
    DYTree->Branch("electron_gsfEta", "vector<float>", &electron_gsfEta);
    DYTree->Branch("electron_gsfPhi", "vector<float>", &electron_gsfPhi);
    DYTree->Branch("electron_gsfCharge", "vector<int>", &electron_gsfCharge);
    DYTree->Branch("electron_scEta", "vector<float>", &electron_scEta);
    DYTree->Branch("electron_scPhi", "vector<float>", &electron_scPhi);
    DYTree->Branch("electron_etaWidth", "vector<float>", &electron_etaWidth);
    DYTree->Branch("electron_phiWidth", "vector<float>", &electron_phiWidth);
    DYTree->Branch("electron_dEtaIn", "vector<float>", &electron_dEtaIn);
    DYTree->Branch("electron_dEtaInSeed", "vector<float>", &electron_dEtaInSeed);
    DYTree->Branch("electron_dPhiIn", "vector<float>", &electron_dPhiIn);
    DYTree->Branch("electron_sigmaIEtaIEta", "vector<float>", &electron_sigmaIEtaIEta);
    DYTree->Branch("electron_Full5x5_SigmaIEtaIEta", "vector<float>", &electron_Full5x5_SigmaIEtaIEta);
    DYTree->Branch("electron_e2x5OverE5x5", "vector<float>", &electron_e2x5OverE5x5);
    DYTree->Branch("electron_e1x5OverE5x5", "vector<float>", &electron_e1x5OverE5x5);
    DYTree->Branch("electron_HoverE", "vector<float>", &electron_HoverE);
    DYTree->Branch("electron_fbrem", "vector<float>", &electron_fbrem);
    DYTree->Branch("electron_eOverP", "vector<float>", &electron_eOverP);
    DYTree->Branch("electron_InvEminusInvP", "vector<float>", &electron_InvEminusInvP);
    DYTree->Branch("electron_dxyVTX", "vector<float>", &electron_dxyVTX);
    DYTree->Branch("electron_dxyerrVTX", "vector<float>", &electron_dxyerrVTX);
    DYTree->Branch("electron_dzVTX", "vector<float>", &electron_dzVTX);
    DYTree->Branch("electron_dzerrVTX", "vector<float>", &electron_dzerrVTX);
    DYTree->Branch("electron_3DIPVTX", "vector<float>", &electron_3DIPVTX);
    DYTree->Branch("electron_3DIPerrVTX", "vector<float>", &electron_3DIPerrVTX);
    DYTree->Branch("electron_dxy", "vector<float>", &electron_dxy);
    DYTree->Branch("electron_sigdxy", "vector<float>", &electron_sigdxy);
    DYTree->Branch("electron_dz", "vector<float>", &electron_dz);
    DYTree->Branch("electron_dxyBS", "vector<float>", &electron_dxyBS);
    DYTree->Branch("electron_dzBS", "vector<float>", &electron_dzBS);
    DYTree->Branch("electron_chIso03", "vector<float>", &electron_chIso03);
    DYTree->Branch("electron_nhIso03", "vector<float>", &electron_nhIso03);
    DYTree->Branch("electron_phIso03", "vector<float>", &electron_phIso03);
    DYTree->Branch("electron_puChIso03", "vector<float>", &electron_puChIso03);
    DYTree->Branch("electron_passConversionVeto", "vector<bool>", &electron_passConversionVeto);
    DYTree->Branch("electron_isGsfCtfScPixChargeConsistent", "vector<bool>", &electron_isGsfCtfScPixChargeConsistent);
    DYTree->Branch("electron_isGsfScPixChargeConsistent", "vector<bool>", &electron_isGsfScPixChargeConsistent);
    DYTree->Branch("electron_isGsfCtfChargeConsistent", "vector<bool>", &electron_isGsfCtfChargeConsistent);
    DYTree->Branch("electron_mHits", "vector<int>", &electron_mHits);
    DYTree->Branch("electron_ecalDriven", "vector<int>", &electron_ecalDriven);
    DYTree->Branch("electron_r9", "vector<float>", &electron_r9);
    DYTree->Branch("electron_scEnergy", "vector<float>", &electron_scEnergy);
    DYTree->Branch("electron_scPreEnergy", "vector<float>", &electron_scPreEnergy);
    DYTree->Branch("electron_scRawEnergy", "vector<float>", &electron_scRawEnergy);
    DYTree->Branch("electron_scEt", "vector<float>", &electron_scEt);
    DYTree->Branch("electron_E15", "vector<float>", &electron_E15);
    DYTree->Branch("electron_E25", "vector<float>", &electron_E25);
    DYTree->Branch("electron_E55", "vector<float>", &electron_E55);
    DYTree->Branch("electron_RelPFIso_dBeta", "vector<float>", &electron_RelPFIso_dBeta);
    DYTree->Branch("electron_RelPFIso_Rho", "vector<float>", &electron_RelPFIso_Rho);
    DYTree->Branch("electron_IDBit", "vector<unsigned int>", &electron_IDBit);
    DYTree->Branch("electron_IDCutBit", "vector<int>", &electron_IDCutBit);
    DYTree->Branch("electron_EnergyUnCorr", "vector<float>", &electron_EnergyUnCorr);
    DYTree->Branch("electron_chMiniIso", "vector<float>", &electron_chMiniIso);
    DYTree->Branch("electron_nhMiniIso", "vector<float>", &electron_nhMiniIso);
    DYTree->Branch("electron_phMiniIso", "vector<float>", &electron_phMiniIso);
    DYTree->Branch("electron_puChMiniIso", "vector<float>", &electron_puChMiniIso);
    DYTree->Branch("electron_trackIso", "vector<float>", &electron_trackIso);
    DYTree->Branch("electron_dr03EcalRecHitSumEt", "vector<float>", &electron_dr03EcalRecHitSumEt);
    DYTree->Branch("electron_dr03HcalDepth1TowerSumEt", "vector<float>", &electron_dr03HcalDepth1TowerSumEt);
    DYTree->Branch("electron_dr03HcalTowerSumEt", "vector<float>", &electron_dr03HcalTowerSumEt);
    DYTree->Branch("electron_dr03TkSumPt", "vector<float>", &electron_dr03TkSumPt);
    DYTree->Branch("electron_ecalPFClusterIso", "vector<float>", &electron_ecalPFClusterIso);
    DYTree->Branch("electron_hcalPFClusterIso", "vector<float>", &electron_hcalPFClusterIso);
    DYTree->Branch("electron_pathbits", "vector<ULong64_t>", &electron_pathbits);
    DYTree->Branch("electron_filterbits", "vector<ULong64_t>", &electron_filterbits);
  }
  
  // -- muon variables -- //
  if( theStoreMuonFlag ){

    DYTree->Branch("muon_PfChargedHadronIsoR04", "vector<float>", &muon_PfChargedHadronIsoR04);
    DYTree->Branch("muon_PfNeutralHadronIsoR04", "vector<float>", &muon_PfNeutralHadronIsoR04);
    DYTree->Branch("muon_PfGammaIsoR04", "vector<float>", &muon_PfGammaIsoR04);
    DYTree->Branch("muon_PFSumPUIsoR04", "vector<float>", &muon_PFSumPUIsoR04);
    DYTree->Branch("muon_PfChargedHadronIsoR03", "vector<float>", &muon_PfChargedHadronIsoR03);
    DYTree->Branch("muon_PfNeutralHadronIsoR03", "vector<float>", &muon_PfNeutralHadronIsoR03);
    DYTree->Branch("muon_PfGammaIsoR03", "vector<float>", &muon_PfGammaIsoR03);
    DYTree->Branch("muon_PFSumPUIsoR03", "vector<float>", &muon_PFSumPUIsoR03);
    DYTree->Branch("muon_TypeBit", "vector<unsigned int>", &muon_TypeBit);
    DYTree->Branch("muon_IDBit", "vector<unsigned int>", &muon_IDBit);
    DYTree->Branch("muon_ishighpt", "vector<bool>", &muon_ishighpt);
    DYTree->Branch("muon_dB", "vector<float>", &muon_dB);
    DYTree->Branch("muon_phi", "vector<float>", &muon_phi);
    DYTree->Branch("muon_eta", "vector<float>", &muon_eta);
    DYTree->Branch("muon_pt", "vector<float>", &muon_pt);
    DYTree->Branch("muon_mass", "vector<float>", &muon_mass);
    DYTree->Branch("muon_trkiso", "vector<float>", &muon_trkiso);
    DYTree->Branch("muon_hcaliso", "vector<float>", &muon_hcaliso);
    DYTree->Branch("muon_ecaliso", "vector<float>", &muon_ecaliso);
    DYTree->Branch("muon_trkisoR05", "vector<float>", &muon_trkisoR05);
    DYTree->Branch("muon_hcalisoR05", "vector<float>", &muon_hcalisoR05);
    DYTree->Branch("muon_ecalisoR05", "vector<float>", &muon_ecalisoR05);
    DYTree->Branch("muon_charge", "vector<int>", &muon_charge);
    DYTree->Branch("muon_nChambers", "vector<int>", &muon_nChambers);
    DYTree->Branch("muon_matchedstations", "vector<int>", &muon_matchedstations);
    DYTree->Branch("muon_stationMask", "vector<int>", &muon_stationMask);
    DYTree->Branch("muon_nSegments", "vector<int>", &muon_nSegments);
    DYTree->Branch("muon_normchi", "vector<float>", &muon_normchi);
    DYTree->Branch("muon_validhits", "vector<int>", &muon_validhits);
    DYTree->Branch("muon_trackerHits", "vector<int>", &muon_trackerHits);
    DYTree->Branch("muon_pixelHits", "vector<int>", &muon_pixelHits);
    DYTree->Branch("muon_validmuonhits", "vector<int>", &muon_validmuonhits);
    DYTree->Branch("muon_trackerLayers", "vector<int>", &muon_trackerLayers);
    DYTree->Branch("muon_qoverp", "vector<float>", &muon_qoverp);
    DYTree->Branch("muon_theta", "vector<float>", &muon_theta);
    DYTree->Branch("muon_lambda", "vector<float>", &muon_lambda);
    DYTree->Branch("muon_dxy", "vector<float>", &muon_dxy);
    DYTree->Branch("muon_d0", "vector<float>", &muon_d0);
    DYTree->Branch("muon_dsz", "vector<float>", &muon_dsz);
    DYTree->Branch("muon_dz", "vector<float>", &muon_dz);
    DYTree->Branch("muon_dxyBS", "vector<float>", &muon_dxyBS);
    DYTree->Branch("muon_dzBS", "vector<float>", &muon_dzBS);
    DYTree->Branch("muon_dszBS", "vector<float>", &muon_dszBS);
    DYTree->Branch("muon_dxyVTX", "vector<float>", &muon_dxyVTX);
    DYTree->Branch("muon_dxyerrVTX", "vector<float>", &muon_dxyerrVTX);
    DYTree->Branch("muon_dzVTX", "vector<float>", &muon_dzVTX);
    DYTree->Branch("muon_dzerrVTX", "vector<float>", &muon_dzerrVTX);
    DYTree->Branch("muon_3DIPVTX", "vector<float>", &muon_3DIPVTX);
    DYTree->Branch("muon_3DIPerrVTX", "vector<float>", &muon_3DIPerrVTX);
    DYTree->Branch("muon_dszVTX", "vector<float>", &muon_dszVTX);
    DYTree->Branch("muon_vx", "vector<float>", &muon_vx);
    DYTree->Branch("muon_vy", "vector<float>", &muon_vy);
    DYTree->Branch("muon_vz", "vector<float>", &muon_vz);
    DYTree->Branch("muon_Best_pt", "vector<float>", &muon_Best_pt);
    DYTree->Branch("muon_Best_ptError", "vector<float>", &muon_Best_ptError);
    DYTree->Branch("muon_Best_eta", "vector<float>", &muon_Best_eta);
    DYTree->Branch("muon_Best_phi", "vector<float>", &muon_Best_phi);
    DYTree->Branch("muon_Inner_pt", "vector<float>", &muon_Inner_pt);
    DYTree->Branch("muon_Inner_ptError", "vector<float>", &muon_Inner_ptError);
    DYTree->Branch("muon_Inner_eta", "vector<float>", &muon_Inner_eta);
    DYTree->Branch("muon_Inner_phi", "vector<float>", &muon_Inner_phi);
    DYTree->Branch("muon_Outer_pt", "vector<float>", &muon_Outer_pt);
    DYTree->Branch("muon_Outer_ptError", "vector<float>", &muon_Outer_ptError);
    DYTree->Branch("muon_Outer_eta", "vector<float>", &muon_Outer_eta);
    DYTree->Branch("muon_Outer_phi", "vector<float>", &muon_Outer_phi);
    DYTree->Branch("muon_GLB_pt", "vector<float>", &muon_GLB_pt);
    DYTree->Branch("muon_GLB_ptError", "vector<float>", &muon_GLB_ptError);
    DYTree->Branch("muon_GLB_eta", "vector<float>", &muon_GLB_eta);
    DYTree->Branch("muon_GLB_phi", "vector<float>", &muon_GLB_phi);
    DYTree->Branch("muon_TuneP_pt", "vector<float>", &muon_TuneP_pt);
    DYTree->Branch("muon_TuneP_ptError", "vector<float>", &muon_TuneP_ptError);
    DYTree->Branch("muon_TuneP_eta", "vector<float>", &muon_TuneP_eta);
    DYTree->Branch("muon_TuneP_phi", "vector<float>", &muon_TuneP_phi);
    DYTree->Branch("muon_TuneP_charge", "vector<float>", &muon_TuneP_charge);
    DYTree->Branch("muon_roch_sf", "vector<float>", &muon_roch_sf);
    DYTree->Branch("muon_roch_sf_up", "vector<float>", &muon_roch_sf_up);
    DYTree->Branch("muon_PfChargedHadronMiniIso", "vector<float>", &muon_PfChargedHadronMiniIso);
    DYTree->Branch("muon_PfNeutralHadronMiniIso", "vector<float>", &muon_PfNeutralHadronMiniIso);
    DYTree->Branch("muon_PfGammaMiniIso", "vector<float>", &muon_PfGammaMiniIso);
    DYTree->Branch("muon_PFSumPUMiniIso", "vector<float>", &muon_PFSumPUMiniIso);
    DYTree->Branch("muon_MVA", "vector<float>", &muon_MVA);
    DYTree->Branch("muon_lowptMVA", "vector<float>", &muon_lowptMVA);
    DYTree->Branch("muon_softMVA", "vector<float>", &muon_softMVA);
    DYTree->Branch("muon_jetPtRatio", "vector<float>", &muon_jetPtRatio);
    DYTree->Branch("muon_jetPtRel", "vector<float>", &muon_jetPtRel);
    DYTree->Branch("muon_chi2LocalPosition", "vector<float>", &muon_chi2LocalPosition);
    DYTree->Branch("muon_trkKink", "vector<float>", &muon_trkKink);
    DYTree->Branch("muon_segmentCompatibility", "vector<float>", &muon_segmentCompatibility);
    DYTree->Branch("muon_trackerValidFraction", "vector<float>", &muon_trackerValidFraction);
    DYTree->Branch("muon_simType", "vector<int>", &muon_simType);
    DYTree->Branch("muon_simExtType", "vector<int>", &muon_simExtType);
    DYTree->Branch("muon_simFlavour", "vector<int>", &muon_simFlavour);
    DYTree->Branch("muon_simHeaviestMotherFlavour", "vector<int>", &muon_simHeaviestMotherFlavour);
    DYTree->Branch("muon_simPdgId", "vector<int>", &muon_simPdgId);
    DYTree->Branch("muon_simMotherPdgId", "vector<int>", &muon_simMotherPdgId);
    DYTree->Branch("muon_simMatchQuality", "vector<float>", &muon_simMatchQuality);
    DYTree->Branch("muon_pathbits", "vector<ULong64_t>", &muon_pathbits);
    DYTree->Branch("muon_filterbits", "vector<ULong64_t>", &muon_filterbits);

  }
  
  // -- LHE info -- //
  if( theStoreLHEFlag ){

    DYTree->Branch("LHE_Px", "vector<float>", &LHE_Px);
    DYTree->Branch("LHE_Py", "vector<float>", &LHE_Py);
    DYTree->Branch("LHE_Pz", "vector<float>", &LHE_Pz);
    DYTree->Branch("LHE_E", "vector<float>", &LHE_E);
    DYTree->Branch("LHE_Status", "vector<int>", &LHE_Status);
    DYTree->Branch("LHE_ID", "vector<int>", &LHE_ID);
    for(auto& it:WeightMap){
      DYTree->Branch("weight_"+it.first,"vector<float>",&weight_[it.first]);
    }
  }
  
  // GEN info
  if( theStoreGENFlag ){
    DYTree->Branch("gen_phi", "vector<float>", &gen_phi);
    DYTree->Branch("gen_eta", "vector<float>", &gen_eta);
    DYTree->Branch("gen_pt", "vector<float>", &gen_pt);
    DYTree->Branch("gen_mass", "vector<float>", &gen_mass);
    DYTree->Branch("gen_charge", "vector<float>", &gen_charge);
    DYTree->Branch("gen_mother_index", "vector<int>", &gen_mother_index);
    DYTree->Branch("gen_status", "vector<int>", &gen_status);
    DYTree->Branch("gen_PID", "vector<int>", &gen_PID);
    DYTree->Branch("gen_isPrompt", "vector<bool>", &gen_isPrompt);
    DYTree->Branch("gen_isPromptFinalState", "vector<bool>", &gen_isPromptFinalState);
    DYTree->Branch("gen_isTauDecayProduct", "vector<bool>", &gen_isTauDecayProduct);
    DYTree->Branch("gen_isPromptTauDecayProduct", "vector<bool>", &gen_isPromptTauDecayProduct);
    DYTree->Branch("gen_isDirectPromptTauDecayProductFinalState", "vector<bool>", &gen_isDirectPromptTauDecayProductFinalState);
    DYTree->Branch("gen_isHardProcess", "vector<bool>", &gen_isHardProcess);
    DYTree->Branch("gen_isLastCopy", "vector<bool>", &gen_isLastCopy);
    DYTree->Branch("gen_isLastCopyBeforeFSR", "vector<bool>", &gen_isLastCopyBeforeFSR);
    DYTree->Branch("gen_isPromptDecayed", "vector<bool>", &gen_isPromptDecayed);
    DYTree->Branch("gen_isDecayedLeptonHadron", "vector<bool>", &gen_isDecayedLeptonHadron);
    DYTree->Branch("gen_fromHardProcessBeforeFSR", "vector<bool>", &gen_fromHardProcessBeforeFSR);
    DYTree->Branch("gen_fromHardProcessDecayed", "vector<bool>", &gen_fromHardProcessDecayed);
    DYTree->Branch("gen_fromHardProcessFinalState", "vector<bool>", &gen_fromHardProcessFinalState);
    DYTree->Branch("gen_isMostlyLikePythia6Status3", "vector<bool>", &gen_isMostlyLikePythia6Status3);
    DYTree->Branch("gen_weight", &gen_weight, "gen_weight/F");
    DYTree->Branch("genWeight_Q", &genWeight_Q, "genWeight_Q/F");
    DYTree->Branch("genWeight_X1", &genWeight_X1, "genWeight_X1/F");
    DYTree->Branch("genWeight_X2", &genWeight_X2, "genWeight_X2/F");
    DYTree->Branch("genWeight_id1", &genWeight_id1, "genWeight_id1/I");
    DYTree->Branch("genWeight_id2", &genWeight_id2, "genWeight_id2/I");
    DYTree->Branch("genWeight_alphaQCD", &genWeight_alphaQCD, "genWeight_alphaQCD/F");
    DYTree->Branch("genWeight_alphaQED", &genWeight_alphaQED, "genWeight_alphaQED/F");
  }
  
  if( theStorePhotonFlag ){
    DYTree->Branch("photon_Energy", "vector<float>", &photon_Energy);
    DYTree->Branch("photon_EnergyUnCorr", "vector<float>", &photon_EnergyUnCorr);
    DYTree->Branch("photon_eta", "vector<float>", &photon_eta);
    DYTree->Branch("photon_phi", "vector<float>", &photon_phi);
    DYTree->Branch("photon_scEta", "vector<float>", &photon_scEta);
    DYTree->Branch("photon_scPhi", "vector<float>", &photon_scPhi);
    DYTree->Branch("photon_HoverE", "vector<float>", &photon_HoverE);
    DYTree->Branch("photon_hasPixelSeed", "vector<bool>", &photon_hasPixelSeed);
    DYTree->Branch("photon_Full5x5_SigmaIEtaIEta", "vector<float>", &photon_Full5x5_SigmaIEtaIEta);
    DYTree->Branch("photon_ChIso", "vector<float>", &photon_ChIso);
    DYTree->Branch("photon_NhIso", "vector<float>", &photon_NhIso);
    DYTree->Branch("photon_PhIso", "vector<float>", &photon_PhIso);
    DYTree->Branch("photon_ChIsoWithEA", "vector<float>", &photon_ChIsoWithEA);
    DYTree->Branch("photon_NhIsoWithEA", "vector<float>", &photon_NhIsoWithEA);
    DYTree->Branch("photon_PhIsoWithEA", "vector<float>", &photon_PhIsoWithEA);
    DYTree->Branch("photon_passMVAID_WP80", "vector<bool>", &photon_passMVAID_WP80);
    DYTree->Branch("photon_passMVAID_WP90", "vector<bool>", &photon_passMVAID_WP90);
    DYTree->Branch("photon_passLooseID", "vector<bool>", &photon_passLooseID);
    DYTree->Branch("photon_passMediumID", "vector<bool>", &photon_passMediumID);
    DYTree->Branch("photon_passTightID", "vector<bool>", &photon_passTightID);
  }
  
  
  
  // Pile-up Reweight
  DYTree->Branch("nPileUp",&nPileUp,"nPileUp/I");
  
  if( theStoreMETFlag ){
    DYTree->Branch("pfMET_pt", &pfMET_pt, "pfMET_pt/F");
    DYTree->Branch("pfMET_phi", &pfMET_phi, "pfMET_phi/F");
    DYTree->Branch("pfMET_SumEt", &pfMET_SumEt, "pfMET_SumEt/F");
    DYTree->Branch("pfMET_Type1_pt", &pfMET_Type1_pt, "pfMET_Type1_pt/F");
    DYTree->Branch("pfMET_Type1_phi", &pfMET_Type1_phi, "pfMET_Type1_phi/F");
    DYTree->Branch("pfMET_Type1_SumEt", &pfMET_Type1_SumEt, "pfMET_Type1_SumEt/F");
    DYTree->Branch("pfMET_Type1_PhiCor_pt", &pfMET_Type1_PhiCor_pt, "pfMET_Type1_PhiCor_pt/F");
    DYTree->Branch("pfMET_Type1_PhiCor_phi", &pfMET_Type1_PhiCor_phi, "pfMET_Type1_PhiCor_phi/F");
    DYTree->Branch("pfMET_Type1_PhiCor_SumEt", &pfMET_Type1_PhiCor_SumEt, "pfMET_Type1_PhiCor_SumEt/F");
    DYTree->Branch("pfMET_pt_shifts", "vector<float>", &pfMET_pt_shifts);
    DYTree->Branch("pfMET_phi_shifts", "vector<float>", &pfMET_phi_shifts);
    DYTree->Branch("pfMET_SumEt_shifts", "vector<float>", &pfMET_SumEt_shifts);
    DYTree->Branch("pfMET_Type1_pt_shifts", "vector<float>", &pfMET_Type1_pt_shifts);
    DYTree->Branch("pfMET_Type1_phi_shifts", "vector<float>", &pfMET_Type1_phi_shifts);
    DYTree->Branch("pfMET_Type1_SumEt_shifts", "vector<float>", &pfMET_Type1_SumEt_shifts);
    DYTree->Branch("pfMET_Type1_PhiCor_pt_shifts", "vector<float>", &pfMET_Type1_PhiCor_pt_shifts);
    DYTree->Branch("pfMET_Type1_PhiCor_phi_shifts", "vector<float>", &pfMET_Type1_PhiCor_phi_shifts);
    DYTree->Branch("pfMET_Type1_PhiCor_SumEt_shifts", "vector<float>", &pfMET_Type1_PhiCor_SumEt_shifts);
  }

  if(theDebugLevel) cout << "[SKFlatMaker::beginJob] finished" << endl;

}

void SKFlatMaker::beginRun(const Run & iRun, const EventSetup & iSetup)
{
  
  if(theDebugLevel) cout << "[SKFlatMaker::beginRun] called" << endl;
  
  vector<string> temp_trigs = {
      "HLT_Mu*", "HLT_Ele*", "HLT_DoubleEle*", "HLT_DoublePhoton*", "HLT_IsoMu*", "HLT_Photon*",
      "HLT_OldMu100_v*", "HLT_TkMu100_v*",
      "HLT_TkMu50_v*",
      "HLT_IsoTkMu24_v*",
      "HLT_TripleMu*",
      "HLT_DiMu*",
      "HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_v*",
      "HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_v*",

/*
      //==== single muon triggers
      "HLT_IsoMu27_v*",
      "HLT_Mu50_v*",
            
      //==== double muon triggers
      "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",
      
      //==== Single Electron
      "HLT_Ele35_WPTight_Gsf_v*",
      
      //==== Double Electron
      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*", // low pt, loose ID & iso, unprescaled 
      "HLT_DoubleEle33_CaloIdL_MW_v*", // loose ID, no isolation
      
      //==== Double Photon
      "HLT_DoublePhoton33_CaloIdL_v*"
*/

  };
  //==== copy this to member variable, HLTName_WildCard
  HLTName_WildCard.clear();
  for(unsigned int i=0; i<temp_trigs.size(); i++ ) HLTName_WildCard.push_back(temp_trigs.at(i));

  bool changedConfig;
  if(!hltConfig_.init(iRun, iSetup, processName, changedConfig)){
    LogError("HLTMuonVal") << "Initialization of HLTConfigProvider failed!!";
    return;
  }

  if(theDebugLevel) cout << "[SKFlatMaker::beginRun] finished" << endl;

}

// ------------ method called once each job just after ending the event loop  ------------ //
void SKFlatMaker::endJob()
{
  if(theDebugLevel) cout << "[SKFlatMaker::endJob] endJob" << endl;
  std::cout <<"[SKFlatMaker::endJob] ++++++++++++++++++++++++++++++++++++++" << std::endl;
  std::cout <<"[SKFlatMaker::endJob] analyzed " << nEvt << " events: " << std::endl;
  std::cout <<"[SKFlatMaker::endJob] ++++++++++++++++++++++++++++++++++++++" << std::endl;
}

///////////////////////////////////////////////////////
// -- makes hlt report and fills it to the ntuple -- //
///////////////////////////////////////////////////////
void SKFlatMaker::hltReport(const edm::Event &iEvent)
{

  if(theDebugLevel) cout << "[SKFlatMaker::hltReport] called" << endl;

  Handle<TriggerResults> trigResult;
  iEvent.getByToken(TriggerToken, trigResult);

  if( !trigResult.failedToGet() ){

    //==== All trigger paths in this event setup
    const edm::TriggerNames trigName = iEvent.triggerNames(*trigResult);

    if(theDebugLevel>=2){
      cout << "[SKFlatMaker::hltReport] trigger names in trigger result (HLT)" << endl;
      cout << "[SKFlatMaker::hltReport] trigName.size() = " << trigName.size() << endl;
      for(int itrig=0; itrig<(int)trigName.size(); itrig++)
        cout << "[SKFlatMaker::hltReport] trigName = " << trigName.triggerName(itrig) << " " << itrig << endl;
    }


    //==== iteration for each Trigger Name WildCards
    //==== e.g., {"HLT_Mu*", "HLT_Ele*"}
    for(unsigned int it_HLTWC = 0; it_HLTWC < HLTName_WildCard.size(); it_HLTWC++ ){
      if(theDebugLevel) cout << "[SKFlatMaker::hltReport] [" << it_HLTWC << "th Input Trigger Name WildCard = " << HLTName_WildCard[it_HLTWC] << "]" << endl;

      //==== find triggers in HLT matched to this WildCard
      std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(trigName.triggerNames(), HLTName_WildCard[it_HLTWC]);

      if( !matches.empty() ){

        //==== iteration for each wildcard-matched triggers
        //==== e.g., {"HLT_Mu8", "HLT_Mu17"}
        BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches){

          //==== Cleaningup
          //==== https://github.com/CMSSNU/SKFlatMaker/issues/9

          if(HLTName_WildCard[it_HLTWC]=="HLT_Mu*"){
            if( (*match).find("HLT_Mu7p5") != std::string::npos ) continue;
            if( (*match).find("Phi") != std::string::npos ) continue;
            if( (*match).find("PFMETNoMu") != std::string::npos ) continue;

          }
          if( (*match).find("Onia") != std::string::npos ) continue;
          if( (*match).find("BTagCSV") != std::string::npos ) continue;
          if( (*match).find("CaloBTagCSV") != std::string::npos ) continue;
          if( (*match).find("Tau") != std::string::npos ) continue;
          if( (*match).find("EBOnly") != std::string::npos ) continue;
          if( (*match).find("R9Id90") != std::string::npos ) continue;
          if( (*match).find("R9Id85") != std::string::npos ) continue;
          if( (*match).find("DisplacedIdL") != std::string::npos ) continue;
          if( (*match).find("HighEta") != std::string::npos ) continue;
          if( (*match).find("EleCleaned") != std::string::npos ) continue;
	  if( (*match).find("_DoubleTrkMu5NoFiltersNoVtx_v") == std::string::npos )
	    if( (*match).find("NoFiltersNoVtx") != std::string::npos ) continue;
          if( (*match).find("WHbb") != std::string::npos ) continue;
	  if( (*match).find("HLT_Ele27_eta2p1_WPTight_Gsf_v") == std::string::npos )
	    if( (*match).find("eta2p1") != std::string::npos ) continue;

          if(theDebugLevel) cout << "[SKFlatMaker::hltReport]   [matched trigger = " << *match << "]" << endl;

          //==== Check if this trigger is fired
          if( trigResult->accept(trigName.triggerIndex(*match)) ){
            HLT_TriggerName.push_back(*match); //save HLT list as a vector
          }

          //==== Save Filter if theStoreHLTObjectFlag
          if(theStoreHLTObjectFlag){
            //==== find modules corresponding to a trigger in HLT configuration
            std::vector<std::string> moduleNames = hltConfig_.moduleLabels( *match );
            if(theDebugLevel){
              cout << "[SKFlatMaker::hltReport]   moduleNames.size() = " << moduleNames.size() << endl;
              for(unsigned int it_md=0; it_md<moduleNames.size(); it_md++){
                if(theDebugLevel) cout << "[SKFlatMaker::hltReport]     " << moduleNames.at(it_md) << endl;
              }
            }
            string this_modulename = "";
            //==== Last module = moduleNames[nmodules-1] is always "hltBoolEnd"
            int nmodules = moduleNames.size();
            if( nmodules-2 >= 0 ){
              this_modulename = moduleNames[nmodules-2];
            }
            HLT_TriggerFilterName.push_back( this_modulename );
          }


        } // END Loop over matches

      } // END matches not empty

    } // END Loop over Trigger Wildcards

    //==== Now, all trigers we are interested in are saved

    if(theStoreHLTObjectFlag){

      //==== HLT Object for lepton matching
      edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerObject;
      iEvent.getByToken(TriggerObjectToken, triggerObject);

      for(pat::TriggerObjectStandAlone obj : *triggerObject){

        obj.unpackPathNames(trigName);
        obj.unpackFilterLabels(iEvent, *trigResult);  //added Suoh

        HLTObject_pt.push_back( obj.pt() );
        HLTObject_eta.push_back( obj.eta() );
        HLTObject_phi.push_back( obj.phi() );

        string FiredFilters = "", FiredPaths = "";

        //==== Path we are interested in are saved in : std::vector<std::string > HLT_TriggerName;
        //==== Filter we are interested in are svaed in : std::vector<std::string > HLT_TriggerFilterName;

        //==== Loop over filters
        //cout << "This HLT Object : " << endl;
        for( size_t i_filter = 0; i_filter < obj.filterLabels().size(); ++i_filter ){
          string this_filter = obj.filterLabels().at(i_filter);
          //cout << "  this_filter = " << this_filter << endl;
          if(std::find( HLT_TriggerFilterName.begin(), HLT_TriggerFilterName.end(), this_filter) != HLT_TriggerFilterName.end() ){
            FiredFilters += this_filter+":";
          }
        } //==== END Filter Loop

        //==== Loop over path
        std::vector<std::string> pathNamesAll  = obj.pathNames(true); // Get path whose last filter is fired
        for( size_t i_filter = 0; i_filter < pathNamesAll.size(); ++i_filter ){
          string this_path = pathNamesAll.at(i_filter);
          //cout << "  this_path = " << this_path << endl;
          if(std::find( HLT_TriggerName.begin(), HLT_TriggerName.end(), this_path ) != HLT_TriggerName.end()){
            FiredPaths += pathNamesAll.at(i_filter)+":";
          }
        } //==== END Filter Loop

        if(theDebugLevel>=2){
          cout << "  ==> FiredFilters = " << FiredFilters << endl;
          cout << "  ==> FiredPaths = " << FiredPaths << endl;
        }

        HLTObject_FiredFilters.push_back( FiredFilters );
        HLTObject_FiredPaths.push_back( FiredPaths );

      } //==== END HLT Object loop

    }

  } // -- end of if( !trigResult.failedToGet() ) -- //

/*
  if( IsData ){
    Handle<TriggerResults> trigResultPAT;
    iEvent.getByToken(TriggerTokenPAT, trigResultPAT);
      
    if( !trigResultPAT.failedToGet() ){
      const edm::TriggerNames trigName = iEvent.triggerNames(*trigResultPAT);

      // cout << "trigger names in trigger result (PAT)" << endl;
      // for(int itrig=0; itrig<(int)trigName.size(); itrig++)
      //   cout << "trigName = " << trigName.triggerName(itrig) << " " << itrig << endl;

      if( trigResultPAT->accept(trigName.triggerIndex("Flag_goodVertices")) ) Flag_goodVertices = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_globalTightHalo2016Filter")) ) Flag_globalTightHalo2016Filter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_globalSuperTightHalo2016Filter")) ) Flag_globalSuperTightHalo2016Filter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_HBHENoiseFilter")) ) Flag_HBHENoiseFilter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_HBHENoiseIsoFilter")) ) Flag_HBHENoiseIsoFilter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_EcalDeadCellTriggerPrimitiveFilter")) ) Flag_EcalDeadCellTriggerPrimitiveFilter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_BadPFMuonFilter")) ) Flag_BadPFMuonFilter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_BadChargedCandidateFilter")) ) Flag_BadChargedCandidateFilter = true;
      if( trigResultPAT->accept(trigName.triggerIndex("Flag_eeBadScFilter")) ) Flag_eeBadScFilter = true;
      //if( trigResultPAT->accept(trigName.triggerIndex("Flag_ecalBadCalibReducedMINIAODFilter")) ) Flag_ecalBadCalibReducedMINIAODFilter = true;
 
    }
  }
*/
  
  //==============
  //==== MiniAOD
  //==============
  
  //save event filter infomation
  if (!iEvent.getByToken(METFilterResultsToken_PAT, METFilterResults)){
    iEvent.getByToken(METFilterResultsToken_RECO, METFilterResults);
  }
  
  const edm::TriggerNames &metNames = iEvent.triggerNames(*METFilterResults);
  
  //*******************************************************************************
  //For Debug printout
  //*******************************************************************************
  //cout << "----------------------------" << endl;
  //for (unsigned int i = 0, n = METFilterResults->size(); i < n; ++i) {
  //  std::cout << "MET Filter " << metNames.triggerName(i).c_str() << "\n";
  //}
    
  edm::Handle< bool > passBadPFMuonDzFilter;
  if(iEvent.getByToken(BadPFMuonDzFilter_token,passBadPFMuonDzFilter)){
    Flag_BadPFMuonDzFilter = (*passBadPFMuonDzFilter);
  }

  for(unsigned int i = 0, n = METFilterResults->size(); i < n; ++i){

    if(strcmp(metNames.triggerName(i).c_str(), "Flag_goodVertices") == 0) Flag_goodVertices = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_globalSuperTightHalo2016Filter") == 0) Flag_globalSuperTightHalo2016Filter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseFilter") == 0) Flag_HBHENoiseFilter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_HBHENoiseIsoFilter") == 0) Flag_HBHENoiseIsoFilter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_EcalDeadCellTriggerPrimitiveFilter") == 0) Flag_EcalDeadCellTriggerPrimitiveFilter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonFilter") == 0) Flag_BadPFMuonFilter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadPFMuonDzFilter") == 0) Flag_BadPFMuonDzFilter = METFilterResults -> accept(i);
    //else if(strcmp(metNames.triggerName(i).c_str(), "Flag_hfNoisyHitsFilter") == 0) Flag_hfNoisyHitsFilter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_BadChargedCandidateFilter") == 0) Flag_BadChargedCandidateFilter = METFilterResults -> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_eeBadScFilter") == 0) Flag_eeBadScFilter = METFilterResults-> accept(i);
    else if(strcmp(metNames.triggerName(i).c_str(), "Flag_ecalBadCalibFilter") == 0) Flag_ecalBadCalibFilter = METFilterResults -> accept(i);
  }
}

///////////////////////////////////
// -- Get Primary vertex info -- //
///////////////////////////////////
void SKFlatMaker::fillPrimaryVertex(const edm::Event &iEvent)
{
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(PrimaryVertexToken, pvHandle);
  const reco::VertexCollection vtx = *(pvHandle.product());
  nPV = pvHandle->size();

  if( vtx.size() > 2 && theDebugLevel > 0) cout << "[SKFlatMaker::fillPrimaryVertex] Reconstructed "<< vtx.size() << " vertices" << endl;
  if (vtx.size() > 0 ){
    PVtrackSize = vtx.front().tracksSize();
    PVchi2 = vtx.front().chi2();
    PVndof = vtx.front().ndof();
    PVnormalizedChi2 = vtx.front().normalizedChi2();
    PVx = vtx.front().x();
    PVy = vtx.front().y();
    PVz = vtx.front().z();
    PVprob = TMath::Prob(PVchi2,(int)PVndof);
   }
  
}

//////////////////////////////
// -- Get Muons info -- //
//////////////////////////////
void SKFlatMaker::fillMuons(const edm::Event &iEvent, const edm::EventSetup& iSetup)
{
  // -- BeamSpot -- //
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(BeamSpotToken, beamSpotHandle);
  reco::BeamSpot beamSpot = (*beamSpotHandle);
  
  // -- Vertex -- //
  math::XYZPoint RefVtx;
  RefVtx.SetXYZ(0, 0, 0);
  
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(PrimaryVertexToken, pvHandle);
  const reco::VertexCollection &vertices = *pvHandle.product();

  // -- What is the purpose of below line? -- //
  for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it){
    RefVtx = it->position();
    break;
  }
  const reco::Vertex &vtx = pvHandle->front();
  
  // muons
  ESHandle<MagneticField> B;
  iSetup.get<IdealMagneticFieldRecord>().get(B);
  
  // -- Call PAT muons -- //
  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByToken(MuonToken, muonHandle);
  using reco::MuonCollection;
  MuonCollection::const_iterator imuon;

  edm::Handle<reco::GenParticleCollection> genParticles;
  if(!IsData){
    iEvent.getByToken(mcLabel_,genParticles);
  }

  //==== HLT Object for lepton matching
  edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerObject;
  iEvent.getByToken(TriggerObjectToken, triggerObject);
  Handle<TriggerResults> trigResult;
  iEvent.getByToken(TriggerToken, trigResult);
  const edm::TriggerNames trigNames = iEvent.triggerNames(*trigResult);

  for( unsigned int i = 0; i != muonHandle->size(); i++ ){
    // cout << "##### Analyze:Start the loop for the muon #####" << endl;
    const pat::Muon imuon = muonHandle->at(i);
 
    //============================================================================================
    //==== Muon selectors (Since 9_4_X) 
    //==== https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Muon_selectors_Since_9_4_X
    //============================================================================================

    //==== DataFormats/MuonReco/interface/Muon.h
    //==== imuon.isStandAloneMuon() : same as type_ & StandAloneMuon
    //==== imuon.isGlobalMuon()  : same as type_ & GlobalMuon;
    //==== imuon.isTrackerMuon() : same as type_ & TrackerMuon;
    //==== imuon.isPFMuon() : same as type_ & PFMuon;

    //==== https://github.com/cms-sw/cmssw/blob/acfa3a98a3433375e7814065546d927851104573/DataFormats/MuonReco/src/MuonSelectors.cc#L1003
    //==== imuon.isTightMuon(vtx) : same as imuon.passed( reco::Muon::CutBasedIdTight );
    //==== imuon.isMediumMuon() : same as imuon.passed( reco::Muon::CutBasedIdMedium );
    //==== imuon.isLooseMuon() : same as imuon.passed( reco::Muon::CutBasedIdLoose );
    //==== imuon.isSoftMuon(vtx) : same as imuon.passed( reco::Muon::SoftCutBasedId );
    //==== imuon.isHighPtMuon(vtx) : same as imuon.passed( reco::Muon::CutBasedIdGlobalHighPt );

    muon_TypeBit.push_back( imuon.type() );
    muon_IDBit.push_back( imuon.selectors() );

    // -- bits 0-1-2-3 = DT stations 1-2-3-4 -- //
    // -- bits 4-5-6-7 = CSC stations 1-2-3-4 -- //
    int _segments = 0;
      
    for( int idet = 1; idet < 4; idet++ ){
      // -- DT (1), CSC (2), RPC (3) -- //
      for( int istation = 1; istation < 5; istation++ ){
        // -- station 1, 2, 3, 4 -- //
        _segments += imuon.numberOfSegments(istation, idet);
      }
    }
    muon_nSegments.push_back( _segments );
      
    // cout << "##### Analyze:Muon Type #####" << endl;
      
      
    // -- reco track information -- //
    reco::TrackRef trackerTrack = imuon.innerTrack();
    reco::TrackRef muonTrack    = imuon.outerTrack();
    reco::TrackRef glbTrack     = imuon.globalTrack();
    
    // cout << "##### Analyze:Muon Tracks #####" << endl;

    muon_validhits.push_back( imuon.numberOfValidHits() );
    muon_chi2LocalPosition.push_back( imuon.combinedQuality().chi2LocalPosition );
    muon_trkKink.push_back( imuon.combinedQuality().trkKink );
    muon_segmentCompatibility.push_back( muon::segmentCompatibility(imuon) );

    //==== Global track
    if( glbTrack.isNonnull() ){

      const reco::HitPattern & glbhit = glbTrack->hitPattern();

      muon_normchi.push_back( glbTrack->normalizedChi2() );
      muon_validmuonhits.push_back( glbhit.numberOfValidMuonHits() );

      muon_qoverp.push_back( glbTrack->qoverp() );
      muon_theta.push_back( glbTrack->theta() );
      muon_lambda.push_back( glbTrack->lambda() );
      muon_dxy.push_back( glbTrack->dxy() );
      muon_d0.push_back( glbTrack->d0() );
      muon_dsz.push_back( glbTrack->dsz() );
      muon_dz.push_back( glbTrack->dz() );
      muon_dxyBS.push_back( glbTrack->dxy(beamSpot.position()) );
      muon_dszBS.push_back( glbTrack->dsz(beamSpot.position()) );
      muon_dzBS.push_back( glbTrack->dz(beamSpot.position()) );
      
      muon_vx.push_back( glbTrack->vx() );
      muon_vy.push_back( glbTrack->vy() );
      muon_vz.push_back( glbTrack->vz() );
      
    }
    //==== No Global track, but innertrack
    else if( trackerTrack.isNonnull() ){

      const reco::HitPattern & trkhit = trackerTrack->hitPattern();

      muon_normchi.push_back( trackerTrack->normalizedChi2()  );
      muon_validmuonhits.push_back( trkhit.numberOfValidMuonHits() );
      
      muon_qoverp.push_back( trackerTrack->qoverp() );
      muon_theta.push_back( trackerTrack->theta() );
      muon_lambda.push_back( trackerTrack->lambda() );
      muon_dxy.push_back( trackerTrack->dxy() );
      muon_d0.push_back( trackerTrack->d0() );
      muon_dsz.push_back( trackerTrack->dsz() );
      muon_dz.push_back( trackerTrack->dz() );
      muon_dxyBS.push_back( trackerTrack->dxy(beamSpot.position()) );
      muon_dszBS.push_back( trackerTrack->dsz(beamSpot.position()) );
      muon_dzBS.push_back( trackerTrack->dz(beamSpot.position()) );
      
      muon_vx.push_back( trackerTrack->vx() );
      muon_vy.push_back( trackerTrack->vy() );
      muon_vz.push_back( trackerTrack->vz() );
      
    }
    //==== No Global, No Tracker -> StandAlone
    else{

      const reco::HitPattern & muonhit = muonTrack->hitPattern();

      muon_normchi.push_back( muonTrack->normalizedChi2()  );
      muon_validmuonhits.push_back( muonhit.numberOfValidMuonHits() );

      muon_qoverp.push_back( muonTrack->qoverp() );
      muon_theta.push_back( muonTrack->theta() );
      muon_lambda.push_back( muonTrack->lambda() );
      muon_dxy.push_back( muonTrack->dxy() );
      muon_d0.push_back( muonTrack->d0() );
      muon_dsz.push_back( muonTrack->dsz() );
      muon_dz.push_back( muonTrack->dz() );
      muon_dxyBS.push_back( muonTrack->dxy(beamSpot.position()) );
      muon_dszBS.push_back( muonTrack->dsz(beamSpot.position()) );
      muon_dzBS.push_back( muonTrack->dz(beamSpot.position()) );

      muon_vx.push_back( muonTrack->vx() );
      muon_vy.push_back( muonTrack->vy() );
      muon_vz.push_back( muonTrack->vz() );

    }

    if( trackerTrack.isNonnull() ){
      const reco::HitPattern & inhit = trackerTrack->hitPattern();
      
      muon_trackerHits.push_back( inhit.numberOfValidTrackerHits() );
      muon_trackerValidFraction.push_back( trackerTrack->validFraction() );
      muon_pixelHits.push_back( inhit.numberOfValidPixelHits() );
      muon_trackerLayers.push_back( inhit.trackerLayersWithMeasurement() );
    }
    else{
      muon_trackerHits.push_back( 0 );
      muon_trackerValidFraction.push_back( 0 );
      muon_pixelHits.push_back( 0 );
      muon_trackerLayers.push_back( 0 );
    }
      
    if( !pvHandle->empty() && !pvHandle->front().isFake() ){
/*
      muon_dxyVTX.push_back( imuon.muonBestTrack()->dxy(vtx.position()) );
      muon_dzVTX.push_back( imuon.muonBestTrack()->dz(vtx.position()) );
      muon_dszVTX.push_back( imuon.muonBestTrack()->dsz(vtx.position()) );
*/

      //==== https://github.com/cms-sw/cmssw/blob/3e032f6f8d88a5a63769f68fe525b4f06d72b55f/PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc#L322-L345
      //==== PAT dB are set with BestTrack

      muon_dxyVTX.push_back( imuon.dB(pat::Muon::PV2D) );
      muon_dxyerrVTX.push_back( imuon.edB(pat::Muon::PV2D) );
      muon_dzVTX.push_back( imuon.dB(pat::Muon::PVDZ) );
      muon_dzerrVTX.push_back( imuon.edB(pat::Muon::PVDZ) );
      muon_3DIPVTX.push_back( imuon.dB(pat::Muon::PV3D) );
      muon_3DIPerrVTX.push_back( imuon.edB(pat::Muon::PV3D) );

      muon_dszVTX.push_back( imuon.muonBestTrack()->dsz(vtx.position()) );

    }
    else{

      muon_dxyVTX.push_back( 9999 );
      muon_dxyerrVTX.push_back( 9999 );
      muon_dzVTX.push_back( 9999 );
      muon_dzerrVTX.push_back( 9999 );
      muon_3DIPVTX.push_back( 9999 );
      muon_3DIPerrVTX.push_back( 9999 );

      muon_dszVTX.push_back( 9999 );
    }
      
    muon_pt.push_back( imuon.pt() );
    muon_mass.push_back( imuon.mass() );
    muon_eta.push_back( imuon.eta() );
    muon_phi.push_back( imuon.phi() );
    
    muon_dB.push_back( imuon.dB() );
      
    // -- Various track informations -- //
    // -- MuonBestTrack -- //
    if( imuon.muonBestTrack().isNonnull() ){
      muon_Best_pt.push_back( imuon.muonBestTrack()->pt() );
      muon_Best_ptError.push_back( imuon.muonBestTrack()->ptError() );
      muon_Best_eta.push_back( imuon.muonBestTrack()->eta() );
      muon_Best_phi.push_back( imuon.muonBestTrack()->phi() );
    }
    else{
      muon_Best_pt.push_back( -999 );
      muon_Best_ptError.push_back( -999 );
      muon_Best_eta.push_back( -999 );
      muon_Best_phi.push_back( -999 );
    }
      
      
    // -- Inner Track -- //
    if( imuon.innerTrack().isNonnull() ){
      muon_Inner_pt.push_back( imuon.innerTrack()->pt() );
      muon_Inner_ptError.push_back( imuon.innerTrack()->ptError() );
      muon_Inner_eta.push_back( imuon.innerTrack()->eta() );
      muon_Inner_phi.push_back( imuon.innerTrack()->phi() );
    }
    else{
      muon_Inner_pt.push_back( -999 );
      muon_Inner_ptError.push_back( -999 );
      muon_Inner_eta.push_back( -999 );
      muon_Inner_phi.push_back( -999 );
    }
      
    // -- Outer Track -- //
    if( imuon.outerTrack().isNonnull() ){
      muon_Outer_pt.push_back( imuon.outerTrack()->pt() );
      muon_Outer_ptError.push_back( imuon.outerTrack()->ptError() );
      muon_Outer_eta.push_back( imuon.outerTrack()->eta() );
      muon_Outer_phi.push_back( imuon.outerTrack()->phi() );
    }
    else{
      muon_Outer_pt.push_back( -999 );
      muon_Outer_ptError.push_back( -999 );
      muon_Outer_eta.push_back( -999 );
      muon_Outer_phi.push_back( -999 );
    }
      
    // -- Global Track -- //
    if( imuon.globalTrack().isNonnull() ){
      muon_GLB_pt.push_back( imuon.globalTrack()->pt() );
      muon_GLB_ptError.push_back( imuon.globalTrack()->ptError() );
      muon_GLB_eta.push_back( imuon.globalTrack()->eta() );
      muon_GLB_phi.push_back( imuon.globalTrack()->phi() );
    }
    else{
      muon_GLB_pt.push_back( -999 );
      muon_GLB_ptError.push_back( -999 );
      muon_GLB_eta.push_back( -999 );
      muon_GLB_phi.push_back( -999 );
    }
      
    // -- tuneP MuonBestTrack -- //
    if( imuon.tunePMuonBestTrack().isNonnull() ){
      muon_TuneP_pt.push_back( imuon.tunePMuonBestTrack()->pt() );
      muon_TuneP_ptError.push_back( imuon.tunePMuonBestTrack()->ptError() );
      muon_TuneP_eta.push_back( imuon.tunePMuonBestTrack()->eta() );
      muon_TuneP_phi.push_back( imuon.tunePMuonBestTrack()->phi() );
      muon_TuneP_charge.push_back( imuon.tunePMuonBestTrack()->charge() );
    }
    else{
      muon_TuneP_pt.push_back( -999 );
      muon_TuneP_ptError.push_back( -999 );
      muon_TuneP_eta.push_back( -999 );
      muon_TuneP_phi.push_back( -999 );
      muon_TuneP_charge.push_back( -999 );
    }
      
    //-- ISOLATIONS GO HERE -- //
    // -- detector based -- //
    muon_trkiso.push_back( imuon.isolationR03().sumPt );
    muon_hcaliso.push_back( imuon.isolationR03().hadEt );
    muon_ecaliso.push_back( imuon.isolationR03().emEt );
    muon_trkisoR05.push_back( imuon.isolationR05().sumPt );
    muon_hcalisoR05.push_back( imuon.isolationR05().hadEt );
    muon_ecalisoR05.push_back( imuon.isolationR05().emEt ); 
    
    // -- pf isolation -- // 
    muon_PfChargedHadronIsoR04.push_back( imuon.pfIsolationR04().sumChargedHadronPt );
    muon_PfNeutralHadronIsoR04.push_back( imuon.pfIsolationR04().sumNeutralHadronEt );
    muon_PfGammaIsoR04.push_back( imuon.pfIsolationR04().sumPhotonEt );
    muon_PFSumPUIsoR04.push_back( imuon.pfIsolationR04().sumPUPt );
    
    muon_PfChargedHadronIsoR03.push_back( imuon.pfIsolationR03().sumChargedHadronPt );
    muon_PfNeutralHadronIsoR03.push_back( imuon.pfIsolationR03().sumNeutralHadronEt );
    muon_PfGammaIsoR03.push_back( imuon.pfIsolationR03().sumPhotonEt );
    muon_PFSumPUIsoR03.push_back( imuon.pfIsolationR03().sumPUPt );

/*
    //==== MiniIso
    PFIsolation this_miniiso = GetMiniIso(pc, imuon.p4(),
                                          miniIsoParams_[0], miniIsoParams_[1], miniIsoParams_[2],
                                          miniIsoParams_[3], miniIsoParams_[4], miniIsoParams_[5],
                                          miniIsoParams_[6], miniIsoParams_[7], miniIsoParams_[8]);

    cout << "=======================" << endl;
    cout << "Pt = " << imuon.pt() << endl;
    cout << "---- R04 ----" << endl;
    cout << "CH = " << imuon.pfIsolationR04().sumChargedHadronPt << endl;
    cout << "NH = " << imuon.pfIsolationR04().sumNeutralHadronEt << endl;
    cout << "Ph = " << imuon.pfIsolationR04().sumPhotonEt << endl;
    cout << "PU = " << imuon.pfIsolationR04().sumPUPt << endl;
    cout << "---- Mini ----" << endl;
    cout << "CH = " << this_miniiso.chargedHadronIso() << endl;
    cout << "HN = " << this_miniiso.neutralHadronIso() << endl;
    cout << "Ph = " << this_miniiso.photonIso() << endl;
    cout << "PU = " << this_miniiso.puChargedHadronIso() << endl;
    cout << "---- MiniAOD MiniIso ----" << endl;
    cout << "CH = " << imuon.miniPFIsolation().chargedHadronIso() << endl;
    cout << "HN = " << imuon.miniPFIsolation().neutralHadronIso() << endl;
    cout << "Ph = " << imuon.miniPFIsolation().photonIso() << endl;
    cout << "PU = " << imuon.miniPFIsolation().puChargedHadronIso() << endl;
    cout << "---- Diff ----" << endl;
    cout << "dCH = " << this_miniiso.chargedHadronIso() - imuon.miniPFIsolation().chargedHadronIso() << endl;
    cout << "dNH = " << this_miniiso.neutralHadronIso() - imuon.miniPFIsolation().neutralHadronIso() << endl;
    cout << "dPh = " << this_miniiso.photonIso() - imuon.miniPFIsolation().photonIso() << endl;
    cout << "dPU = " << this_miniiso.puChargedHadronIso() - imuon.miniPFIsolation().puChargedHadronIso() << endl;

    muon_PfChargedHadronMiniIso.push_back( this_miniiso.chargedHadronIso() );
    muon_PfNeutralHadronMiniIso.push_back( this_miniiso.neutralHadronIso() );
    muon_PfGammaMiniIso.push_back( this_miniiso.photonIso() );
    muon_PFSumPUMiniIso.push_back( this_miniiso.puChargedHadronIso() );
*/

    muon_PfChargedHadronMiniIso.push_back( imuon.miniPFIsolation().chargedHadronIso() );
    muon_PfNeutralHadronMiniIso.push_back( imuon.miniPFIsolation().neutralHadronIso() );
    muon_PfGammaMiniIso.push_back( imuon.miniPFIsolation().photonIso() );
    muon_PFSumPUMiniIso.push_back( imuon.miniPFIsolation().puChargedHadronIso() );

    // -- Else -- //
    muon_charge.push_back( imuon.charge() );
    muon_nChambers.push_back( imuon.numberOfChambers() ); // -- # of chambers -- //
    muon_matchedstations.push_back( imuon.numberOfMatchedStations() ); // -- # of chambers with matched segments -- //
    muon_stationMask.push_back( imuon.stationMask() ); // -- bit map of stations with matched segments -- //

    //==== Muon mva
    //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#Lepton_MVA
    muon_MVA.push_back( imuon.mvaValue() );
    //muon_lowptMVA.push_back( imuon.lowptMvaValue() ); // TODO not supported in CMSSW_10_2_10
    muon_softMVA.push_back( imuon.softMvaValue() );

    //==== Muon miniiso variable
    muon_jetPtRatio.push_back( imuon.jetPtRatio() );
    muon_jetPtRel.push_back( imuon.jetPtRel() );

    //==== Muon sim-matching variable
    //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2#SimHit_Muon_matching_since_CMSSW
    muon_simType.push_back( imuon.simType() );
    muon_simExtType.push_back( imuon.simExtType() );
    muon_simFlavour.push_back( imuon.simFlavour() );
    muon_simHeaviestMotherFlavour.push_back( imuon.simHeaviestMotherFlavour() );
    muon_simPdgId.push_back( imuon.simPdgId() );
    muon_simMotherPdgId.push_back( imuon.simMotherPdgId() );
    //muon_simMatchQuality.push_back( imuon.simMatchQuality() );  // TODO not supported in CMSSW_10_2_10

    muon_ishighpt.push_back(muon::isHighPtMuon(imuon, vtx));
    
    //==== Rochestor

    double this_roccor = 1.;
    double this_roccor_err = 0.;

    //==== https://indico.cern.ch/event/460734/contributions/1131590/attachments/1183866/1715221/rochcor_run2_MuonPOG_110815.pdf
    //==== Recommended pt range is pt>20 GeV
    //==== JH asked to the authors, and they said pt 10~20 GeV MIGHT be okay..
    //==== HOWEVER, few GeVs of muon shuold not use this correction,
    //==== because the authors did not consider the radiations of low pt muons inside detectors
    //==== So for now, if pt<10 GeV, no correction applied
    if(imuon.pt()>10.){

      //==== Data
      if(IsData){
        this_roccor = rc.kScaleDT(imuon.charge(), imuon.pt(), imuon.eta(), imuon.phi(), 0, 0); //data
        this_roccor_err = rc.kScaleDTerror(imuon.charge(), imuon.pt(), imuon.eta(), imuon.phi());
      }
      //==== MC
      else{

        gRandom->SetSeed( evtNum ); // to make the seed always same
        double u = gRandom->Rndm();

        int n_trackerLayersWithMeasurement = 0;
        if(trackerTrack.isNonnull()) n_trackerLayersWithMeasurement = imuon.innerTrack()->hitPattern().trackerLayersWithMeasurement();


        //==== Using GenPt

        //cout << "[This Reco Muon] pt = " << imuon.pt() << ", eta = " << imuon.eta() << ", phi = " << imuon.phi() << endl;
        //cout << "isTight = " << imuon.isTightMuon(vtx) << endl;
        //cout << "isMediumMuon = " << imuon.isMediumMuon() << endl;
        //cout << "isLoose = " << imuon.isLooseMuon() << endl;
        //cout << "dxy = " << muon_dxyVTX.at(i) << endl;
        //cout << "dz = " << muon_dzVTX.at(i) << endl;

        int counter=0;
        double this_genpt=-999.;
        double mindr = 0.2;
        for( reco::GenParticleCollection::const_iterator genptl = genParticles->begin(); genptl != genParticles->end(); ++genptl, ++counter) {

          //int idx = -1;
          //for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ){
          //  if( genptl->mother()==&(*mit) ){
          //    idx = std::distance(genParticles->begin(),mit);
          //    break;
          //  }
          //}
          //cout << counter << "\t" << genptl->pdgId() << "\t" << idx << "\t" << "gen pt = " << genptl->pt() << ", eta = " << genptl->eta() << ", phi = " << genptl->phi() << endl;

          if( fabs(genptl->pdgId()) != 13 ) continue;
          if( genptl->status() != 1 ) continue;

          double this_dr = reco::deltaR( imuon.eta(), imuon.phi(), genptl->eta(), genptl->phi() );
          if( this_dr < mindr ){
            this_genpt = genptl->pt();
            mindr = this_dr;
          }
        }
        //cout << this_genpt << endl;

        if(this_genpt>0){
          this_roccor     = rc.kSpreadMC     (imuon.charge(), imuon.pt(), imuon.eta(), imuon.phi(), this_genpt, 0, 0);
          this_roccor_err = rc.kSpreadMCerror(imuon.charge(), imuon.pt(), imuon.eta(), imuon.phi(), this_genpt);
        }
        else{
          this_roccor     = rc.kSmearMC     (imuon.charge(), imuon.pt(), imuon.eta(), imuon.phi(), n_trackerLayersWithMeasurement, u, 0, 0);
          this_roccor_err = rc.kSmearMCerror(imuon.charge(), imuon.pt(), imuon.eta(), imuon.phi(), n_trackerLayersWithMeasurement, u);
        }

      }

    } // END IF pt>10 GeV
/*
    cout << "--------------------" << endl;
    cout << "pt = " << imuon.pt() << endl;
    cout << "this_roccor = " << this_roccor << endl;
    cout << "this_roccor_err = " << this_roccor_err << endl;
*/

/*
    //==== Debugging for large correction muons
    if(this_roccor<-20){

      cout << "#### RECO ####" << endl;
      cout << "this_roccor = " << this_roccor << endl;
      cout << "pt = " << imuon.pt() << endl;
      cout << "eta = " << imuon.eta() << endl;
      cout << "phi = " << imuon.phi() << endl;
      cout << "muon_normchi = " << muon_normchi.at(i) << endl;
      cout << "isLooseMuon = " << imuon.isLooseMuon() << endl;
      cout << "#### GEN ####" << endl;
      cout << "Index\tPID\tMotherIndex\tpt\teta\tphi" << endl;
      int counter = 0;
      for( reco::GenParticleCollection::const_iterator genptl = genParticles->begin(); genptl != genParticles->end(); ++genptl, ++counter) {

        int idx = -1;
        for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ){
          if( genptl->mother()==&(*mit) ){
            idx = std::distance(genParticles->begin(),mit);
            break;
          }
        }

        cout << counter << "\t" << genptl->pdgId() << "\t" << idx << "\t" << "gen pt = " << genptl->pt() << ", eta = " << genptl->eta() << ", phi = " << genptl->phi() << endl;

      }

    }
*/

    muon_roch_sf.push_back( this_roccor );
    muon_roch_sf_up.push_back( this_roccor+this_roccor_err );

    // -- HLT object matching -- //
    ULong64_t pathbits=0;
    ULong64_t filterbits=0;
    for(pat::TriggerObjectStandAlone obj : *triggerObject){
      if(deltaR(obj,imuon)<0.3){
	obj.unpackPathNames(trigNames);
	obj.unpackFilterLabels(iEvent, *trigResult);
	if(obj.path("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<0;
	if(obj.path("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*",true,true)) pathbits|=ULong64_t(1)<<1;
	if(obj.path("HLT_IsoMu24_v*",true,true)) pathbits|=ULong64_t(1)<<2;
	if(obj.path("HLT_IsoMu27_v*",true,true)) pathbits|=ULong64_t(1)<<3;
	if(obj.path("HLT_IsoTkMu24_v*",true,true)) pathbits|=ULong64_t(1)<<4;
	if(obj.path("HLT_Mu17_Mu8_SameSign_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<5;
	if(obj.path("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*",true,true)) pathbits|=ULong64_t(1)<<6;
	if(obj.path("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*",true,true)) pathbits|=ULong64_t(1)<<7;
	if(obj.path("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<8;
	if(obj.path("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",true,true)) pathbits|=ULong64_t(1)<<9;
	if(obj.path("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<10;
	if(obj.path("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",true,true)) pathbits|=ULong64_t(1)<<11;
	if(obj.path("HLT_Mu17_TrkIsoVVL_v*",true,true)) pathbits|=ULong64_t(1)<<12;
	if(obj.path("HLT_Mu17_v*",true,true)) pathbits|=ULong64_t(1)<<13;
	if(obj.path("HLT_Mu18_Mu9_SameSign_v*",true,true)) pathbits|=ULong64_t(1)<<14;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<15;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<16;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<17;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<18;
	if(obj.path("HLT_Mu50_v*",true,true)) pathbits|=ULong64_t(1)<<19;
	if(obj.path("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*",true,true)) pathbits|=ULong64_t(1)<<20;
	if(obj.path("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<21;
	if(obj.path("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<22;
	if(obj.path("HLT_Mu8_TrkIsoVVL_v*",true,true)) pathbits|=ULong64_t(1)<<23;
	if(obj.path("HLT_Mu8_v*",true,true)) pathbits|=ULong64_t(1)<<24;
	if(obj.path("HLT_OldMu100_v*",true,true)) pathbits|=ULong64_t(1)<<25;
	if(obj.path("HLT_TkMu100_v*",true,true)) pathbits|=ULong64_t(1)<<26;
	if(obj.path("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<27;
	if(obj.path("HLT_TkMu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",true,true)) pathbits|=ULong64_t(1)<<28;
	if(obj.path("HLT_TkMu50_v*",true,true)) pathbits|=ULong64_t(1)<<29;
	if(obj.path("HLT_TripleMu_10_5_5_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<30;
	if(obj.path("HLT_TripleMu_12_10_5_v*",true,true)) pathbits|=ULong64_t(1)<<31;
	if(obj.path("HLT_TrkMu12_DoubleTrkMu5NoFiltersNoVtx_v*",true,true)) pathbits|=ULong64_t(1)<<32;
	if(obj.path("HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx_v*",true,true)) pathbits|=ULong64_t(1)<<33;
	if(obj.filter("hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9")) filterbits|=ULong64_t(1)<<0;
	if(obj.filter("hltDiMuon178Mass3p8Filtered")) filterbits|=ULong64_t(1)<<1;
	if(obj.filter("hltDiMuon178Mass8Filtered")) filterbits|=ULong64_t(1)<<2;
	if(obj.filter("hltDiMuon178RelTrkIsoFiltered0p4")) filterbits|=ULong64_t(1)<<3;
	if(obj.filter("hltDiMuon178RelTrkIsoFiltered0p4DzFiltered0p2")) filterbits|=ULong64_t(1)<<4;
	if(obj.filter("hltDiMuonGlb17Glb8DzFiltered0p2")) filterbits|=ULong64_t(1)<<5;
	if(obj.filter("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4")) filterbits|=ULong64_t(1)<<6;
	if(obj.filter("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2")) filterbits|=ULong64_t(1)<<7;
	if(obj.filter("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4")) filterbits|=ULong64_t(1)<<8;
	if(obj.filter("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2")) filterbits|=ULong64_t(1)<<9;
	if(obj.filter("hltDiMuonGlbFiltered17TrkFiltered8")) filterbits|=ULong64_t(1)<<10;
	if(obj.filter("hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4")) filterbits|=ULong64_t(1)<<11;
	if(obj.filter("hltDiMuonTrk17Trk8RelTrkIsoFiltered0p4DzFiltered0p2")) filterbits|=ULong64_t(1)<<12;
	if(obj.filter("hltDiTkMuonTkFiltered17TkFiltered8")) filterbits|=ULong64_t(1)<<13;
	if(obj.filter("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105")) filterbits|=ULong64_t(1)<<14;
	if(obj.filter("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5")) filterbits|=ULong64_t(1)<<15;
	if(obj.filter("hltL3crIsoL1sMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p09")) filterbits|=ULong64_t(1)<<16;
	if(obj.filter("hltL3crIsoL1sMu22Or25L1f0L2f10QL3f27QL3trkIsoFiltered0p07")) filterbits|=ULong64_t(1)<<17;
	if(obj.filter("hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07")) filterbits|=ULong64_t(1)<<18;
	if(obj.filter("hltL3fL1DoubleMu155fFiltered17")) filterbits|=ULong64_t(1)<<19;
	if(obj.filter("hltL3fL1DoubleMu155fPreFiltered8")) filterbits|=ULong64_t(1)<<20;
	if(obj.filter("hltL3fL1DoubleMu157fFiltered18")) filterbits|=ULong64_t(1)<<21;
	if(obj.filter("hltL3fL1DoubleMu157fFiltered9")) filterbits|=ULong64_t(1)<<22;
	if(obj.filter("hltL3fL1DoubleMu7EG7f0Filtered9")) filterbits|=ULong64_t(1)<<23;
	if(obj.filter("hltL3fL1Mu6DoubleEG10f0Filtered8")) filterbits|=ULong64_t(1)<<24;
	if(obj.filter("hltL3fL1TripleMu553f0Filtered10105")) filterbits|=ULong64_t(1)<<25;
	if(obj.filter("hltL3fL1TripleMu553f0Filtered1055")) filterbits|=ULong64_t(1)<<26;
	if(obj.filter("hltL3fL1TripleMu553f0PreFiltered555")) filterbits|=ULong64_t(1)<<27;
	if(obj.filter("hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17")) filterbits|=ULong64_t(1)<<28;
	if(obj.filter("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17")) filterbits|=ULong64_t(1)<<29;
	if(obj.filter("hltL3fL1sDoubleMu114TkFiltered17Q")) filterbits|=ULong64_t(1)<<30;
	if(obj.filter("hltL3fL1sMu10lqL1f0L2f10L3Filtered17")) filterbits|=ULong64_t(1)<<31;
	if(obj.filter("hltL3fL1sMu15DQlqL1f0L2f10L3Filtered17")) filterbits|=ULong64_t(1)<<32;
	if(obj.filter("hltL3fL1sMu1lqL1f0L2f10L3Filtered17TkIsoFiltered0p4")) filterbits|=ULong64_t(1)<<33;
	if(obj.filter("hltL3fL1sMu22L1f0Tkf24QL3trkIsoFiltered0p09")) filterbits|=ULong64_t(1)<<34;
	if(obj.filter("hltL3fL1sMu5L1f0L2f5L3Filtered8")) filterbits|=ULong64_t(1)<<35;
	if(obj.filter("hltL3fL1sMu5L1f0L2f5L3Filtered8TkIsoFiltered0p4")) filterbits|=ULong64_t(1)<<36;
	if(obj.filter("hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8")) filterbits|=ULong64_t(1)<<37;
	if(obj.filter("hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8")) filterbits|=ULong64_t(1)<<38;
	if(obj.filter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter")) filterbits|=ULong64_t(1)<<39;
	if(obj.filter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")) filterbits|=ULong64_t(1)<<40;
	if(obj.filter("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23")) filterbits|=ULong64_t(1)<<41;
	if(obj.filter("hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8")) filterbits|=ULong64_t(1)<<42;
	if(obj.filter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter")) filterbits|=ULong64_t(1)<<43;
	if(obj.filter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8")) filterbits|=ULong64_t(1)<<44;
	if(obj.filter("hltMu9Ele9DZFilter")) filterbits|=ULong64_t(1)<<45;
	if(obj.filter("hltTripleTrkMuFiltered5NoVtx")) filterbits|=ULong64_t(1)<<46;
      }
    }
    muon_pathbits.push_back(pathbits);
    muon_filterbits.push_back(filterbits);
  } // -- End of imuon iteration -- //
  
}

//////////////////////////////
// -- Get Electrons info -- //
//////////////////////////////
void SKFlatMaker::fillElectrons(const edm::Event &iEvent, const edm::EventSetup& iSetup)
{
  
  if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] called" << endl;
  // -- BeamSpot -- //
  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(BeamSpotToken, beamSpotHandle);
  reco::BeamSpot beamSpot = (*beamSpotHandle);
  
  // -- Primary vertex -- //
  edm::Handle<reco::VertexCollection> pvHandle;
  iEvent.getByToken(PrimaryVertexToken, pvHandle);
  const reco::VertexCollection vtx = *(pvHandle.product());
  
  // -- electron -- //
  edm::Handle< edm::View<pat::Electron> > ElecHandle;
  iEvent.getByToken(ElectronToken, ElecHandle);

  edm::Handle< std::vector<reco::Conversion> > conversions;
  iEvent.getByToken(ConversionsToken, conversions);
  
  // -- B-field for vertex variables (ee, emu) -- //
  ESHandle<MagneticField> B;
  iSetup.get<IdealMagneticFieldRecord>().get(B);
  
  // -- muon for emu vertex -- //
  edm::Handle< std::vector<pat::Muon> > muonHandle;
  iEvent.getByToken(MuonToken, muonHandle);

  //==== HLT Object for lepton matching
  edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerObject;
  iEvent.getByToken(TriggerObjectToken, triggerObject);
  Handle<TriggerResults> trigResult;
  iEvent.getByToken(TriggerToken, trigResult);
  const edm::TriggerNames trigNames = iEvent.triggerNames(*trigResult);

  edm::FileInPath eaConstantsFile(electron_EA_NHandPh_file);
  EffectiveAreas effectiveAreas(eaConstantsFile.fullPath());

  if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] for ElecHandle starts, ElecHandle->size() : " << ElecHandle->size() << endl;
  
  for(int i=0; i< (int)ElecHandle->size(); i++){
    const auto el = ElecHandle->ptrAt(i);
    
    electron_MVAIso.push_back( el -> userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values") );
    electron_MVANoIso.push_back( el -> userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values") );

    if(el->hasUserFloat("ecalTrkEnergyPostCorr")){

      //==== UnCorrected
      electron_EnergyUnCorr.push_back( el->userFloat("ecalTrkEnergyPreCorr") );

      electron_Energy.push_back( el->userFloat("ecalTrkEnergyPostCorr") );
      electron_Energy_Scale_Up.push_back( el->userFloat("energyScaleUp") );
      electron_Energy_Scale_Down.push_back( el->userFloat("energyScaleDown") );
      electron_Energy_Smear_Up.push_back( el->userFloat("energySigmaUp") );
      electron_Energy_Smear_Down.push_back( el->userFloat("energySigmaDown") );

    }
    else{

      //==== UnCorrected
      electron_EnergyUnCorr.push_back( el->energy() );

      electron_Energy.push_back( el->energy() );
      electron_Energy_Scale_Up.push_back( el->energy() );
      electron_Energy_Scale_Down.push_back( el->energy() );
      electron_Energy_Smear_Up.push_back( el->energy() );
      electron_Energy_Smear_Down.push_back( el->energy() );

    }
    electron_eta.push_back( el->eta() );
    electron_phi.push_back( el->phi() );

/*
    //==== Debugging lines for egamma correction
    cout << "==== Electron ====" << endl;
    cout << "el->energy() = " << el->energy() << endl;
    cout << "el->userFloat(\"ecalTrkEnergyPreCorr\") = " << el->userFloat("ecalTrkEnergyPreCorr") << endl;
    cout << "el->userFloat(\"ecalTrkEnergyPostCorr\") = " << el->userFloat("ecalTrkEnergyPostCorr") << endl;
*/

    electron_charge.push_back( el->charge() );
    electron_fbrem.push_back( el->fbrem() );
    electron_eOverP.push_back( el->eSuperClusterOverP() );
    electron_ecalDriven.push_back( el->ecalDrivenSeed() );
    if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] basic info. input" << endl;
    
    
    // -- Information from SuperCluster -- //
    electron_scEnergy.push_back( el->superCluster()->energy() );
    electron_scPreEnergy.push_back( el->superCluster()->preshowerEnergy() );
    electron_scRawEnergy.push_back( el->superCluster()->rawEnergy() );
    electron_scEta.push_back( el->superCluster()->eta() );
    electron_scPhi.push_back( el->superCluster()->phi() );
    double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
    double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());
    electron_scEt.push_back( el->superCluster()->energy()*(Rt/R) );
    
    if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] supercluster" << endl;

    
    // -- Information from ECAL  -- //
    electron_sigmaIEtaIEta.push_back( el->full5x5_sigmaIetaIeta() );
    electron_Full5x5_SigmaIEtaIEta.push_back( el->full5x5_sigmaIetaIeta() );
    electron_E15.push_back( el->e1x5() );
    electron_E25.push_back( el->e2x5Max() );
    electron_E55.push_back( el->e5x5() );

    //==== HEEP
    const double e5x5 = el->full5x5_e5x5();
    const double e2x5OverE5x5 = e5x5!=0 ? el->full5x5_e2x5Max()/e5x5 : 0;
    const double e1x5OverE5x5 = e5x5!=0 ? el->full5x5_e1x5()/e5x5 : 0;
    electron_e2x5OverE5x5.push_back( e2x5OverE5x5 );
    electron_e1x5OverE5x5.push_back( e1x5OverE5x5 );

    // electron_HoverE.push_back( el->hcalOverEcal() );
    electron_HoverE.push_back( el->hadronicOverEm() ); // -- https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleHadronicOverEMCut.cc#L40 -- //
    electron_etaWidth.push_back( el->superCluster()->etaWidth() );
    electron_phiWidth.push_back( el->superCluster()->phiWidth() );
    electron_r9.push_back( el->r9() );
    
    if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] ECAL" << endl;
    
    
    // -- Information from ECAL & Track -- //
    electron_dEtaIn.push_back( el->deltaEtaSuperClusterTrackAtVtx() );
    electron_dPhiIn.push_back( el->deltaPhiSuperClusterTrackAtVtx() );
    
    // -- https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDEtaInSeedCut.cc#L30-L33 -- //
    electron_dEtaInSeed.push_back( el->superCluster().isNonnull() && el->superCluster()->seed().isNonnull() ?
el->deltaEtaSuperClusterTrackAtVtx() - el->superCluster()->eta() + el->superCluster()->seed()->eta() : std::numeric_limits<float>::max() );
    
    // -- |1/E-1/p| = |1/E - EoverPinner/E| is computed below. The if protects against ecalEnergy == inf or zero -- //
    if( el->ecalEnergy() == 0 ) electron_InvEminusInvP.push_back( 1e30 );
    else if(  !std::isfinite( el->ecalEnergy() )  ) electron_InvEminusInvP.push_back( 1e30 );
    else electron_InvEminusInvP.push_back( fabs( 1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() ) );
    
    if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] ECAL & track" << endl;
    
    // -- Isolation -- //
    double pfCharged = el->pfIsolationVariables().sumChargedHadronPt;
    double pfNeutral = el->pfIsolationVariables().sumNeutralHadronEt;
    double pfPhoton = el->pfIsolationVariables().sumPhotonEt;
    double pfChargedFromPU = el->pfIsolationVariables().sumPUPt;
    electron_chIso03.push_back( pfCharged );
    electron_nhIso03.push_back( pfNeutral );
    electron_phIso03.push_back( pfPhoton );
    electron_puChIso03.push_back( pfChargedFromPU );
    electron_RelPFIso_dBeta.push_back( (pfCharged + max<float>( 0.0, pfNeutral + pfPhoton - 0.5 * pfChargedFromPU))/(el->pt()) );
    
    float abseta = fabs(el->superCluster()->eta());
    float eA = effectiveAreas.getEffectiveArea(abseta);
    electron_RelPFIso_Rho.push_back( (pfCharged + max<float>( 0.0, pfNeutral + pfPhoton - Rho * eA))/(el->pt()) );

    //==== MiniIso
/*
    cout << "=======================" << endl;
    cout << "Pt = " << el->pt() << endl;
    cout << "isEE = " << el->isEE() << endl;
    cout << "---- R03 ----" << endl;
    cout << "CH = " << el->pfIsolationVariables().sumChargedHadronPt << endl;
    cout << "NH = " << el->pfIsolationVariables().sumNeutralHadronEt << endl;
    cout << "Ph = " << el->pfIsolationVariables().sumPhotonEt << endl;
    cout << "PU = " << el->pfIsolationVariables().sumPUPt << endl;
    cout << "---- MiniAOD MiniIso ----" << endl;
    cout << "CH = " << el->miniPFIsolation().chargedHadronIso() << endl;
    cout << "HN = " << el->miniPFIsolation().neutralHadronIso() << endl;
    cout << "Ph = " << el->miniPFIsolation().photonIso() << endl;
    cout << "PU = " << el->miniPFIsolation().puChargedHadronIso() << endl;
*/

    electron_chMiniIso.push_back( el->miniPFIsolation().chargedHadronIso() );
    electron_nhMiniIso.push_back( el->miniPFIsolation().neutralHadronIso() );
    electron_phMiniIso.push_back( el->miniPFIsolation().photonIso() );
    electron_puChMiniIso.push_back( el->miniPFIsolation().puChargedHadronIso() );

    // For UL heepTrkPtIso is not saved as userFloat.
    // https://github.com/cms-egamma/EgammaPostRecoTools/blob/3257bda882e7b34be4adbe2cd1e83dac8b533897/python/EgammaPostRecoTools.py#L452
    // https://hypernews.cern.ch/HyperNews/CMS/get/egamma/2268.html
    //electron_trackIso.push_back( el->userFloat("heepTrkPtIso") );
    electron_trackIso.push_back( el->dr03TkSumPtHEEP() );

    electron_dr03EcalRecHitSumEt.push_back( el->dr03EcalRecHitSumEt() );
    electron_dr03HcalDepth1TowerSumEt.push_back( el->dr03HcalDepth1TowerSumEt() );
    electron_dr03HcalTowerSumEt.push_back( el->dr03HcalTowerSumEt() );
    electron_dr03TkSumPt.push_back( el->dr03TkSumPt() );
    electron_ecalPFClusterIso.push_back( el->ecalPFClusterIso() );
    electron_hcalPFClusterIso.push_back( el->hcalPFClusterIso() );

    // cout << "##### fillElectrons: Before elecTrk #####" << endl;
    
    if(theDebugLevel) cout << "[SKFlatMaker::fillElectrons] Isolation" << endl;
    
    // -- Track - Impact Parameter, Conversion rejection, Converted -- //
    reco::GsfTrackRef elecTrk = el->gsfTrack();

    electron_passConversionVeto.push_back( el->passConversionVeto() );
    electron_isGsfCtfScPixChargeConsistent.push_back( el->isGsfCtfScPixChargeConsistent() );
    electron_isGsfScPixChargeConsistent.push_back( el->isGsfScPixChargeConsistent() );
    electron_isGsfCtfChargeConsistent.push_back( el->isGsfCtfChargeConsistent() );

    //==== https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleMissingHitsCut.cc#L34-L41
    constexpr auto missingHitType = reco::HitPattern::MISSING_INNER_HITS;
    electron_mHits.push_back( elecTrk->hitPattern().numberOfLostHits(missingHitType) );

    electron_sigdxy.push_back( elecTrk->dxy() / elecTrk->dxyError() );
    electron_dxy.push_back( elecTrk->dxy() );
    electron_dz.push_back( elecTrk->dz() );

    electron_dxyBS.push_back( elecTrk->dxy(beamSpot.position()) );
    electron_dzBS.push_back( elecTrk->dz(beamSpot.position()) );

    if( !pvHandle->empty() && !pvHandle->front().isFake() ){

/*
      electron_dxyVTX.push_back( elecTrk->dxy(vtx.position()) );
      electron_dzVTX.push_back( elecTrk->dz(vtx.position()) );
*/

      electron_dxyVTX.push_back( el->dB(pat::Electron::PV2D) );
      electron_dxyerrVTX.push_back( el->edB(pat::Electron::PV2D) );
      electron_dzVTX.push_back( el->dB(pat::Electron::PVDZ) );
      electron_dzerrVTX.push_back( el->edB(pat::Electron::PVDZ) );
      electron_3DIPVTX.push_back( el->dB(pat::Electron::PV3D) );
      electron_3DIPerrVTX.push_back( el->edB(pat::Electron::PV3D) );

    }
    else{

      electron_dxyVTX.push_back( -999 );
      electron_dxyerrVTX.push_back( -999 );
      electron_dzVTX.push_back( -999 );
      electron_dzerrVTX.push_back( -999 );
      electron_3DIPVTX.push_back( -999 );
      electron_3DIPerrVTX.push_back( -999 );

    }

    if( elecTrk.isNonnull() ){
      electron_gsfpt.push_back( elecTrk->pt() );
      electron_gsfEta.push_back( elecTrk->eta() );
      electron_gsfPhi.push_back( elecTrk->phi() );
      electron_gsfCharge.push_back( elecTrk->charge() );
    }
    else{
      electron_gsfpt.push_back( -999 );
      electron_gsfEta.push_back( -999 );
      electron_gsfPhi.push_back( -999 );
      electron_gsfCharge.push_back( -999 );
    }
      
    //==== ID Booleans

    //==== Cut Based
    //==== 94x-V2 : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Recipe_for_regular_users_formats
    //==== MVA Based
    //==== Fall17 : https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#VID_based_recipe_provides_pass_f
    //==== HEEP ID
    //==== V7 : https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2#Configuring_and_Running_VID_in_a

    unsigned int IDBit = 0;
    for(unsigned int it_ID=0; it_ID<electron_IDtoSave.size(); it_ID++){
      if(el->electronID(electron_IDtoSave.at(it_ID))){
        IDBit |= (1 << it_ID);
      }
      else{
        IDBit &= ~(1 << it_ID);
      }

      //==== bits for each cuts
      //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideCMSPhysicsObjectSchoolAACHEN2019EGamma, Exercise 2
      int this_idcutbid = -1;
      if( el->hasUserInt(electron_IDtoSave.at(it_ID)) ){
        this_idcutbid = el->userInt(electron_IDtoSave.at(it_ID));
      }
      electron_IDCutBit.push_back( this_idcutbid );
    }
    electron_IDBit.push_back( IDBit );
/*
    cout << "Electron ID Bit = " << IDBit << endl;
    cout << "--> CB veto = " << el->electronID("cutBasedElectronID-Fall17-94X-V2-veto") << endl;
    cout << "--> CB loose = " << el->electronID("cutBasedElectronID-Fall17-94X-V2-loose") << endl;
    cout << "--> CB medium = " << el->electronID("cutBasedElectronID-Fall17-94X-V2-medium") << endl;
    cout << "--> CB tight = " << el->electronID("cutBasedElectronID-Fall17-94X-V2-tight") << endl;
    cout << "--> MVA noIso WP80 = " << el->electronID("mvaEleID-Fall17-noIso-V1-wp80") << endl;
    cout << "--> MVA noIso WP90 = " << el->electronID("mvaEleID-Fall17-noIso-V1-wp90") << endl;
    cout << "--> MVA Iso WP80 = " << el->electronID("mvaEleID-Fall17-iso-V1-wp80") << endl;
    cout << "--> MVA Iso WP90 = " << el->electronID("mvaEleID-Fall17-iso-V1-wp90") << endl;
    cout << "--> Heepv70 = " << el->electronID("heepElectronID-HEEPV70") << endl;
*/
/*
    //==== checking pogboolean and cut-by-hand
    if(isPassLoose){

      double this_hovere = el->hadronicOverEm();
      double corrE = el->userFloat("ecalTrkEnergyPostCorr");
      double uncorrE = el->userFloat("ecalTrkEnergyPreCorr");
      double scE = el->superCluster()->energy();

      if(abseta<1.479){

        double this_loosecut_corrE   = 0.05 + 1.12/corrE   + 0.0368*Rho/corrE;
        double this_loosecut_uncorrE = 0.05 + 1.12/uncorrE + 0.0368*Rho/uncorrE;

        if( !(this_hovere < this_loosecut_uncorrE) ){
          cout << "---- |scEta| = " << abseta << endl;
          cout << "hadronicOverEm() = " << this_hovere << endl;
          cout << "this_loosecut_corrE = " << this_loosecut_corrE << endl;
          cout << "this_loosecut_uncorrE = " << this_loosecut_uncorrE << endl;
        }

      }
      else{

        double this_loosecut_corrE   = 0.0414 + 0.5/corrE   + 0.201*Rho/corrE;
        double this_loosecut_uncorrE = 0.0414 + 0.5/uncorrE + 0.201*Rho/uncorrE;
        double this_loosecut_scE     = 0.0414 + 0.5/scE   + 0.201*Rho/scE;

        if( !(this_hovere < this_loosecut_uncorrE) ){
          cout << "fabs(el->superCluster()->eta()) = " << fabs(el->superCluster()->eta()) << endl;
          cout << "fabs(el->eta()) = " << fabs(el->eta()) << endl;
          cout << "el->electronID(\"cutBasedElectronID-Fall17-94X-V1-loose\") = " << el -> electronID("cutBasedElectronID-Fall17-94X-V1-loose") << endl;
          cout << "el->hadronicOverEm() = " << el->hadronicOverEm() << endl;
          cout << "el->userFloat(\"ecalTrkEnergyPreCorr\") = " << el->userFloat("ecalTrkEnergyPreCorr") << endl;
          cout << "el->userFloat(\"ecalTrkEnergyPostCorr\") = " << el->userFloat("ecalTrkEnergyPostCorr") << endl;
          cout << "el->superCluster()->energy() = " << el->superCluster()->energy() << endl;
          cout << "Rho = " << Rho << endl;
          cout << "this_loosecut_corrE = " << this_loosecut_corrE << endl;
          cout << "this_loosecut_uncorrE = " << this_loosecut_uncorrE << endl;
          cout << "this_loosecut_scE = " << this_loosecut_scE << endl;
        }

      }
    }
*/
    // -- HLT object matching -- //
    ULong64_t pathbits=0;
    ULong64_t filterbits=0;
    for(pat::TriggerObjectStandAlone obj : *triggerObject){
      if(deltaR(obj,*el)<0.3){
	obj.unpackPathNames(trigNames);
	obj.unpackFilterLabels(iEvent, *trigResult);
	if(obj.path("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<0;
	if(obj.path("HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v*",true,true)) pathbits|=ULong64_t(1)<<1;
	if(obj.path("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",true,true)) pathbits|=ULong64_t(1)<<2;
	if(obj.path("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*",true,true)) pathbits|=ULong64_t(1)<<3;
	if(obj.path("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",true,true)) pathbits|=ULong64_t(1)<<4;
	if(obj.path("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v*",true,true)) pathbits|=ULong64_t(1)<<5;
	if(obj.path("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",true,true)) pathbits|=ULong64_t(1)<<6;
	if(obj.path("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<7;
	if(obj.path("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<8;
	if(obj.path("HLT_Ele27_WPTight_Gsf_v*",true,true)) pathbits|=ULong64_t(1)<<9;
	if(obj.path("HLT_Ele27_eta2p1_WPTight_Gsf_v*",true,true)) pathbits|=ULong64_t(1)<<10;
	if(obj.path("HLT_Ele28_WPTight_Gsf_v*",true,true)) pathbits|=ULong64_t(1)<<11;
	if(obj.path("HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*",true,true)) pathbits|=ULong64_t(1)<<12;
	if(obj.path("HLT_Ele32_WPTight_Gsf_v*",true,true)) pathbits|=ULong64_t(1)<<13;
	if(obj.path("HLT_Ele35_WPTight_Gsf_v*",true,true)) pathbits|=ULong64_t(1)<<14;
	if(obj.path("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30_v*",true,true)) pathbits|=ULong64_t(1)<<15;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<16;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<17;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<18;
	if(obj.path("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<19;
	if(obj.path("HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v*",true,true)) pathbits|=ULong64_t(1)<<20;
	if(obj.path("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",true,true)) pathbits|=ULong64_t(1)<<21;
	if(obj.path("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",true,true)) pathbits|=ULong64_t(1)<<22;
	if(obj.filter("hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter")) filterbits|=ULong64_t(1)<<0;
	if(obj.filter("hltEGL1SingleEGOrFilter")) filterbits|=ULong64_t(1)<<1;
	if(obj.filter("hltEle12CaloIdLTrackIdLIsoVLTrackIsoFilter")) filterbits|=ULong64_t(1)<<2;
	if(obj.filter("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter")) filterbits|=ULong64_t(1)<<3;
	if(obj.filter("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter")) filterbits|=ULong64_t(1)<<4;
	if(obj.filter("hltEle17CaloIdLTrackIdLIsoVLTrackIsoFilter")) filterbits|=ULong64_t(1)<<5;
	if(obj.filter("hltEle17CaloIdMGsfTrackIdMDphiFilter")) filterbits|=ULong64_t(1)<<6;
	if(obj.filter("hltEle23CaloIdLTrackIdLIsoVLJet30TrackIsoFilter")) filterbits|=ULong64_t(1)<<7;
	if(obj.filter("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter")) filterbits|=ULong64_t(1)<<8;
	if(obj.filter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter")) filterbits|=ULong64_t(1)<<9;
	if(obj.filter("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter")) filterbits|=ULong64_t(1)<<10;
	if(obj.filter("hltEle27WPTightGsfTrackIsoFilter")) filterbits|=ULong64_t(1)<<11;
	if(obj.filter("hltEle27erWPTightGsfTrackIsoFilter")) filterbits|=ULong64_t(1)<<12;
	if(obj.filter("hltEle28WPTightGsfTrackIsoFilter")) filterbits|=ULong64_t(1)<<13;
	if(obj.filter("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter")) filterbits|=ULong64_t(1)<<14;
	if(obj.filter("hltEle32WPTightGsfTrackIsoFilter")) filterbits|=ULong64_t(1)<<15;
	if(obj.filter("hltEle35noerWPTightGsfTrackIsoFilter")) filterbits|=ULong64_t(1)<<16;
	if(obj.filter("hltEle8CaloIdLTrackIdLIsoVLTrackIsoFilter")) filterbits|=ULong64_t(1)<<17;
	if(obj.filter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLDZFilter")) filterbits|=ULong64_t(1)<<18;
	if(obj.filter("hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")) filterbits|=ULong64_t(1)<<19;
	if(obj.filter("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")) filterbits|=ULong64_t(1)<<20;
	if(obj.filter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLDZFilter")) filterbits|=ULong64_t(1)<<21;
	if(obj.filter("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter")) filterbits|=ULong64_t(1)<<22;
	if(obj.filter("hltMu9Ele9DZFilter")) filterbits|=ULong64_t(1)<<23;
      }
    }
    electron_pathbits.push_back(pathbits);
    electron_filterbits.push_back(filterbits);

  } // -- end of for(int i=0; i< (int)ElecHandle->size(); i++): 1st electron iteration -- //
  
  // cout << "##### End of fillElectrons #####" << endl;
}

////////////////////////
// -- Get LHE info -- //
////////////////////////
void SKFlatMaker::fillLHEInfo(const edm::Event &iEvent)
{
  Handle<LHEEventProduct> LHEInfo;
  iEvent.getByToken(LHEEventProductToken, LHEInfo);
  if(!LHEInfo.isValid()) iEvent.getByToken(LHEEventProductSourceToken, LHEInfo);
  if(!LHEInfo.isValid()) return;

  //==== LHE object
  const lhef::HEPEUP& lheEvent = LHEInfo->hepeup();
  std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
  for(size_t idxParticle = 0; idxParticle<lheParticles.size(); ++idxParticle){

    LHE_Px.push_back( lheParticles[idxParticle][0] );
    LHE_Py.push_back( lheParticles[idxParticle][1] );
    LHE_Pz.push_back( lheParticles[idxParticle][2] );
    LHE_E.push_back( lheParticles[idxParticle][3] );

    LHE_Status.push_back( lheEvent.ISTUP[idxParticle] );
    LHE_ID.push_back( lheEvent.IDUP[idxParticle] );

  }

  // -- PDf weights for theoretical uncertainties: scale, PDF replica and alphaS variation -- //
  // -- ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW -- //
  int nWeight = (int)LHEInfo->weights().size();

  //==== https://indico.cern.ch/event/690726/contributions/2890915/attachments/1598277/2532794/2017_MC.pdf
  //==== Check which convention

  //==== Save id-weight as a map
  map<int,double> map_id_to_weight;
  for(int i=0; i<nWeight; i++){
    int this_id = stoi(LHEInfo->weights()[i].id.c_str());
    map_id_to_weight[this_id] = LHEInfo->weights()[i].wgt;
    if(theDebugLevel) cout << "[SKFlatMaker::fillLHEInfo] map_id_to_weight["<<this_id<<"] = " << map_id_to_weight[this_id] << endl;
  }

  double nominalweight= LHEInfo->originalXWGTUP();
  if(WeightMap["Scale"].size()) nominalweight=map_id_to_weight[WeightMap["Scale"].at(0)];
  if(theDebugLevel){
    if(!WeightMap["Scale"].size()) cout << "[SKFlatMaker::fillLHEInfo] central = LHEInfo->originalXWGTUP()" << endl;
    else cout << "[SKFlatMaker::fillLHEInfo] central = " << WeightMap["Scale"].at(0) << endl;
    cout << "[SKFlatMaker::fillLHEInfo] nominalweight = " << nominalweight << endl;
  }

  //==============================
  //==== weights
  //==============================

  for(auto [wname,ids]:WeightMap){
    if(wname=="AlphaS"){
      for(unsigned int i=0;i<ids.size();i++){
	int id=ids[i];
	if(theDebugLevel) cout << "[SKFlatMaker::fillLHEInfo] weights for "+wname+"; adding id = " << id << endl;
	weight_[wname].push_back( 1 + (map_id_to_weight[id]-nominalweight)/nominalweight * PDFAlphaSScaleValue_.at(i));
      }      
    }      
    else if(wname=="PSSyst") continue;
    else{
      for(int id:ids){
	if(theDebugLevel) cout << "[SKFlatMaker::fillLHEInfo] weights for "+wname+"; adding id = " << id << endl;
	weight_[wname].push_back( map_id_to_weight[id]/nominalweight );
      }
    }
  }

  if(theDebugLevel){
    cout << "[SKFlatMaker::fillLHEInfo] LHEInfo->originalXWGTUP() = " << LHEInfo->originalXWGTUP() << endl;
    for(auto [wname,weights]:weight_){
      if(wname=="PSSyst") continue;
      cout << "[SKFlatMaker::fillLHEInfo] [weight_["+wname+"]]" << endl;
      for(unsigned int i=0;i<weights.size();i++) cout << "[SKFlatMaker::fillLHEInfo] " << weights.at(i) << endl;
    }
  }
}

////////////////////////
// -- Get GEN info -- //
////////////////////////
void SKFlatMaker::fillGENInfo(const edm::Event &iEvent)
{
  //cout << "fill pdf info" << endl;
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(mcLabel_,genParticles);
  
  
  int counter=0;
  for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it , ++counter) {      

    if(!theKeepAllGen && counter > 30) continue;
    
    gen_PID.push_back( it->pdgId() );
    gen_pt.push_back( it->pt() );
    gen_mass.push_back( it->mass() );
    gen_charge.push_back( it->charge() );
    gen_eta.push_back( it->eta() );
    gen_phi.push_back( it->phi() );
    gen_status.push_back( it->status() );
    
    //Flags (Ref: https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf)
    gen_isPrompt.push_back( it->statusFlags().isPrompt() ); //not from hadron, muon or tau decay 
    gen_isPromptFinalState.push_back( it->isPromptFinalState() ); //isPrompt && final state (status==1)
    gen_isTauDecayProduct.push_back( it->statusFlags().isTauDecayProduct() ); //is directly or indirectly from a tau decay
    gen_isPromptTauDecayProduct.push_back( it->statusFlags().isPromptTauDecayProduct() ); //is directly or indirectly from a tau decay, where the tau did not come from a hadron decay
    gen_isDirectPromptTauDecayProductFinalState.push_back( it->isDirectPromptTauDecayProductFinalState() ); // is the direct decay product from a tau decay (ie no intermediate hadron), where the tau did not come from a hadron decay && final state
    gen_isHardProcess.push_back( it->isHardProcess() );
    gen_isLastCopy.push_back( it->isLastCopy() );
    gen_isLastCopyBeforeFSR.push_back( it->isLastCopyBeforeFSR() );
    gen_isPromptDecayed.push_back( it->isPromptDecayed() );
    gen_isDecayedLeptonHadron.push_back( it->statusFlags().isDecayedLeptonHadron() );
    gen_fromHardProcessBeforeFSR.push_back( it->fromHardProcessBeforeFSR() );
    gen_fromHardProcessDecayed.push_back( it->fromHardProcessDecayed() );
    gen_fromHardProcessFinalState.push_back( it->fromHardProcessFinalState() );
    gen_isMostlyLikePythia6Status3.push_back( it->fromHardProcessBeforeFSR() );
    
    int idx = -1;
    for( reco::GenParticleCollection::const_iterator mit = genParticles->begin(); mit != genParticles->end(); ++mit ){
      if( it->mother()==&(*mit) ){
        idx = std::distance(genParticles->begin(),mit);
        break;
      }
    }
    
    gen_mother_index.push_back( idx );
    
  }
   
  edm::Handle<GenEventInfoProduct> genEvtInfo;
  iEvent.getByToken(GenEventInfoToken, genEvtInfo);
  gen_weight = genEvtInfo->weight();

  //PS Syst ref.: https://twiki.cern.ch/twiki/bin/viewauth/CMS/HowToPDF#Parton_shower_weights
  for(auto [wname,ids]:WeightMap){
    if(wname!="PSSyst") continue;
    for(int id:ids){
      if(theDebugLevel){
        if(id >(int) genEvtInfo->weights().size()-1){ printf("Error: PSSyst id:%d > PSweight vec. size(%d)-1\n", id, (int) genEvtInfo->weights().size()); continue; }
        printf("[SKFlatMaker::fillGENInfo] weights for %s, id = %d, weight: %.3e\n", wname.Data(), id, genEvtInfo->weights().at(id)/gen_weight);
      }
      weight_[wname].push_back( genEvtInfo->weights().at(id)/gen_weight );
    }
  }
  genWeight_Q = genEvtInfo->pdf()->scalePDF;
  genWeight_X1 = genEvtInfo->pdf()->x.first;
  genWeight_X2 = genEvtInfo->pdf()->x.second;
  genWeight_id1 = genEvtInfo->pdf()->id.first;
  genWeight_id2 = genEvtInfo->pdf()->id.second;
  genWeight_alphaQCD = genEvtInfo->alphaQCD();
  genWeight_alphaQED = genEvtInfo->alphaQED();
  
}


/////////////////////////
// Get Photons info -- // 
/////////////////////////
void SKFlatMaker::fillPhotons(const edm::Event &iEvent)
{
  
  edm::Handle< edm::View<pat::Photon> > PhotonHandle;
  iEvent.getByToken(PhotonToken, PhotonHandle);
  
  EffectiveAreas photon_EA_CH( photon_EA_CH_file.fullPath() );
  EffectiveAreas photon_EA_HN( photon_EA_HN_file.fullPath() );
  EffectiveAreas photon_EA_Ph( photon_EA_Ph_file.fullPath() );

  for(size_t i=0; i< PhotonHandle->size(); ++i){
    const auto pho = PhotonHandle->ptrAt(i);

    if( pho -> hasUserFloat( "ecalEnergyPostCorr" ) ){
      photon_Energy.push_back( pho -> userFloat("ecalEnergyPostCorr") );
      photon_EnergyUnCorr.push_back( pho -> userFloat("ecalEnergyPreCorr") );
    }
    else{
      photon_Energy.push_back( pho -> energy() );
      photon_EnergyUnCorr.push_back( pho -> energy() );
    }
    photon_eta.push_back( pho->eta() );
    photon_phi.push_back( pho->phi() );
    
    photon_scEta.push_back( pho->superCluster()->eta() );
    photon_scPhi.push_back( pho->superCluster()->phi() );
    
    photon_HoverE.push_back( pho->hadTowOverEm() );
    photon_hasPixelSeed.push_back( pho->hasPixelSeed() );
    //photon_Full5x5_SigmaIEtaIEta.push_back( (*full5x5SigmaIEtaIEtaMap)[ pho ] );
    photon_Full5x5_SigmaIEtaIEta.push_back( pho -> full5x5_sigmaIetaIeta() );
    
    
    float chIso = pho -> chargedHadronIso();
    float nhIso = pho -> neutralHadronIso();
    float phIso = pho -> photonIso();
    
    photon_ChIso.push_back( chIso );
    photon_NhIso.push_back( nhIso );
    photon_PhIso.push_back( phIso );
    
    float abseta = fabs( pho->superCluster()->eta());
    photon_ChIsoWithEA.push_back( std::max( float(0.0), chIso - Rho*photon_EA_CH.getEffectiveArea(abseta) ) );
    photon_NhIsoWithEA.push_back( std::max( float(0.0), nhIso - Rho*photon_EA_HN.getEffectiveArea(abseta) ) );
    photon_PhIsoWithEA.push_back( std::max( float(0.0), phIso - Rho*photon_EA_Ph.getEffectiveArea(abseta) ) );
    
    bool isPassLoose  = pho -> photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
    bool isPassMedium  = pho -> photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
    bool isPassTight  = pho  -> photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
    bool isPassMVA_WP80 = pho -> photonID("mvaPhoID-RunIIFall17-v2-wp80");
    bool isPassMVA_WP90 = pho -> photonID("mvaPhoID-RunIIFall17-v2-wp90");
  
    photon_passMVAID_WP80.push_back( isPassMVA_WP80 );
    photon_passMVAID_WP90.push_back( isPassMVA_WP90 );
    photon_passLooseID.push_back( isPassLoose );
    photon_passMediumID.push_back( isPassMedium );
    photon_passTightID.push_back( isPassTight );




  }
  
}


/////////////////////////
// -- Get METs info -- // 
/////////////////////////
void SKFlatMaker::fillMET(const edm::Event &iEvent)
{
  edm::Handle< std::vector<pat::MET> > metHandle;
  iEvent.getByToken(MetToken,metHandle);
  
  if( (metHandle->size() > 1) && (theDebugLevel > 0)) cout << "[SKFlatMaker::fillMET] # of METs = " << metHandle->size() << endl;
  
  pfMET_pt = metHandle->front().uncorPt();
  pfMET_phi = metHandle->front().uncorPhi();
  pfMET_SumEt = metHandle->front().uncorSumEt();
  
  pfMET_Type1_pt = metHandle->front().pt();
  pfMET_Type1_phi = metHandle->front().phi();
  pfMET_Type1_SumEt = metHandle->front().sumEt();
  
  pfMET_Type1_PhiCor_pt = metHandle->front().corPt(pat::MET::Type1XY);
  pfMET_Type1_PhiCor_phi = metHandle->front().corPhi(pat::MET::Type1XY);
  pfMET_Type1_PhiCor_SumEt = metHandle->front().corSumEt(pat::MET::Type1XY);


  //cout << "[MET Check]" << endl;
  //cout << "Uncorrected\t" << pfMET_pt << "\t" << pfMET_phi << endl;
  //cout << "Cor(Raw)\t\t" <<  metHandle->front().corPt(pat::MET::Raw) << "\t" << metHandle->front().corPhi(pat::MET::Raw) << endl;
  //cout << "Default(Type1)\t"     << pfMET_Type1_pt << "\t" << pfMET_Type1_phi << endl;
  //cout << "Cor(Type1)\t\t\t" <<  metHandle->front().corPt(pat::MET::Type1) << "\t" << metHandle->front().corPhi(pat::MET::Type1) << endl;


  //==== Uncertainties
  //==== https://github.com/cms-sw/cmssw/blob/4dbb008c8f4473dc9beb26171d7dde863880d02e/DataFormats/PatCandidates/interface/MET.h#L151-L157

  for(int i=0; i<pat::MET::METUncertaintySize; i++){

    pfMET_pt_shifts.push_back( metHandle->front().shiftedPt(pat::MET::METUncertainty(i), pat::MET::Raw) );
    pfMET_phi_shifts.push_back( metHandle->front().shiftedPhi(pat::MET::METUncertainty(i), pat::MET::Raw) );
    pfMET_SumEt_shifts.push_back( metHandle->front().shiftedSumEt(pat::MET::METUncertainty(i), pat::MET::Raw) );

    pfMET_Type1_pt_shifts.push_back( metHandle->front().shiftedPt(pat::MET::METUncertainty(i), pat::MET::Type1) );
    pfMET_Type1_phi_shifts.push_back( metHandle->front().shiftedPhi(pat::MET::METUncertainty(i), pat::MET::Type1) );
    pfMET_Type1_SumEt_shifts.push_back( metHandle->front().shiftedSumEt(pat::MET::METUncertainty(i), pat::MET::Type1) );

    pfMET_Type1_PhiCor_pt_shifts.push_back( metHandle->front().shiftedPt(pat::MET::METUncertainty(i), pat::MET::Type1XY) );
    pfMET_Type1_PhiCor_phi_shifts.push_back( metHandle->front().shiftedPhi(pat::MET::METUncertainty(i), pat::MET::Type1XY) );
    pfMET_Type1_PhiCor_SumEt_shifts.push_back( metHandle->front().shiftedSumEt(pat::MET::METUncertainty(i), pat::MET::Type1XY) );

  }


}

/////////////////////////
// -- Get Jets info -- // 
/////////////////////////
void SKFlatMaker::fillJet(const edm::Event &iEvent)
{

  // edm::Handle<edm::View<pat::Jet> > jetHandle;
  edm::Handle< std::vector<pat::Jet> > jetHandle;
  iEvent.getByToken(JetToken,jetHandle);
  
  if( jetHandle->size() > 0 && theDebugLevel > 0) 
    cout << "[SKFlatMaker::fillJet] # of Jets = " << jetHandle->size() << endl;
  
  if(jetHandle->size() == 0) return;

  const JetCollection& ForJER_jets = *jetHandle;
  unsigned int runNum_uint = static_cast <unsigned int> (iEvent.id().run());
  unsigned int lumiNum_uint = static_cast <unsigned int> (iEvent.id().luminosityBlock());
  unsigned int evNum_uint = static_cast <unsigned int> (iEvent.id().event());
  unsigned int jet0eta = uint32_t(ForJER_jets.empty() ? 0 : ForJER_jets[0].eta()/0.01);
  int m_nomVar=1;
  std::uint32_t seed = jet0eta + m_nomVar + (lumiNum_uint<<10) + (runNum_uint<<20) + evNum_uint;
  m_random_generator.seed(seed);

  for (vector<pat::Jet>::const_iterator jets_iter = jetHandle->begin(); jets_iter != jetHandle->end(); ++jets_iter){

    jet_pt.push_back( jets_iter->pt() );
    jet_eta.push_back( jets_iter->eta() );
    jet_phi.push_back( jets_iter->phi() );
    jet_charge.push_back( jets_iter->jetCharge() );
    jet_area.push_back( jets_iter->jetArea() );
    jet_partonFlavour.push_back( jets_iter->partonFlavour() );
    jet_hadronFlavour.push_back( jets_iter->hadronFlavour() );

    //=========================================================================
    //==== Taggers
    //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    //=========================================================================

    //=== CSVv2

    jet_CSVv2.push_back( jets_iter->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

    //==== DeepCSV (B)
/*
    cout << "==============================" << endl;
    cout << "---- BvsAll ----" << endl;
    cout << "methodA = " << jets_iter->bDiscriminator("pfDeepCSVJetTags:probb")+jets_iter->bDiscriminator("pfDeepCSVJetTags:probbb") << endl;
    cout << "methodB = " << jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") << endl;
*/
    //==== methodA
    //jet_DeepCSV.push_back( jets_iter->bDiscriminator("pfDeepCSVJetTags:probb")+jets_iter->bDiscriminator("pfDeepCSVJetTags:probbb") );
    //==== methodB
    jet_DeepCSV.push_back( jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") );

    //==== DeepFlavour

    jet_DeepFlavour_b.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probb"));
    jet_DeepFlavour_bb.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    jet_DeepFlavour_lepb.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    jet_DeepFlavour_c.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probc"));
    jet_DeepFlavour_uds.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probuds"));
    jet_DeepFlavour_g.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probg"));

    //==== Old Charm Tagger

    jet_CvsL.push_back( jets_iter->bDiscriminator("pfCombinedCvsLJetTags") );
    jet_CvsB.push_back( jets_iter->bDiscriminator("pfCombinedCvsBJetTags") );

    //==== DeepCSV charm tagger

/*
    cout << "---- CvsL ----" << endl;
    cout << "methodA = " << deepcharm_c/(deepcharm_c+deepcharm_udsg) << endl;
    cout << "methodB = " << jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL") << endl;
    cout << "---- CvsB ----" << endl;
    cout << "methodA = " << deepcharm_c/(deepcharm_c+deepcharm_b+deepcharm_bb) << endl;
    cout << "methodB = " << jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB") << endl;
*/
    //==== methodA : it gives value = 0.5 or 0.333 when these are -1.. 
    //double deepcharm_c = jets_iter->bDiscriminator("pfDeepCSVJetTags:probc");
    //double deepcharm_udsg = jets_iter->bDiscriminator("pfDeepCSVJetTags:probudsg");
    //double deepcharm_b = jets_iter->bDiscriminator("pfDeepCSVJetTags:probb");
    //double deepcharm_bb = jets_iter->bDiscriminator("pfDeepCSVJetTags:probbb");
    //jet_DeepCvsL.push_back( deepcharm_c/(deepcharm_c+deepcharm_udsg) );
    //jet_DeepCvsB.push_back( deepcharm_c/(deepcharm_c+deepcharm_b+deepcharm_bb) );
    //==== methodB
    jet_DeepCvsL.push_back( jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL") );
    jet_DeepCvsB.push_back( jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB") );


    jet_chargedHadronEnergyFraction.push_back( jets_iter->chargedHadronEnergyFraction() );
    jet_neutralHadronEnergyFraction.push_back( jets_iter->neutralHadronEnergyFraction() );
    jet_neutralEmEnergyFraction.push_back( jets_iter->neutralEmEnergyFraction() );
    jet_chargedEmEnergyFraction.push_back( jets_iter->chargedEmEnergyFraction() );
    jet_muonEnergyFraction.push_back( jets_iter->muonEnergyFraction() );
    jet_chargedMultiplicity.push_back( jets_iter->chargedMultiplicity() );
    jet_neutralMultiplicity.push_back( jets_iter->neutralMultiplicity() );

    //==== JetID
    //==== https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
    double NHF  = jets_iter->neutralHadronEnergyFraction();
    double NEMF = jets_iter->neutralEmEnergyFraction();
    double CHF  = jets_iter->chargedHadronEnergyFraction();
    double MUF  = jets_iter->muonEnergyFraction();
    double CEMF = jets_iter->chargedEmEnergyFraction();
    double NumConst = jets_iter->chargedMultiplicity()+jets_iter->neutralMultiplicity();
    double NumNeutralParticle =jets_iter->neutralMultiplicity();
    double CHM      = jets_iter->chargedMultiplicity();
    double eta = jets_iter->eta();

    bool tightJetID        = false;
    bool tightLepVetoJetID = false;

    if(DataYear==2016){
      if(fabs(eta)<=2.4){
	tightJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
	tightLepVetoJetID = tightJetID && CEMF<0.8 && MUF <0.8; 
      }else if(fabs(eta)<=2.7){
	tightJetID = ( NEMF<0.99 && NHF < 0.9 );
	tightLepVetoJetID = tightJetID;
      }else if(fabs(eta)<=3.0){
	tightJetID = ( NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1 );
	tightLepVetoJetID = tightJetID;
      }else{
	tightJetID = ( NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 );
	tightLepVetoJetID = tightJetID;
      }
    }
    else if(DataYear==2017||DataYear==2018){
      if(fabs(eta)<=2.6){
	tightJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
	tightLepVetoJetID = tightJetID && CEMF<0.8 && MUF <0.8; 
      }else if(fabs(eta)<=2.7){
	tightJetID = ( CHM>0 && NEMF<0.99 && NHF < 0.9 );
	tightLepVetoJetID = tightJetID && CEMF<0.8 && MUF <0.8;
      }else if(fabs(eta)<=3.0){
	tightJetID = ( NEMF>0.01 && NEMF<0.99 && NumNeutralParticle>1 );
	tightLepVetoJetID = tightJetID;
      }else{
	tightJetID = ( NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 );
	tightLepVetoJetID = tightJetID;
      }
    }
    else{
      cout << "DataYear is not set : DataYear = " << DataYear << endl;
      exit(EXIT_FAILURE);
    }
      
    jet_tightJetID.push_back(tightJetID);
    jet_tightLepVetoJetID.push_back(tightLepVetoJetID);

    int partonPdgId = jets_iter->genParton() ? jets_iter->genParton()->pdgId() : 0;
    jet_partonPdgId.push_back( partonPdgId );

    if( jets_iter->hasUserFloat("vtxNtracks") ) jet_vtxNtracks.push_back( jets_iter->userFloat("vtxNtracks") );
    jet_m.push_back( jets_iter->mass() );
    jet_energy.push_back( jets_iter->energy() );
    if( jets_iter->hasUserFloat("pileupJetId:fullDiscriminant") ) jet_PileupJetId.push_back( jets_iter->userFloat("pileupJetId:fullDiscriminant") );

    //==========
    //==== JEC
    //==========

/*
    //==== Printing JEC levels
    cout << "[AK4]" << endl;
    cout << "currentJECSet = " << jets_iter->currentJECSet() << endl;
    cout << "currentJECLevel = " << jets_iter->currentJECLevel() << endl;
    cout << "availableJECLevels() : " << endl;
    const std::vector<std::string> aaa = jets_iter->availableJECLevels(jets_iter->currentJECSet());
    for(unsigned int z=0; z< aaa.size(); z++){
      cout << "  " << aaa.at(z) << endl;
    }
*/

    //==== https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties
    jet_jecUnc->setJetEta( jets_iter->eta() );
    jet_jecUnc->setJetPt( jets_iter->pt() ); // here you must use the CORRECTED jet pt
    double unc = jet_jecUnc->getUncertainty(true);
    jet_shiftedEnUp.push_back( 1.+unc );
    jet_shiftedEnDown.push_back( 1.-unc ); // I found jet_jecUnc->getUncertainty(true) = jet_jecUnc->getUncertainty(false)

    //==== For cross check
    //==== methodB
    //jet_jecUnc_methodB->setJetEta( jets_iter->eta() );
    //jet_jecUnc_methodB->setJetPt( jets_iter->pt() ); 
    //double unc_methodB = jet_jecUnc_methodB->getUncertainty(true);
    //cout << "jec unc methodA = " << unc << endl;
    //cout << "jec unc methodB = " << unc_methodB << endl << endl;

    //=======================================================================================================================
    //==== Save JEC variables for Lepton-aware JEC
    //==== https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/PhysicsTools/NanoAOD/plugins/LeptonJetVarProducer.cc#L154-L184
    //=======================================================================================================================

    //==== Save L1 Correction
    //==== auto L1Correctedjet = jets_iter->correctedP4("L1FastJet"); : "JET corrected up to L1FastJet"
    //==== jets_iter->jecFactor("L1FastJet") = L1Correctedjet.pt()/jets_iter->pt()
    //==== i.e., "FINAL JET" times "jecFactor(L1FastJet)" = "JET corrected up to L1FastJet"
    //==== So, to get L1 correction alone,
    //==== jecFactor("L1FastJet")/jecFactor("Uncorrected")

    jet_JECL1FastJet.push_back( jets_iter->jecFactor("L1FastJet")/jets_iter->jecFactor("Uncorrected") );

    //==== Save Full Correction
    //==== auto uncorjet =  jets_iter->correctedP4("Uncorrected");
    //==== jets_iter->jecFactor("Uncorrected") = uncorjet.pt()/jets_iter->pt()
    auto uncorjet =  jets_iter->correctedP4("Uncorrected");
    //cout << "uncorjet.pt()/jets_iter->pt() = " << uncorjet.pt()/jets_iter->pt() << endl;
    //cout << ( jets_iter->jecFactor("Uncorrected") ) / ( uncorjet.pt()/jets_iter->pt() ) << endl;
    jet_JECFull.push_back( jets_iter->pt()/uncorjet.pt() );

/*
    //==== debug for JEC updator
    cout << "[AK4]" << endl;
    cout << "jets_iter->pt() = " << jets_iter->pt() << endl;
    cout << "uncorjet.pt() = " << uncorjet.pt() << endl;
    cout << "ratio = " << jets_iter->pt()/uncorjet.pt() << endl;
    cout << "1/ratio = " << uncorjet.pt()/jets_iter->pt() << endl;
    cout << "unc = " << unc << endl;
    cout << "jet_JECL1FastJet = " << jets_iter->jecFactor("L1FastJet")/jets_iter->jecFactor("Uncorrected") << endl;
    cout << "jet_JECFull = " << jets_iter->pt()/uncorjet.pt() << endl;
*/

    if(!IsData){

      //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
      //==== https://github.com/cms-sw/cmssw/blob/CMSSW_9_4_X/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h
      //==== https://github.com/cms-sw/cmssw/blob/aa14a75177e9f010b3b8d0af0440598735c49311/PhysicsTools/PatUtils/python/patPFMETCorrections_cff.py

      //==========
      //==== JER
      //==========

      double jer = jet_resolution.getResolution({{JME::Binning::JetPt, jets_iter->pt()}, {JME::Binning::JetEta, jets_iter->eta()}, {JME::Binning::Rho, Rho}});
      double jer_sf = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jets_iter->pt()},{JME::Binning::JetEta, jets_iter->eta()}}, Variation::NOMINAL);
      double jer_sf_UP = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jets_iter->pt()},{JME::Binning::JetEta, jets_iter->eta()}}, Variation::UP);
      double jer_sf_DOWN = jet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jets_iter->pt()},{JME::Binning::JetEta, jets_iter->eta()}}, Variation::DOWN);

      if(theDebugLevel){

        cout << "#### Jet ####" << endl;
        cout << "pt = " << jets_iter->pt() << endl;
        cout << "eta = " << jets_iter->eta() << endl;
        cout << "Rho = " << Rho << endl;
        cout << "jer = " << jer << endl;
        cout << "jer_sf = " << jer_sf << endl;
        cout << "jer_sf_UP = " << jer_sf_UP << endl;

      }

      const reco::GenJet* genJet = nullptr;
      genJet = match(*jets_iter, jets_iter->pt() * jer, "AK4");
      double smearFactor = 1., smearFactor_UP = 1., smearFactor_DOWN = 1.;
      if(genJet){

        //==== Case 1: we have a "good" gen jet matched to the reco jet

        double dPt = jets_iter->pt() - genJet->pt();
        smearFactor = 1 + (jer_sf - 1.) * dPt / jets_iter->pt();
        smearFactor_UP = 1 + (jer_sf_UP - 1.) * dPt / jets_iter->pt();
        smearFactor_DOWN = 1 + (jer_sf_DOWN - 1.) * dPt / jets_iter->pt();

      }
      else if (jer_sf > 1) {

        //==== Case 2: we don't have a gen jet. Smear jet pt using a random gaussian variation

        double sigma = jer * std::sqrt(std::max(jer_sf*jer_sf-1,0.));
        std::normal_distribution<> d(0, sigma);
        smearFactor = 1. + d(m_random_generator);

        double sigma_UP = jer * std::sqrt(std::max(jer_sf_UP*jer_sf_UP-1,0.));
        std::normal_distribution<> d_UP(0, sigma_UP);
        smearFactor_UP = 1. + d_UP(m_random_generator);

        double sigma_DOWN = jer * std::sqrt(std::max(jer_sf_DOWN*jer_sf_DOWN-1,0.));
        std::normal_distribution<> d_DOWN(0, sigma_DOWN);
        smearFactor_DOWN = 1. + d_DOWN(m_random_generator);  

      }
      else{

      }

      //==== Negative or too small smearFactor. We would change direction of the jet
      //==== and this is not what we want.
      //==== Recompute the smearing factor in order to have jet.energy() == MIN_JET_ENERGY
      double newSmearFactor = MIN_JET_ENERGY / jets_iter->energy();

      smearFactor = std::max(smearFactor,newSmearFactor);
      smearFactor_UP = std::max(smearFactor_UP,newSmearFactor);
      smearFactor_DOWN = std::max(smearFactor_DOWN,newSmearFactor);

      if(theDebugLevel){

        cout << "--->" << endl;
        cout << "smearFactor = " << smearFactor << endl;
        cout << "smearFactor_UP = " << smearFactor_UP << endl;
        cout << "smearFactor_DOWN = " << smearFactor_DOWN << endl;

      }

      jet_smearedRes.push_back(smearFactor);
      jet_smearedResUp.push_back(smearFactor_UP);
      jet_smearedResDown.push_back(smearFactor_DOWN);

    }

    //cout << "QGTagger:qgLikelihood = " << jets_iter->userFloat("QGTagger:qgLikelihood") << endl;


  } 
  
}

/////////////////////////////
// -- Get FatJets info -- // 
////////////////////////////
void SKFlatMaker::fillFatJet(const edm::Event &iEvent)
{

  edm::Handle< std::vector<pat::Jet> > jetHandle;
  iEvent.getByToken(FatJetToken,jetHandle);

  if( jetHandle->size() > 0 && theDebugLevel > 0)
    cout << "[SKFlatMaker::fillFatJet] # of FatJets = " << jetHandle->size() << endl;

  if(jetHandle->size() == 0) return;

  const JetCollection& ForJER_jets = *jetHandle;
  unsigned int runNum_uint = static_cast <unsigned int> (iEvent.id().run());
  unsigned int lumiNum_uint = static_cast <unsigned int> (iEvent.id().luminosityBlock());
  unsigned int evNum_uint = static_cast <unsigned int> (iEvent.id().event());
  unsigned int jet0eta = uint32_t(ForJER_jets.empty() ? 0 : ForJER_jets[0].eta()/0.01);
  int m_nomVar=1;
  std::uint32_t seed = jet0eta + m_nomVar + (lumiNum_uint<<10) + (runNum_uint<<20) + evNum_uint;
  m_random_generator.seed(seed);

  for (vector<pat::Jet>::const_iterator jets_iter = jetHandle->begin(); jets_iter != jetHandle->end(); ++jets_iter){

    //==== https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1785/1/1.html
    //==== FatJet with pt 30~170 GeV are only for SM jet analysis
    //==== We can apply pt cut here

    //if(jets_iter->pt()<=170.) continue;
    if(! (jets_iter->hasPFSpecific()) ) continue; // With the new JEC, new pt > 170 is not safe anymore

    //==== https://hypernews.cern.ch/HyperNews/CMS/get/jet-algorithms/443/2.html
    fatjet_pt.push_back( jets_iter->pt() );
    fatjet_eta.push_back( jets_iter->eta() );
    fatjet_phi.push_back( jets_iter->phi() );
    fatjet_charge.push_back( jets_iter->jetCharge() );
    fatjet_area.push_back( jets_iter->jetArea() );
    fatjet_partonFlavour.push_back( jets_iter->partonFlavour() );
    fatjet_hadronFlavour.push_back( jets_iter->hadronFlavour() );

    //=========================================================================
    //==== Taggers
    //==== https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    //=========================================================================

    //=== CSVv2

    fatjet_CSVv2.push_back( jets_iter->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

    //==== DeepCSV (B)
/*
    cout << "==============================" << endl;
    cout << "---- BvsAll ----" << endl;
    cout << "methodA = " << jets_iter->bDiscriminator("pfDeepCSVJetTags:probb")+jets_iter->bDiscriminator("pfDeepCSVJetTags:probbb") << endl;
    cout << "methodB = " << jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") << endl;
*/
    //==== methodA
    //fatjet_DeepCSV.push_back( jets_iter->bDiscriminator("pfDeepCSVJetTags:probb")+jets_iter->bDiscriminator("pfDeepCSVJetTags:probbb") );
    //==== methodB
    fatjet_DeepCSV.push_back( jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:BvsAll") );

    //==== DeepFlavour

    fatjet_DeepFlavour_b.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probb"));
    fatjet_DeepFlavour_bb.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probbb"));
    fatjet_DeepFlavour_lepb.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:problepb"));
    fatjet_DeepFlavour_c.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probc"));
    fatjet_DeepFlavour_uds.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probuds"));
    fatjet_DeepFlavour_g.push_back( jets_iter->bDiscriminator("pfDeepFlavourJetTags:probg"));

    //==== Old Charm Tagger

    fatjet_CvsL.push_back( jets_iter->bDiscriminator("pfCombinedCvsLJetTags") );
    fatjet_CvsB.push_back( jets_iter->bDiscriminator("pfCombinedCvsBJetTags") );

    //==== DeepCSV charm tagger

/*
    cout << "---- CvsL ----" << endl;
    cout << "methodA = " << deepcharm_c/(deepcharm_c+deepcharm_udsg) << endl;
    cout << "methodB = " << jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL") << endl;
    cout << "---- CvsB ----" << endl;
    cout << "methodA = " << deepcharm_c/(deepcharm_c+deepcharm_b+deepcharm_bb) << endl;
    cout << "methodB = " << jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB") << endl;
*/
    //==== methodA
    //double deepcharm_c = jets_iter->bDiscriminator("pfDeepCSVJetTags:probc");
    //double deepcharm_udsg = jets_iter->bDiscriminator("pfDeepCSVJetTags:probudsg");
    //double deepcharm_b = jets_iter->bDiscriminator("pfDeepCSVJetTags:probb");
    //double deepcharm_bb = jets_iter->bDiscriminator("pfDeepCSVJetTags:probbb");
    //fatjet_DeepCvsL.push_back( deepcharm_c/(deepcharm_c+deepcharm_udsg) );
    //fatjet_DeepCvsB.push_back( deepcharm_c/(deepcharm_c+deepcharm_b+deepcharm_bb) );
    //==== methodB
    fatjet_DeepCvsL.push_back( jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsL") );
    fatjet_DeepCvsB.push_back( jets_iter->bDiscriminator("pfDeepCSVDiscriminatorsJetTags:CvsB") );


    fatjet_chargedHadronEnergyFraction.push_back( jets_iter->chargedHadronEnergyFraction() );
    fatjet_neutralHadronEnergyFraction.push_back( jets_iter->neutralHadronEnergyFraction() );
    fatjet_neutralEmEnergyFraction.push_back( jets_iter->neutralEmEnergyFraction() );
    fatjet_chargedEmEnergyFraction.push_back( jets_iter->chargedEmEnergyFraction() );
    fatjet_muonEnergyFraction.push_back( jets_iter->muonEnergyFraction() );
    fatjet_chargedMultiplicity.push_back( jets_iter->chargedMultiplicity() );
    fatjet_neutralMultiplicity.push_back( jets_iter->neutralMultiplicity() );

    //==== JetID
    //==== https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL
    // For AK8 jets, the corresponding (CHS or PUPPI) AK4 jet ID should be used. 
    double NHF  = jets_iter->neutralHadronEnergyFraction();
    double NEMF = jets_iter->neutralEmEnergyFraction();
    double CHF  = jets_iter->chargedHadronEnergyFraction();
    double MUF  = jets_iter->muonEnergyFraction();
    double CEMF = jets_iter->chargedEmEnergyFraction();
    double NumConst = jets_iter->chargedMultiplicity()+jets_iter->neutralMultiplicity();
    double NumNeutralParticle =jets_iter->neutralMultiplicity();
    double CHM      = jets_iter->chargedMultiplicity();
    double eta = jets_iter->eta();

    bool tightJetID        = false;
    bool tightLepVetoJetID = false;

    if(DataYear==2016){
      if(fabs(eta)<=2.4){
	tightJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
	tightLepVetoJetID = tightJetID && CEMF<0.8 && MUF <0.8; 
      }else if(fabs(eta)<=2.7){
	tightJetID = ( NEMF<0.99 && NHF < 0.9 );
	tightLepVetoJetID = tightJetID;
      }else if(fabs(eta)<=3.0){
	tightJetID = ( NEMF>0.0 && NEMF<0.99 && NHF<0.9 && NumNeutralParticle>1 );
	tightLepVetoJetID = tightJetID;
      }else{
	tightJetID = ( NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 );
	tightLepVetoJetID = tightJetID;
      }
    }
    else if(DataYear==2017||DataYear==2018){
      if(fabs(eta)<=2.6){
	tightJetID = ( CHM>0 && CHF>0 && NumConst>1 && NEMF<0.9 && NHF < 0.9 ); 
	tightLepVetoJetID = tightJetID && CEMF<0.8 && MUF <0.8; 
      }else if(fabs(eta)<=2.7){
	tightJetID = ( CHM>0 && NEMF<0.99 && NHF < 0.9 );
	tightLepVetoJetID = tightJetID && CEMF<0.8 && MUF <0.8;
      }else if(fabs(eta)<=3.0){
	tightJetID = ( NEMF>0.01 && NEMF<0.99 && NumNeutralParticle>1 );
	tightLepVetoJetID = tightJetID;
      }else{
	tightJetID = ( NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 );
	tightLepVetoJetID = tightJetID;
      }
    }
    else{
      cout << "DataYear is not set : DataYear = " << DataYear << endl;
      exit(EXIT_FAILURE);
    }

    fatjet_tightJetID.push_back(tightJetID);
    fatjet_tightLepVetoJetID.push_back(tightLepVetoJetID);

    int partonPdgId = jets_iter->genParton() ? jets_iter->genParton()->pdgId() : 0;
    fatjet_partonPdgId.push_back( partonPdgId );

    if( jets_iter->hasUserFloat("vtxNtracks") ) fatjet_vtxNtracks.push_back( jets_iter->userFloat("vtxNtracks") );
    fatjet_m.push_back( jets_iter->mass() );
    fatjet_energy.push_back( jets_iter->energy() );

    fatjet_puppi_tau1.push_back( jets_iter->userFloat("NjettinessAK8Puppi:tau1") );
    fatjet_puppi_tau2.push_back( jets_iter->userFloat("NjettinessAK8Puppi:tau2") );
    fatjet_puppi_tau3.push_back( jets_iter->userFloat("NjettinessAK8Puppi:tau3") );
    fatjet_puppi_tau4.push_back( jets_iter->userFloat("NjettinessAK8Puppi:tau4") );

    //==== https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorUncertainties

    //==== This was checked by jalmond and jskim. This uncorrected softdrop mass gives same result as using userfloat.
    TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
    auto const & sdSubjetsPuppi = jets_iter->subjets("SoftDropPuppi");
    for ( auto const & it : sdSubjetsPuppi ) {
      puppi_softdrop_subjet.SetPtEtaPhiM(it->correctedP4("Uncorrected").pt(),it->correctedP4("Uncorrected").eta(),it->correctedP4("Uncorrected").phi(),it->correctedP4("Uncorrected").mass());
      puppi_softdrop+=puppi_softdrop_subjet;
    }

    fatjet_softdropmass.push_back( puppi_softdrop.M() );
    
    //==== Printing JEC levels
    //cout << "[AK8]" << endl;
    //cout << "currentJECSet = " << jets_iter->currentJECSet() << endl;
    //cout << "currentJECLevel = " << jets_iter->currentJECLevel() << endl;
    //cout << "availableJECLevels() : " << endl;
    //const std::vector<std::string> aaa = jets_iter->availableJECLevels(jets_iter->currentJECSet());
    //for(unsigned int z=0; z< aaa.size(); z++){
    //  cout << "  " << aaa.at(z) << endl;
    // }


    fatjet_jecUnc->setJetEta( jets_iter->eta() );
    fatjet_jecUnc->setJetPt( jets_iter->pt() ); // here you must use the CORRECTED jet pt
    double unc = fatjet_jecUnc->getUncertainty(true);
    fatjet_shiftedEnUp.push_back( 1.+unc );
    fatjet_shiftedEnDown.push_back( 1.-unc ); // I found fatjet_jecUnc->getUncertainty(true) = fatjet_jecUnc->getUncertainty(false)

/*
    cout << "[AK8]" << endl;
    cout << "jets_iter->pt() = " << jets_iter->pt() << endl;
*/

    if(!IsData){

      //==========
      //==== JER
      //==========

      double jer = fatjet_resolution.getResolution({{JME::Binning::JetPt, jets_iter->pt()}, {JME::Binning::JetEta, jets_iter->eta()}, {JME::Binning::Rho, Rho}});
      double jer_sf = fatjet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jets_iter->pt()},{JME::Binning::JetEta, jets_iter->eta()}}, Variation::NOMINAL);
      double jer_sf_UP = fatjet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jets_iter->pt()},{JME::Binning::JetEta, jets_iter->eta()}}, Variation::UP);
      double jer_sf_DOWN = fatjet_resolution_sf.getScaleFactor({{JME::Binning::JetPt, jets_iter->pt()},{JME::Binning::JetEta, jets_iter->eta()}}, Variation::DOWN);

      if(theDebugLevel){

        cout << "#### FatJet ####" << endl;
        cout << "pt = " << jets_iter->pt() << endl;
        cout << "eta = " << jets_iter->eta() << endl;
        cout << "Rho = " << Rho << endl;
        cout << "jer = " << jer << endl;
        cout << "jer_sf = " << jer_sf << endl;
        cout << "jer_sf_UP = " << jer_sf_UP << endl;

      }

      const reco::GenJet* genJet = nullptr;
      genJet = match(*jets_iter, jets_iter->pt() * jer, "AK8");
      double smearFactor = 1., smearFactor_UP = 1., smearFactor_DOWN = 1.;
      if(genJet){

        //==== Case 1: we have a "good" gen jet matched to the reco jet

        double dPt = jets_iter->pt() - genJet->pt();
        smearFactor = 1 + (jer_sf - 1.) * dPt / jets_iter->pt();
        smearFactor_UP = 1 + (jer_sf_UP - 1.) * dPt / jets_iter->pt();
        smearFactor_DOWN = 1 + (jer_sf_DOWN - 1.) * dPt / jets_iter->pt();

      }
      else if (jer_sf > 1) {

        //==== Case 2: we don't have a gen jet. Smear jet pt using a random gaussian variation

        double sigma = jer * std::sqrt(std::max(jer_sf*jer_sf-1,0.));
        std::normal_distribution<> d(0, sigma);
        smearFactor = 1. + d(m_random_generator);

        double sigma_UP = jer * std::sqrt(std::max(jer_sf_UP*jer_sf_UP-1,0.));
        std::normal_distribution<> d_UP(0, sigma_UP);
        smearFactor_UP = 1. + d_UP(m_random_generator);

        double sigma_DOWN = jer * std::sqrt(std::max(jer_sf_DOWN*jer_sf_DOWN-1,0.));
        std::normal_distribution<> d_DOWN(0, sigma_DOWN);
        smearFactor_DOWN = 1. + d_DOWN(m_random_generator);

      }
      else{

      }

      //==== Negative or too small smearFactor. We would change direction of the jet
      //==== and this is not what we want.
      //==== Recompute the smearing factor in order to have jet.energy() == MIN_JET_ENERGY
      double newSmearFactor = MIN_JET_ENERGY / jets_iter->energy();

      smearFactor = std::max(smearFactor,newSmearFactor);
      smearFactor_UP = std::max(smearFactor_UP,newSmearFactor);
      smearFactor_DOWN = std::max(smearFactor_DOWN,newSmearFactor);

      if(theDebugLevel){

        cout << "--->" << endl;
        cout << "smearFactor = " << smearFactor << endl;
        cout << "smearFactor_UP = " << smearFactor_UP << endl;
        cout << "smearFactor_DOWN = " << smearFactor_DOWN << endl;

      }

      fatjet_smearedRes.push_back(smearFactor);
      fatjet_smearedResUp.push_back(smearFactor_UP);
      fatjet_smearedResDown.push_back(smearFactor_DOWN);

    }

    //==== LSF variables

    std::vector<reco::CandidatePtr> pfConstituents = jets_iter->getJetConstituents();
    std::vector<fastjet::PseudoJet> lClusterParticles;
    for(unsigned int ic=0; ic<pfConstituents.size(); ic++) {
      reco::CandidatePtr pfcand = pfConstituents[ic];
      fastjet::PseudoJet   pPart(pfcand->px(),pfcand->py(),pfcand->pz(),pfcand->energy());
      lClusterParticles.emplace_back(pPart);
    }
    std::sort(lClusterParticles.begin(),lClusterParticles.end(),JetTools::orderPseudoJet);

    float lepCPt(-100), lepCEta(-100), lepCPhi(-100);
    float lepCId(0);
    float this_lsf(0);

    if(JetTools::leptons((*jets_iter),3)> 0 && JetTools::leptons(*jets_iter,7)<0.8){
      lepCPt = JetTools::leptons(*jets_iter,3);
      lepCEta = JetTools::leptons(*jets_iter,5);
      lepCPhi = JetTools::leptons(*jets_iter,6);
      lepCId = JetTools::leptons(*jets_iter,4);
      std::vector<fastjet::PseudoJet> vSubC_3;
      this_lsf = JetTools::lsf(lClusterParticles, vSubC_3, lepCPt, lepCEta, lepCPhi, lepCId, 2.0, 3);

    }

    fatjet_LSF.push_back(this_lsf);
    fatjet_LSFlep_PID.push_back( lepCId );
    fatjet_LSFlep_Pt.push_back( lepCPt );
    fatjet_LSFlep_Eta.push_back( lepCEta );
    fatjet_LSFlep_Phi.push_back( lepCPhi );

  }

}

void SKFlatMaker::endRun(const Run & iRun, const EventSetup & iSetup)
{
  if(theDebugLevel) cout << "[SKFlatMaker::endRun] called" << endl;

  // -- only when LHE information is available (ex> aMC@NLO, Powheg) case. Samples generated by pythia8 doesn't work! -- //
  if( this->theStoreLHEFlag ){
    // -- LHE information -- //
    // -- ref: https://twiki.cern.ch/twiki/bin/viewauth/CMS/LHEReaderCMSSW#Retrieving_the_weights -- //
    edm::Handle<LHERunInfoProduct> LHERunInfo;
    iRun.getByToken(LHERunInfoProductToken, LHERunInfo);
    if(!LHERunInfo.isValid()) iRun.getByToken(LHERunInfoProductSourceToken, LHERunInfo);
    if(!LHERunInfo.isValid()) return;

    cout << "[SKFlatMaker::endRun] ##### Information about PDF weights #####" << endl;
    LHERunInfoProduct myLHERunInfoProduct = *(LHERunInfo.product());
    typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
    for(headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      std::cout << "[SKFlatMaker::endRun,tag] " << iter->tag() << std::endl;
      std::vector<std::string> lines = iter->lines();
      for (unsigned int iLine = 0; iLine<lines.size(); iLine++){
        std::cout << "[SKFlatMaker::endRun,lines] " << lines.at(iLine);
      }
    }
    cout << "[SKFlatMaker::endRun] ##### End of information about PDF weights #####" << endl;
  }
}
void SKFlatMaker::beginLuminosityBlock(const edm::LuminosityBlock & iLumi, const EventSetup & iSetup){
  if(printGenWeightNames){
    edm::Handle<GenLumiInfoHeader> genLumiInfoHeader;
    iLumi.getByToken(GenLumiInfoHeaderToken, genLumiInfoHeader);
    if(genLumiInfoHeader.isValid()){
      for(unsigned int i=0,n=genLumiInfoHeader->weightNames().size();i<n;i++){
	std::cout << "[SKFlatMaker::beginLuminosityBlock] <weight> id="<<i<<" name='"<<genLumiInfoHeader->weightNames().at(i)<<"' </weight>"<<std::endl;
      }
    }
    printGenWeightNames=false;
  }
}

template<class T>
const reco::GenJet* SKFlatMaker::match(const T& jet, double resolution, TString whichjet){

  edm::Handle<reco::GenJetCollection> this_genjets;

  double m_dR_max = 0.4/2.;
  if(whichjet=="AK4"){
    this_genjets = m_genJets;
    m_dR_max = 0.4/2.;
  }
  else if(whichjet=="AK8"){
    this_genjets = m_genFatJets;
    m_dR_max = 0.8/2.;
  }
  else{
    cout << "[const reco::GenJet* SKFlatMaker::match] whichjet = " << whichjet << endl;
  }

  const reco::GenJetCollection& genJets = *this_genjets;

  //==== Try to find a gen jet matching
  //==== dR < m_dR_max
  //==== dPt < m_dPt_max_factor * resolution

  double m_dPt_max_factor = 3.;

  double min_dR = std::numeric_limits<double>::infinity();
  const reco::GenJet* matched_genJet = nullptr;

  for (const auto& genJet: genJets) {
      double dR = deltaR(genJet, jet);

      if (dR > min_dR)
          continue;

      if (dR < m_dR_max ) {
          double dPt = std::abs(genJet.pt() - jet.pt());
          if (dPt > m_dPt_max_factor * resolution)
              continue;

          min_dR = dR;
          matched_genJet = &genJet;
      }
  }

  return matched_genJet;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SKFlatMaker);
