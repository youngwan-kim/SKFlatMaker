#ifndef SKFlatMaker_H
#define SKFlatMaker_H

////////////////////////////////
// -- system include files -- //
////////////////////////////////
#include <memory>
#include <iostream>
#include <random>

//////////////////////
// -- FrameWorks -- //
//////////////////////
#include "FWCore/Framework/interface/Frameworkfwd.h" // -- Forward declarations of types in the EDM(Event Data Model).
#include "FWCore/Framework/interface/EDAnalyzer.h" // -- EDAnalyzer is the base class for all analyzer "modules".
#include "FWCore/Framework/interface/Event.h" // -- This is the primary interface for accessing EDProducts from a single collision and inserting new derived products.
#include "FWCore/Framework/interface/MakerMacros.h" // -- including temporary until a better solution can be found (leads to more physical coupling than is probably necessary)
#include "FWCore/ParameterSet/interface/ParameterSet.h" // -- Declaration for ParameterSet(parameter set) and related types 
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h" // -- Foward declarations of variables for ParameterSet
#include "FWCore/Framework/interface/ESHandle.h" // -- No description
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RegexMatch.h"


/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/MuonReco/interface/MuonFwd.h" // -- Forward declarations of muon variables
#include "DataFormats/MuonReco/interface/Muon.h" // -- A reconstructed Muon. (tracker alone, muon detector alone, combined muon plus tracker)
#include "DataFormats/PatCandidates/interface/Muon.h" // -- Analysis-level muon class, pat::Muon implements the analysis-level muon class within the 'pat' namespace. PAT(Physics Analysis Toolkit)
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h" // -- Variables shows how muon is cosmic-like
#include "SKFlatMaker/SKFlatMaker/interface/RoccoR.h"

/////////////////////////
// -- For Electrons -- //
/////////////////////////
#include "DataFormats/PatCandidates/interface/Photon.h" // -- Analysis-level Photon class, pat::Photon implements the analysis-level photon class within the 'pat' namespace
#include "DataFormats/PatCandidates/interface/Electron.h" // -- Analysis-level electron class,  pat::Electron implements the analysis-level electron class within the 'pat' namespace
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h" // -- It seems as a function refer effective areas of ECAL
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h" // -- Analysis-level MET class, pat::MET implements an analysis-level missing energy class as a 4-vector within the 'pat' namespace
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h" // --  MET made from Particle Flow candidates, with various fraction functions ex)photonEtFraction()
#include "DataFormats/METReco/interface/PFMETCollection.h"

////////////////////
// -- For Jets -- //
////////////////////
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h" // -- Analysis-level calorimeter jet class, Jet implements the analysis-level calorimeter jet class within the 'pat' namespace.
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
#include "SKFlatMaker/SKFlatMaker/interface/JetTools.h"

////////////////////////////
// -- For GenParticles -- //
////////////////////////////
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h" // -- LHE info like PDF
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h" // -- ??
#include "DataFormats/HepMCCandidate/interface/GenParticle.h" // -- 
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" // -- ???
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" // -- contains information related to the details of the pileup simulation for a given event, ex) Nvtx, z of PU, sum(pt) of vtxs


//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h" // -- Forward definitions for tracker variables
#include "DataFormats/TrackReco/interface/Track.h" // -- reconstructed tracks that are stored in the AOD and RECO. also contains a reference to more detailed information(RECO)
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h" // -- A composite Candidate  with error matrix and other vertex fix information. 
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h" // -- Foward declaration of VertexCompositeCandidate class's variables
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h" // -- Least-squares vertex fitter implemented in the Kalman Filter formalism
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h" // -- 
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h" // -- Helper class to build TransientTrack from the persistent Track.
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

////////////////////
// -- Triggers -- //
////////////////////
#include "FWCore/Common/interface/TriggerNames.h" // -- Used to access the names and indices of the triggers corresponding to a particular TriggerResults object
#include "FWCore/Common/interface/TriggerResultsByName.h" // --  Class which provides methods to access trigger results
#include "DataFormats/Common/interface/TriggerResults.h" // -- The trigger path results are maintained here as a sequence of entries, one per trigger path
#include "DataFormats/HLTReco/interface/TriggerEvent.h" // -- The single EDProduct to be saved for each event (AOD case) describing the (HLT) trigger table
#include "DataFormats/HLTReco/interface/TriggerObject.h" // --  A single trigger object (e.g., an isolated muon, or MET) described by its 4-momentum and physics type
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h" // -- Analysis-level trigger object class (stand-alone). (within the 'pat' namespace.)
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" // -- This class provides access routines to get hold of the HLT Configuration


////////////////
// -- Else -- //
////////////////
#include "DataFormats/Math/interface/LorentzVector.h" // -- (pt, eta, phi, E) OR (pt, eta, phi, M)
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h" // -- Base class for surfaces and volumes positioned in global 3D space.
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h" // -- The Muon DetUnit geometry.
#include "DataFormats/Common/interface/RefToBaseVector.h" // -- 
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h" // -- Foward definition of NamedCompositeCandidate class
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h" // -- A Candidate composed of daughters. The daughters are owned by the composite candidate.
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h" // -- Foward definition of CompositeCandidate class
#include "DataFormats/Candidate/interface/CompositeCandidate.h" // A Candidate composed of daughters. The daughters are owned by the composite candidate.
#include "DataFormats/Common/interface/View.h" // -- Provide access to the collected elements contained by any WrapperBase that is a sequence.
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" // -- 
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h" // -- Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" // -- contains information related to the details of the pileup simulation for a given event
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"


////////////////
// -- ROOT -- //
////////////////
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TMath.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/foreach.hpp>
#include <TRandom.h>

using namespace std;
using namespace pat;
using namespace edm;

namespace reco { class CandCommonVertexFitterBase; class VertexCompositeCandidate; class CandCommonVertexFitter; }

class SKFlatMaker : public edm::EDAnalyzer
{
 public:
  explicit SKFlatMaker(const edm::ParameterSet&);
  ~SKFlatMaker();
  
 private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(const edm::Run &, const edm::EventSetup & );
  virtual void endRun(const edm::Run &, const edm::EventSetup & );
  
  virtual void fillPrimaryVertex(const edm::Event &iEvent);  // fill primary vertex information
  virtual void fillMET(const edm::Event &iEvent);            // fill MET information
  virtual void fillPhotons(const edm::Event &iEvent);
  virtual void fillMuons(const edm::Event &iEvent, const edm::EventSetup& iSetup);
  virtual void fillElectrons(const edm::Event &iEvent, const edm::EventSetup& iSetup);
  virtual void fillJet(const edm::Event &iEvent);            // fill jet and b-tagging information
  virtual void fillFatJet(const edm::Event &iEvent);            // fill jet and b-tagging information
  virtual void hltReport(const edm::Event &iEvent);          // fill list of triggers fired in an event
  virtual void fillLHEInfo(const edm::Event &iEvent);
  virtual void fillGENInfo(const edm::Event &iEvent);            // fill MET information

  int DataYear;
  int theDebugLevel;                   // 0 no prints, 1 some, 2 lots
  std::string processName;
  std::string theElectronID;
  
  HLTConfigProvider hltConfig_;
  
  // -- Tokens (for 76X) -- //

  edm::EDGetTokenT< bool > BadPFMuonDzFilter_token;

  edm::EDGetTokenT< std::vector<pat::Muon> >            MuonToken;
  edm::EDGetTokenT< edm::View<pat::Electron> >          ElectronToken;
  edm::EDGetTokenT< edm::View<pat::Photon> >            PhotonToken;
  edm::EDGetTokenT< std::vector<pat::Jet> >             JetToken;
  edm::EDGetTokenT< reco::GenJetCollection >            genJetToken;
  edm::EDGetTokenT< std::vector<pat::Jet> >             FatJetToken;
  edm::EDGetTokenT< reco::GenJetCollection >            genFatJetToken;
  edm::EDGetTokenT< std::vector<pat::MET> >             MetToken;

  edm::EDGetTokenT< LHEEventProduct >               LHEEventProductToken;
  edm::EDGetTokenT< LHEEventProduct >               LHEEventProductSourceToken;
  edm::EDGetTokenT< LHERunInfoProduct >             LHERunInfoProductToken;
  edm::EDGetTokenT< LHERunInfoProduct >             LHERunInfoProductSourceToken;
  edm::EDGetTokenT< reco::GenParticleCollection>    mcLabel_;

  edm::EDGetTokenT< edm::TriggerResults >          METFilterResultsToken_PAT;
  edm::EDGetTokenT< edm::TriggerResults >          METFilterResultsToken_RECO;
  
  edm::EDGetTokenT< double >                          RhoToken;
  edm::EDGetTokenT< std::vector<reco::Conversion> >   ConversionsToken;
  edm::EDGetTokenT< std::vector< reco::GsfTrack > >   GsfTrackToken;
  
  edm::EDGetTokenT< edm::TriggerResults >                          TriggerToken;
  edm::EDGetTokenT< edm::TriggerResults >                          TriggerTokenPAT;
  edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> >    TriggerObjectToken;
  
  edm::EDGetTokenT< GenEventInfoProduct >                GenEventInfoToken;
  edm::EDGetTokenT< reco::BeamSpot >                     BeamSpotToken;
  edm::EDGetTokenT< reco::VertexCollection >             PrimaryVertexToken;
  edm::EDGetTokenT< edm::View<reco::Track> >             TrackToken;
  edm::EDGetTokenT< std::vector< PileupSummaryInfo > >   PileUpInfoToken;

  edm::Handle<edm::TriggerResults> METFilterResults;

  edm::EDGetTokenT< double > prefweight_token;
  edm::EDGetTokenT< double > prefweightup_token;
  edm::EDGetTokenT< double > prefweightdown_token;
  
  edm::FileInPath electron_EA_NHandPh_file;
  edm::FileInPath photon_EA_CH_file;
  edm::FileInPath photon_EA_HN_file;
  edm::FileInPath photon_EA_Ph_file;
  
  // -- Store flags -- // 
  bool theStorePriVtxFlag;                // Yes or No to store primary vertex
  bool theStoreJetFlag;                // Yes or No to store Jet
  bool theStoreFatJetFlag;                // Yes or No to store FatJet
  bool theStoreMETFlag;                // Yes or No to store MET 
  bool theStoreHLTReportFlag;             // Yes or No to store HLT reuslts (list of triggers fired)
  bool theStoreHLTObjectFlag;
  bool theStoreMuonFlag;
  bool theStoreElectronFlag;
  bool theStoreLHEFlag;
  bool theStoreGENFlag;
  bool theStorePhotonFlag;
  bool theStoreL1PrefireFlag;
  bool theKeepAllGen;
  bool IsData;

  // double theCrossSection;
  // double theFilterEfficiency;
  // double theTotalNevents;
  double theBDiscriminant;
  // double theIntegratedLumi;
  double theBDiscriminant_alg1;
  double theBDiscriminant_alg2;
  double theBDiscriminant_alg3;
  int theNHighQualLeptons;
  
  edm::ESHandle<TransientTrackBuilder> theTTBuilder;

  std::vector<std::string > HLTName_WildCard;
  
  bool Flag_goodVertices;
  bool Flag_globalSuperTightHalo2016Filter;
  bool Flag_HBHENoiseFilter;
  bool Flag_HBHENoiseIsoFilter;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool Flag_BadPFMuonFilter;
  bool Flag_BadPFMuonDzFilter;
  //bool Flag_hfNoisyHitsFilter; //in MiniAODv2
  bool Flag_BadChargedCandidateFilter;
  bool Flag_eeBadScFilter;
  bool Flag_ecalBadCalibFilter;
  
  unsigned int nPileUp;
  
  std::vector< double > PDFAlphaSScaleValue_;
  map<TString,vector<int>> WeightMap;

  TTree *DYTree;
  
  int nEvt;
  
  // Invariant Mass distribution of SS(OS) di-muon events GG (2 global)
  // GlbTree
  int runNum;
  unsigned long long evtNum;
  int lumiBlock;
  double sumEt;
  double photonEt;
  double chargedHadronEt;
  double neutralHadronEt;
  float Rho;
  
  int Nelectrons;
  
  // PV
  int nPV;
  int PVtrackSize;
  float PVchi2;
  float PVndof;
  float PVnormalizedChi2;
  float PVx;
  float PVy;
  float PVz;
  double PVprob;

  //==== L1 Prefire reweights

  float L1PrefireReweight_Central;
  float L1PrefireReweight_Up;
  float L1PrefireReweight_Down;

  //==== trigger object

  vector<string> HLT_TriggerName;
  vector<string> HLT_TriggerFilterName;
  vector<float> HLTObject_pt;
  vector<float> HLTObject_eta;
  vector<float> HLTObject_phi;
  vector<string> HLTObject_FiredFilters;
  vector<string> HLTObject_FiredPaths;

  //==== Jet

  vector<float> jet_pt;
  vector<float> jet_eta;
  vector<float> jet_phi;
  vector<float> jet_charge;
  vector<float> jet_area;
  vector<int> jet_partonFlavour;
  vector<int> jet_hadronFlavour;
  vector<float> jet_CSVv2;
  vector<float> jet_DeepCSV;
  vector<float> jet_CvsL;
  vector<float> jet_CvsB;
  vector<float> jet_DeepFlavour_b;
  vector<float> jet_DeepFlavour_bb;
  vector<float> jet_DeepFlavour_lepb;
  vector<float> jet_DeepFlavour_c;
  vector<float> jet_DeepFlavour_uds;
  vector<float> jet_DeepFlavour_g;
  vector<float> jet_DeepCvsL;
  vector<float> jet_DeepCvsB;
  vector<float> jet_chargedHadronEnergyFraction;
  vector<float> jet_neutralHadronEnergyFraction;
  vector<float> jet_neutralEmEnergyFraction;
  vector<float> jet_chargedEmEnergyFraction;
  vector<float> jet_muonEnergyFraction;
  vector<int> jet_chargedMultiplicity;
  vector<int> jet_neutralMultiplicity;
  vector<bool> jet_tightJetID;
  vector<bool> jet_tightLepVetoJetID;
  vector<int> jet_partonPdgId;
  vector<int> jet_vtxNtracks;
  vector<float> jet_m;
  vector<float> jet_energy;
  vector<float> jet_PileupJetId;
  vector<float> jet_shiftedEnUp;
  vector<float> jet_shiftedEnDown;
  vector<float> jet_smearedRes;
  vector<float> jet_smearedResUp;
  vector<float> jet_smearedResDown;
  vector<float> jet_JECL1FastJet;
  vector<float> jet_JECFull;

  //==== JEC
  JetCorrectionUncertainty *jet_jecUnc;
  //JetCorrectionUncertainty *jet_jecUnc_methodB; // for cross check
  std::string jet_payloadName_;
  JetCorrectionUncertainty *fatjet_jecUnc;
  std::string fatjet_payloadName_;

  //==== JER

  static constexpr const double MIN_JET_ENERGY = 1e-2;
  std::mt19937 m_random_generator;
  template<class T>
  const reco::GenJet* match(const T& jet, double resolution, TString whichjet);

  JME::JetResolution jet_resolution;
  JME::JetResolutionScaleFactor jet_resolution_sf;
  edm::Handle<reco::GenJetCollection> m_genJets;

  JME::JetResolution fatjet_resolution;
  JME::JetResolutionScaleFactor fatjet_resolution_sf;
  edm::Handle<reco::GenJetCollection> m_genFatJets;

  //==== FatJet

  vector<float> fatjet_pt;
  vector<float> fatjet_eta;
  vector<float> fatjet_phi;
  vector<float> fatjet_charge;
  vector<float> fatjet_area;
  vector<int> fatjet_partonFlavour;
  vector<int> fatjet_hadronFlavour;
  vector<float> fatjet_CSVv2;
  vector<float> fatjet_DeepCSV;
  vector<float> fatjet_CvsL;
  vector<float> fatjet_CvsB;
  vector<float> fatjet_DeepFlavour_b;
  vector<float> fatjet_DeepFlavour_bb;
  vector<float> fatjet_DeepFlavour_lepb;
  vector<float> fatjet_DeepFlavour_c;
  vector<float> fatjet_DeepFlavour_uds;
  vector<float> fatjet_DeepFlavour_g;
  vector<float> fatjet_DeepCvsL;
  vector<float> fatjet_DeepCvsB;
  vector<bool> fatjet_tightJetID;
  vector<bool> fatjet_tightLepVetoJetID;
  vector<int> fatjet_partonPdgId;
  vector<int> fatjet_vtxNtracks;
  vector<float> fatjet_m;
  vector<float> fatjet_energy;
  vector<float> fatjet_puppi_tau1;
  vector<float> fatjet_puppi_tau2;
  vector<float> fatjet_puppi_tau3;
  vector<float> fatjet_puppi_tau4;
  vector<float> fatjet_softdropmass;
  vector<float> fatjet_chargedHadronEnergyFraction;
  vector<float> fatjet_neutralHadronEnergyFraction;
  vector<float> fatjet_neutralEmEnergyFraction;
  vector<float> fatjet_chargedEmEnergyFraction;
  vector<float> fatjet_muonEnergyFraction;
  vector<int> fatjet_chargedMultiplicity;
  vector<int> fatjet_neutralMultiplicity;
  vector<float> fatjet_shiftedEnUp;
  vector<float> fatjet_shiftedEnDown;
  vector<float> fatjet_smearedRes;
  vector<float> fatjet_smearedResUp;
  vector<float> fatjet_smearedResDown;
  vector<float> fatjet_LSF;
  vector<float> fatjet_LSFlep_PID;
  vector<float> fatjet_LSFlep_Pt;
  vector<float> fatjet_LSFlep_Eta;
  vector<float> fatjet_LSFlep_Phi;

  //==== Electron

  vector<std::string> electron_IDtoSave;
  vector<float> electron_MVAIso;
  vector<float> electron_MVANoIso;
  vector<float> electron_Energy;
  vector<float> electron_Energy_Scale_Up;
  vector<float> electron_Energy_Scale_Down;
  vector<float> electron_Energy_Smear_Up;
  vector<float> electron_Energy_Smear_Down;
  vector<float> electron_eta;
  vector<float> electron_phi;
  vector<int> electron_charge;
  vector<float> electron_gsfpt;
  vector<float> electron_gsfEta;
  vector<float> electron_gsfPhi;
  vector<int> electron_gsfCharge;
  vector<float> electron_scEta;
  vector<float> electron_scPhi;
  vector<float> electron_etaWidth;
  vector<float> electron_phiWidth;
  vector<float> electron_dEtaIn;
  vector<float> electron_dEtaInSeed;
  vector<float> electron_dPhiIn;
  vector<float> electron_sigmaIEtaIEta;
  vector<float> electron_Full5x5_SigmaIEtaIEta;
  vector<float> electron_e2x5OverE5x5;
  vector<float> electron_e1x5OverE5x5;
  vector<float> electron_HoverE;
  vector<float> electron_fbrem;
  vector<float> electron_eOverP;
  vector<float> electron_InvEminusInvP;
  vector<float> electron_dxyVTX;
  vector<float> electron_dxyerrVTX;
  vector<float> electron_dzVTX;
  vector<float> electron_dzerrVTX;
  vector<float> electron_3DIPVTX;
  vector<float> electron_3DIPerrVTX;
  vector<float> electron_dxy;
  vector<float> electron_sigdxy;
  vector<float> electron_dz;
  vector<float> electron_dxyBS;
  vector<float> electron_dzBS;
  vector<float> electron_chIso03;
  vector<float> electron_nhIso03;
  vector<float> electron_phIso03;
  vector<float> electron_puChIso03;
  vector<bool> electron_passConversionVeto;
  vector<bool> electron_isGsfCtfScPixChargeConsistent;
  vector<bool> electron_isGsfScPixChargeConsistent;
  vector<bool> electron_isGsfCtfChargeConsistent;
  vector<int> electron_mHits;
  vector<int> electron_ecalDriven;
  vector<float> electron_r9;
  vector<float> electron_scEnergy;
  vector<float> electron_scPreEnergy;
  vector<float> electron_scRawEnergy;
  vector<float> electron_scEt;
  vector<float> electron_E15;
  vector<float> electron_E25;
  vector<float> electron_E55;
  vector<float> electron_RelPFIso_dBeta;
  vector<float> electron_RelPFIso_Rho;
  vector<unsigned int> electron_IDBit;
  vector<int> electron_IDCutBit;
  vector<float> electron_EnergyUnCorr;
  vector<float> electron_chMiniIso;
  vector<float> electron_nhMiniIso;
  vector<float> electron_phMiniIso;
  vector<float> electron_puChMiniIso;
  vector<float> electron_trackIso;
  vector<float> electron_dr03EcalRecHitSumEt;
  vector<float> electron_dr03HcalDepth1TowerSumEt;
  vector<float> electron_dr03HcalTowerSumEt;
  vector<float> electron_dr03TkSumPt;
  vector<float> electron_ecalPFClusterIso;
  vector<float> electron_hcalPFClusterIso;
  vector<ULong64_t> electron_pathbits;
  vector<ULong64_t> electron_filterbits;

  //==== Muon

  vector<float> muon_PfChargedHadronIsoR04;
  vector<float> muon_PfNeutralHadronIsoR04;
  vector<float> muon_PfGammaIsoR04;
  vector<float> muon_PFSumPUIsoR04;
  vector<float> muon_PfChargedHadronIsoR03;
  vector<float> muon_PfNeutralHadronIsoR03;
  vector<float> muon_PfGammaIsoR03;
  vector<float> muon_PFSumPUIsoR03;
  vector<unsigned int> muon_TypeBit;
  vector<unsigned int> muon_IDBit;
  vector<bool> muon_ishighpt;
  vector<float> muon_dB;
  vector<float> muon_phi;
  vector<float> muon_eta;
  vector<float> muon_pt;
  vector<float> muon_mass;
  vector<float> muon_trkiso;
  vector<float> muon_hcaliso;
  vector<float> muon_ecaliso;
  vector<float> muon_trkisoR05;
  vector<float> muon_hcalisoR05;
  vector<float> muon_ecalisoR05;
  vector<int> muon_charge;
  vector<int> muon_nChambers;
  vector<int> muon_matchedstations;
  vector<int> muon_stationMask;
  vector<int> muon_nSegments;
  vector<float> muon_normchi;
  vector<int> muon_validhits;
  vector<int> muon_trackerHits;
  vector<int> muon_pixelHits;
  vector<int> muon_validmuonhits;
  vector<int> muon_trackerLayers;
  vector<int> muon_trackerHitsGLB;
  vector<int> muon_trackerLayersGLB;
  vector<int> muon_pixelHitsGLB;
  vector<float> muon_qoverp;
  vector<float> muon_theta;
  vector<float> muon_lambda;
  vector<float> muon_dxy;
  vector<float> muon_d0;
  vector<float> muon_dsz;
  vector<float> muon_dz;
  vector<float> muon_dxyBS;
  vector<float> muon_dzBS;
  vector<float> muon_dszBS;
  vector<float> muon_dxyVTX;
  vector<float> muon_dxyerrVTX;
  vector<float> muon_dzVTX;
  vector<float> muon_dzerrVTX;
  vector<float> muon_3DIPVTX;
  vector<float> muon_3DIPerrVTX;
  vector<float> muon_dszVTX;
  vector<float> muon_vx;
  vector<float> muon_vy;
  vector<float> muon_vz;
  vector<float> muon_Best_pt;
  vector<float> muon_Best_ptError;
  vector<float> muon_Best_eta;
  vector<float> muon_Best_phi;
  vector<float> muon_Inner_pt;
  vector<float> muon_Inner_ptError;
  vector<float> muon_Inner_eta;
  vector<float> muon_Inner_phi;
  vector<float> muon_Outer_pt;
  vector<float> muon_Outer_ptError;
  vector<float> muon_Outer_eta;
  vector<float> muon_Outer_phi;
  vector<float> muon_GLB_pt;
  vector<float> muon_GLB_ptError;
  vector<float> muon_GLB_eta;
  vector<float> muon_GLB_phi;
  vector<float> muon_TuneP_pt;
  vector<float> muon_TuneP_ptError;
  vector<float> muon_TuneP_eta;
  vector<float> muon_TuneP_phi;
  vector<float> muon_TuneP_charge;
  vector<float> muon_roch_sf;
  vector<float> muon_roch_sf_up;
  vector<float> muon_PfChargedHadronMiniIso;
  vector<float> muon_PfNeutralHadronMiniIso;
  vector<float> muon_PfGammaMiniIso;
  vector<float> muon_PFSumPUMiniIso;
  vector<float> muon_MVA;
  vector<float> muon_lowptMVA;
  vector<float> muon_softMVA;
  vector<float> muon_jetPtRatio;
  vector<float> muon_jetPtRel;
  vector<int> muon_simType;
  vector<int> muon_simExtType;
  vector<int> muon_simFlavour;
  vector<int> muon_simHeaviestMotherFlavour;
  vector<int> muon_simPdgId;
  vector<int> muon_simMotherPdgId;
  vector<float> muon_simMatchQuality;
  vector<ULong64_t> muon_pathbits;
  vector<ULong64_t> muon_filterbits;

  //==== Rochestor correction
  RoccoR rc;

  //==== LHE

  vector<float> LHE_Px;
  vector<float> LHE_Py;
  vector<float> LHE_Pz;
  vector<float> LHE_E;
  vector<int> LHE_Status;
  vector<int> LHE_ID;
  map<TString,vector<float>> weight_;
 
  //==== GEN

  vector<float> gen_phi;
  vector<float> gen_eta;
  vector<float> gen_pt;
  vector<float> gen_mass;
  vector<float> gen_charge;
  vector<int> gen_mother_index;
  vector<int> gen_status;
  vector<int> gen_PID;
  vector<bool> gen_isPrompt;
  vector<bool> gen_isPromptFinalState;
  vector<bool> gen_isTauDecayProduct;
  vector<bool> gen_isPromptTauDecayProduct;
  vector<bool> gen_isDirectPromptTauDecayProductFinalState;
  vector<bool> gen_isHardProcess;
  vector<bool> gen_isLastCopy;
  vector<bool> gen_isLastCopyBeforeFSR;
  vector<bool> gen_isPromptDecayed;
  vector<bool> gen_isDecayedLeptonHadron;
  vector<bool> gen_fromHardProcessBeforeFSR;
  vector<bool> gen_fromHardProcessDecayed;
  vector<bool> gen_fromHardProcessFinalState;
  vector<bool> gen_isMostlyLikePythia6Status3;
  float gen_weight;
  float genWeight_Q;
  float genWeight_X1;
  float genWeight_X2;
  int genWeight_id1;
  int genWeight_id2;
  float genWeight_alphaQCD;
  float genWeight_alphaQED;
  
  //==== Photon information
  vector<float> photon_Energy;
  vector<float> photon_EnergyUnCorr;
  vector<float> photon_eta;
  vector<float> photon_phi;
  vector<float> photon_scEta;
  vector<float> photon_scPhi;
  vector<float> photon_HoverE;
  vector<bool> photon_hasPixelSeed;
  vector<float> photon_Full5x5_SigmaIEtaIEta;
  vector<float> photon_ChIso;
  vector<float> photon_NhIso;
  vector<float> photon_PhIso;
  vector<float> photon_ChIsoWithEA;
  vector<float> photon_NhIsoWithEA;
  vector<float> photon_PhIsoWithEA;
  vector<bool> photon_passMVAID_WP80;
  vector<bool> photon_passMVAID_WP90;
  vector<bool> photon_passLooseID;
  vector<bool> photon_passMediumID;
  vector<bool> photon_passTightID;

  // Effective area constants for all isolation types
  // EffectiveAreas effAreaChHadrons_;
  // EffectiveAreas effAreaNeuHadrons_;
  // EffectiveAreas effAreaPhotons_;

  //==== MET
  float pfMET_pt;
  float pfMET_phi;
  float pfMET_SumEt;
  float pfMET_Type1_pt;
  float pfMET_Type1_phi;
  float pfMET_Type1_SumEt;
  float pfMET_Type1_PhiCor_pt;
  float pfMET_Type1_PhiCor_phi;
  float pfMET_Type1_PhiCor_SumEt;
  vector<float> pfMET_pt_shifts;
  vector<float> pfMET_phi_shifts;
  vector<float> pfMET_SumEt_shifts;
  vector<float> pfMET_Type1_pt_shifts;
  vector<float> pfMET_Type1_phi_shifts;
  vector<float> pfMET_Type1_SumEt_shifts;
  vector<float> pfMET_Type1_PhiCor_pt_shifts;
  vector<float> pfMET_Type1_PhiCor_phi_shifts;
  vector<float> pfMET_Type1_PhiCor_SumEt_shifts;

};
#endif
