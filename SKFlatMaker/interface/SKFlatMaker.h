#ifndef SKFlatMaker_H
#define SKFlatMaker_H

////////////////////////////////
// -- system include files -- //
////////////////////////////////
#include <memory>
#include <iostream>

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
#include "DataFormats/PatCandidates/interface/Jet.h" // -- Analysis-level calorimeter jet class, Jet implements the analysis-level calorimeter jet class within the 'pat' namespace.

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
  virtual void hltReport(const edm::Event &iEvent);          // fill list of triggers fired in an event
  virtual void fillLHEInfo(const edm::Event &iEvent);
  virtual void fillGENInfo(const edm::Event &iEvent);            // fill MET information
  virtual void fillTT(const edm::Event&);
  
  bool reorder(double &a, double &b)
  {
    return a > b;
  }
  
  int theDebugLevel;                   // 0 no prints, 1 some, 2 lots
  std::string processName;
  std::string theElectronID;
  
  HLTConfigProvider hltConfig_;
  
  // -- Tokens (for 76X) -- //
  edm::EDGetTokenT< std::vector<pat::Muon> >             MuonToken;
  edm::EDGetTokenT< edm::View<pat::Electron> >                         ElectronToken;
  edm::EDGetTokenT< edm::View<pat::Electron> >                                          UnCorrElectronToken;
  edm::EDGetTokenT< edm::View<pat::Photon> >                   PhotonToken;
  edm::EDGetTokenT< edm::View<pat::Photon> >                                            UnCorrPhotonToken;
  edm::EDGetTokenT< std::vector<pat::Jet> >             JetToken;
  edm::EDGetTokenT< std::vector<pat::MET> >             MetToken;
  //edm::EDGetTokenT<pat::METCollection>                                           MetToken;

  edm::EDGetTokenT< LHEEventProduct >               LHEEventProductToken;
  edm::EDGetTokenT< LHERunInfoProduct >              LHERunInfoProductToken;
  edm::EDGetTokenT< reco::GenParticleCollection>                                        mcLabel_;
  
  edm::EDGetTokenT< edm::TriggerResults >                                               METFilterResultsToken_PAT;
  edm::EDGetTokenT< edm::TriggerResults >                                               METFilterResultsToken_RECO;
  
  edm::EDGetTokenT< double >                      RhoToken;
  edm::EDGetTokenT< edm::ValueMap<float> >                                              mvaIsoValuesMapToken;
  edm::EDGetTokenT< edm::ValueMap<float> >                                              mvaNoIsoValuesMapToken;
  edm::EDGetTokenT< std::vector<reco::Conversion> > 		           		ConversionsToken;
  edm::EDGetTokenT< std::vector< reco::GsfTrack > > 			         	GsfTrackToken;
  
  edm::EDGetTokenT< edm::TriggerResults > 						TriggerToken;
  edm::EDGetTokenT< edm::TriggerResults > 						TriggerTokenPAT;
  edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > 	                TriggerObjectToken;
  
  edm::EDGetTokenT< GenEventInfoProduct > 						GenEventInfoToken;
  edm::EDGetTokenT< reco::BeamSpot > 						       	BeamSpotToken;
  edm::EDGetTokenT< reco::VertexCollection > 						PrimaryVertexToken;
  edm::EDGetTokenT< edm::View<reco::Track> > 						TrackToken;
  edm::EDGetTokenT< std::vector< PileupSummaryInfo > > 	                   		PileUpInfoToken;
  
  
  //edm::Handle<bool> ifilterbadChCand;
  //edm::Handle<bool> ifilterbadPFMuon;
  edm::Handle<edm::TriggerResults> METFilterResults;
  
  // // -- Photon information -- //
  // edm::InputTag thePhotonLabel;
  // edm::InputTag full5x5SigmaIEtaIEtaMapLabel; 
  // edm::InputTag phoChargedIsolationLabel; 
  // edm::InputTag phoNeutralHadronIsolationLabel; 
  // edm::InputTag phoPhotonIsolationLabel;     
  edm::FileInPath effAreaChHadronsFile;
  edm::FileInPath effAreaNeuHadronsFile;
  edm::FileInPath effAreaPhotonsFile;  
  
    
  // -- Store flags -- // 
  bool theStorePriVtxFlag;                // Yes or No to store primary vertex
  bool theStoreJetFlag;                // Yes or No to store Jet
  bool theStoreMETFlag;                // Yes or No to store MET 
  bool theStoreHLTReportFlag;             // Yes or No to store HLT reuslts (list of triggers fired)
  bool theStoreMuonFlag;
  bool theStoreElectronFlag;
  bool theStoreLHEFlag;
  bool theStoreGENFlag;
  bool theStorePhotonFlag;
  bool theStoreTTFlag;
  bool theKeepAllGen;
  bool isMC;  //turn gen on and off
  bool theApplyFilter;
  int theFilterType;
  
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

  std::vector<std::string > MuonHLT;
  std::vector<int > MuonHLTPS;
  std::vector<std::string > ListHLT;
  std::vector<int > ListHLTPS;
  std::vector<std::string > trigModuleNames;
  std::vector<std::string > trigModuleNames_preFil;
  
  bool Flag_goodVertices;
  bool Flag_globalTightHalo2016Filter;
  bool Flag_HBHENoiseFilter;
  bool Flag_HBHENoiseIsoFilter;
  bool Flag_EcalDeadCellTriggerPrimitiveFilter;
  bool Flag_BadPFMuonFilter;
  bool Flag_BadChargedCandidateFilter;
  bool Flag_eeBadScFilter;
  bool Flag_ecalBadCalibFilter;
  
  // Pile-up Reweight
  // edm::LumiReWeighting LumiWeights_;
  // reweight::PoissonMeanShifter PShiftUp_;
  // reweight::PoissonMeanShifter PShiftDown_;
  // edm::LumiReWeighting LumiWeightsMuonPhys_;
  // reweight::PoissonMeanShifter PShiftUpMuonPhys_;
  // reweight::PoissonMeanShifter PShiftDownMuonPhys_;
  
  // std::vector<double> PileUpRD_;
  // std::vector<double> PileUpRDMuonPhys_;
  // std::vector<double> PileUpMC_;
  
  unsigned int nPileUp;
  double pileUpReweightIn;
  double pileUpReweight;
  double pileUpReweightPlus;
  double pileUpReweightMinus;
  
  double pileUpReweightInMuonPhys;
  double pileUpReweightMuonPhys;
  double pileUpReweightPlusMuonPhys;
  double pileUpReweightMinusMuonPhys;
  
  TTree *DYTree;
  
  int nEvt;
  
  //FILL a tree
  static const int MPSIZE = 2000;
  
  // Invariant Mass distribution of SS(OS) di-muon events GG (2 global)
  // GlbTree
  int runNum;
  unsigned long long evtNum;
  int lumiBlock;
  int nMuon;
  double PUweight;
  double sumEt;
  double photonEt;
  double chargedHadronEt;
  double neutralHadronEt;
  
  // double MET_sumEt;
  // double MET_pt;
  // double MET_px;
  // double MET_py;
  // double MET_phi;
  // double pfMET_sumEt;
  // double pfMET_pt;
  // double pfMET_px;
  // double pfMET_py;
  // double pfMET_phi;
  int Njets;
  int Nelectrons;
  int Nmuons;
  int Nbtagged;
  int NbtaggedCloseMuon;
  
  // -- PDf weights -- //
  std::vector< double > PDFWeights;
  
  // -- Flags in re-miniAOD -- //
  bool Flag_duplicateMuons;
  bool Flag_badMuons;
  bool Flag_noBadMuons;
  
  // PV
  int nVertices;
  int PVtrackSize;
  double PVchi2;
  double PVndof;
  double PVnormalizedChi2;
  double PVx;
  double PVy;
  double PVz;
  double PVprob;
  
  // trigger object
  int _HLT_ntrig;
  int _HLT_trigType[MPSIZE];
  int _HLT_trigFired[MPSIZE];
  std::vector<std::string> _HLT_trigName;
  std::vector<int> _HLT_trigPS;
  double _HLT_trigPt[MPSIZE];
  double _HLT_trigEta[MPSIZE];
  double _HLT_trigPhi[MPSIZE];
  
  // Jet
  double JETbDiscriminant[MPSIZE];
  double JETcharge[MPSIZE];
  int JETflavour[MPSIZE];
  int JETntracks[MPSIZE];
  double JETpt[MPSIZE];
  double JETeta[MPSIZE];
  double JETphi[MPSIZE];
  
  int Nbtagged_alg1;
  int NbtaggedCloseMuon_alg1;
  int Nbtagged_alg2;
  int NbtaggedCloseMuon_alg2;
  int Nbtagged_alg3;
  int NbtaggedCloseMuon_alg3;
  double JETbDiscriminant_alg1[MPSIZE];
  double JETbDiscriminant_alg2[MPSIZE];
  double JETbDiscriminant_alg3[MPSIZE];
  
  double Jet_pT[MPSIZE]; 
  double Jet_eta[MPSIZE]; 
  double Jet_phi[MPSIZE]; 
  double Jet_Charge[MPSIZE]; 
  double Jet_area[MPSIZE];
  double Jet_rho[MPSIZE];
  int Jet_Flavor[MPSIZE]; 
  int Jet_Hadron[MPSIZE];
  double Jet_bTag[MPSIZE]; 
  double Jet_CHfrac[MPSIZE]; 
  double Jet_NHfrac[MPSIZE]; 
  double Jet_NHEMfrac[MPSIZE]; 
  double Jet_CHEMfrac[MPSIZE]; 
  int Jet_CHmulti[MPSIZE]; 
  int Jet_NHmulti[MPSIZE];
  
  
  // PAT Electron
  double Electron_MVAIso[MPSIZE];
  double Electron_MVANoIso[MPSIZE];
  double Electron_et[MPSIZE];
  double Electron_caloEnergy[MPSIZE];
  double Electron_Energy[MPSIZE];
  double Electron_pT[MPSIZE];
  double Electron_Px[MPSIZE];
  double Electron_Py[MPSIZE];
  double Electron_Pz[MPSIZE];
  double Electron_eta[MPSIZE];
  double Electron_phi[MPSIZE];
  int Electron_charge[MPSIZE];
  double Electron_gsfpT[MPSIZE];
  double Electron_gsfPx[MPSIZE];
  double Electron_gsfPy[MPSIZE];
  double Electron_gsfPz[MPSIZE];
  double Electron_gsfEta[MPSIZE];
  double Electron_gsfPhi[MPSIZE];
  int Electron_gsfCharge[MPSIZE];
  double Electron_etaSC[MPSIZE];
  double Electron_phiSC[MPSIZE];
  double Electron_etaWidth[MPSIZE];
  double Electron_phiWidth[MPSIZE];
  double Electron_dEtaIn[MPSIZE];
  double Electron_dEtaInSeed[MPSIZE];
  double Electron_dPhiIn[MPSIZE];
  double Electron_sigmaIEtaIEta[MPSIZE];
  double Electron_Full5x5_SigmaIEtaIEta[MPSIZE];
  double Electron_HoverE[MPSIZE];
  double Electron_fbrem[MPSIZE];
  double Electron_eOverP[MPSIZE];
  double Electron_energyEC[MPSIZE];
  double Electron_Pnorm[MPSIZE];
  double Electron_InvEminusInvP[MPSIZE];
  double Electron_dxyVTX[MPSIZE];
  double Electron_dzVTX[MPSIZE];
  double Electron_dxy[MPSIZE];
  double Electron_sigdxy[MPSIZE];
  double Electron_dz[MPSIZE];
  double Electrron_ip3D[MPSIZE];
  double Electrron_sigip3D[MPSIZE];
  double Electron_dxyBS[MPSIZE];
  double Electron_dzBS[MPSIZE];
  double Electron_AEff03[MPSIZE];
  double Electron_chIso03[MPSIZE];
  double Electron_chIso04[MPSIZE];
  double Electron_nhIso03[MPSIZE];
  double Electron_nhIso04[MPSIZE];
  double Electron_phIso03[MPSIZE];
  double Electron_phIso04[MPSIZE];
  double Electron_pcIso03[MPSIZE];
  double Electron_pcIso04[MPSIZE];
  double Electron_relIsoCom03[MPSIZE];
  double Electron_relIsoCom04[MPSIZE];
  double Electron_relIsoBeta03[MPSIZE];
  double Electron_relIsoBeta04[MPSIZE];
  double Electron_relIsoRho03[MPSIZE];
  bool Electron_hasConversion[MPSIZE];
  int Electron_mHits[MPSIZE];
  
  int Electron_crack[MPSIZE];
  int Electron_ecalDriven[MPSIZE];
  double Electron_isoEMHADDepth1[MPSIZE];
  double Electron_25over55[MPSIZE];
  double Electron_15over55[MPSIZE];
  double Electron_isoHADDepth2[MPSIZE];
  double Electron_isoPtTrks[MPSIZE];
  double Electron_modIsoEMHADDepth1[MPSIZE];
  double Electron_modIsoPtTrks[MPSIZE];
  double Electron_modIsoEMHADDepth1Orig[MPSIZE];
  double Electron_modIsoPtTrksOrig[MPSIZE];
  double Electron_ambGsf0Pt[MPSIZE];
  double Electron_ambGsf0Eta[MPSIZE];
  double Electron_ambGsf0Phi[MPSIZE];
  double Electron_ambGsf0Charge[MPSIZE];
  double Electron_ambGsf1Pt[MPSIZE];
  double Electron_ambGsf1Eta[MPSIZE];
  double Electron_ambGsf1Phi[MPSIZE];
  double Electron_ambGsf1Charge[MPSIZE];
  double Electron_ambGsf2Pt[MPSIZE];
  double Electron_ambGsf2Eta[MPSIZE];
  double Electron_ambGsf2Phi[MPSIZE];
  double Electron_ambGsf2Charge[MPSIZE];
  double Electron_ambGsf3Pt[MPSIZE];
  double Electron_ambGsf3Eta[MPSIZE];
  double Electron_ambGsf3Phi[MPSIZE];
  double Electron_ambGsf3Charge[MPSIZE];
  
  double Electron_r9[MPSIZE];
  double Electron_EnergySC[MPSIZE];
  double Electron_preEnergySC[MPSIZE];
  double Electron_rawEnergySC[MPSIZE];
  double Electron_etSC[MPSIZE];
  double Electron_E15[MPSIZE];
  double Electron_E25[MPSIZE];
  double Electron_E55[MPSIZE];
  double Electron_ChIso03FromPU[MPSIZE];
  double Electron_RelPFIso_dBeta[MPSIZE];
  double Electron_RelPFIso_Rho[MPSIZE];
  bool Electron_passConvVeto[MPSIZE];
  bool Electron_passVetoID[MPSIZE];
  bool Electron_passLooseID[MPSIZE];
  bool Electron_passMediumID[MPSIZE];
  bool Electron_passTightID[MPSIZE];
  bool Electron_passMVAID_noIso_WP80[MPSIZE];
  bool Electron_passMVAID_noIso_WP90[MPSIZE];
  bool Electron_passMVAID_iso_WP80[MPSIZE];
  bool Electron_passMVAID_iso_WP90[MPSIZE];
  bool Electron_passMVAID_WP80[MPSIZE];
  bool Electron_passMVAID_WP90[MPSIZE];
  bool Electron_passHEEPID[MPSIZE];
  
  std::vector<double> vtxTrkDiEChi2;
  std::vector<double> vtxTrkDiEProb;
  std::vector<double> vtxTrkDiENdof;
  std::vector<double> vtxTrkDiE1Pt;
  std::vector<double> vtxTrkDiE2Pt;
  
  // -- emu vertex -- //
  std::vector<double> vtxTrkEMuChi2;
  std::vector<double> vtxTrkEMuProb;
  std::vector<double> vtxTrkEMuNdof;
  std::vector<double> vtxTrkEMu1Pt;
  std::vector<double> vtxTrkEMu2Pt;
  std::vector<double> vtxTrkEMuChi2_TuneP;
  std::vector<double> vtxTrkEMuProb_TuneP;
  std::vector<double> vtxTrkEMuNdof_TuneP;
  std::vector<double> vtxTrkEMu1Pt_TuneP;
  std::vector<double> vtxTrkEMu2Pt_TuneP;
  
  // -- Un-corrected electrons -- //
  int nUnCorrElectron;
  double Electron_pTUnCorr[MPSIZE];
  double Electron_etaUnCorr[MPSIZE];
  double Electron_phiUnCorr[MPSIZE];
  double Electron_PxUnCorr[MPSIZE];
  double Electron_PyUnCorr[MPSIZE];
  double Electron_PzUnCorr[MPSIZE];
  double Electron_EnergyUnCorr[MPSIZE];
  double Electron_EnergySCUnCorr[MPSIZE];
  double Electron_etaSCUnCorr[MPSIZE];
  double Electron_phiSCUnCorr[MPSIZE];
  double Electron_etSCUnCorr[MPSIZE];
  
  // Pat Muon
  //pf isolations
  double Muon_PfChargedHadronIsoR05[MPSIZE];
  double Muon_PfNeutralHadronIsoR05[MPSIZE];
  double Muon_PfGammaIsoR05[MPSIZE];
  double Muon_PfChargedHadronIsoR04[MPSIZE];
  double Muon_PfNeutralHadronIsoR04[MPSIZE];
  double Muon_PfGammaIsoR04[MPSIZE];
  double Muon_PFSumPUIsoR04[MPSIZE];
  double Muon_PfChargedHadronIsoR03[MPSIZE];
  double Muon_PfNeutralHadronIsoR03[MPSIZE];
  double Muon_PfGammaIsoR03[MPSIZE];
  double Muon_PFSumPUIsoR03[MPSIZE];
  
  int isPFmuon[MPSIZE];
  int isGLBmuon[MPSIZE];
  int isTRKmuon[MPSIZE];
  int isSTAmuon[MPSIZE];
  
  int Muon_muonType[MPSIZE];
  int Muon_nTrig[MPSIZE];
  int Muon_triggerObjectType[MPSIZE];
  int Muon_filterName[MPSIZE];
  double Muon_dB[MPSIZE];
  double Muon_phi[MPSIZE];
  double Muon_eta[MPSIZE];
  double Muon_pT[MPSIZE];
  double Muon_cktpT[MPSIZE];
  double Muon_cktPx[MPSIZE];
  double Muon_cktPy[MPSIZE];
  double Muon_cktPz[MPSIZE];
  double Muon_cktpTError[MPSIZE];
  double Muon_Px[MPSIZE];
  double Muon_Py[MPSIZE];
  double Muon_Pz[MPSIZE];
  double Muon_sumtrkpt[MPSIZE];
  double Muon_trkiso[MPSIZE];
  double Muon_hcaliso[MPSIZE];
  double Muon_ecaliso[MPSIZE];
  double Muon_trkisoR05[MPSIZE];
  double Muon_hcalisoR05[MPSIZE];
  double Muon_ecalisoR05[MPSIZE];
  int Muon_charge[MPSIZE];
  int Muon_nChambers[MPSIZE];
  int Muon_nMatches[MPSIZE];
  int Muon_nMatchesRPCLayers[MPSIZE];
  int Muon_stationMask[MPSIZE];
  int Muon_nSegments[MPSIZE];
  double Muon_chi2dof[MPSIZE];
  int Muon_nhits[MPSIZE];
  int Muon_trackerHits[MPSIZE];
  int Muon_pixelHits[MPSIZE];
  int Muon_muonHits[MPSIZE];
  int Muon_trackerLayers[MPSIZE];
  int Muon_trackerHitsGLB[MPSIZE];
  int Muon_trackerLayersGLB[MPSIZE];
  int Muon_pixelHitsGLB[MPSIZE];
  double Muon_qoverp[MPSIZE];
  double Muon_theta[MPSIZE];
  double Muon_lambda[MPSIZE];
  double Muon_dxy[MPSIZE];
  double Muon_d0[MPSIZE];
  double Muon_dsz[MPSIZE];
  double Muon_dz[MPSIZE];
  double Muon_dxyBS[MPSIZE];
  double Muon_dzBS[MPSIZE];
  double Muon_dszBS[MPSIZE];
  double Muon_dxyVTX[MPSIZE];
  double Muon_dzVTX[MPSIZE];
  double Muon_dszVTX[MPSIZE];
  double Muon_dxycktVTX[MPSIZE];
  double Muon_dzcktVTX[MPSIZE];
  double Muon_dszcktVTX[MPSIZE];
  double Muon_vx[MPSIZE];
  double Muon_vy[MPSIZE];
  double Muon_vz[MPSIZE];
  std::vector<double> CosAngle;
  std::vector<double> vtxTrkChi2;
  std::vector<double> vtxTrkProb;
  std::vector<double> vtxTrkNdof;
  std::vector<double> vtxTrkCkt1Pt;
  std::vector<double> vtxTrkCkt2Pt;
  
  std::vector<double> CosAngle_TuneP;
  std::vector<double> vtxTrk1Pt_TuneP;
  std::vector<double> vtxTrk2Pt_TuneP;
  std::vector<double> vtxTrkChi2_TuneP;
  std::vector<double> vtxTrkNdof_TuneP;
  std::vector<double> vtxTrkProb_TuneP;
  
  //Various track informations
  //MuonBestTrack
  double Muon_Best_pT[MPSIZE];
  double Muon_Best_pTError[MPSIZE];
  double Muon_Best_Px[MPSIZE];
  double Muon_Best_Py[MPSIZE];
  double Muon_Best_Pz[MPSIZE];
  double Muon_Best_eta[MPSIZE];
  double Muon_Best_phi[MPSIZE];
  //Inner Track
  double Muon_Inner_pT[MPSIZE];
  double Muon_Inner_pTError[MPSIZE];
  double Muon_Inner_Px[MPSIZE];
  double Muon_Inner_Py[MPSIZE];
  double Muon_Inner_Pz[MPSIZE];
  double Muon_Inner_eta[MPSIZE];
  double Muon_Inner_phi[MPSIZE];
  //Outer Track
  double Muon_Outer_pT[MPSIZE];
  double Muon_Outer_pTError[MPSIZE];
  double Muon_Outer_Px[MPSIZE];
  double Muon_Outer_Py[MPSIZE];
  double Muon_Outer_Pz[MPSIZE];
  double Muon_Outer_eta[MPSIZE];
  double Muon_Outer_phi[MPSIZE];
  //Global Track
  double Muon_GLB_pT[MPSIZE];
  double Muon_GLB_pTError[MPSIZE];
  double Muon_GLB_Px[MPSIZE];
  double Muon_GLB_Py[MPSIZE];
  double Muon_GLB_Pz[MPSIZE];
  double Muon_GLB_eta[MPSIZE];
  double Muon_GLB_phi[MPSIZE];
  
  //tuneP MuonBestTrack
  double Muon_TuneP_pT[MPSIZE];
  double Muon_TuneP_pTError[MPSIZE];
  double Muon_TuneP_Px[MPSIZE];
  double Muon_TuneP_Py[MPSIZE];
  double Muon_TuneP_Pz[MPSIZE];
  double Muon_TuneP_eta[MPSIZE];
  double Muon_TuneP_phi[MPSIZE];
  
  // LHE
  int nLHEParticle;
  double LHELepton_Px[MPSIZE];
  double LHELepton_Py[MPSIZE];
  double LHELepton_Pz[MPSIZE];
  double LHELepton_E[MPSIZE];
  int LHELepton_ID[MPSIZE];
  int LHELepton_status[MPSIZE];
  
  // GEN
  int GENnPair;
  double GENLepton_phi[MPSIZE];
  double GENLepton_eta[MPSIZE];
  double GENLepton_pT[MPSIZE];
  double GENLepton_Px[MPSIZE];
  double GENLepton_Py[MPSIZE];
  double GENLepton_Pz[MPSIZE];
  double GENLepton_E[MPSIZE];
  int GENLepton_mother[MPSIZE];
  double GENLepton_mother_pT[MPSIZE];
  int GENLepton_mother_index[MPSIZE];
  int GENLepton_charge[MPSIZE];
  int GENLepton_status[MPSIZE];
  int GENLepton_ID[MPSIZE];
  int GENLepton_isPrompt[MPSIZE];
  int GENLepton_isPromptFinalState[MPSIZE];
  int GENLepton_isTauDecayProduct[MPSIZE];
  int GENLepton_isPromptTauDecayProduct[MPSIZE];
  int GENLepton_isDirectPromptTauDecayProductFinalState[MPSIZE];
  int GENLepton_isHardProcess[MPSIZE];
  int GENLepton_isLastCopy[MPSIZE];
  int GENLepton_isLastCopyBeforeFSR[MPSIZE];
  int GENLepton_isPromptDecayed[MPSIZE];
  int GENLepton_isDecayedLeptonHadron[MPSIZE];
  int GENLepton_fromHardProcessBeforeFSR[MPSIZE];
  int GENLepton_fromHardProcessDecayed[MPSIZE];
  int GENLepton_fromHardProcessFinalState[MPSIZE];
  int GENLepton_isMostlyLikePythia6Status3[MPSIZE];
  double GENEvt_weight;
  double GENEvt_QScale;
  double GENEvt_x1;
  double GENEvt_x2;
  int GENEvt_id1;
  int GENEvt_id2;
  double GENEvt_alphaQCD;
  double GENEvt_alphaQED;
  
  
  // -- Photon information -- //
  int nPhotons;
  double Photon_pT[MPSIZE];
  double Photon_eta[MPSIZE];
  double Photon_phi[MPSIZE];
  double Photon_etaSC[MPSIZE];
  double Photon_phiSC[MPSIZE];
  double Photon_HoverE[MPSIZE];
  int Photon_hasPixelSeed[MPSIZE];
  double Photon_Full5x5_SigmaIEtaIEta[MPSIZE];
  double Photon_ChIso[MPSIZE];
  double Photon_NhIso[MPSIZE];
  double Photon_PhIso[MPSIZE];
  double Photon_ChIsoWithEA[MPSIZE];
  double Photon_NhIsoWithEA[MPSIZE];
  double Photon_PhIsoWithEA[MPSIZE];
  bool Photon_passMVAID_WP80[MPSIZE];
  bool Photon_passMVAID_WP90[MPSIZE];
  bool Photon_passLooseID[MPSIZE];
  bool Photon_passMediumID[MPSIZE];
  bool Photon_passTightID[MPSIZE];

  int nUnCorrPhoton;
  double Photon_pTUnCorr[MPSIZE];
  double Photon_etaUnCorr[MPSIZE];
  double Photon_phiUnCorr[MPSIZE];
  
  // Effective area constants for all isolation types
  // EffectiveAreas effAreaChHadrons_;
  // EffectiveAreas effAreaNeuHadrons_;
  // EffectiveAreas effAreaPhotons_;
  
  int NTT;
  double TTrack_dxy[MPSIZE];
  double TTrack_dxyErr[MPSIZE];
  double TTrack_d0[MPSIZE];
  double TTrack_d0Err[MPSIZE];
  double TTrack_dsz[MPSIZE];
  double TTrack_dszErr[MPSIZE];
  double TTrack_dz[MPSIZE];
  double TTrack_dzErr[MPSIZE];
  double TTrack_dxyBS[MPSIZE];
  double TTrack_dszBS[MPSIZE];
  double TTrack_dzBS[MPSIZE];
  double TTrack_pT[MPSIZE];
  double TTrack_Px[MPSIZE];
  double TTrack_Py[MPSIZE];
  double TTrack_Pz[MPSIZE];
  double TTrack_eta[MPSIZE];
  double TTrack_phi[MPSIZE];
  double TTrack_charge[MPSIZE];
  
  // -- MET -- //
  double pfMET_pT;
  double pfMET_phi; 
  double pfMET_Px; 
  double pfMET_Py; 
  double pfMET_SumEt;
  
  double pfMET_Type1_pT;
  double pfMET_Type1_phi; 
  double pfMET_Type1_Px; 
  double pfMET_Type1_Py; 
  double pfMET_Type1_SumEt;
  
  double pfMET_Type1_PhiCor_pT;
  double pfMET_Type1_PhiCor_phi; 
  double pfMET_Type1_PhiCor_Px; 
  double pfMET_Type1_PhiCor_Py; 
  double pfMET_Type1_PhiCor_SumEt;
};
#endif
