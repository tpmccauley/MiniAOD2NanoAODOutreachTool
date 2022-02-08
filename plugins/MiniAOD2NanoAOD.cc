// -*- C++ -*-
//
// Package:    MiniAOD2NanoAOD
// Class:      MiniAOD2NanoAOD

#include <memory>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "math.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "FWCore/ParameterSet/interface/Registry.h"

#include "DataFormats/JetReco/interface/PileupJetIdentifier.h"
#include "DataFormats/Common/interface/ValueMap.h"

const static std::vector<std::string> interestingTriggers = {
    "HLT_IsoMu24_eta2p1",
    "HLT_IsoMu24",
    "HLT_IsoMu17_eta2p1_LooseIsoPFTau20",
};

template <typename T>
void subtractInvisible(T g, reco::Candidate::LorentzVector& p4) {
  auto daughters = (*g).daughterRefVector();
  for (auto d = daughters.begin(); d != daughters.end(); d++) {
    const auto pdgId = (*d)->pdgId();
    if (std::abs(pdgId) == 12 || std::abs(pdgId) == 14 ||
        std::abs(pdgId) == 16 || std::abs(pdgId) == 18) {
      p4 = p4 - (*d)->p4();
    }
    subtractInvisible(*d, p4);
  }
}

template <typename T>
int findBestVisibleMatch(T& gens, reco::Candidate::LorentzVector& p4) {
  float minDeltaR = 999.0;
  int idx = -1;
  for (auto g = gens.begin(); g != gens.end(); g++) {
    auto tmp_p4 = g->p4();
    subtractInvisible(g, tmp_p4);
    const auto tmp = deltaR(tmp_p4, p4);
    if (tmp < minDeltaR) {
      minDeltaR = tmp;
      idx = g - gens.begin();
    }
  }
  return idx;
}

template <typename T>
int findBestMatch(T& gens, reco::Candidate::LorentzVector& p4) {
  float minDeltaR = 999.0;
  int idx = -1;
  for (auto g = gens.begin(); g != gens.end(); g++) {
    const auto tmp = deltaR(g->p4(), p4);
    if (tmp < minDeltaR) {
      minDeltaR = tmp;
      idx = g - gens.begin();
    }
  }
  return idx;
}

class MiniAOD2NanoAOD : public edm::EDAnalyzer {
public:
  explicit MiniAOD2NanoAOD(const edm::ParameterSet &);
  ~MiniAOD2NanoAOD();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  bool providesGoodLumisection(const edm::Event &iEvent);

  edm::EDGetTokenT<edm::TriggerResults> triggerToken;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken;
  edm::EDGetTokenT<std::vector<pat::Muon> > muonToken;
  edm::EDGetTokenT<std::vector<pat::Electron> > electronToken;

  //edm::EDGetTokenT<reco::PFTauCollection> tauToken;  
  //edm::EDGetTokenT<reco::PFTauDiscriminator> tauDiscriminatorTokens[12]; 

  edm::EDGetTokenT<std::vector<pat::Photon> > photonToken;
  edm::EDGetTokenT<std::vector<pat::MET> > metToken;
  edm::EDGetTokenT<edm::View<pat::Jet> > jetToken;
  
  TTree *tree;

  // Event information
  Int_t value_run;
  UInt_t value_lumi_block;
  ULong64_t value_event;

  // Trigger
  const static int max_trig = 1000;
  bool value_trig[max_trig];

  // Vertices
  int value_ve_n;
  float value_ve_x;
  float value_ve_y;
  float value_ve_z;

  // Muons
  const static int max_mu = 1000;
  UInt_t value_mu_n;
  float value_mu_pt[max_mu];
  float value_mu_eta[max_mu];
  float value_mu_phi[max_mu];
  float value_mu_mass[max_mu];
  int value_mu_charge[max_mu];
  float value_mu_pfreliso03all[max_mu];
  float value_mu_pfreliso04all[max_mu];
  bool value_mu_tightid[max_mu];
  bool value_mu_softid[max_mu];
  float value_mu_dxy[max_mu];
  float value_mu_dxyErr[max_mu];
  float value_mu_dz[max_mu];
  float value_mu_dzErr[max_mu];
  int value_mu_genpartidx[max_mu];
  int value_mu_jetidx[max_mu];

  // Electrons
  const static int max_el = 1000;
  UInt_t value_el_n;
  float value_el_pt[max_el];
  float value_el_eta[max_el];
  float value_el_phi[max_el];
  float value_el_mass[max_el];
  int value_el_charge[max_el];
  float value_el_pfreliso03all[max_el];
  float value_el_dxy[max_el];
  float value_el_dxyErr[max_el];
  float value_el_dz[max_el];
  float value_el_dzErr[max_el];
  bool value_el_cutbasedid[max_el];
  bool value_el_pfid[max_el];
  int value_el_genpartidx[max_el];
  int value_el_jetidx[max_el];
  /*
  // Taus
  const static int max_tau = 1000;
  UInt_t value_tau_n;
  float value_tau_pt[max_tau];
  float value_tau_eta[max_tau];
  float value_tau_phi[max_tau];
  float value_tau_mass[max_tau];
  int value_tau_charge[max_tau];
  int value_tau_decaymode[max_tau];
  float value_tau_reliso_all[max_tau];
  int value_tau_genpartidx[max_tau];
  int value_tau_jetidx[max_tau];
  bool value_tau_iddecaymode[max_tau];
  float value_tau_idisoraw[max_tau];
  bool value_tau_idisovloose[max_tau];
  bool value_tau_idisoloose[max_tau];
  bool value_tau_idisomedium[max_tau];
  bool value_tau_idisotight[max_tau];
  bool value_tau_idantieleloose[max_tau];
  bool value_tau_idantielemedium[max_tau];
  bool value_tau_idantieletight[max_tau];
  bool value_tau_idantimuloose[max_tau];
  bool value_tau_idantimumedium[max_tau];
  bool value_tau_idantimutight[max_tau];
  */

  // Photons
  const static int max_ph = 1000;
  UInt_t value_ph_n;
  float value_ph_pt[max_ph];
  float value_ph_eta[max_ph];
  float value_ph_phi[max_ph];
  float value_ph_mass[max_ph];
  int value_ph_charge[max_ph];
  float value_ph_pfreliso03all[max_ph];
  int value_ph_genpartidx[max_ph];
  int value_ph_jetidx[max_ph];
 
  // MET
  float value_met_pt;
  float value_met_phi;
  float value_met_sumet;
  float value_met_significance;
  float value_met_covxx;
  float value_met_covxy;
  float value_met_covyy;

  // Jets
  const static int max_jet = 1000;
  UInt_t value_jet_n;
  float value_jet_pt[max_jet];
  float value_jet_eta[max_jet];
  float value_jet_phi[max_jet];
  float value_jet_mass[max_jet];
  int value_jet_puid[max_jet];
  float value_jet_puid_disc[max_jet];
  float value_jet_btag[max_jet];

};

MiniAOD2NanoAOD::MiniAOD2NanoAOD(const edm::ParameterSet &iConfig)
{
  edm::Service<TFileService> fs;

  tree = fs->make<TTree>("Events", "Events");

  // Event information
  tree->Branch("run", &value_run);
  tree->Branch("luminosityBlock", &value_lumi_block);
  tree->Branch("event", &value_event);

  // Trigger
  triggerToken = consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "HLT"));

  for(size_t i = 0; i < interestingTriggers.size(); i++) {
    tree->Branch(interestingTriggers[i].c_str(), value_trig + i, (interestingTriggers[i] + "/O").c_str());
  }
  
  // Vertices
  vertexToken = consumes<reco::VertexCollection>(edm::InputTag("offlineSlimmedPrimaryVertices"));
  
  tree->Branch("PV_npvs", &value_ve_n, "PV_npvs/I");
  tree->Branch("PV_x", &value_ve_x, "PV_x/F");
  tree->Branch("PV_y", &value_ve_y, "PV_y/F");
  tree->Branch("PV_z", &value_ve_z, "PV_z/F");

  // Muons
  muonToken = consumes<std::vector<pat::Muon> >(edm::InputTag("slimmedMuons"));
  
  tree->Branch("nMuon", &value_mu_n, "nMuon/i");
  tree->Branch("Muon_pt", value_mu_pt, "Muon_pt[nMuon]/F");
  tree->Branch("Muon_eta", value_mu_eta, "Muon_eta[nMuon]/F");
  tree->Branch("Muon_phi", value_mu_phi, "Muon_phi[nMuon]/F");
  tree->Branch("Muon_mass", value_mu_mass, "Muon_mass[nMuon]/F");
  tree->Branch("Muon_charge", value_mu_charge, "Muon_charge[nMuon]/I");
  tree->Branch("Muon_pfRelIso03_all", value_mu_pfreliso03all, "Muon_pfRelIso03_all[nMuon]/F");
  tree->Branch("Muon_pfRelIso04_all", value_mu_pfreliso04all, "Muon_pfRelIso04_all[nMuon]/F");
  tree->Branch("Muon_tightId", value_mu_tightid, "Muon_tightId[nMuon]/O");
  tree->Branch("Muon_softId", value_mu_softid, "Muon_softId[nMuon]/O");
  tree->Branch("Muon_dxy", value_mu_dxy, "Muon_dxy[nMuon]/F");
  tree->Branch("Muon_dxyErr", value_mu_dxyErr, "Muon_dxyErr[nMuon]/F");
  tree->Branch("Muon_dz", value_mu_dz, "Muon_dz[nMuon]/F");
  tree->Branch("Muon_dzErr", value_mu_dzErr, "Muon_dzErr[nMuon]/F");
  tree->Branch("Muon_jetIdx", value_mu_jetidx, "Muon_jetIdx[nMuon]/I");
  tree->Branch("Muon_genPartIdx", value_mu_genpartidx, "Muon_genPartIdx[nMuon]/I");

  // Electrons
  electronToken = consumes<std::vector<pat::Electron> >(edm::InputTag("slimmedElectrons"));
  
  tree->Branch("nElectron", &value_el_n, "nElectron/i");
  tree->Branch("Electron_pt", value_el_pt, "Electron_pt[nElectron]/F");
  tree->Branch("Electron_eta", value_el_eta, "Electron_eta[nElectron]/F");
  tree->Branch("Electron_phi", value_el_phi, "Electron_phi[nElectron]/F");
  tree->Branch("Electron_mass", value_el_mass, "Electron_mass[nElectron]/F");
  tree->Branch("Electron_charge", value_el_charge, "Electron_charge[nElectron]/I");
  tree->Branch("Electron_pfRelIso03_all", value_el_pfreliso03all, "Electron_pfRelIso03_all[nElectron]/F");
  tree->Branch("Electron_dxy", value_el_dxy, "Electron_dxy[nElectron]/F");
  tree->Branch("Electron_dxyErr", value_el_dxyErr, "Electron_dxyErr[nElectron]/F");
  tree->Branch("Electron_dz", value_el_dz, "Electron_dz[nElectron]/F");
  tree->Branch("Electron_dzErr", value_el_dzErr, "Electron_dzErr[nElectron]/F");
  tree->Branch("Electron_cutBasedId", value_el_cutbasedid, "Electron_cutBasedId[nElectron]/O");
  tree->Branch("Electron_pfId", value_el_pfid, "Electron_pfId[nElectron]/O");
  tree->Branch("Electron_jetIdx", value_el_jetidx, "Electron_jetIdx[nElectron]/I");
  tree->Branch("Electron_genPartIdx", value_el_genpartidx, "Electron_genPartIdx[nElectron]/I");

  // Taus
  /*
  tauToken = consumes<reco::PFTauCollection>(edm::InputTag("hpsPFTauProducer"));

  tauDiscriminatorTokens[0] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByDecayModeFinding"));
  tauDiscriminatorTokens[1] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr"));
  tauDiscriminatorTokens[2] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByVLooseCombinedIsolationDBSumPtCorr"));
  tauDiscriminatorTokens[3] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr"));
  tauDiscriminatorTokens[4] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr"));
  tauDiscriminatorTokens[5] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr"));
  tauDiscriminatorTokens[6] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByLooseElectronRejection"));
  tauDiscriminatorTokens[7] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByMediumElectronRejection"));
  tauDiscriminatorTokens[8] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByTightElectronRejection"));
  tauDiscriminatorTokens[9] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByLooseMuonRejection"));
  tauDiscriminatorTokens[10] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByMediumMuonRejection"));
  tauDiscriminatorTokens[11] = consumes<reco::PFTauDiscriminator>(edm::InputTag("hpsPFTauDiscriminationByTightMuonRejection"));

  tree->Branch("nTau", &value_tau_n, "nTau/i");
  tree->Branch("Tau_pt", value_tau_pt, "Tau_pt[nTau]/F");
  tree->Branch("Tau_eta", value_tau_eta, "Tau_eta[nTau]/F");
  tree->Branch("Tau_phi", value_tau_phi, "Tau_phi[nTau]/F");
  tree->Branch("Tau_mass", value_tau_mass, "Tau_mass[nTau]/F");
  tree->Branch("Tau_charge", value_tau_charge, "Tau_charge[nTau]/I");
  tree->Branch("Tau_decayMode", value_tau_decaymode, "Tau_decayMode[nTau]/I");
  tree->Branch("Tau_relIso_all", value_tau_reliso_all, "Tau_relIso_all[nTau]/F");
  tree->Branch("Tau_jetIdx", value_tau_jetidx, "Tau_jetIdx[nTau]/I");
  tree->Branch("Tau_genPartIdx", value_tau_genpartidx, "Tau_genPartIdx[nTau]/I");
  tree->Branch("Tau_idDecayMode", value_tau_iddecaymode, "Tau_idDecayMode[nTau]/O");
  tree->Branch("Tau_idIsoRaw", value_tau_idisoraw, "Tau_idIsoRaw[nTau]/F");
  tree->Branch("Tau_idIsoVLoose", value_tau_idisovloose, "Tau_idIsoVLoose[nTau]/O");
  tree->Branch("Tau_idIsoLoose", value_tau_idisoloose, "Tau_idIsoLoose[nTau]/O");
  tree->Branch("Tau_idIsoMedium", value_tau_idisomedium, "Tau_idIsoMedium[nTau]/O");
  tree->Branch("Tau_idIsoTight", value_tau_idisotight, "Tau_idIsoTight[nTau]/O");
  tree->Branch("Tau_idAntiEleLoose", value_tau_idantieleloose, "Tau_idAntiEleLoose[nTau]/O");
  tree->Branch("Tau_idAntiEleMedium", value_tau_idantielemedium, "Tau_idAntiEleMedium[nTau]/O");
  tree->Branch("Tau_idAntiEleTight", value_tau_idantieletight, "Tau_idAntiEleTight[nTau]/O");
  tree->Branch("Tau_idAntiMuLoose", value_tau_idantimuloose, "Tau_idAntiMuLoose[nTau]/O");
  tree->Branch("Tau_idAntiMuMedium", value_tau_idantimumedium, "Tau_idAntiMuMedium[nTau]/O");
  tree->Branch("Tau_idAntiMuTight", value_tau_idantimutight, "Tau_idAntiMuTight[nTau]/O");
  */

  // Photons
  photonToken = consumes<std::vector<pat::Photon> >(edm::InputTag("slimmedPhotons"));

  tree->Branch("nPhoton", &value_ph_n, "nPhoton/i");
  tree->Branch("Photon_pt", value_ph_pt, "Photon_pt[nPhoton]/F");
  tree->Branch("Photon_eta", value_ph_eta, "Photon_eta[nPhoton]/F");
  tree->Branch("Photon_phi", value_ph_phi, "Photon_phi[nPhoton]/F");
  tree->Branch("Photon_mass", value_ph_mass, "Photon_mass[nPhoton]/F");
  tree->Branch("Photon_charge", value_ph_charge, "Photon_charge[nPhoton]/I");
  tree->Branch("Photon_pfRelIso03_all", value_ph_pfreliso03all, "Photon_pfRelIso03_all[nPhoton]/F");
  tree->Branch("Photon_jetIdx", value_ph_jetidx, "Photon_jetIdx[nPhoton]/I");
  tree->Branch("Photon_genPartIdx", value_ph_genpartidx, "Photon_genPartIdx[nPhoton]/I");

  // MET
  metToken = consumes<std::vector<pat::MET> >(edm::InputTag("slimmedMETs"));

  tree->Branch("MET_pt", &value_met_pt, "MET_pt/F");
  tree->Branch("MET_phi", &value_met_phi, "MET_phi/F");
  tree->Branch("MET_sumet", &value_met_sumet, "MET_sumet/F");
  tree->Branch("MET_significance", &value_met_significance, "MET_significance/F");
  tree->Branch("MET_CovXX", &value_met_covxx, "MET_CovXX/F");
  tree->Branch("MET_CovXY", &value_met_covxy, "MET_CovXY/F");
  tree->Branch("MET_CovYY", &value_met_covyy, "MET_CovYY/F");

  // Jets
  jetToken = consumes<edm::View<pat::Jet> >(edm::InputTag("slimmedJets"));

  tree->Branch("nJet", &value_jet_n, "nJet/i");
  tree->Branch("Jet_pt", value_jet_pt, "Jet_pt[nJet]/F");
  tree->Branch("Jet_eta", value_jet_eta, "Jet_eta[nJet]/F");
  tree->Branch("Jet_phi", value_jet_phi, "Jet_phi[nJet]/F");
  tree->Branch("Jet_mass", value_jet_mass, "Jet_mass[nJet]/F");
}

MiniAOD2NanoAOD::~MiniAOD2NanoAOD() {}

void MiniAOD2NanoAOD::analyze(const edm::Event &iEvent,
                          const edm::EventSetup &iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  // Event information
  value_run = iEvent.run();
  value_lumi_block = iEvent.luminosityBlock();
  value_event = iEvent.id().event();

  // Trigger results
  Handle<TriggerResults> trigger;
  iEvent.getByToken(triggerToken, trigger);

  auto psetRegistry = edm::pset::Registry::instance();
  auto triggerParams = psetRegistry->getMapped(trigger->parameterSetID());
  TriggerNames triggerNames(*triggerParams);
  TriggerResultsByName triggerByName(&(*trigger), &triggerNames);
  for (size_t i = 0; i < interestingTriggers.size(); i++) {
    value_trig[i] = false;
  }
  const auto names = triggerByName.triggerNames();
  for (size_t i = 0; i < names.size(); i++) {
    const auto name = names[i];
    for (size_t j = 0; j < interestingTriggers.size(); j++) {
      const auto interest = interestingTriggers[j];
      if (name.find(interest) == 0) {
        const auto substr = name.substr(interest.length(), 2);
        if (substr.compare("_v") == 0) {
          const auto status = triggerByName.state(name);
          if (status == 1) {
            value_trig[j] = true;
            break;
          }
        }
      }
    }
  }

 
  // Vertex
  Handle<VertexCollection> vertices;
  iEvent.getByToken(vertexToken, vertices);

  value_ve_n = vertices->size();
  value_ve_x = vertices->begin()->x();
  value_ve_y = vertices->begin()->y();
  value_ve_z = vertices->begin()->z();
  math::XYZPoint pv(vertices->begin()->position());
  
  // Muons
  Handle<std::vector<pat::Muon> > muons;
  iEvent.getByToken(muonToken, muons);

  value_mu_n = 0;
  const float mu_min_pt = 3;
 
  for (auto it = muons->begin(); it != muons->end(); it++) {
    if (it->pt() > mu_min_pt) {
    
      value_mu_pt[value_mu_n] = it->pt();
      value_mu_eta[value_mu_n] = it->eta();
      value_mu_phi[value_mu_n] = it->phi();
      value_mu_charge[value_mu_n] = it->charge();
      value_mu_mass[value_mu_n] = it->mass();
      if (it->isPFMuon() && it->isPFIsolationValid()) {
        auto iso03 = it->pfIsolationR03();
        value_mu_pfreliso03all[value_mu_n] =
            (iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt)/it->pt();
        auto iso04 = it->pfIsolationR04();
        value_mu_pfreliso04all[value_mu_n] =
            (iso04.sumChargedHadronPt + iso04.sumNeutralHadronEt + iso04.sumPhotonEt)/it->pt();
      } else {
        value_mu_pfreliso03all[value_mu_n] = -999;
        value_mu_pfreliso04all[value_mu_n] = -999;
      }
      value_mu_tightid[value_mu_n] = muon::isTightMuon(*it, *vertices->begin());
      value_mu_softid[value_mu_n] = muon::isSoftMuon(*it, *vertices->begin());
      auto trk = it->globalTrack();
      if (trk.isNonnull()) {
        value_mu_dxy[value_mu_n] = trk->dxy(pv);
        value_mu_dz[value_mu_n] = trk->dz(pv);
        value_mu_dxyErr[value_mu_n] = trk->d0Error();
        value_mu_dzErr[value_mu_n] = trk->dzError();
      } else {
        value_mu_dxy[value_mu_n] = -999;
        value_mu_dxyErr[value_mu_n] = -999;
        value_mu_dz[value_mu_n] = -999;
        value_mu_dzErr[value_mu_n] = -999;
      }
      value_mu_genpartidx[value_mu_n] = -1;
      value_mu_jetidx[value_mu_n] = -1;
      value_mu_n++;
    }
  }

  // Electrons
  Handle<std::vector<pat::Electron> > electrons;
  iEvent.getByToken(electronToken, electrons);

  value_el_n = 0;
  const float el_min_pt = 5;

  for (auto it = electrons->begin(); it != electrons->end(); it++) {
    if (it->pt() > el_min_pt) {
     
      value_el_pt[value_el_n] = it->pt();
      value_el_eta[value_el_n] = it->eta();
      value_el_phi[value_el_n] = it->phi();
      value_el_charge[value_el_n] = it->charge();
      value_el_mass[value_el_n] = it->mass();

      /*
	need to fix
      value_el_cutbasedid[value_el_n] = it->passingCutBasedPreselection();
      value_el_pfid[value_el_n] = it->passingPflowPreselection();

      if (it->passingPflowPreselection()) {
        auto iso03 = it->pfIsolationVariables();
        value_el_pfreliso03all[value_el_n] =
            (iso03.sumChargedHadronPt + iso03.sumNeutralHadronEt + iso03.sumPhotonEt)/it->pt();
      } else {
        value_el_pfreliso03all[value_el_n] = -999;
      }
      */

      auto trk = it->gsfTrack();
      value_el_dxy[value_el_n] = trk->dxy(pv);
      value_el_dz[value_el_n] = trk->dz(pv);
      value_el_dxyErr[value_el_n] = trk->d0Error();
      value_el_dzErr[value_el_n] = trk->dzError();
      value_el_jetidx[value_el_n] = -1;
      value_el_genpartidx[value_el_n] = -1;
      value_el_n++;
    }
  }

  /*
  // Taus
  // References for Tau collections and IDs:
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID#53X
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/NutShellRecipeFor5312AndNewer

  Handle<PFTauCollection> taus;
  iEvent.getByToken(tauToken, taus);

  Handle<PFTauDiscriminator> tausLooseIso, tausVLooseIso, tausMediumIso, tausTightIso,
                             tausDecayMode, tausLooseEleRej, tausMediumEleRej,
                             tausTightEleRej, tausLooseMuonRej, tausMediumMuonRej,
                             tausTightMuonRej, tausRawIso;

  iEvent.getByToken(tauDiscriminatorTokens[0], tausDecayMode);
  iEvent.getByToken(tauDiscriminatorTokens[1], tausRawIso);
  iEvent.getByToken(tauDiscriminatorTokens[2], tausVLooseIso);
  iEvent.getByToken(tauDiscriminatorTokens[3], tausLooseIso);
  iEvent.getByToken(tauDiscriminatorTokens[4], tausMediumIso);
  iEvent.getByToken(tauDiscriminatorTokens[5], tausTightIso);
  iEvent.getByToken(tauDiscriminatorTokens[6], tausLooseEleRej);
  iEvent.getByToken(tauDiscriminatorTokens[7], tausMediumEleRej);
  iEvent.getByToken(tauDiscriminatorTokens[8], tausTightEleRej);
  iEvent.getByToken(tauDiscriminatorTokens[9], tausLooseMuonRej);
  iEvent.getByToken(tauDiscriminatorTokens[10], tausMediumMuonRej);
  iEvent.getByToken(tauDiscriminatorTokens[11], tausTightMuonRej);

  const float tau_min_pt = 15;
  value_tau_n = 0;
  std::vector<PFTau> selectedTaus;
  for (auto it = taus->begin(); it != taus->end(); it++) {
    if (it->pt() > tau_min_pt) {
      selectedTaus.emplace_back(*it);
      value_tau_pt[value_tau_n] = it->pt();
      value_tau_eta[value_tau_n] = it->eta();
      value_tau_phi[value_tau_n] = it->phi();
      value_tau_charge[value_tau_n] = it->charge();
      value_tau_mass[value_tau_n] = it->mass();
      value_tau_decaymode[value_tau_n] = it->decayMode();
      // Discriminators
      const auto idx = it - taus->begin();
      value_tau_iddecaymode[value_tau_n] = tausDecayMode->operator[](idx).second;
      value_tau_idisoraw[value_tau_n] = tausRawIso->operator[](idx).second;
      value_tau_idisovloose[value_tau_n] = tausVLooseIso->operator[](idx).second;
      value_tau_idisoloose[value_tau_n] = tausLooseIso->operator[](idx).second;
      value_tau_idisomedium[value_tau_n] = tausMediumIso->operator[](idx).second;
      value_tau_idisotight[value_tau_n] = tausTightIso->operator[](idx).second;
      value_tau_idantieleloose[value_tau_n] = tausLooseEleRej->operator[](idx).second;
      value_tau_idantielemedium[value_tau_n] = tausMediumEleRej->operator[](idx).second;
      value_tau_idantieletight[value_tau_n] = tausTightEleRej->operator[](idx).second;
      value_tau_idantimuloose[value_tau_n] = tausLooseMuonRej->operator[](idx).second;
      value_tau_idantimumedium[value_tau_n] = tausMediumMuonRej->operator[](idx).second;
      value_tau_idantimutight[value_tau_n] = tausTightMuonRej->operator[](idx).second;

      value_tau_reliso_all[value_tau_n] = (it->isolationPFChargedHadrCandsPtSum() + it->isolationPFGammaCandsEtSum()) / it->pt();
      value_tau_jetidx[value_tau_n] = -1;
      value_tau_genpartidx[value_tau_n] = -1;
      value_tau_n++;
    }
  }
  */

  // Photons
  Handle<std::vector<pat::Photon> > photons;
  iEvent.getByToken(photonToken, photons);

  value_ph_n = 0;
  const float ph_min_pt = 5;
 
  for (auto it = photons->begin(); it != photons->end(); it++) {
    if (it->pt() > ph_min_pt) {
      
      value_ph_pt[value_ph_n] = it->pt();
      value_ph_eta[value_ph_n] = it->eta();
      value_ph_phi[value_ph_n] = it->phi();
      value_ph_charge[value_ph_n] = it->charge();
      value_ph_mass[value_ph_n] = it->mass();
      value_ph_pfreliso03all[value_ph_n] = it->ecalRecHitSumEtConeDR03() / it->pt();
      value_ph_jetidx[value_ph_n] = -1;
      value_ph_genpartidx[value_ph_n] = -1;
      value_ph_n++;
    }
  }
 
  // MET
  Handle<std::vector<pat::MET> > met;
  iEvent.getByToken(metToken, met);

  value_met_pt = met->begin()->pt();
  value_met_phi = met->begin()->phi();
  value_met_sumet = met->begin()->sumEt();
  value_met_significance = met->begin()->significance();
  auto cov = met->begin()->getSignificanceMatrix();
  value_met_covxx = cov[0][0];
  value_met_covxy = cov[0][1];
  value_met_covyy = cov[1][1];
 
  Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(jetToken, jets);

  const float jet_min_pt = 15;
  value_jet_n = 0;
 
  for ( size_t i=0; i<jets->size(); ++i ) {

    auto jet = jets->refAt(i);

    if ( jet->pt() > jet_min_pt ) {
            
      value_jet_pt[value_jet_n] = jet->pt();
      value_jet_eta[value_jet_n] = jet->eta();
      value_jet_phi[value_jet_n] = jet->phi();
      value_jet_mass[value_jet_n] = jet->mass();

      value_jet_n++;

    }

  }
   
  // Fill event
  tree->Fill();
}

void MiniAOD2NanoAOD::beginJob() {}

void MiniAOD2NanoAOD::endJob() {}

DEFINE_FWK_MODULE(MiniAOD2NanoAOD);
