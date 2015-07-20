#ifndef NERO_FATJETS_H
#define NERO_FATJETS_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BareFatJets.hpp"
#include "RecoBTag/SecondaryVertex/interface/TemplatedSimpleSecondaryVertexComputer.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "RecoBTag/SecondaryVertex/interface/TrackKinematics.h"
#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/Njettiness.hh"


class NeroFatJets : virtual public NeroCollection,
    virtual public BareFatJets
{
    public:
      typedef reco::CandIPTagInfo IPTagInfo;
      typedef reco::VertexCompositePtrCandidate Vertex;
      typedef reco::TemplatedSecondaryVertexTagInfo<IPTagInfo,Vertex> SVTagInfo;
      typedef typename IPTagInfo::input_container Tracks;
      typedef typename IPTagInfo::input_container::value_type TrackRef;
      NeroFatJets(double r0);
      ~NeroFatJets();
      int analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);
      virtual inline string name(){return "NeroFatJets";};
      void setBranchAddresses(TTree *);
      void defineBranches(TTree *);
      void clear();

      bool doSubstructure = 0;

      // --- specific fuctions
      // --- Handle
      edm::Handle<pat::JetCollection> handle;
      edm::Handle<reco::VertexCollection> primaryVertices;

      // --- Token
      edm::EDGetTokenT<pat::JetCollection> token;

      void SetPrefix(std::string n)       { prefix = n;         }

      void SetSubjetsName (std::string n) { SubjetsName = n;    }
      std::string GetSubjetsName()        { return SubjetsName; }

      void SetR0 (float n)                { R0 = n;             }
      float GetR0()                       { return R0;          }

      void runBTagging(const pat::Jet* jet);
      void vertexKinematicsAndCharge(const Vertex & vertex, reco::TrackKinematics & vertexKinematics, Int_t & charge);
      void recalcNsubjettiness(const pat::Jet &jet,
                                            const SVTagInfo &svTagInfo,
                                            float & tau1,
                                            float & tau2,
                                            std::vector<fastjet::PseudoJet> &currentAxes);
      void setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);
      void setTracksPVBase(const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);
      void computeTrkVars(const IPTagInfo * ipTagInfo, const SVTagInfo * svTagInfo, reco::TrackKinematics & allKinematics);
      void computeSVVars(const IPTagInfo * ipTagInfo, const SVTagInfo * svTagInfo);
      void setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);
      void setTracksPVBase(const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);
      void setTracksSV (const TrackRef & trackRef, const SVTagInfo * svTagInfo, int & isFromSV, int & iSV, float & SVweight);


    protected:
      const GenericMVAJetTagComputer *computer;
      std::string SubjetsName = "SubJets";
      float R0                = -1;
      fastjet::contrib::Njettiness njettiness;
      std::string prefix = "fatjet";

        /////////////////////////////////
       //       btagging vars         //
      /////////////////////////////////
      std::vector<float> *tau1IVF;                       // per fatjet
      std::vector<float> *tau2IVF;
      std::vector<unsigned int> *nSV;
      std::vector<float> *zRatio;
      std::vector<float> *tauDot;
      std::vector<int> *nTracks;

      std::vector<std::vector<float>*> *svMass;          // per secondary vertex in fatjet (but no more than 4 saved)
      std::vector<std::vector<float>*> *svEnergyRatio;
      std::vector<std::vector<float>*> *svPt;

      std::vector<vector<float>*> *track_p;       // per track_ in fatjet
      std::vector<vector<float>*> *track_eta;
      std::vector<vector<float>*> *track_phi;
      std::vector<vector<float>*> *track_ptRel;
      std::vector<vector<float>*> *track_pPar;
      std::vector<vector<float>*> *track_etaRel;
      std::vector<vector<float>*> *track_deltaR;
      std::vector<vector<float>*> *track_ptRatio;
      std::vector<vector<float>*> *track_pParRatio;
      std::vector<vector<float>*> *track_sip2dVal;
      std::vector<vector<float>*> *track_sip2dSig;
      std::vector<vector<float>*> *track_sip3dVal;
      std::vector<vector<float>*> *track_sip3dSig;
      std::vector<vector<float>*> *track_decayLenVal;
      std::vector<vector<float>*> *track_decayLenSig;
      std::vector<vector<float>*> *track_jetDistVal;
      std::vector<vector<float>*> *track_jetDistSig;
      std::vector<vector<float>*> *track_chi2;
      std::vector<vector<int>*> *track_nTotalHits;
      std::vector<vector<int>*> *track_nPixelHits;
      std::vector<vector<float>*> *track_dxy;
      std::vector<vector<float>*> *track_dz;
      std::vector<vector<float>*> *track_IP2D;
      std::vector<vector<float>*> *track_IP2Dsig;
      std::vector<vector<float>*> *track_IP;
      std::vector<vector<float>*> *track_IPsig;
      std::vector<vector<float>*> *track_IP2Derr;
      std::vector<vector<float>*> *track_IPerr;
      std::vector<vector<float>*> *track_prob;
      std::vector<vector<float>*> *track_PVWeight;
      std::vector<vector<float>*> *track_SVWeight;
      std::vector<vector<float>*> *track_charge;
      std::vector<vector<int>*> *track_nHitStrip;
      std::vector<vector<int>*> *track_nHitTIB;
      std::vector<vector<int>*> *track_nHitTID;
      std::vector<vector<int>*> *track_nHitTOB;
      std::vector<vector<int>*> *track_nHitTEC;
      std::vector<vector<int>*> *track_nHitPXB;
      std::vector<vector<int>*> *track_nHitPXF;
      std::vector<vector<int>*> *track_isHitL1;
      std::vector<vector<int>*> *track_PV;
      std::vector<vector<int>*> *track_fromSV;
      std::vector<vector<int>*> *track_SV;


      //grooming vars
      std::vector<TLorentzVector> *prunedMomentum;
      std::vector<TLorentzVector> *trimmedMomentum;
      std::vector<TLorentzVector> *sdMomentum;
      std::vector<float> *prunedJEC0;
      std::vector<float> *trimmedJEC0;
      std::vector<float> *sdJEC0;

};

oid NeroFatJets::clear() {
  BareFatJets::clear();
  tau1IVF->clear();
  tau2IVF->clear();
  nSV->clear();
  zRatio->clear();
  tauDot->clear();
  svMass->clear();
  svEnergyRatio->clear();
  svPt->clear();
  prunedMomentum->clear();
  trimmedMomentum->clear();
  sdMomentum->clear();
  prunedJEC0->clear();
  trimmedJEC0->clear();
  sdJEC0->clear();
  nTracks->clear();
  track_momentum->clear();
  track_eta->clear();
  track_phi->clear();
  track_ptRel->clear();
  track_pPar->clear();
  track_etaRel->clear();
  track_deltaR->clear();
  track_ptRatio->clear();
  track_pParRatio->clear();
  track_sip2dVal->clear();
  track_sip2dSig->clear();
  track_sip3dVal->clear();
  track_sip3dSig->clear();
  track_decayLenVal->clear();
  track_decayLenSig->clear();
  track_jetDistVal->clear();
  track_jetDistSig->clear();
  track_chi2->clear();
  track_nTotalHits->clear();
  track_nPixelHits->clear();
  track_dxy->clear();
  track_dz->clear();
  track_IP2D->clear();
  track_IP2Dsig->clear();
  track_IP->clear();
  track_IPsig->clear();
  track_IP2Derr->clear();
  track_IPerr->clear();
  track_prob->clear();
  track_PVWeight->clear();
  track_SVWeight->clear();
  track_charge->clear();
  track_nHitStrip->clear();
  track_nHitTIB->clear();
  track_nHitTID->clear();
  track_nHitTOB->clear();
  track_nHitTEC->clear();
  track_nHitPXB->clear();
  track_nHitPXF->clear();
  track_isHitL1->clear();
  track_PV->clear();
  track_fromSV->clear();
  track_SV->clear();

}

void NeroFatJets::defineBranches(TTree * t){
  BareFatJets::defineBranches(t,prefix);
  tau1IVF = new vector<float>;
  tau2IVF = new vector<float>;
  nSV = new vector<unsigned int>;
  zRatio = new vector<float>;
  tauDot = new vector<float>;
  svMass = new vector<vector<float>*>;
  svEnergyRatio = new vector<vector<float>*>;
  svPt = new vector<vector<float>*>;
  prunedMomentum = new vector<TLorentzVector>;
  trimmedMomentum = new vector<TLorentzVector>;
  sdMomentum = new vector<TLorentzVector>;
  prunedJEC0 = new vector<float>;
  trimmedJEC0 = new vector<float>;
  sdJEC0 = new vector<float>;
  nTracks = new vector<int>;
  track_momentum = new vector<vector<float>*>;
  track_eta = new vector<vector<float>*>;
  track_phi = new vector<vector<float>*>;
  track_ptRel = new vector<vector<float>*>;
  track_pPar = new vector<vector<float>*>;
  track_etaRel = new vector<vector<float>*>;
  track_deltaR = new vector<vector<float>*>;
  track_ptRatio = new vector<vector<float>*>;
  track_pParRatio = new vector<vector<float>*>;
  track_sip2dVal = new vector<vector<float>*>;
  track_sip2dSig = new vector<vector<float>*>;
  track_sip3dVal = new vector<vector<float>*>;
  track_sip3dSig = new vector<vector<float>*>;
  track_decayLenVal = new vector<vector<float>*>;
  track_decayLenSig = new vector<vector<float>*>;
  track_jetDistVal = new vector<vector<float>*>;
  track_jetDistSig = new vector<vector<float>*>;
  track_chi2 = new vector<vector<float>*>;
  track_nTotalHits = new vector<vector<int>*>;
  track_nPixelHits = new vector<vector<int>*>;
  track_dxy = new vector<vector<float>*>;
  track_dz = new vector<vector<float>*>;
  track_IP2D = new vector<vector<float>*>;
  track_IP2Dsig = new vector<vector<float>*>;
  track_IP = new vector<vector<float>*>;
  track_IPsig = new vector<vector<float>*>;
  track_IP2Derr = new vector<vector<float>*>;
  track_IPerr = new vector<vector<float>*>;
  track_prob = new vector<vector<float>*>;
  track_PVWeight = new vector<vector<float>*>;
  track_SVWeight = new vector<vector<float>*>;
  track_charge = new vector<vector<float>*>;
  track_nHitStrip = new vector<vector<int>*>;
  track_nHitTIB = new vector<vector<int>*>;
  track_nHitTID = new vector<vector<int>*>;
  track_nHitTOB = new vector<vector<int>*>;
  track_nHitTEC = new vector<vector<int>*>;
  track_nHitPXB = new vector<vector<int>*>;
  track_nHitPXF = new vector<vector<int>*>;
  track_isHitL1 = new vector<vector<int>*>;
  track_PV = new vector<vector<int>*>;
  track_fromSV = new vector<vector<int>*>;
  track_SV = new vector<vector<int>*>;

  t->Branch((prefix+string("_tau1IVF")).c_str(),"vector<float>",&tau1IVF);
  t->Branch((prefix+string("_tau2IVF")).c_str(),"vector<float>",&tau2IVF);
  t->Branch((prefix+string("_nSV")).c_str(),"vector<unsigned int>",&nSV);
  t->Branch((prefix+string("_zRatio")).c_str(),"vector<float>",&zRatio);
  t->Branch((prefix+string("_tauDot")).c_str(),"vector<float>",&tauDot);
  t->Branch((prefix+string("_svMass")).c_str(),"vector<vector<float>*>",&svMass);
  t->Branch((prefix+string("_svEnergyRatio")).c_str(),"vector<vector<float>*>",&svEnergyRatio);
  t->Branch((prefix+string("_svPt")).c_str(),"vector<vector<float>*>",&svPt);
  t->Branch((prefix+string("_prunedMomentum")).c_str(),"vector<TLorentzVector>",&prunedMomentum);
  t->Branch((prefix+string("_trimmedMomentum")).c_str(),"vector<TLorentzVector>",&trimmedMomentum);
  t->Branch((prefix+string("_sdMomentum")).c_str(),"vector<TLorentzVector>",&sdMomentum);
  t->Branch((prefix+string("_prunedJEC0")).c_str(),"vector<float>",&prunedJEC0);
  t->Branch((prefix+string("_trimmedJEC0")).c_str(),"vector<float>",&trimmedJEC0);
  t->Branch((prefix+string("_sdJEC0")).c_str(),"vector<float>",&sdJEC0);
  t->Branch((prefix+string("_nTracks")).c_str(),"vector<int>",&nTracks);
  t->Branch((prefix+string("_track_momentum")).c_str(),"vector<vector<float>*>",&track_momentum);
  t->Branch((prefix+string("_track_eta")).c_str(),"vector<vector<float>*>",&track_eta);
  t->Branch((prefix+string("_track_phi")).c_str(),"vector<vector<float>*>",&track_phi);
  t->Branch((prefix+string("_track_ptRel")).c_str(),"vector<vector<float>*>",&track_ptRel);
  t->Branch((prefix+string("_track_pPar")).c_str(),"vector<vector<float>*>",&track_pPar);
  t->Branch((prefix+string("_track_etaRel")).c_str(),"vector<vector<float>*>",&track_etaRel);
  t->Branch((prefix+string("_track_deltaR")).c_str(),"vector<vector<float>*>",&track_deltaR);
  t->Branch((prefix+string("_track_ptRatio")).c_str(),"vector<vector<float>*>",&track_ptRatio);
  t->Branch((prefix+string("_track_pParRatio")).c_str(),"vector<vector<float>*>",&track_pParRatio);
  t->Branch((prefix+string("_track_sip2dVal")).c_str(),"vector<vector<float>*>",&track_sip2dVal);
  t->Branch((prefix+string("_track_sip2dSig")).c_str(),"vector<vector<float>*>",&track_sip2dSig);
  t->Branch((prefix+string("_track_sip3dVal")).c_str(),"vector<vector<float>*>",&track_sip3dVal);
  t->Branch((prefix+string("_track_sip3dSig")).c_str(),"vector<vector<float>*>",&track_sip3dSig);
  t->Branch((prefix+string("_track_decayLenVal")).c_str(),"vector<vector<float>*>",&track_decayLenVal);
  t->Branch((prefix+string("_track_decayLenSig")).c_str(),"vector<vector<float>*>",&track_decayLenSig);
  t->Branch((prefix+string("_track_jetDistVal")).c_str(),"vector<vector<float>*>",&track_jetDistVal);
  t->Branch((prefix+string("_track_jetDistSig")).c_str(),"vector<vector<float>*>",&track_jetDistSig);
  t->Branch((prefix+string("_track_chi2")).c_str(),"vector<vector<float>*>",&track_chi2);
  t->Branch((prefix+string("_track_nTotalHits")).c_str(),"vector<vector<int>*>",&track_nTotalHits);
  t->Branch((prefix+string("_track_nPixelHits")).c_str(),"vector<vector<int>*>",&track_nPixelHits);
  t->Branch((prefix+string("_track_dxy")).c_str(),"vector<vector<float>*>",&track_dxy);
  t->Branch((prefix+string("_track_dz")).c_str(),"vector<vector<float>*>",&track_dz);
  t->Branch((prefix+string("_track_IP2D")).c_str(),"vector<vector<float>*>",&track_IP2D);
  t->Branch((prefix+string("_track_IP2Dsig")).c_str(),"vector<vector<float>*>",&track_IP2Dsig);
  t->Branch((prefix+string("_track_IP")).c_str(),"vector<vector<float>*>",&track_IP);
  t->Branch((prefix+string("_track_IPsig")).c_str(),"vector<vector<float>*>",&track_IPsig);
  t->Branch((prefix+string("_track_IP2Derr")).c_str(),"vector<vector<float>*>",&track_IP2Derr);
  t->Branch((prefix+string("_track_IPerr")).c_str(),"vector<vector<float>*>",&track_IPerr);
  t->Branch((prefix+string("_track_prob")).c_str(),"vector<vector<float>*>",&track_prob);
  t->Branch((prefix+string("_track_PVWeight")).c_str(),"vector<vector<float>*>",&track_PVWeight);
  t->Branch((prefix+string("_track_SVWeight")).c_str(),"vector<vector<float>*>",&track_SVWeight);
  t->Branch((prefix+string("_track_charge")).c_str(),"vector<vector<float>*>",&track_charge);
  t->Branch((prefix+string("_track_nHitStrip")).c_str(),"vector<vector<int>*>",&track_nHitStrip);
  t->Branch((prefix+string("_track_nHitTIB")).c_str(),"vector<vector<int>*>",&track_nHitTIB);
  t->Branch((prefix+string("_track_nHitTID")).c_str(),"vector<vector<int>*>",&track_nHitTID);
  t->Branch((prefix+string("_track_nHitTOB")).c_str(),"vector<vector<int>*>",&track_nHitTOB);
  t->Branch((prefix+string("_track_nHitTEC")).c_str(),"vector<vector<int>*>",&track_nHitTEC);
  t->Branch((prefix+string("_track_nHitPXB")).c_str(),"vector<vector<int>*>",&track_nHitPXB);
  t->Branch((prefix+string("_track_nHitPXF")).c_str(),"vector<vector<int>*>",&track_nHitPXF);
  t->Branch((prefix+string("_track_isHitL1")).c_str(),"vector<vector<int>*>",&track_isHitL1);
  t->Branch((prefix+string("_track_PV")).c_str(),"vector<vector<int>*>",&track_PV);
  t->Branch((prefix+string("_track_fromSV")).c_str(),"vector<vector<int>*>",&track_fromSV);
  t->Branch((prefix+string("_track_SV")).c_str(),"vector<vector<int>*>",&track_SV);

}

void NeroFatJets::setBranchAddresses(TTree *t){
  //
  BareFatJets::setBranchAddresses(t,prefix);
  tau1IVF = new vector<float>;
  tau2IVF = new vector<float>;
  nSV = new vector<unsigned int>;
  zRatio = new vector<float>;
  tauDot = new vector<float>;
  svMass = new vector<vector<float>*>;
  svEnergyRatio = new vector<vector<float>*>;
  svPt = new vector<vector<float>*>;
  prunedMomentum = new vector<TLorentzVector>;
  trimmedMomentum = new vector<TLorentzVector>;
  sdMomentum = new vector<TLorentzVector>;
  prunedJEC0 = new vector<float>;
  trimmedJEC0 = new vector<float>;
  sdJEC0 = new vector<float>;
  nTracks = new vector<int>;
  track_momentum = new vector<vector<float>*>;
  track_eta = new vector<vector<float>*>;
  track_phi = new vector<vector<float>*>;
  track_ptRel = new vector<vector<float>*>;
  track_pPar = new vector<vector<float>*>;
  track_etaRel = new vector<vector<float>*>;
  track_deltaR = new vector<vector<float>*>;
  track_ptRatio = new vector<vector<float>*>;
  track_pParRatio = new vector<vector<float>*>;
  track_sip2dVal = new vector<vector<float>*>;
  track_sip2dSig = new vector<vector<float>*>;
  track_sip3dVal = new vector<vector<float>*>;
  track_sip3dSig = new vector<vector<float>*>;
  track_decayLenVal = new vector<vector<float>*>;
  track_decayLenSig = new vector<vector<float>*>;
  track_jetDistVal = new vector<vector<float>*>;
  track_jetDistSig = new vector<vector<float>*>;
  track_chi2 = new vector<vector<float>*>;
  track_nTotalHits = new vector<vector<int>*>;
  track_nPixelHits = new vector<vector<int>*>;
  track_dxy = new vector<vector<float>*>;
  track_dz = new vector<vector<float>*>;
  track_IP2D = new vector<vector<float>*>;
  track_IP2Dsig = new vector<vector<float>*>;
  track_IP = new vector<vector<float>*>;
  track_IPsig = new vector<vector<float>*>;
  track_IP2Derr = new vector<vector<float>*>;
  track_IPerr = new vector<vector<float>*>;
  track_prob = new vector<vector<float>*>;
  track_PVWeight = new vector<vector<float>*>;
  track_SVWeight = new vector<vector<float>*>;
  track_charge = new vector<vector<float>*>;
  track_nHitStrip = new vector<vector<int>*>;
  track_nHitTIB = new vector<vector<int>*>;
  track_nHitTID = new vector<vector<int>*>;
  track_nHitTOB = new vector<vector<int>*>;
  track_nHitTEC = new vector<vector<int>*>;
  track_nHitPXB = new vector<vector<int>*>;
  track_nHitPXF = new vector<vector<int>*>;
  track_isHitL1 = new vector<vector<int>*>;
  track_PV = new vector<vector<int>*>;
  track_fromSV = new vector<vector<int>*>;
  track_SV = new vector<vector<int>*>;

  t->SetBranchAddress((prefix+string("_tau1IVF")).c_str(),&tau1IVF);
  t->SetBranchAddress((prefix+string("_tau2IVF")).c_str(),&tau2IVF);
  t->SetBranchAddress((prefix+string("_nSV")).c_str(),&nSV);
  t->SetBranchAddress((prefix+string("_zRatio")).c_str(),&zRatio);
  t->SetBranchAddress((prefix+string("_tauDot")).c_str(),&tauDot);
  t->SetBranchAddress((prefix+string("_svMass")).c_str(),&svMass);
  t->SetBranchAddress((prefix+string("_svEnergyRatio")).c_str(),&svEnergyRatio);
  t->SetBranchAddress((prefix+string("_svPt")).c_str(),&svPt);
  t->SetBranchAddress((prefix+string("_prunedMomentum")).c_str(),&prunedMomentum);
  t->SetBranchAddress((prefix+string("_trimmedMomentum")).c_str(),&trimmedMomentum);
  t->SetBranchAddress((prefix+string("_sdMomentum")).c_str(),&sdMomentum);
  t->SetBranchAddress((prefix+string("_prunedJEC0")).c_str(),&prunedJEC0);
  t->SetBranchAddress((prefix+string("_trimmedJEC0")).c_str(),&trimmedJEC0);
  t->SetBranchAddress((prefix+string("_sdJEC0")).c_str(),&sdJEC0);
  t->SetBranchAddress((prefix+string("_nTracks")).c_str(),&nTracks);
  t->SetBranchAddress((prefix+string("_track_momentum")).c_str(),&track_momentum);
  t->SetBranchAddress((prefix+string("_track_eta")).c_str(),&track_eta);
  t->SetBranchAddress((prefix+string("_track_phi")).c_str(),&track_phi);
  t->SetBranchAddress((prefix+string("_track_ptRel")).c_str(),&track_ptRel);
  t->SetBranchAddress((prefix+string("_track_pPar")).c_str(),&track_pPar);
  t->SetBranchAddress((prefix+string("_track_etaRel")).c_str(),&track_etaRel);
  t->SetBranchAddress((prefix+string("_track_deltaR")).c_str(),&track_deltaR);
  t->SetBranchAddress((prefix+string("_track_ptRatio")).c_str(),&track_ptRatio);
  t->SetBranchAddress((prefix+string("_track_pParRatio")).c_str(),&track_pParRatio);
  t->SetBranchAddress((prefix+string("_track_sip2dVal")).c_str(),&track_sip2dVal);
  t->SetBranchAddress((prefix+string("_track_sip2dSig")).c_str(),&track_sip2dSig);
  t->SetBranchAddress((prefix+string("_track_sip3dVal")).c_str(),&track_sip3dVal);
  t->SetBranchAddress((prefix+string("_track_sip3dSig")).c_str(),&track_sip3dSig);
  t->SetBranchAddress((prefix+string("_track_decayLenVal")).c_str(),&track_decayLenVal);
  t->SetBranchAddress((prefix+string("_track_decayLenSig")).c_str(),&track_decayLenSig);
  t->SetBranchAddress((prefix+string("_track_jetDistVal")).c_str(),&track_jetDistVal);
  t->SetBranchAddress((prefix+string("_track_jetDistSig")).c_str(),&track_jetDistSig);
  t->SetBranchAddress((prefix+string("_track_chi2")).c_str(),&track_chi2);
  t->SetBranchAddress((prefix+string("_track_nTotalHits")).c_str(),&track_nTotalHits);
  t->SetBranchAddress((prefix+string("_track_nPixelHits")).c_str(),&track_nPixelHits);
  t->SetBranchAddress((prefix+string("_track_dxy")).c_str(),&track_dxy);
  t->SetBranchAddress((prefix+string("_track_dz")).c_str(),&track_dz);
  t->SetBranchAddress((prefix+string("_track_IP2D")).c_str(),&track_IP2D);
  t->SetBranchAddress((prefix+string("_track_IP2Dsig")).c_str(),&track_IP2Dsig);
  t->SetBranchAddress((prefix+string("_track_IP")).c_str(),&track_IP);
  t->SetBranchAddress((prefix+string("_track_IPsig")).c_str(),&track_IPsig);
  t->SetBranchAddress((prefix+string("_track_IP2Derr")).c_str(),&track_IP2Derr);
  t->SetBranchAddress((prefix+string("_track_IPerr")).c_str(),&track_IPerr);
  t->SetBranchAddress((prefix+string("_track_prob")).c_str(),&track_prob);
  t->SetBranchAddress((prefix+string("_track_PVWeight")).c_str(),&track_PVWeight);
  t->SetBranchAddress((prefix+string("_track_SVWeight")).c_str(),&track_charge);
  t->SetBranchAddress((prefix+string("_track_charge")).c_str(),&track_SVWeight);
  t->SetBranchAddress((prefix+string("_track_nHitStrip")).c_str(),&track_nHitStrip);
  t->SetBranchAddress((prefix+string("_track_nHitTIB")).c_str(),&track_nHitTIB);
  t->SetBranchAddress((prefix+string("_track_nHitTID")).c_str(),&track_nHitTID);
  t->SetBranchAddress((prefix+string("_track_nHitTOB")).c_str(),&track_nHitTOB);
  t->SetBranchAddress((prefix+string("_track_nHitTEC")).c_str(),&track_nHitTEC);
  t->SetBranchAddress((prefix+string("_track_nHitPXB")).c_str(),&track_nHitPXB);
  t->SetBranchAddress((prefix+string("_track_nHitPXF")).c_str(),&track_nHitPXF);
  t->SetBranchAddress((prefix+string("_track_isHitL1")).c_str(),&track_isHitL1);
  t->SetBranchAddress((prefix+string("_track_PV")).c_str(),&track_PV);
  t->SetBranchAddress((prefix+string("_track_fromSV")).c_str(),&track_fromSV);
  t->SetBranchAddress((prefix+string("_track_SV")).c_str(),&track_SV);

}

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
