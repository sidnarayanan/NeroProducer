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

      std::vector<vector<float>*> *trackMomentum;       // per track in fatjet
      std::vector<vector<float>*> *trackEta;
      std::vector<vector<float>*> *trackPhi;
      std::vector<vector<float>*> *trackPtRel;
      std::vector<vector<float>*> *trackPPar;
      std::vector<vector<float>*> *trackEtaRel;
      std::vector<vector<float>*> *trackDeltaR;
      std::vector<vector<float>*> *trackPtRatio;
      std::vector<vector<float>*> *trackPParRatio;
      std::vector<vector<float>*> *trackSip2dVal;
      std::vector<vector<float>*> *trackSip2dSig;
      std::vector<vector<float>*> *trackSip3dVal;
      std::vector<vector<float>*> *trackSip3dSig;
      std::vector<vector<float>*> *trackDecayLenVal;
      std::vector<vector<float>*> *trackDecayLenSig;
      std::vector<vector<float>*> *trackJetDistVal;
      std::vector<vector<float>*> *trackJetDistSig;
      std::vector<vector<float>*> *trackChi2;
      std::vector<vector<int>*> *trackNTotalHits;
      std::vector<vector<int>*> *trackNPixelHits;


      //grooming vars
      std::vector<float> *prunedPt;
      std::vector<float> *prunedEta;
      std::vector<float> *prunedPhi;
      std::vector<float> *prunedJEC0;
      std::vector<float> *trimmedPt;
      std::vector<float> *trimmedEta;
      std::vector<float> *trimmedPhi;
      std::vector<float> *trimmedJEC0;
      std::vector<float> *sdPt;
      std::vector<float> *sdEta;
      std::vector<float> *sdPhi;
      std::vector<float> *sdJEC0;

};


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
