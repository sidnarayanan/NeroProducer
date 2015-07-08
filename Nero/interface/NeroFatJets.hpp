#ifndef NERO_FATJETS_H
#define NERO_FATJETS_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BareFatJets.hpp"
#include "RecoBTag/SecondaryVertex/interface/TemplatedSimpleSecondaryVertexComputer.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
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
      int analyze(const edm::Event& iEvent);
      virtual inline string name(){return "NeroFatJets";};
      void setBranchAddresses(TTree *);
      void defineBranches(TTree *);
      void clear();


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

      void doBTagging(const pat::Jet* jet);
      void vertexKinematicsAndCharge(const Vertex & vertex, reco::TrackKinematics & vertexKinematics, Int_t & charge);
      void recalcNsubjettiness(const pat::Jet &jet,
                                            const SVTagInfo &svTagInfo,
                                            float & tau1,
                                            float & tau2,
                                            std::vector<fastjet::PseudoJet> &currentAxes);
      void setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);
      void setTracksPVBase(const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight);


    protected:
      std::string SubjetsName = "SubJets";
      float R0                = -1;
      fastjet::contrib::Njettiness njettiness;
      std::string prefix = "fatjet";
      

      std::vector<float> *tau1IVF;
      std::vector<float> *tau2IVF;
      std::vector<unsigned int> *nSV;
      std::vector<float> *zRatio;
      std::vector<float> *tauDot;
      std::vector<float> *svMass0;
      std::vector<float> *svEnergyRatio0;
      std::vector<float> *svEnergyRatio1;
      std::vector<float> *svPt0;



};


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
