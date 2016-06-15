#ifndef NERO_PUPPIJETS_H
#define NERO_PUPPIJETS_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BarePuppiJets.hpp"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

//
#include "NeroProducer/Nero/interface/NeroPF.hpp"


class NeroPuppiJets : virtual public NeroCollection, virtual public BarePuppiJets
{
    public:
        NeroPuppiJets();
        ~NeroPuppiJets();
        int analyze(const edm::Event& iEvent, const edm::EventSetup&iSetup);
        int analyze(const edm::Event& iEvent){return 2;} // never called
        virtual inline string name(){return "NeroPuppiJets";};
        void init() override;

        // --- specific fuctions
        static bool JetId(const pat::Jet &, string id);

        // --- Handle
        edm::Handle<pat::JetCollection> handle;	
        edm::Handle<double> rho_handle;

        // --- Token
        edm::EDGetTokenT<pat::JetCollection> token;
        edm::EDGetTokenT<double> rho_token;

        // --- configuration
        float mMinPt;
        int   mMinNjets;
        float mMinEta;
        string mMinId;
        bool reclustered=false; // means we need to apply JEC
        FactorizedJetCorrector *mMCJetCorrector;   // needed for puppi fat jets
        FactorizedJetCorrector *mDataJetCorrector; 

        // extra info
        NeroPF *pf;

};


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
