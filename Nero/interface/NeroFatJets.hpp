#ifndef NERO_FATJETS_H
#define NERO_FATJETS_H

#include "NeroProducer/Nero/interface/NeroCollection.hpp"
#include "NeroProducer/Core/interface/BareFatJets.hpp"


class NeroFatJets : virtual public NeroCollection,
    virtual public BareFatJets
{
    public:
        NeroFatJets();
        ~NeroFatJets();
        int analyze(const edm::Event& iEvent);
        virtual inline string name(){return "NeroFatJets";};

        // --- specific fuctions
        // --- Handle
        edm::Handle<pat::JetCollection> handle;

        // --- Token
        edm::EDGetTokenT<pat::JetCollection> token;

        void SetSubjetsName (std::string n) { fSubjetsName = n;    }
        std::string SubjetsName()           { return fSubjetsName; }

      protected:
        std::string fSubjetsName;



};


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
