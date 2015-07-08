#ifndef BARE_FATJETS_H
#define BARE_FATJETS_H

#include "NeroProducer/Core/interface/BareCollection.hpp"
#include "NeroProducer/Core/interface/BareP4.hpp"


class BareFatJets : virtual public BareP4
{
    public:
        BareFatJets();
        ~BareFatJets();
        virtual void clear();
        virtual void defineBranches(TTree *t, std::string prefix="fatjet");
        virtual void setBranchAddresses(TTree* t, std::string prefix="fatjet");
        virtual inline string name(){return "BareFatJets";};

        // -- variables
        //TClonesArray  *p4;
        vector<float> *rawPt;
        vector<int>   *flavour;
        vector<float> *tau1;
        vector<float> *tau2;
        vector<float> *tau3;

        vector<float> *trimmedMass;
        vector<float> *prunedMass;
        vector<float> *filteredMass;
        vector<float> *softdropMass;

        TClonesArray  *subjets;
        vector<int>   *hasSubjet;
        vector<float> *subjetBtag;

};

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
