#ifndef BARE_TRIGGER_H
#define BARE_TRIGGER_H


#include "NeroProducer/Core/interface/BareCollection.hpp"


class BareTrigger: virtual public BareCollection
{
    public:
        BareTrigger();
        ~BareTrigger();
        void clear();
        void defineBranches(TTree*);
        virtual void setBranchAddresses(TTree*);
        virtual inline string name(){ return "BareTrigger";}

        //variables
        vector<int>  *triggerFired; // bool vector are not supported in TTree
        vector<float>  *triggerPrescale; // bool vector are not supported in TTree

        //configuration
        vector<string> *triggerNames;

        vector<int>  *triggerLeps;
        vector<int>  *triggerJets;
        vector<int>  *triggerTaus;
        vector<int>  *triggerPhotons;

};

#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
