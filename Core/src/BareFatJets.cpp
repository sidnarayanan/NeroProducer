#include "NeroProducer/Core/interface/BareFatJets.hpp"


BareFatJets::BareFatJets(){
    p4 = NULL;
    rawPt = NULL;
    flavour = NULL;
    tau1 = NULL;
    tau2 = NULL;
    tau3 = NULL;
    trimmedMass = NULL;
    prunedMass = NULL;
    filteredMass = NULL;
    softdropMass = NULL;
    subjets = NULL;
    hasSubjet = NULL;
    subjetBtag = NULL;

}

BareFatJets::~BareFatJets(){
}

void BareFatJets::clear(){
    // This function clear all the internal storage and init it to an arbitrary value
    BareP4::clear();
    p4 -> Clear();
    rawPt -> clear();
    flavour -> clear();
    tau1 -> clear();
    tau2 -> clear();
    tau3 -> clear();
    trimmedMass -> clear();
    prunedMass -> clear();
    filteredMass -> clear();
    softdropMass -> clear();
    subjets->Clear();
    subjetBtag ->clear();
    hasSubjet->clear();
}

void BareFatJets::defineBranches(TTree *t, std::string prefix){
    //
    BareP4::defineBranches(t, prefix );
    //

    rawPt = new vector<float>;
    flavour = new vector<int>;
    tau1 = new vector<float>;
    tau2 = new vector<float>;
    tau3 = new vector<float>;
    trimmedMass = new vector<float>;
    prunedMass = new vector<float>;
    filteredMass = new vector<float>;
    softdropMass = new vector<float>;
    subjets = new TClonesArray("TLorentzVector", 20);
    hasSubjet = new vector<int>;
    subjetBtag = new vector<float>;
    t->Branch((prefix+string("_rawPt")).c_str(),"vector<float>",&rawPt);
    t->Branch((prefix+string("_flavour")).c_str(),"vector<int>",&flavour);
    t->Branch((prefix+string("_tau1")).c_str(),"vector<float>",&tau1);
    t->Branch((prefix+string("_tau2")).c_str(),"vector<float>",&tau2);
    t->Branch((prefix+string("_tau3")).c_str(),"vector<float>",&tau3);
    t->Branch((prefix+string("_trimmedMass")).c_str(),"vector<float>",&trimmedMass);
    t->Branch((prefix+string("_prunedMass")).c_str(),"vector<float>",&prunedMass);
    t->Branch((prefix+string("_filteredMass")).c_str(),"vector<float>",&filteredMass);
    t->Branch((prefix+string("_softdropMass")).c_str(),"vector<float>",&softdropMass);
    t->Branch((prefix+string("_subjets")).c_str(),"TClonesArray",&subjets, 128000, 0);
    t->Branch((prefix+string("_hasSubjet")).c_str(),"vector<int>",&hasSubjet);
    t->Branch((prefix+string("_subjetBtag")).c_str(),"vector<float>",&subjetBtag);

}

void BareFatJets::setBranchAddresses(TTree *t,string prefix){
    //
    BareP4::setBranchAddresses(t,prefix.c_str());

    rawPt = new vector<float>;
    // -- Jet Flavour by PAT
    flavour = new vector<int>;
    //
    tau1 = new vector<float>;
    tau2 = new vector<float>;
    tau3 = new vector<float>;

    trimmedMass = new vector<float>;
    prunedMass = new vector<float>;
    filteredMass = new vector<float>;
    softdropMass = new vector<float>;
    subjets = new TClonesArray("TLorentzVector", 20);
    hasSubjet =  new vector<int>;
    subjetBtag =  new vector<float>;

    t->SetBranchAddress((prefix+string("_rawPt")).c_str(),&rawPt);
    t->SetBranchAddress((prefix+string("_flavour")).c_str(),&flavour);
    t->SetBranchAddress((prefix+string("_tau1")).c_str(),&tau1);
    t->SetBranchAddress((prefix+string("_tau2")).c_str(),&tau2);
    t->SetBranchAddress((prefix+string("_tau3")).c_str(),&tau3);
    t->SetBranchAddress((prefix+string("_trimmedMass")).c_str(),&trimmedMass);
    t->SetBranchAddress((prefix+string("_prunedMass")).c_str(),&prunedMass);
    t->SetBranchAddress((prefix+string("_filteredMass")).c_str(),&filteredMass);
    t->SetBranchAddress((prefix+string("_softdropMass")).c_str(),&softdropMass);
    t->SetBranchAddress((prefix+string("_subjets")).c_str(),&subjets);
    t->SetBranchAddress((prefix+string("_hasSubjet")).c_str(),&hasSubjet);
    t->SetBranchAddress((prefix+string("_subjetBtag")).c_str(),&subjetBtag);

}
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
