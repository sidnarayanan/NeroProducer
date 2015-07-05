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

void BareFatJets::defineBranches(TTree *t){
    //
    BareP4::defineBranches(t, "fatjet" );
    //
    rawPt = new vector<float>;
    t->Branch("fatjetRawPt","vector<float>",&rawPt);
    // -- Jet Flavour by PAT
    flavour = new vector<int>;
    t->Branch("fatjetFlavour","vector<int>",&flavour);
    //
    tau1 = new vector<float>;
    t->Branch("fatjetTau1","vector<float>",&tau1);
    tau2 = new vector<float>;
    t->Branch("fatjetTau2","vector<float>",&tau2);
    tau3 = new vector<float>;
    t->Branch("fatjetTau3","vector<float>",&tau3);

    //
    trimmedMass = new vector<float>;
    t->Branch("fatjetTrimmedMass","vector<float>",&trimmedMass);
    prunedMass = new vector<float>;
    t->Branch("fatjetPrunedMass","vector<float>",&prunedMass);
    filteredMass = new vector<float>;
    t->Branch("fatjetFilteredMass","vector<float>",&filteredMass);
    softdropMass = new vector<float>;
    t->Branch("fatjetSoftdropMass","vector<float>",&softdropMass);

    subjets = new TClonesArray("TLorentzVector", 20);
    t->Branch("subjets","TClonesArray", &subjets, 128000, 0);
    hasSubjet =  new vector<int>;
    t->Branch("hasSubjet","vector<int>",&hasSubjet);
    subjetBtag =  new vector<float>;
    t->Branch("subjetBtag","vector<float>",&subjetBtag);

}

void BareFatJets::setBranchAddresses(TTree *t){
    //
    BareP4::setBranchAddresses(t,"fatjet");

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

    t->SetBranchAddress("fatjetRawPt"	,&rawPt);
    t->SetBranchAddress("fatjetFlavour" ,&flavour);
    t->SetBranchAddress("fatjetTau1"	,&tau1);
    t->SetBranchAddress("fatjetTau2"	,&tau2);
    t->SetBranchAddress("fatjetTau3"	,&tau3);

    t->SetBranchAddress("fatjetTrimmedMass"	,&trimmedMass);
    t->SetBranchAddress("fatjetPrunedMass"	,&prunedMass);
    t->SetBranchAddress("fatjetFilteredMass"	,&filteredMass);
    t->SetBranchAddress("fatjetSoftdropMass"	,&softdropMass);

    subjets = new TClonesArray("TLorentzVector", 20);
    t->SetBranchAddress("subjets"	,&subjets);
    t->SetBranchAddress("hasSubjet",&hasSubjet);
    t->SetBranchAddress("subjetBtag",&subjetBtag);
}
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
