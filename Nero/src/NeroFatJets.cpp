#include "NeroProducer/Nero/interface/NeroFatJets.hpp"
#include "NeroProducer/Nero/interface/NeroJets.hpp" // JetId
#include "NeroProducer/Nero/interface/Nero.hpp"

NeroFatJets::NeroFatJets() : BareFatJets(){
}

NeroFatJets::~NeroFatJets(){
}

int NeroFatJets::analyze(const edm::Event& iEvent){

    if ( mOnlyMc  ) return 0;

    // maybe handle should be taken before
    iEvent.getByToken(token, handle);

    int ijetRef = -1;
    int nsubjet = 0;
    for (const pat::Jet& j : *handle)
    {
        ijetRef++;
        if (j.pt() < 100 ) continue;

        // JET ID
        if ( !NeroJets::JetId(j,"loose") ) continue;

        // GET  ValueMaps

        // Fill output object
        //p4 -> AddLast(new TLorentzVector(j.px(), j.py(), j.pz(), j.energy())  );
        new ( (*p4)[p4->GetEntriesFast()]) TLorentzVector(j.px(), j.py(), j.pz(), j.energy());

        rawPt -> push_back (j.pt()*j.jecFactor("Uncorrected"));
        flavour -> push_back( j.partonFlavour() );
        tau1 -> push_back(j.userFloat("Njettiness:tau1"));
        tau2 -> push_back(j.userFloat("Njettiness:tau2"));
        tau3 -> push_back(j.userFloat("Njettiness:tau3"));

        trimmedMass ->push_back(j.userFloat("PFJetsCHSTrimmedMass"));
        prunedMass  ->push_back(j.userFloat("PFJetsCHSPrunedMass"));
        filteredMass->push_back(j.userFloat("PFJetsCHSFilteredMass"));
        softdropMass->push_back(j.userFloat("PFJetsCHSSoftDropMass"));
        hasSubjet->push_back(j.hasSubjets(fSubjetsName));

        auto Subjets = j.subjets(fSubjetsName);
        for ( auto const & i : Subjets ) {
            new ( (*subjets)[nsubjet]) TLorentzVector(i->px(), i->py(), i->pz(), i->energy());
            subjetBtag->push_back(i->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            nsubjet++;
        }

    }
    return 0;
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
