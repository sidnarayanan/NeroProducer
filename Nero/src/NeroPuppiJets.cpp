#include "NeroProducer/Nero/interface/NeroPuppiJets.hpp"
#include "NeroProducer/Nero/interface/Nero.hpp"
#include "NeroProducer/Nero/interface/NeroJets.hpp" // JetId
#include "NeroProducer/Core/interface/BareFunctions.hpp"
#include <cstdlib>
#include <string>

//JES


NeroPuppiJets::NeroPuppiJets() : 
        NeroCollection(),
        BarePuppiJets(),
        mMCJetCorrector(0),
        mDataJetCorrector(0)
{
    mMinPt = 25.;
    mMinNjets = 0;
    mMinEta = 4.7;
    mMinId = "noid";
}

NeroPuppiJets::~NeroPuppiJets(){
  BareFunctions::Delete(mMCJetCorrector);
  BareFunctions::Delete(mDataJetCorrector);
}

void NeroPuppiJets::init()
{
   BarePuppiJets::init();
   std::string jecDir = "jec/";

   std::vector<JetCorrectorParameters> mcParams;
   mcParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_MC_L1FastJet_AK4PFPuppi.txt"));
   mcParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_MC_L2Relative_AK4PFPuppi.txt"));
   mcParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_MC_L3Absolute_AK4PFPuppi.txt"));
   mcParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_MC_L2L3Residual_AK4PFPuppi.txt"));
   mMCJetCorrector = new FactorizedJetCorrector(mcParams);
 
   std::vector<JetCorrectorParameters> dataParams;
   dataParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_DATA_L1FastJet_AK4PFPuppi.txt"));
   dataParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_DATA_L2Relative_AK4PFPuppi.txt"));
   dataParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_DATA_L3Absolute_AK4PFPuppi.txt"));
   dataParams.push_back(JetCorrectorParameters(jecDir + "Fall15_25nsV2_DATA_L2L3Residual_AK4PFPuppi.txt"));
   mDataJetCorrector = new FactorizedJetCorrector(dataParams);
}


int NeroPuppiJets::analyze(const edm::Event& iEvent, const edm::EventSetup &iSetup){

    if ( mOnlyMc  ) return 0;

    // maybe handle should be taken before
    iEvent.getByToken(token, handle);
    iEvent.getByToken(rho_token,rho_handle);

    if ( not handle.isValid() ) cout<<"[NeroPuppiJets]::[analyze]::[ERROR] handle is not valid"<<endl;
    assert(rho_handle.isValid());

    FactorizedJetCorrector *corrector = ( iEvent.isRealData() ) ? mDataJetCorrector : mMCJetCorrector;
    for (const pat::Jet& j : *handle)
    {

        if ( fabs(j.eta()) > mMinEta ) continue;
        if ( !JetId(j,mMinId) ) continue;

        // do JEC now
        double this_pt = j.pt(), this_rawpt=0, jecFactor=1;
        if (reclustered) {
            this_rawpt = this_pt;
            if (fabs(j.eta())<5.191) {
              corrector->setJetPt(j.pt());
              corrector->setJetEta(j.eta());
              corrector->setJetPhi(j.phi());
              corrector->setJetE(j.energy());
              corrector->setRho(*rho_handle);
              corrector->setJetA(j.jetArea());
              corrector->setJetEMF(-99.0);
              jecFactor = corrector->getCorrection();
              this_pt *= jecFactor;
            }
        } else {
            this_rawpt = j.pt()*j.jecFactor("Uncorrected");
        }

        if (this_pt < mMinPt || this_rawpt < mMinPt) continue;

        new ( (*p4)[p4->GetEntriesFast()]) TLorentzVector(j.px()*jecFactor, j.py()*jecFactor, j.pz()*jecFactor, j.energy()*jecFactor);
        rawPt  -> push_back (this_rawpt);

        float charge =  0.;
        float charge_den = 0.;

        for( size_t idx =0; idx < j.numberOfDaughters() ; ++idx)
        {
            pat::PackedCandidate *cand = ( pat::PackedCandidate* ) j.daughter(idx) ; 
            if (cand->charge() !=0 ) {  
                charge     += cand->charge() * cand->puppiWeight() * ( j.px()*cand->px() + j.py()*cand->py() + j.pz()*cand->pz()  ) ;
                charge_den +=                  cand->puppiWeight() * ( j.px()*cand->px() + j.py()*cand->py() + j.pz()*cand->pz()  ) ;
            }
        }

        if ( charge_den == 0 ) { charge=0.0 ; charge_den =1.0;}  //  guard, if no jet id

        // Fill output object   
        bDiscr -> push_back( j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

        unsigned bits=0;
        bits |=  (1 * JetBaseline);
        bits |= JetId(j,"monojet") * mjId;
        bits |= JetId(j,"monojetloose") * mjIdLoose;
        bits |= JetId(j,"monojet2015") * mjId2015;
        bits |= JetId(j,"loose") * JetLoose;
        bits |= JetId(j,"tight") * JetTight;

        selBits -> push_back( bits);
        Q -> push_back (charge/charge_den);
    }

    if ( int(rawPt -> size()) < mMinNjets ) return 1;

    return 0;
}

bool NeroPuppiJets::JetId(const pat::Jet &j, std::string id)
{
    bool jetid = false;

    //float NHF    = j.neutralHadronEnergyFraction();
    //float NEMF   = j.neutralEmEnergyFraction();
    //float CHF    = j.chargedHadronEnergyFraction();
    //float MUF    = j.muonEnergyFraction();
    //float CEMF   = j.chargedEmEnergyFraction();
    //int NumConst = j.chargedMultiplicity()+j.neutralMultiplicity();
    //int CHM      = j.chargedMultiplicity();
    //int NumNeutralParticle =j.neutralMultiplicity(); 
    //float eta = j.eta();

    jetid= NeroJets::JetId(j,id);

    return jetid;
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
