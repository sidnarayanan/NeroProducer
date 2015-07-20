// system include files
#include <memory>

#include "NeroProducer/Nero/interface/NeroFatJets.hpp"
#include "NeroProducer/Nero/interface/NeroJets.hpp" // JetId
#include "NeroProducer/Nero/interface/Nero.hpp"

#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/CandIPTagInfo.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/CandSoftLeptonTagInfo.h"
#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/GenericMVAJetTagComputerWrapper.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputer.h"
#include "RecoBTau/JetTagComputer/interface/JetTagComputerRecord.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"


NeroFatJets::NeroFatJets(double r0) : BareFatJets(),
                             R0(r0),
                             njettiness(fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1.0,R0)),
                             tau1IVF(0),
                             tau2IVF(0),
                             nSV(0),
                             zRatio(0),
                             tauDot(0),
                             svMass0(0),
                             svEnergyRatio0(0),
                             svEnergyRatio1(0),
                             svPt0(0),
                             computer(0)
{
}

NeroFatJets::~NeroFatJets(){
}

void NeroFatJets::clear() {
  BareFatJets::clear();
  tau1IVF->clear();
  tau2IVF->clear();
  nSV->clear();
  zRatio->clear();
  tauDot->clear();
  svMass->clear();
  svEnergyRatio->clear();
  svPt->clear();
  prunedPt->clear();
  prunedEta->clear();
  prunedPhi->clear();
  prunedJEC0->clear();
  trimmedPt->clear();
  trimmedEta->clear();
  trimmedPhi->clear();
  trimmedJEC0->clear();
  sdPt->clear();
  sdEta->clear();
  sdPhi->clear();
  sdJEC0->clear();
  nTracks->clear();
  trackMomentum->clear();
  trackEta->clear();
  trackPhi->clear();
  trackPtRel->clear();
  trackPPar->clear();
  trackEtaRel->clear();
  trackDeltaR->clear();
  trackPtRatio->clear();
  trackPParRatio->clear();
  trackSip2dVal->clear();
  trackSip2dSig->clear();
  trackSip3dVal->clear();
  trackSip3dSig->clear();
  trackDecayLenVal->clear();
  trackDecayLenSig->clear();
  trackJetDistVal->clear();
  trackJetDistSig->clear();
  trackChi2->clear();
  trackNTotalHits->clear();
  trackNPixelHits->clear();

}

void NeroFatJets::defineBranches(TTree * t){
  BareFatJets::defineBranches(t,prefix);
  tau1IVF = new vector<float>;
  tau2IVF = new vector<float>;
  nSV = new vector<unsigned int>;
  zRatio = new vector<float>;
  tauDot = new vector<float>;
  svMass = new vector<vector<float>*>;
  svEnergyRatio = new vector<vector<float>*>;
  svPt = new vector<vector<float>*>;
  prunedPt = new vector<float>;
  prunedEta = new vector<float>;
  prunedPhi = new vector<float>;
  prunedJEC0 = new vector<float>;
  trimmedPt = new vector<float>;
  trimmedEta = new vector<float>;
  trimmedPhi = new vector<float>;
  trimmedJEC0 = new vector<float>;
  sdPt = new vector<float>;
  sdEta = new vector<float>;
  sdPhi = new vector<float>;
  sdJEC0 = new vector<float>;
  nTracks = new vector<int>;
  trackMomentum = new vector<vector<float>*>;
  trackEta = new vector<vector<float>*>;
  trackPhi = new vector<vector<float>*>;
  trackPtRel = new vector<vector<float>*>;
  trackPPar = new vector<vector<float>*>;
  trackEtaRel = new vector<vector<float>*>;
  trackDeltaR = new vector<vector<float>*>;
  trackPtRatio = new vector<vector<float>*>;
  trackPParRatio = new vector<vector<float>*>;
  trackSip2dVal = new vector<vector<float>*>;
  trackSip2dSig = new vector<vector<float>*>;
  trackSip3dVal = new vector<vector<float>*>;
  trackSip3dSig = new vector<vector<float>*>;
  trackDecayLenVal = new vector<vector<float>*>;
  trackDecayLenSig = new vector<vector<float>*>;
  trackJetDistVal = new vector<vector<float>*>;
  trackJetDistSig = new vector<vector<float>*>;
  trackChi2 = new vector<vector<float>*>;
  trackNTotalHits = new vector<vector<int>*>;
  trackNPixelHits = new vector<vector<int>*>;

  t->Branch((prefix+string("_tau1IVF")).c_str(),"vector<float>",&tau1IVF);
  t->Branch((prefix+string("_tau2IVF")).c_str(),"vector<float>",&tau2IVF);
  t->Branch((prefix+string("_nSV")).c_str(),"vector<unsigned int>",&nSV);
  t->Branch((prefix+string("_zRatio")).c_str(),"vector<float>",&zRatio);
  t->Branch((prefix+string("_tauDot")).c_str(),"vector<float>",&tauDot);
  t->Branch((prefix+string("_svMass")).c_str(),"vector<vector<float>*>",&svMass);
  t->Branch((prefix+string("_svEnergyRatio")).c_str(),"vector<vector<float>*>",&svEnergyRatio);
  t->Branch((prefix+string("_svPt")).c_str(),"vector<vector<float>*>",&svPt);
  t->Branch((prefix+string("_prunedPt")).c_str(),"vector<float>",&prunedPt);
  t->Branch((prefix+string("_prunedEta")).c_str(),"vector<float>",&prunedEta);
  t->Branch((prefix+string("_prunedPhi")).c_str(),"vector<float>",&prunedPhi);
  t->Branch((prefix+string("_prunedJEC0")).c_str(),"vector<float>",&prunedJEC0);
  t->Branch((prefix+string("_trimmedPt")).c_str(),"vector<float>",&trimmedPt);
  t->Branch((prefix+string("_trimmedEta")).c_str(),"vector<float>",&trimmedEta);
  t->Branch((prefix+string("_trimmedPhi")).c_str(),"vector<float>",&trimmedPhi);
  t->Branch((prefix+string("_trimmedJEC0")).c_str(),"vector<float>",&trimmedJEC0);
  t->Branch((prefix+string("_sdPt")).c_str(),"vector<float>",&sdPt);
  t->Branch((prefix+string("_sdEta")).c_str(),"vector<float>",&sdEta);
  t->Branch((prefix+string("_sdPhi")).c_str(),"vector<float>",&sdPhi);
  t->Branch((prefix+string("_sdJEC0")).c_str(),"vector<float>",&sdJEC0);
  t->Branch((prefix+string("_nTracks")).c_str(),"vector<int>",&nTracks);
  t->Branch((prefix+string("_trackMomentum")).c_str(),"vector<vector<float>*>",&trackMomentum);
  t->Branch((prefix+string("_trackEta")).c_str(),"vector<vector<float>*>",&trackEta);
  t->Branch((prefix+string("_trackPhi")).c_str(),"vector<vector<float>*>",&trackPhi);
  t->Branch((prefix+string("_trackPtRel")).c_str(),"vector<vector<float>*>",&trackPtRel);
  t->Branch((prefix+string("_trackPPar")).c_str(),"vector<vector<float>*>",&trackPPar);
  t->Branch((prefix+string("_trackEtaRel")).c_str(),"vector<vector<float>*>",&trackEtaRel);
  t->Branch((prefix+string("_trackDeltaR")).c_str(),"vector<vector<float>*>",&trackDeltaR);
  t->Branch((prefix+string("_trackPtRatio")).c_str(),"vector<vector<float>*>",&trackPtRatio);
  t->Branch((prefix+string("_trackPParRatio")).c_str(),"vector<vector<float>*>",&trackPParRatio);
  t->Branch((prefix+string("_trackSip2dVal")).c_str(),"vector<vector<float>*>",&trackSip2dVal);
  t->Branch((prefix+string("_trackSip2dSig")).c_str(),"vector<vector<float>*>",&trackSip2dSig);
  t->Branch((prefix+string("_trackSip3dVal")).c_str(),"vector<vector<float>*>",&trackSip3dVal);
  t->Branch((prefix+string("_trackSip3dSig")).c_str(),"vector<vector<float>*>",&trackSip3dSig);
  t->Branch((prefix+string("_trackDecayLenVal")).c_str(),"vector<vector<float>*>",&trackDecayLenVal);
  t->Branch((prefix+string("_trackDecayLenSig")).c_str(),"vector<vector<float>*>",&trackDecayLenSig);
  t->Branch((prefix+string("_trackJetDistVal")).c_str(),"vector<vector<float>*>",&trackJetDistVal);
  t->Branch((prefix+string("_trackJetDistSig")).c_str(),"vector<vector<float>*>",&trackJetDistSig);
  t->Branch((prefix+string("_trackChi2")).c_str(),"vector<vector<float>*>",&trackChi2);
  t->Branch((prefix+string("_trackNTotalHits")).c_str(),"vector<vector<int>*>",&trackNTotalHits);
  t->Branch((prefix+string("_trackNPixelHits")).c_str(),"vector<vector<int>*>",&trackNPixelHits);

}

void NeroFatJets::setBranchAddresses(TTree *t){
  //
  BareFatJets::setBranchAddresses(t,prefix);
  tau1IVF = new vector<float>;
  tau2IVF = new vector<float>;
  nSV = new vector<unsigned int>;
  zRatio = new vector<float>;
  tauDot = new vector<float>;
  svMass = new vector<vector<float>*>;
  svEnergyRatio = new vector<vector<float>*>;
  svPt = new vector<vector<float>*>;
  prunedPt = new vector<float>;
  prunedEta = new vector<float>;
  prunedPhi = new vector<float>;
  prunedJEC0 = new vector<float>;
  trimmedPt = new vector<float>;
  trimmedEta = new vector<float>;
  trimmedPhi = new vector<float>;
  trimmedJEC0 = new vector<float>;
  sdPt = new vector<float>;
  sdEta = new vector<float>;
  sdPhi = new vector<float>;
  sdJEC0 = new vector<float>;
  nTracks = new vector<int>;
  trackMomentum = new vector<vector<float>*>;
  trackEta = new vector<vector<float>*>;
  trackPhi = new vector<vector<float>*>;
  trackPtRel = new vector<vector<float>*>;
  trackPPar = new vector<vector<float>*>;
  trackEtaRel = new vector<vector<float>*>;
  trackDeltaR = new vector<vector<float>*>;
  trackPtRatio = new vector<vector<float>*>;
  trackPParRatio = new vector<vector<float>*>;
  trackSip2dVal = new vector<vector<float>*>;
  trackSip2dSig = new vector<vector<float>*>;
  trackSip3dVal = new vector<vector<float>*>;
  trackSip3dSig = new vector<vector<float>*>;
  trackDecayLenVal = new vector<vector<float>*>;
  trackDecayLenSig = new vector<vector<float>*>;
  trackJetDistVal = new vector<vector<float>*>;
  trackJetDistSig = new vector<vector<float>*>;
  trackChi2 = new vector<vector<float>*>;
  trackNTotalHits = new vector<vector<int>*>;
  trackNPixelHits = new vector<vector<int>*>;

  t->SetBranchAddress((prefix+string("_tau1IVF")).c_str(),&tau1IVF);
  t->SetBranchAddress((prefix+string("_tau2IVF")).c_str(),&tau2IVF);
  t->SetBranchAddress((prefix+string("_nSV")).c_str(),&nSV);
  t->SetBranchAddress((prefix+string("_zRatio")).c_str(),&zRatio);
  t->SetBranchAddress((prefix+string("_tauDot")).c_str(),&tauDot);
  t->SetBranchAddress((prefix+string("_svMass")).c_str(),&svMass);
  t->SetBranchAddress((prefix+string("_svEnergyRatio")).c_str(),&svEnergyRatio);
  t->SetBranchAddress((prefix+string("_svPt")).c_str(),&svPt);
  t->SetBranchAddress((prefix+string("_prunedPt")).c_str(),&prunedPt);
  t->SetBranchAddress((prefix+string("_prunedEta")).c_str(),&prunedEta);
  t->SetBranchAddress((prefix+string("_prunedPhi")).c_str(),&prunedPhi);
  t->SetBranchAddress((prefix+string("_prunedJEC0")).c_str(),&prunedJEC0);
  t->SetBranchAddress((prefix+string("_trimmedPt")).c_str(),&trimmedPt);
  t->SetBranchAddress((prefix+string("_trimmedEta")).c_str(),&trimmedEta);
  t->SetBranchAddress((prefix+string("_trimmedPhi")).c_str(),&trimmedPhi);
  t->SetBranchAddress((prefix+string("_trimmedJEC0")).c_str(),&trimmedJEC0);
  t->SetBranchAddress((prefix+string("_sdPt")).c_str(),&sdPt);
  t->SetBranchAddress((prefix+string("_sdEta")).c_str(),&sdEta);
  t->SetBranchAddress((prefix+string("_sdPhi")).c_str(),&sdPhi);
  t->SetBranchAddress((prefix+string("_sdJEC0")).c_str(),&sdJEC0);
  t->SetBranchAddress((prefix+string("_nTracks")).c_str(),&nTracks);
  t->SetBranchAddress((prefix+string("_trackMomentum")).c_str(),&trackMomentum);
  t->SetBranchAddress((prefix+string("_trackEta")).c_str(),&trackEta);
  t->SetBranchAddress((prefix+string("_trackPhi")).c_str(),&trackPhi);
  t->SetBranchAddress((prefix+string("_trackPtRel")).c_str(),&trackPtRel);
  t->SetBranchAddress((prefix+string("_trackPPar")).c_str(),&trackPPar);
  t->SetBranchAddress((prefix+string("_trackEtaRel")).c_str(),&trackEtaRel);
  t->SetBranchAddress((prefix+string("_trackDeltaR")).c_str(),&trackDeltaR);
  t->SetBranchAddress((prefix+string("_trackPtRatio")).c_str(),&trackPtRatio);
  t->SetBranchAddress((prefix+string("_trackPParRatio")).c_str(),&trackPParRatio);
  t->SetBranchAddress((prefix+string("_trackSip2dVal")).c_str(),&trackSip2dVal);
  t->SetBranchAddress((prefix+string("_trackSip2dSig")).c_str(),&trackSip2dSig);
  t->SetBranchAddress((prefix+string("_trackSip3dVal")).c_str(),&trackSip3dVal);
  t->SetBranchAddress((prefix+string("_trackSip3dSig")).c_str(),&trackSip3dSig);
  t->SetBranchAddress((prefix+string("_trackDecayLenVal")).c_str(),&trackDecayLenVal);
  t->SetBranchAddress((prefix+string("_trackDecayLenSig")).c_str(),&trackDecayLenSig);
  t->SetBranchAddress((prefix+string("_trackJetDistVal")).c_str(),&trackJetDistVal);
  t->SetBranchAddress((prefix+string("_trackJetDistSig")).c_str(),&trackJetDistSig);
  t->SetBranchAddress((prefix+string("_trackChi2")).c_str(),&trackChi2);
  t->SetBranchAddress((prefix+string("_trackNTotalHits")).c_str(),&trackNTotalHits);
  t->SetBranchAddress((prefix+string("_trackNPixelHits")).c_str(),&trackNPixelHits);

}


int NeroFatJets::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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
        hasSubjet->push_back(j.hasSubjets(SubjetsName));


        // B Tagging
        if (doSubstructure) {
          runBTagging(&j);

          // get groomed jet infos
          TString sR0 = TString::Format("%i",int(R0*10));
          TString sPt(":Pt");  TString sEta(":Eta");  TString sPhi(":Phi"); TString sJEC0(":jecFactor0");
          TString sPruned = TString("Pruned") + sR0;
          TString sTrimmed = TString("Trimmed") + sR0;
          TString sSD = TString("SoftDrop") + sR0;

          prunedPt->push_back(j.userFloat( (sPruned+sPt).Data() ));
          prunedEta->push_back(j.userFloat( (sPruned+sEta).Data() ));
          prunedPhi->push_back(j.userFloat( (sPruned+sPhi).Data() ));
          prunedJEC0->push_back(j.userFloat( (sPruned+sJEC0).Data() ));
          trimmedPt->push_back(j.userFloat( (sTrimmed+sPt).Data() ));
          trimmedEta->push_back(j.userFloat( (sTrimmed+sEta).Data() ));
          trimmedPhi->push_back(j.userFloat( (sTrimmed+sPhi).Data() ));
          trimmedJEC0->push_back(j.userFloat( (sTrimmed+sJEC0).Data() ));
          sdPt->push_back(j.userFloat( (sSD+sPt).Data() ));
          sdEta->push_back(j.userFloat( (sSD+sEta).Data() ));
          sdPhi->push_back(j.userFloat( (sSD+sPhi).Data() ));
          sdJEC0->push_back(j.userFloat( (sSD+sJEC0).Data() ));
        }

        auto Subjets = j.subjets(SubjetsName);
        for ( auto const & i : Subjets ) {
            new ( (*subjets)[nsubjet]) TLorentzVector(i->px(), i->py(), i->pz(), i->energy());
            subjetBtag->push_back(i->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            nsubjet++;
        }

    }
    return 0;
}

void NeroFatJets::runBTagging(const pat::Jet* jet) {
  const IPTagInfo * ipTagInfo = jet->tagInfoCandIP("pfImpactParameter");
  const SVTagInfo * svTagInfo = jet->tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");

  std::vector<fastjet::PseudoJet> currentAxes;
  float tau1IVF_tmp = tau1->back();
  float tau2IVF_tmp = tau2->back();
  recalcNsubjettiness(*jet,*svTagInfo,tau1IVF_tmp,tau2IVF_tmp,currentAxes);
  tau1IVF->push_back(tau1IVF_tmp);
  tau2IVF->push_back(tau2IVF_tmp);

  const Tracks & selectedTracks( ipTagInfo->selectedTracks() );
  reco::TrackKinematics allKinematics;
  unsigned int trackSize = selectedTracks.size();
  for (unsigned int itt=0; itt < trackSize; ++itt) {
    const reco::Track & ptrack = *(reco::btag::toTrack(selectedTracks[itt]));
    const TrackRef ptrackRef = selectedTracks[itt];
    int trackPV;
    float trackPVWeight;
    setTracksPV(ptrackRef,primaryVertices,trackPV, trackPVWeight);
    if (!trackPV && trackPVWeight > 0)
      allKinematics.add(ptrack,trackPVWeight);
  }
  math::XYZTLorentzVector allSum =  allKinematics.weightedVectorSum() ; //allKinematics.vectorSum()

  std::map<double, unsigned int> VTXmass;
  unsigned int nSV_tmp = svTagInfo->nVertices();
  nSV->push_back(nSV_tmp);
  unsigned int nMaxVertices = 128; // probably not more than this
  float svtxMass[nMaxVertices];
  float svtxEnergyRatio[nMaxVertices];
  float maxSVDeltaRToJet = R0-0.1+(R0-0.8)*0.1/0.7;
  edm::RefToBase<reco::Jet> rJet = ipTagInfo->jet();
  math::XYZVector jetDir = rJet->momentum().Unit();
  for (unsigned int vtx = 0; vtx < nSV_tmp; ++vtx) {

    const Vertex &vertex = svTagInfo->secondaryVertex(vtx);
    svtxMass[vtx] = vertex.p4().mass();

    Int_t totcharge=0;
    reco::TrackKinematics vertexKinematics;
    // get the vertex kinematics and charge
    vertexKinematicsAndCharge(vertex, vertexKinematics, totcharge);

    math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
    GlobalVector flightDir = svTagInfo->flightDirection(vtx);

    // JetInfo[iJetColl].SV_deltaR_jet[JetInfo[iJetColl].nSV]     = ( reco::deltaR(flightDir, jetDir) );
    // JetInfo[iJetColl].SV_deltaR_sum_jet[JetInfo[iJetColl].nSV] = ( reco::deltaR(vertexSum, jetDir) );
    // JetInfo[iJetColl].SV_deltaR_sum_dir[JetInfo[iJetColl].nSV] = ( reco::deltaR(vertexSum, flightDir) );

    // Line::PositionType pos(GlobalPoint(position(vertex).x(),position(vertex).y(),position(vertex).z()));
    // Line trackline(pos,flightDir);
    // // get the Jet  line
    // Line::PositionType pos2(GlobalPoint(pv->x(),pv->y(),pv->z()));
    // Line::DirectionType dir2(GlobalVector(jetDir.x(),jetDir.y(),jetDir.z()));
    // Line jetline(pos2,dir2);
    // // now compute the distance between the two lines
    // JetInfo[iJetColl].SV_vtxDistJetAxis[JetInfo[iJetColl].nSV] = (jetline.distance(trackline)).mag();

    svtxEnergyRatio[vtx] = vertexSum.E() / allSum.E();
    if (reco::deltaR2(flightDir, jetDir)<(maxSVDeltaRToJet*maxSVDeltaRToJet))
      VTXmass[svtxMass[vtx]]=vtx;
  } //// if secondary vertices present
  int cont=0;
  GlobalVector flightDir0, flightDir1;
  reco::Candidate::LorentzVector svP4_0 , svP4_1;
  vector<float>* svEnergyRatio_tmp = new vector<float>;
  vector<float>* svMass_tmp = new vector<float>;
  vector<float>* svPt_tmp = new vector<float>;
  for ( std::map<double, unsigned int>::reverse_iterator iVtx=VTXmass.rbegin(); iVtx!=VTXmass.rend(); ++iVtx)
  {
    ++cont;
    const Vertex &vertex = svTagInfo->secondaryVertex(iVtx->second);
    float SV_EnergyRatio = svtxEnergyRatio[iVtx->second];
    svEnergyRatio_tmp->push_back(SV_EnergyRatio);
    svMass_tmp->push_back(vertex.p4().mass);
    svPt_tmp->push_back(vertex.p4().pt());
    if (cont==1)
    {
      flightDir0 = svTagInfo->flightDirection(iVtx->second);
      svP4_0= vertex.p4();
      float tauDot_tmp;
      if (reco::deltaR2(flightDir0,currentAxes[1])<reco::deltaR2(flightDir0,currentAxes[0]))
        tauDot_tmp = (currentAxes[1].px()*flightDir0.x()+currentAxes[1].py()*flightDir0.y()+currentAxes[1].pz()*flightDir0.z())/(sqrt(currentAxes[1].modp2())*flightDir0.mag());
      else
        tauDot_tmp = (currentAxes[0].px()*flightDir0.x()+currentAxes[0].py()*flightDir0.y()+currentAxes[0].pz()*flightDir0.z())/(sqrt(currentAxes[0].modp2())*flightDir0.mag());
      tauDot->push_back(tauDot_tmp);
    }
    if (cont==2)
    {
      svEnergyRatio1->push_back(SV_EnergyRatio);
      flightDir1 = svTagInfo->flightDirection(iVtx->second);
      svP4_1 = vertex.p4();
      zRatio->push_back(reco::deltaR(flightDir0,flightDir1)*(svPt0->back())/(svP4_0+svP4_1).mass());
    }
    if (cont==4)
      break;
  }
  if (cont<2) {
    zRatio->push_back(-1.);
  }
  if (cont<1) {
    tauDot->push_back(-1);
  }
  svEnergyRatio->push_back(svEnergyRatio_tmp);
  svMass->push_back(svMass_tmp);
  svPt->push_back(svPt_tmp);

  // TRACK INFO
  reco::TaggingVariableList ipVars = ipTagInfo->taggingVariables();
  reco::TaggingVariableList svVars = svTagInfo->taggingVariables();
  nTracks->push_back(ipTagInfo->selectedTracks().size());
  storeTrkVars(ipVars, svVars);


}

void NeroFatJets::storeTrkVars(reco::TaggingVariableList &ipVars, reco::TaggingVariableList &svVars) {
  trackMomentum->push_back( &(ipVars.getList(reco::btau::trackMomentum,false) );
  trackEta->push_back( &(ipVars.getList(reco::btau::trackEta,false) );
  trackPhi->push_back( &(ipVars.getList(reco::btau::trackPhi,false) );
  trackPtRel->push_back( &(ipVars.getList(reco::btau::trackPtRel,false) );
  trackPPar->push_back( &(ipVars.getList(reco::btau::trackPPar,false) );
  trackEtaRel->push_back( &(ipVars.getList(reco::btau::trackEtaRel,false) );
  trackDeltaR->push_back( &(ipVars.getList(reco::btau::trackDeltaR,false) );
  trackPtRatio->push_back( &(ipVars.getList(reco::btau::trackPtRatio,false) );
  trackPParRatio->push_back( &(ipVars.getList(reco::btau::trackPParRatio,false) );
  trackSip2dVal->push_back( &(ipVars.getList(reco::btau::trackSip2dVal,false) );
  trackSip2dSig->push_back( &(ipVars.getList(reco::btau::trackSip2dSig,false) );
  trackSip3dVal->push_back( &(ipVars.getList(reco::btau::trackSip3dVal,false) );
  trackSip3dSig->push_back( &(ipVars.getList(reco::btau::trackSip3dSig,false) );
  trackDecayLenVal->push_back( &(ipVars.getList(reco::btau::trackDecayLenVal,false) );
  trackDecayLenSig->push_back( &(ipVars.getList(reco::btau::trackDecayLenSig,false) );
  trackJetDistVal->push_back( &(ipVars.getList(reco::btau::trackJetDistVal,false) );
  trackJetDistSig->push_back( &(ipVars.getList(reco::btau::trackJetDistSig,false) );
  trackChi2->push_back( &(ipVars.getList(reco::btau::trackChi2,false) );
  trackNTotalHits->push_back( &(ipVars.getList(reco::btau::trackNTotalHits,false) );
  trackNPixelHits->push_back( &(ipVars.getList(reco::btau::trackNPixelHits,false) );


}

void NeroFatJets::setTracksPVBase(const reco::TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight)
{
  iPV = -1;
  PVweight = 0.;

  const reco::TrackBaseRef trackBaseRef( trackRef );

  typedef reco::VertexCollection::const_iterator IV;
  typedef reco::Vertex::trackRef_iterator IT;

  for(IV iv=pvHandle->begin(); iv!=pvHandle->end(); ++iv)
  {
    const reco::Vertex & vtx = *iv;
    // loop over tracks in vertices
    for(IT it=vtx.tracks_begin(); it!=vtx.tracks_end(); ++it)
    {
      const reco::TrackBaseRef & baseRef = *it;
      // one of the tracks in the vertex is the same as the track considered in the function
      if( baseRef == trackBaseRef )
      {
        float w = vtx.trackWeight(baseRef);
        // select the vertex for which the track has the highest weight
        if( w > PVweight )
        {
          PVweight = w;
          iPV = ( iv - pvHandle->begin() );
          break;
        }
      }
    }
  }
}

void NeroFatJets::setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight)
{
  iPV = -1;
  PVweight = 0.;

  const pat::PackedCandidate * pcand = dynamic_cast<const pat::PackedCandidate *>(trackRef.get());

  if(pcand) // MiniAOD case
  {
    if( pcand->fromPV() == pat::PackedCandidate::PVUsedInFit )
    {
      iPV = 0;
      PVweight = 1.;
    }
  }
  else
  {
    const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(trackRef.get());

    setTracksPVBase(pfcand->trackRef(), pvHandle, iPV, PVweight);
  }
}

void NeroFatJets::vertexKinematicsAndCharge(const Vertex & vertex, reco::TrackKinematics & vertexKinematics, Int_t & charge)
{
  const std::vector<reco::CandidatePtr> & tracks = vertex.daughterPtrVector();

  for(std::vector<reco::CandidatePtr>::const_iterator track = tracks.begin(); track != tracks.end(); ++track) {
    const reco::Track& mytrack = *(*track)->bestTrack();
    vertexKinematics.add(mytrack, 1.0);
    charge+=mytrack.charge();
  }
}

void NeroFatJets::recalcNsubjettiness(const pat::Jet &jet,
                                      const SVTagInfo &svTagInfo,
                                      float & tau1,
                                      float & tau2,
                                      std::vector<fastjet::PseudoJet> &currentAxes) {
  // recalculate nsubjettiness removing the bjet
  std::vector<fastjet::PseudoJet> fjParticles;
  std::vector<reco::CandidatePtr> svDaughters;

  // loop over IVF vertices and push them in the vector of FastJet constituents and also collect their daughters
  for(size_t i=0; i<svTagInfo.nVertices(); ++i)
  {
    const reco::VertexCompositePtrCandidate & vtx = svTagInfo.secondaryVertex(i);

    fjParticles.push_back( fastjet::PseudoJet( vtx.px(), vtx.py(), vtx.pz(), vtx.energy() ) );

    const std::vector<reco::CandidatePtr> & daughters = vtx.daughterPtrVector();
    svDaughters.insert(svDaughters.end(), daughters.begin(), daughters.end());
  }

  // loop over jet constituents and select those that are not daughters of IVF vertices
  std::vector<reco::CandidatePtr> constituentsOther;
  for(const reco::CandidatePtr & daughter : jet.daughterPtrVector())
  {
    if (std::find(svDaughters.begin(), svDaughters.end(), daughter) == svDaughters.end())
      constituentsOther.push_back( daughter );
  }

  // loop over jet constituents that are not daughters of IVF vertices and push them in the vector of FastJet constituents
  for(const reco::CandidatePtr & constit : constituentsOther)
  {
    if ( constit.isNonnull() && constit.isAvailable() )
      fjParticles.push_back( fastjet::PseudoJet( constit->px(), constit->py(), constit->pz(), constit->energy() ) );
    else
      edm::LogWarning("MissingJetConstituent") << "Jet constituent required for N-subjettiness computation is missing!";
  }

  // re-calculate N-subjettiness
  tau1 = njettiness.getTau(1, fjParticles);
  tau2 = njettiness.getTau(2, fjParticles);
  currentAxes = njettiness.currentAxes();
}

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
