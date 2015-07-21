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
                             njettiness(fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1.0,R0))
{
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

        trimmedMass ->push_back(j.userFloat("PFJetsCHSTrimmedMass"));
        prunedMass  ->push_back(j.userFloat("PFJetsCHSPrunedMass"));
        filteredMass->push_back(j.userFloat("PFJetsCHSFilteredMass"));
        softdropMass->push_back(j.userFloat("PFJetsCHSSoftDropMass"));
        hasSubjet->push_back(j.hasSubjets(SubjetsName));

        if (doSubstructure) {
          // B Tagging
          runBTagging(&j);

          // get groomed jet infos
          TString sR0 = TString::Format("%i",int(R0*10));
          TString sPt(":Pt");  TString sEta(":Eta");  TString sPhi(":Phi"); TString sJEC0(":jecFactor0"); TString sMass(":Mass");
          TString sPruned = TString("Pruned") + sR0;
          TString sTrimmed = TString("Trimmed") + sR0;
          TString sSD = TString("SoftDrop") + sR0;
          TString sNjettiness = TString("Njettiness") + sR0;
          TLorentzVector p_tmp;
          p_tmp.SetPtEtaPhiM(j.userFloat((sPruned+sPt).Data()),
                              j.userFloat((sPruned+sEta).Data()),
                              j.userFloat((sPruned+sPhi).Data()),
                              j.userFloat((sPruned+sMass).Data()));
          prunedMomentum->push_back(p_tmp);
          p_tmp.SetPtEtaPhiM(j.userFloat((sTrimmed+sPt).Data()),
                              j.userFloat((sTrimmed+sEta).Data()),
                              j.userFloat((sTrimmed+sPhi).Data()),
                              j.userFloat((sTrimmed+sMass).Data()));
          trimmedMomentum->push_back(p_tmp);
          p_tmp.SetPtEtaPhiM(j.userFloat((sSD+sPt).Data()),
                              j.userFloat((sSD+sEta).Data()),
                              j.userFloat((sSD+sPhi).Data()),
                              j.userFloat((sSD+sMass).Data()));
          sdMomentum->push_back(p_tmp);
          prunedJEC0->push_back(j.userFloat( (sPruned+sJEC0).Data() ));
          trimmedJEC0->push_back(j.userFloat( (sTrimmed+sJEC0).Data() ));
          sdJEC0->push_back(j.userFloat( (sSD+sJEC0).Data() ));
          tau1 -> push_back(j.userFloat( (sNjettiness+":tau1").Data() ));
          tau2 -> push_back(j.userFloat( (sNjettiness+":tau2").Data() ));
          tau3 -> push_back(j.userFloat( (sNjettiness+":tau3").Data() ));
          // add tau4 for ttbar rejection?
        } else {
          tau1 -> push_back(j.userFloat("Njettiness:tau1"));
          tau2 -> push_back(j.userFloat("Njettiness:tau2"));
          tau3 -> push_back(j.userFloat("Njettiness:tau3"));
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

  // TRACK INFO
  reco::TrackKinematics allKinematics;
  computeTrkVars(ipTagInfo, svTagInfo, allKinematics);
  math::XYZTLorentzVector allSum =  allKinematics.weightedVectorSum() ; //allKinematics.vectorSum()

  std::map<double, unsigned int> VTXmass;
  unsigned int nSV_tmp = svTagInfo->nVertices();
  nSV->push_back(nSV_tmp);
  float maxSVDeltaRToJet = R0-0.1+(R0-0.8)*0.1/0.7;
  edm::RefToBase<reco::Jet> rJet = ipTagInfo->jet();
  math::XYZVector jetDir = rJet->momentum().Unit();
  for (unsigned int vtx = 0; vtx < nSV_tmp; ++vtx) {
    const Vertex &vertex = svTagInfo->secondaryVertex(vtx);
    double svtxMass = vertex.p4().mass();
    GlobalVector flightDir = svTagInfo->flightDirection(vtx);
    if (reco::deltaR2(flightDir, jetDir)<(maxSVDeltaRToJet*maxSVDeltaRToJet))
      VTXmass[svtxMass]=vtx;
  } //// if secondary vertices present
  computeSVVars(VTXmass, svTagInfo, allSum.E());
}

template<typename T>
void NeroFatJets::fillTaggingVar(reco::TaggingVariableList const & varList, std::vector<vector<T>*>* target, reco::btau::TaggingVariableName kVar) {
  std::vector<T> * tmp_vec = new vector<T>;
  std::vector<T> vals = varList.getList(kVar,false);
  std::copy(vals.begin(), vals.end(), back_inserter(*tmp_vec));
  target->push_back(tmp_vec);
}

void NeroFatJets::computeTrkVars(const IPTagInfo * ipTagInfo, const SVTagInfo * svTagInfo, reco::TrackKinematics & allKinematics) {
  const reco::TaggingVariableList ipVars = ipTagInfo->taggingVariables();
  nTracks->push_back(ipTagInfo->selectedTracks().size());
  fillTaggingVar(ipVars,track_momentum,reco::btau::trackMomentum);
  fillTaggingVar(ipVars,track_eta,reco::btau::trackEta);
  fillTaggingVar(ipVars,track_phi,reco::btau::trackPhi);
  fillTaggingVar(ipVars,track_ptRel,reco::btau::trackPtRel);
  fillTaggingVar(ipVars,track_pPar,reco::btau::trackPPar);
  fillTaggingVar(ipVars,track_etaRel,reco::btau::trackEtaRel);
  fillTaggingVar(ipVars,track_deltaR,reco::btau::trackDeltaR);
  fillTaggingVar(ipVars,track_ptRatio,reco::btau::trackPtRatio);
  fillTaggingVar(ipVars,track_pParRatio,reco::btau::trackPParRatio);
  fillTaggingVar(ipVars,track_sip2dVal,reco::btau::trackSip2dVal);
  fillTaggingVar(ipVars,track_sip2dSig,reco::btau::trackSip2dSig);
  fillTaggingVar(ipVars,track_sip3dVal,reco::btau::trackSip3dVal);
  fillTaggingVar(ipVars,track_sip3dSig,reco::btau::trackSip3dSig);
  fillTaggingVar(ipVars,track_decayLenVal,reco::btau::trackDecayLenVal);
  fillTaggingVar(ipVars,track_decayLenSig,reco::btau::trackDecayLenSig);
  fillTaggingVar(ipVars,track_jetDistVal,reco::btau::trackJetDistVal);
  fillTaggingVar(ipVars,track_jetDistSig,reco::btau::trackJetDistSig);
  fillTaggingVar(ipVars,track_chi2,reco::btau::trackChi2);
  fillTaggingVar(ipVars,track_nTotalHits,reco::btau::trackNTotalHits);
  fillTaggingVar(ipVars,track_nPixelHits,reco::btau::trackNPixelHits);

  std::vector<float> tmp_dxy;
  std::vector<float> tmp_dz;
  std::vector<float> tmp_IP2D ;
  std::vector<float> tmp_IP2Dsig ;
  std::vector<float> tmp_IP ;
  std::vector<float> tmp_IPsig ;
  std::vector<float> tmp_IP2Derr ;
  std::vector<float> tmp_IPerr ;
  std::vector<float> tmp_prob ;
  std::vector<float> tmp_PVWeight ;
  std::vector<float> tmp_SVWeight ;
  std::vector<int> tmp_nHitStrip ;
  std::vector<int> tmp_nHitTIB ;
  std::vector<int> tmp_nHitTID ;
  std::vector<int> tmp_nHitTOB ;
  std::vector<int> tmp_nHitTEC ;
  std::vector<int> tmp_nHitPXB ;
  std::vector<int> tmp_nHitPXF ;
  std::vector<int> tmp_isHitL1 ;
  std::vector<int> tmp_PV ;
  std::vector<int> tmp_fromSV ;
  std::vector<int> tmp_SV ;
  const Tracks & selectedTracks( ipTagInfo->selectedTracks() );
  unsigned int trackSize = selectedTracks.size();
  const reco::Vertex *pv = &(*primaryVertices->begin());
  for (unsigned int itt=0; itt < trackSize; ++itt) {
    const reco::Track & ptrack = *(reco::btag::toTrack(selectedTracks[itt]));
    const TrackRef ptrackRef = selectedTracks[itt];

    tmp_dxy.push_back(ptrack.dxy(pv->position()));
    tmp_dz.push_back(ptrack.dz(pv->position()));
    tmp_IP2D.push_back(ipTagInfo->impactParameterData()[itt].ip2d.value());
    tmp_IP2Dsig.push_back(ipTagInfo->impactParameterData()[itt].ip2d.significance());
    tmp_IP.push_back(ipTagInfo->impactParameterData()[itt].ip3d.value());
    tmp_IPsig.push_back(ipTagInfo->impactParameterData()[itt].ip3d.significance());
    tmp_IP2Derr.push_back(ipTagInfo->impactParameterData()[itt].ip2d.error());
    tmp_IPerr.push_back(ipTagInfo->impactParameterData()[itt].ip3d.error());
    tmp_prob.push_back(ipTagInfo->probabilities(0)[itt]);

    tmp_nHitStrip.push_back(ptrack.hitPattern().numberOfValidStripHits());
    tmp_nHitTIB.push_back(ptrack.hitPattern().numberOfValidStripTIBHits());
    tmp_nHitTID.push_back(ptrack.hitPattern().numberOfValidStripTIDHits());
    tmp_nHitTOB.push_back(ptrack.hitPattern().numberOfValidStripTOBHits());
    tmp_nHitTEC.push_back(ptrack.hitPattern().numberOfValidStripTECHits());
    tmp_nHitPXB.push_back(ptrack.hitPattern().numberOfValidPixelBarrelHits());
    tmp_nHitPXF.push_back(ptrack.hitPattern().numberOfValidPixelEndcapHits());
    tmp_isHitL1.push_back(ptrack.hitPattern().hasValidHitInFirstPixelBarrel());
    int PV, fromSV, SV;
    float PVWeight, SVWeight;
    setTracksPV(ptrackRef, primaryVertices, PV, PVWeight);
    tmp_PV.push_back(PV);
    tmp_PVWeight.push_back(PVWeight);
    if (!PV && PVWeight > 0)
      allKinematics.add(ptrack,PVWeight);
    if (svTagInfo) {
      setTracksSV(ptrackRef,svTagInfo,fromSV, SV, SVWeight);
    } else {
      fromSV = 0;
      SV = -1;
      SVWeight = 0;
    }
    tmp_fromSV.push_back(fromSV);
    tmp_SV.push_back(SV);
    tmp_SVWeight.push_back(SVWeight);
  }
  track_dxy->push_back(&tmp_dxy);
  track_dz->push_back(&tmp_dz);
  track_IP2D->push_back(&tmp_IP2D);
  track_IP2Dsig->push_back(&tmp_IP2Dsig);
  track_IP->push_back(&tmp_IP);
  track_IPsig->push_back(&tmp_IPsig);
  track_IP2Derr->push_back(&tmp_IP2Derr);
  track_IPerr->push_back(&tmp_IPerr);
  track_prob->push_back(&tmp_prob);
  track_PVWeight->push_back(&tmp_PVWeight);
  track_SVWeight->push_back(&tmp_SVWeight);
  track_nHitStrip->push_back(&tmp_nHitStrip);
  track_nHitTIB->push_back(&tmp_nHitTIB);
  track_nHitTID->push_back(&tmp_nHitTID);
  track_nHitTOB->push_back(&tmp_nHitTOB);
  track_nHitTEC->push_back(&tmp_nHitTEC);
  track_nHitPXB->push_back(&tmp_nHitPXB);
  track_nHitPXF->push_back(&tmp_nHitPXF);
  track_isHitL1->push_back(&tmp_isHitL1);
  track_PV->push_back(&tmp_PV);
  track_fromSV->push_back(&tmp_fromSV);
  track_SV->push_back(&tmp_SV);
}

void NeroFatJets::computeSVVars(const std::map<double, unsigned int> VTXmass, const SVTagInfo * svTagInfo, float allEnergy) {
  std::vector<float> tmp_mass;
  std::vector<float> tmp_pt;
  std::vector<float> tmp_eta;
  std::vector<float> tmp_phi;
  std::vector<float> tmp_energyRatio;
  std::vector<int> tmp_charge;
  std::vector<float> tmp_dirX;
  std::vector<float> tmp_dirY;
  std::vector<float> tmp_dirZ;
  std::vector<float> tmp_deltaRJet;
  std::vector<float> tmp_deltaRSumJet;
  std::vector<float> tmp_deltaRSumDir;
  std::vector<float> tmp_vtxDistJetAxis;
  std::vector<float> tmp_flight;
  std::vector<float> tmp_flightErr;
  std::vector<float> tmp_flight2D;
  std::vector<float> tmp_flight2DErr;
  std::vector<float> tmp_nTrk;
  GlobalVector flightDir0, flightDir1;
  reco::Candidate::LorentzVector svP4_0 , svP4_1;
  const reco::Vertex *pv = &(*primaryVertices->begin());
  int cont=0;
  // fill FatJet::fSVData in order of decreasing secondary vertex mass
  for ( std::map<double, unsigned int>::reverse_iterator iVtx=VTXmass.rbegin(); iVtx!=VTXmass.rend(); ++iVtx) {
    ++cont;
    unsigned int idx = iVtx->second;
    const Vertex &vertex = svTagInfo->secondaryVertex(idx);
    int totcharge=0;
    reco::TrackKinematics vertexKinematics;
    vertexKinematicsAndCharge(vertex, vertexKinematics, totcharge);
    math::XYZTLorentzVector vertexSum = vertexKinematics.weightedVectorSum();
    tmp_energyRatio.push_back(vertexSum.E() / allEnergy);
    tmp_charge.push_back(totcharge);

    //svx kinematics
    tmp_mass.push_back(iVtx->first);
    tmp_pt.push_back(vertex.p4().pt());
    tmp_eta.push_back(vertex.p4().eta());
    tmp_phi.push_back(vertex.p4().phi());

    GlobalVector flightDir = svTagInfo->flightDirection(idx);
    tmp_dirX.push_back(flightDir.x());
    tmp_dirY.push_back(flightDir.y());
    tmp_dirZ.push_back(flightDir.z());
    tmp_deltaRJet.push_back(reco::deltaR(flightDir, jetDir));
    tmp_deltaRSumJet.push_back(reco::deltaR(vertexSum, jetDir));
    tmp_deltaRSumDir.push_back(reco::deltaR(vertexSum, flightDir));

    Line::PositionType pos(GlobalPoint(vertex.position().x(),vertex.position().y(),vertex.position().z()));
    Line trackline(pos,flightDir);
    Line::PositionType pos2(GlobalPoint(pv->x(),pv->y(),pv->z()));
    Line::DirectionType dir2(GlobalVector(jetDir.x(),jetDir.y(),jetDir.z()));
    Line jetline(pos2,dir2);
    tmp_vtxDistJetAxis.push_back((jetline.distance(trackline)).mag());

    tmp_flight.push_back(svTagInfo->flightDistance(idx).value());
    tmp_flightErr.push_back(svTagInfo->flightDistance(idx).error());
    tmp_flight2D.push_back(svTagInfo->flightDistance(idx, true).value());
    tmp_flight2DErr.push_back(svTagInfo->flightDistance(idx, true).error());
    tmp_nTrk.push_back(svTagInfo->secondaryVertex(idx).numberOfSourceCandidatePtrs());

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
      flightDir1 = svTagInfo->flightDirection(iVtx->second);
      svP4_1 = vertex.p4();
      zRatio->push_back(reco::deltaR(flightDir0,flightDir1)*(svP4_0.pt())/(svP4_0+svP4_1).mass());
    }
  }
  sv_mass->push_back(&tmp_mass);
  sv_pt->push_back(&tmp_pt);
  sv_eta->push_back(&tmp_eta);
  sv_phi->push_back(&tmp_phi);
  sv_energyRatio->push_back(&tmp_energyRatio);
  sv_charge->push_back(&tmp_charge);
  sv_dirX->push_back(&tmp_dirX);
  sv_dirY->push_back(&tmp_dirY);
  sv_dirZ->push_back(&tmp_dirZ);
  sv_deltaRJet->push_back(&tmp_deltaRJet);
  sv_deltaRSumJet->push_back(&tmp_deltaRSumJet);
  sv_deltaRSumDir->push_back(&tmp_deltaRSumDir);
  sv_vtxDistJetAxis->push_back(&tmp_vtxDistJetAxis);
  sv_flight->push_back(&tmp_flight);
  sv_flightErr->push_back(&tmp_flightErr);
  sv_flight2D->push_back(&tmp_flight2D);
  sv_flight2DErr->push_back(&tmp_flight2DErr);
  sv_nTrk->push_back(&tmp_nTrk);

  if (cont<2) zRatio->push_back(-1.);
  if (cont<1) tauDot->push_back(-1);

  reco::TaggingVariableList svVars = svTagInfo->taggingVariables();
  

}

void NeroFatJets::setTracksSV (const TrackRef & trackRef, const SVTagInfo * svTagInfo, int & isFromSV, int & iSV, float & SVweight)
{
  isFromSV = 0;
  iSV = -1;
  SVweight = 0.;
  typedef std::vector<reco::CandidatePtr>::const_iterator IT;
  size_t nSV = svTagInfo->nVertices();
  for(size_t iv=0; iv<nSV; ++iv)  {
    const Vertex & vtx = svTagInfo->secondaryVertex(iv);
    // one of the tracks in the vertex is the same as the track considered in the function
    const std::vector<reco::CandidatePtr> & tracks = vtx.daughterPtrVector();
    if( std::find(tracks.begin(),tracks.end(),trackRef) != tracks.end() )    {
      SVweight = 1.;
      isFromSV = 1;
      iSV = iv;
    // select the first vertex for which the track is used in the fit
    }
    // (reco::VertexCompositePtrCandidate does not store track weights so can't select the vertex for which the track has the highest weight)
    if(iSV>=0)
      break;
  }
}

void NeroFatJets::setTracksPV(const TrackRef & trackRef, const edm::Handle<reco::VertexCollection> & pvHandle, int & iPV, float & PVweight)
{
  iPV = -1;
  PVweight = 0.;
  const reco::PFCandidate * pfcand = dynamic_cast<const reco::PFCandidate *>(trackRef.get());
  setTracksPVBase(pfcand->trackRef(), pvHandle, iPV, PVweight);
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
