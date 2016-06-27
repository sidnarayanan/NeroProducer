import FWCore.ParameterSet.Config as cms
from NeroProducer.Nero.Nero_cfi import nero

print " ------- LOADING MONO TOP CONFIGURATION -------- "

nero.info = cms.string("NeroMonotop")

nero.triggerNames = cms.vstring([
                                 'HLT_PFMET170_NoiseCleaned_v',                      # MET
                                 'HLT_PFMET170_JetIdCleaned_v',
                                 'HLT_PFMET170_HBHECleaned_v',
                                 'HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v',
                                 'HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v',
                                 'HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v',
                                 'HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v',
                                 'HLT_IsoMu18_v',                                    # MUON
                                 'HLT_IsoMu20_v',
                                 'HLT_IsoMu22_v',
                                 'HLT_IsoMu24',
                                 'HLT_IsoMu27_v',
                                 'HLT_IsoTkMu18_v',
                                 'HLT_IsoTkMu24_v',
                                 'HLT_Ele25_eta2p1_WPTight_Gsf_v',                    # ELECTRON
                                 'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
                                 'HLT_Ele27_WPTight_Gsf_v',
                                 'HLT_Ele35_WPLoose_Gsf_v',
                                 'HLT_Ele105_CaloIdVT_GsfTrkIdT_v', 
                                 'HLT_ECALHT800_v',                                    # ELECTRON+PHOTON
                                 'HLT_Photon175_v',                                    # PHOTON
                                 'HLT_Photon165_HE10_v',
                                 'HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v', 
                                 'HLT_Photon135_PFMET100_v',
                                 'HLT_Photon300_NoHE_v',
                               ])
nero.minJetPt         = cms.double (15.)
nero.minJetEta        = cms.double (4.7)
nero.minJetId         = cms.string ('noid')
nero.minPuppiJetPt    = cms.double (15.)
nero.minPuppiJetEta   = cms.double (4.7)
nero.minPuppiJetId    = cms.string ('noid')
nero.minCA15PuppiPt   = cms.double (150.)
nero.minCA15PuppiEta  = cms.double (2.5)
nero.minCA15PuppiId   = cms.string ('noid')
nero.puppiCA15        = cms.InputTag('packedPatJetsPFCA15Puppi')
nero.minTauId         = cms.string ('decayModeFinding')
nero.maxTauIso        = cms.double (5)
nero.extendTau        = cms.bool(True)
nero.extendMet        = cms.bool(True)
nero.doReclustering   = cms.bool(True)
nero.doCA15           = cms.bool(True)
woof=True
if woof:
  nero.doPuppi          = cms.bool(True)
  nero.metsPuppi        = cms.InputTag("type1PuppiMET")
  nero.metsPuppiUncorrected = cms.InputTag("pfMETPuppi")
  nero.puppijets        = cms.InputTag('patJetsPFAK4Puppi')
