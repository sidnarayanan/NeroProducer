import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_4_1/RelValProdTTbar_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/0A9E2CED-C9EC-E411-A8E4-003048FFCBA8.root'
        )
)

process.skimmedMETs  = cms.EDFilter("PATMETSelector",
                                    src = cms.InputTag("slimmedMETs"),
                                    cut = cms.string("energy > 200"),
                                    filter = cms.bool(True)
                                    )
process.filter = cms.Path(process.skimmedMETs)

##Output File
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.outPut = cms.OutputModule("PoolOutputModule",
                                  SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('filter')),
                                  outputCommands = cms.untracked.vstring(
        "keep *",
        "drop *_*_*_SKIM",
        "keep patJets_selectedPatJetsAK1p5PFCHS_*_*",
        "keep *_skimmedMETs_*_SKIM"),
                                  fileName = cms.untracked.string("skim.root")
                                  )

process.endpath = cms.EndPath(process.outPut)
#process.outPutPath = cms.EndPath(process.outPut)

#################################################
## Remake jets
#################################################

## Filter out neutrinos from packed GenParticles
process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))

## Define GenJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')
process.ak1p5GenJetsNoNu = ak4GenJets.clone(rParam = cms.double(1.5), src = 'packedGenParticlesForJetsNoNu')

## Select charged hadron subtracted packed PF candidates
process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

## Define PFJetsCHS
process.ak4PFJetsCHS = ak4PFJets.clone(src = 'pfCHS', doAreaFastjet = True)
process.ak1p5PFJetsCHS = ak4PFJets.clone(rParam = cms.double(1.5), src = 'pfCHS', doAreaFastjet = True)

## SoftDrop fat GenJets (two jet collections are produced, fat jets and subjets)                                                                                                   
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.ak1p5GenJetsNoNuSoftDrop = ak4GenJets.clone(
    rParam = cms.double(1.5),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0   = cms.double(1.5),
    useExplicitGhosts = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)


#################################################
## Remake PAT jets
#################################################

from PhysicsTools.PatAlgos.tools.jetTools import *

## b-tag discriminators
bTagDiscriminators = [
    'pfCombinedInclusiveSecondaryVertexV2BJetTags'
]

addJetCollection(
    process,
    labelName = 'AK1p5PFCHS',
    jetSource = cms.InputTag('ak1p5PFJetsCHS'),
    pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    svSource = cms.InputTag('slimmedSecondaryVertices'),
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None'),
    genJetCollection = cms.InputTag('ak1p5GenJetsNoNu'),
    genParticles = cms.InputTag('prunedGenParticles'),
    algo = 'AK',
    rParam = 1.5
)

getattr(process,'selectedPatJetsAK1p5PFCHS').cut = cms.string('pt > 200')

process.p = cms.Path(process.selectedPatJetsAK1p5PFCHS)

from PhysicsTools.PatAlgos.tools.pfTools import *
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

process.options = cms.untracked.PSet( 
        wantSummary = cms.untracked.bool(True),
        allowUnscheduled = cms.untracked.bool(True)
)
