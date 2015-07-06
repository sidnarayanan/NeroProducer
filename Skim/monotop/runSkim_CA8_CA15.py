import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

process = cms.Process("SKIM")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(10000)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# process.source = cms.Source("PoolSource",
#     fileNames = cms.untracked.vstring(
#         'file:/afs/cern.ch/work/z/zdemirag/work/MonoJet/CMSSW_7_4_1/src/MiniAodAnalyzer/MiniAodAnalyzer/80CF5456-B9EC-E411-93DA-002618FDA248.root'
#         )
# )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_7_4_1/RelValTTbar_13/MINIAODSIM/PU50ns_MCRUN2_74_V8_gensim_740pre7-v1/00000/7EC72BA9-44EC-E411-9DF3-0025905A60B0.root'
    )
)

process.skimmedMETs  = cms.EDFilter("PATMETSelector",
                                    src = cms.InputTag("slimmedMETs"),
                                    cut = cms.string("energy > 100"),
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
        "keep patJets_selectedPatJetsPFCHS8_*_*",
        "keep patJets_packedPatJetsPFCHS8_*_*",
        "keep patJets_selectedPatJetsPFCHS15_*_*",
        "keep patJets_packedPatJetsPFCHS15_*_*",
        "keep patJets_selectedPatJetsSoftDropPFCHS8_*_*",
        "keep patJets_selectedPatJetsSoftDropSubjetsPFCHS8_*_*",
        "keep patJets_selectedPatJetsSoftDropPFCHSPacked8_*_*",
        "keep patJets_selectedPatJetsPrunedPFCHS8_*_*",
        "keep patJets_selectedPatJetsPrunedSubjetsPFCHS8_*_*",
        "keep patJets_selectedPatJetsPrunedPFCHSPacked8_*_*",
        "keep patJets_selectedPatJetsSoftDropPFCHS15_*_*",
        "keep patJets_selectedPatJetsSoftDropSubjetsPFCHS158_*_*",
        "keep patJets_selectedPatJetsSoftDropPFCHSPacked15_*_*",
        "keep patJets_selectedPatJetsPrunedPFCHS15_*_*",
        "keep patJets_selectedPatJetsPrunedSubjetsPFCHS15_*_*",
        "keep patJets_selectedPatJetsPrunedPFCHSPacked15_*_*",
        "keep *_skimmedMETs_*_SKIM"
        ),
                                  fileName = cms.untracked.string("skim.root")
                                  )

process.endpath = cms.EndPath(process.outPut)
#process.outPutPath = cms.EndPath(process.outPut)

#################################################
## Initialize jet stuff
#################################################
## Jet energy corrections
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')

algoLabel = 'CA'
jetAlgo = 'CambridgeAachen'
jetRadii = [0.8,1.5]

## Postfix
postfix = "PFlow"
genParticles = 'prunedGenParticles'
jetSource = 'ak4PFJets'
genJetCollection = 'ak4GenJetsNoNu'
pfCandidates = 'packedPFCandidates'
pvSource = 'offlineSlimmedPrimaryVertices'
svSource = 'slimmedSecondaryVertices'
muSource = 'slimmedMuons'
elSource = 'slimmedElectrons'

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")

# from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
## Select isolated collections
process.selectedMuons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedMuons"), cut = cms.string('''abs(eta)<2.5 && pt>10. &&
   (pfIsolationR04().sumChargedHadronPt+
max(0.,pfIsolationR04().sumNeutralHadronEt+
pfIsolationR04().sumPhotonEt-
0.50*pfIsolationR04().sumPUPt))/pt < 0.20 &&
(isPFMuon && (isGlobalMuon || isTrackerMuon) )'''))
process.selectedElectrons = cms.EDFilter("CandPtrSelector", src = cms.InputTag("slimmedElectrons"), cut = cms.string('''abs(eta)<2.5 && pt>20. &&
gsfTrack.isAvailable() &&
gsfTrack.hitPattern().numberOfLostHits(\'MISSING_INNER_HITS\') < 2 &&
(pfIsolationVariables().sumChargedHadronPt+
max(0.,pfIsolationVariables().sumNeutralHadronEt+
pfIsolationVariables().sumPhotonEt-
0.5*pfIsolationVariables().sumPUPt))/pt < 0.15'''))

## Do projections
process.pfCHS = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedPFCandidates"), cut = cms.string("fromPV"))
process.pfNoMuonCHS =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfCHS"), veto = cms.InputTag("selectedMuons"))
process.pfNoElectronsCHS = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuonCHS"), veto = cms.InputTag("selectedElectrons"))

process.pfNoMuon =  cms.EDProducer("CandPtrProjector", src = cms.InputTag("packedPFCandidates"), veto = cms.InputTag("selectedMuons"))
process.pfNoElectrons = cms.EDProducer("CandPtrProjector", src = cms.InputTag("pfNoMuon"), veto = cms.InputTag("selectedElectrons"))

process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
process.ak4GenJetsNoNu = ak4GenJets.clone(src = 'packedGenParticlesForJetsNoNu')

process.ak4PFJets = ak4PFJets.clone(src = 'pfNoElectronsCHS', doAreaFastjet = True)

## Load standard PAT objects (here we only need PAT muons but the framework will figure out what it needs to run using the unscheduled mode)
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")


bTagInfos = [
    'pfImpactParameterTagInfos'
   ,'pfSecondaryVertexTagInfos'
   ,'pfInclusiveSecondaryVertexFinderTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
]
## b-tag discriminators
bTagDiscriminators = [
     'pfJetBProbabilityBJetTags'
    ,'pfJetProbabilityBJetTags'
    ,'pfPositiveOnlyJetBProbabilityBJetTags'
    ,'pfPositiveOnlyJetProbabilityBJetTags'
    ,'pfNegativeOnlyJetBProbabilityBJetTags'
    ,'pfNegativeOnlyJetProbabilityBJetTags'
    ,'pfTrackCountingHighPurBJetTags'
    ,'pfTrackCountingHighEffBJetTags'
    ,'pfNegativeTrackCountingHighPurBJetTags'
    ,'pfNegativeTrackCountingHighEffBJetTags'
    ,'pfSimpleSecondaryVertexHighEffBJetTags'
    ,'pfSimpleSecondaryVertexHighPurBJetTags'
    ,'pfNegativeSimpleSecondaryVertexHighEffBJetTags'
    ,'pfNegativeSimpleSecondaryVertexHighPurBJetTags'
    ,'pfCombinedSecondaryVertexV2BJetTags'
    ,'pfPositiveCombinedSecondaryVertexV2BJetTags'
    ,'pfNegativeCombinedSecondaryVertexV2BJetTags'
    ,'pfCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'softPFMuonBJetTags'
    ,'positiveSoftPFMuonBJetTags'
    ,'negativeSoftPFMuonBJetTags'
    ,'softPFElectronBJetTags'
    ,'positiveSoftPFElectronBJetTags'
    ,'negativeSoftPFElectronBJetTags'
]


from PhysicsTools.PatAlgos.tools.jetTools import *
## Switch the default jet collection (done in order to use the above specified b-tag infos and discriminators)
switchJetCollection(
    process,
    jetSource = cms.InputTag(jetSource),
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag(genJetCollection),
    genParticles = cms.InputTag(genParticles),
    postfix = postfix
)



#################################################
## Remake jets
#################################################

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.genJetsNoNu8 = ak4GenJets.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(.8),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu"))
)
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.PFJetsCHS8 = ak4PFJets.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(.8),
    src = (getattr(process,"ak4PFJets").src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(250)
)
## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPruned8 = ak4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(.8),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu")),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.PFJetsCHSPruned8 = ak4PFJetsPruned.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(.8),
    src = getattr(process,"ak4PFJets").src,
    srcPVs = getattr(process,"ak4PFJets").srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(250)
)
## SoftDrop fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
process.genJetsNoNuSoftDrop8 = ak4GenJets.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(.8),
    R0 = cms.double(.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsSoftDrop_cfi import ak4PFJetsSoftDrop
process.PFJetsCHSSoftDrop8 = ak4PFJetsSoftDrop.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(.8),
    R0 = cms.double(.8),
    src = getattr(process,"ak4PFJets").src,
    srcPVs = getattr(process,"ak4PFJets").srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(250)
)

##

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.genJetsNoNu15 = ak4GenJets.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(1.5),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu"))
)
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.PFJetsCHS15 = ak4PFJets.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(1.5),
    src = (getattr(process,"ak4PFJets").src),
    srcPVs = (getattr(process,"ak4PFJets").srcPVs),
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(250)
)
## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPruned15 = ak4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(1.5),
    src = (cms.InputTag("packedGenParticlesForJetsNoNu")),
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.PFJetsCHSPruned15 = ak4PFJetsPruned.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(1.5),
    src = getattr(process,"ak4PFJets").src,
    srcPVs = getattr(process,"ak4PFJets").srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(250)
)
## SoftDrop fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
process.genJetsNoNuSoftDrop15 = ak4GenJets.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(1.5),
    R0 = cms.double(1.5),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"),
    useSoftDrop = cms.bool(True),
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsSoftDrop_cfi import ak4PFJetsSoftDrop
process.PFJetsCHSSoftDrop15 = ak4PFJetsSoftDrop.clone(
    jetAlgorithm = cms.string(jetAlgo),
    rParam = cms.double(1.5),
    R0 = cms.double(1.5),
    src = getattr(process,"ak4PFJets").src,
    srcPVs = getattr(process,"ak4PFJets").srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(250)
)


##########################################################
## MAKE PAT JETS
##########################################################

addJetCollection(
    process,
    labelName='PFCHS8',
    jetSource=cms.InputTag('PFJetsCHS8'),
    algo=algoLabel,           # needed for jet flavor clustering
    rParam=.8, # needed for jet flavor clustering
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK8,
    genJetCollection = cms.InputTag('genJetsNoNu8'),
    genParticles = cms.InputTag(genParticles),
    # postfix = postfix
)
getattr(process,'selectedPatJetsPFCHS8').cut = cms.string("abs(eta) < " + str(2.5))

addJetCollection(
    process,
    jetSource=cms.InputTag('PFJetsCHSSoftDrop8'),
    labelName='SoftDropPFCHS8',
    algo=algoLabel,
    btagInfos = ['None'],
    btagDiscriminators = ['None'],
    jetCorrections=jetCorrectionsAK8,
    genJetCollection = cms.InputTag('genJetsNoNu8'),
    genParticles = cms.InputTag(genParticles),
    getJetMCFlavour = False, # jet flavor disabled
    # postfix = postfix
)
addJetCollection(
    process,
    labelName='SoftDropSubjetsPFCHS8',
    jetSource=cms.InputTag('PFJetsCHSSoftDrop8','SubJets'),
    algo=algoLabel,           # needed for subjet flavor clustering
    rParam=.8, # needed for subjet flavor clustering
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag('genJetsNoNuSoftDrop8','SubJets'),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = True,  # needed for subjet b tagging
    svClustering = True, # needed for subjet b tagging
    fatJets = cms.InputTag('PFJetsCHS8'),              # needed for subjet flavor clustering
    groomedFatJets = cms.InputTag('PFJetsCHSSoftDrop8'), # needed for subjet flavor clustering
    runIVF = False,
    # postfix = postfix
)
## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsSoftDropPFCHSPacked8 = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsSoftDropPFCHS8"),
    subjetSrc=cms.InputTag("selectedPatJetsSoftDropSubjetsPFCHS8")
)

addJetCollection(
    process,
    labelName='PrunedPFCHS8',
    jetSource=cms.InputTag('PFJetsCHSPruned8'),
    algo=algoLabel,
    btagInfos = ['None'],
    btagDiscriminators = ['None'],
    jetCorrections=jetCorrectionsAK8,
    genJetCollection = cms.InputTag('genJetsNoNu8'),
    genParticles = cms.InputTag(genParticles),
    getJetMCFlavour = False, # jet flavor disabled
    # postfix = postfix
)
addJetCollection(
    process,
    labelName='PrunedSubjetsPFCHS8',
    jetSource=cms.InputTag('PFJetsCHSPruned8','SubJets'),
    algo=algoLabel,           # needed for subjet flavor clustering
    rParam=.8, # needed for subjet flavor clustering
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag('genJetsNoNuPruned8','SubJets'),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = True,  # needed for subjet b tagging
    svClustering = True, # needed for subjet b tagging
    fatJets = cms.InputTag('PFJetsCHS8'),              # needed for subjet flavor clustering
    groomedFatJets = cms.InputTag('PFJetsCHSPruned8'), # needed for subjet flavor clustering
    runIVF = False,
    # postfix = postfix
)

## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsPrunedPFCHSPacked8 = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsPrunedPFCHS8"),
    subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPFCHS8")
)

## Pack fat jets with subjets
process.packedPatJetsPFCHS8 = cms.EDProducer("JetSubstructurePacker",
        jetSrc = cms.InputTag('selectedPatJetsPFCHS8'),
        distMax = cms.double(.8),
        algoTags = cms.VInputTag(),
        algoLabels = cms.vstring(),
        fixDaughters = cms.bool(False)
)
process.packedPatJetsPFCHS8.algoTags.append( cms.InputTag('selectedPatJetsSoftDropPFCHSPacked8') )
process.packedPatJetsPFCHS8.algoLabels.append( 'SoftDrop' )
process.packedPatJetsPFCHS8.algoTags.append( cms.InputTag('selectedPatJetsPrunedPFCHSPacked8') )
process.packedPatJetsPFCHS8.algoLabels.append( 'Pruned' )

##

addJetCollection(
    process,
    labelName='PFCHS15',
    jetSource=cms.InputTag('PFJetsCHS15'),
    algo=algoLabel,           # needed for jet flavor clustering
    rParam=1.5, # needed for jet flavor clustering
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK8,
    genJetCollection = cms.InputTag('genJetsNoNu15'),
    genParticles = cms.InputTag(genParticles),
    # postfix = postfix
)
getattr(process,'selectedPatJetsPFCHS15').cut = cms.string("abs(eta) < " + str(2.5))

addJetCollection(
    process,
    labelName='SoftDropPFCHS15',
    jetSource=cms.InputTag('PFJetsCHSSoftDrop15'),
    algo=algoLabel,
    btagInfos = ['None'],
    btagDiscriminators = ['None'],
    jetCorrections=jetCorrectionsAK8,
    genJetCollection = cms.InputTag('genJetsNoNu15'),
    genParticles = cms.InputTag(genParticles),
    getJetMCFlavour = False, # jet flavor disabled
    # postfix = postfix
)
addJetCollection(
    process,
    labelName='SoftDropSubjetsPFCHS15',
    jetSource=cms.InputTag('PFJetsCHSSoftDrop15','SubJets'),
    algo=algoLabel,           # needed for subjet flavor clustering
    rParam=1.5, # needed for subjet flavor clustering
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag('genJetsNoNuSoftDrop15','SubJets'),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = True,  # needed for subjet b tagging
    svClustering = True, # needed for subjet b tagging
    fatJets = cms.InputTag('PFJetsCHS15'),              # needed for subjet flavor clustering
    groomedFatJets = cms.InputTag('PFJetsCHSSoftDrop15'), # needed for subjet flavor clustering
    runIVF = False,
    # postfix = postfix
)
## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsSoftDropPFCHSPacked15 = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsSoftDropPFCHS15"),
    subjetSrc=cms.InputTag("selectedPatJetsSoftDropSubjetsPFCHS15")
)

addJetCollection(
    process,
    labelName='PrunedPFCHS15',
    jetSource=cms.InputTag('PFJetsCHSPruned15'),
    algo=algoLabel,
    btagInfos = ['None'],
    btagDiscriminators = ['None'],
    jetCorrections=jetCorrectionsAK8,
    genJetCollection = cms.InputTag('genJetsNoNu15'),
    genParticles = cms.InputTag(genParticles),
    getJetMCFlavour = False, # jet flavor disabled
    # postfix = postfix
)
addJetCollection(
    process,
    labelName='PrunedSubjetsPFCHS15',
    jetSource=cms.InputTag('PFJetsCHSPruned15','SubJets'),
    algo=algoLabel,           # needed for subjet flavor clustering
    rParam=1.5, # needed for subjet flavor clustering
    pfCandidates = cms.InputTag(pfCandidates),
    pvSource = cms.InputTag(pvSource),
    svSource = cms.InputTag(svSource),
    muSource = cms.InputTag(muSource),
    elSource = cms.InputTag(elSource),
    btagInfos = bTagInfos,
    btagDiscriminators = bTagDiscriminators,
    jetCorrections = jetCorrectionsAK4,
    genJetCollection = cms.InputTag('genJetsNoNuPruned15','SubJets'),
    genParticles = cms.InputTag(genParticles),
    explicitJTA = True,  # needed for subjet b tagging
    svClustering = True, # needed for subjet b tagging
    fatJets = cms.InputTag('PFJetsCHS15'),              # needed for subjet flavor clustering
    groomedFatJets = cms.InputTag('PFJetsCHSPruned15'), # needed for subjet flavor clustering
    runIVF = False,
    # postfix = postfix
)

## Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsPrunedPFCHSPacked15 = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsPrunedPFCHS15"),
    subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPFCHS15")
)

## Pack fat jets with subjets
process.packedPatJetsPFCHS15 = cms.EDProducer("JetSubstructurePacker",
        jetSrc = cms.InputTag('selectedPatJetsPFCHS15'),
        distMax = cms.double(1.5),
        algoTags = cms.VInputTag(),
        algoLabels = cms.vstring(),
        fixDaughters = cms.bool(False)
)
process.packedPatJetsPFCHS15.algoTags.append( cms.InputTag('selectedPatJetsSoftDropPFCHSPacked15') )
process.packedPatJetsPFCHS15.algoLabels.append( 'SoftDrop' )
process.packedPatJetsPFCHS15.algoTags.append( cms.InputTag('selectedPatJetsPrunedPFCHSPacked15') )
process.packedPatJetsPFCHS15.algoLabels.append( 'Pruned' )


##############################################
## RUN
##############################################

process.p = cms.Path(
                     process.selectedPatJetsPFCHS8*
                     process.selectedPatJetsSoftDropPFCHS8*
                     process.selectedPatJetsSoftDropSubjetsPFCHS8*
                     process.selectedPatJetsSoftDropPFCHSPacked8*
                     process.selectedPatJetsPrunedPFCHS8*
                     process.selectedPatJetsPrunedSubjetsPFCHS8*
                     process.selectedPatJetsPrunedPFCHSPacked8*
                     process.packedPatJetsPFCHS8*
                     process.selectedPatJetsPFCHS15*
                     process.selectedPatJetsSoftDropPFCHS15*
                     process.selectedPatJetsSoftDropSubjetsPFCHS15*
                     process.selectedPatJetsSoftDropPFCHSPacked15*
                     process.selectedPatJetsPrunedPFCHS15*
                     process.selectedPatJetsPrunedSubjetsPFCHS15*
                     process.selectedPatJetsPrunedPFCHSPacked15*
                     process.packedPatJetsPFCHS15
                     )

from PhysicsTools.PatAlgos.tools.pfTools import *
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False),
        allowUnscheduled = cms.untracked.bool(True)
)
