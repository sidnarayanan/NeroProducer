import FWCore.ParameterSet.Config as cms

process = cms.Process("nero")

process.load("FWCore.MessageService.MessageLogger_cfi")
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# fileList = ['/store/relval/CMSSW_7_4_1/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/80CF5456-B9EC-E411-93DA-002618FDA248.root']
fileList = ['/store/relval/CMSSW_7_4_1/RelValProdTTbar_13/MINIAODSIM/MCRUN2_74_V9_gensim71X-v1/00000/0A9E2CED-C9EC-E411-A8E4-003048FFCBA8.root']
# fileList = ['file:/afs/cern.ch/user/s/snarayan/work/private/CMSSW_7_4_5/src/NeroProducer/Skim/monotop/skim.root']


### do not remove the line below!
###FILELIST###

process.source = cms.Source("PoolSource",
    	fileNames = cms.untracked.vstring(fileList)
    )

# ---- define the output file -------------------------------------------
process.TFileService = cms.Service("TFileService",
			closeFileFast = cms.untracked.bool(True),
			fileName = cms.string("NeroNtuples.root"),
                )
process.dump=cms.EDAnalyzer('EventContentAnalyzer')

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


# ------------------------QG-----------------------------------------------
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets            = cms.InputTag("slimmedJets")    # Could be reco::PFJetCollection or pat::JetCollection (both AOD and miniAOD)
process.QGTagger.srcVertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
process.QGTagger.useQualityCuts = cms.bool(False)
##----------------GLOBAL TAG ---------------------------
# used by photon id and jets
# process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#data
#process.GlobalTag.globaltag = 'GR_R_52_V9::All'
#mc https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions#Global_Tags_for_Run2_MC_Producti
process.GlobalTag.globaltag = 'MCRUN2_74_V9::All'

#-----------------------ELECTRON ID-------------------------------
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
useAOD=False

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_PHYS14_PU20bx25_V2_cff']
#                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV51_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

### PHOTONS
switchOnVIDPhotonIdProducer(process, dataFormat) ### PHOTON
pho_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']

for idmod in pho_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
##ISO
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")
process.load("RecoEgamma/ElectronIdentification/ElectronIDValueMapProducer_cfi")

## SKIM INFO
process.load('NeroProducer.Skim.infoProducerSequence_cff')

process.load('NeroProducer.Nero.Nero_cfi')
process.load('NeroProducer.Nero.NeroMonotop_cfi')
#process.load('NeroProducer.Nero.NeroChargedHiggs_cfi')

from PhysicsTools.PatAlgos.tools.pfTools import *
adaptPVs(process, pvCollection=cms.InputTag('offlineSlimmedPrimaryVertices'))

## Add full JetFlavourInfo and TagInfos to PAT jets
for r0 in ['8','15']:
    for m in ['patJetsPFCHS'+r0, 'patJetsSoftDropSubjetsPFCHS'+r0, 'patJetsPrunedSubjetsPFCHS'+r0]:
        if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
            print "Switching 'addTagInfos' for " + m + " to 'True'"
            setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
        if hasattr(process,m):
            print "Switching 'addJetFlavourInfo' for " + m + " to 'True'"
            setattr( getattr(process,m), 'addJetFlavourInfo', cms.bool(True) )

#------------------------------------------------------
process.p = cms.Path(
                # process.dump*
		        process.infoProducerSequence *
                # process.selectedPatJetsPFCHS8*
                # process.selectedPatJetsSoftDropPFCHS8*
                # process.selectedPatJetsSoftDropSubjetsPFCHS8*
                # process.selectedPatJetsSoftDropPFCHSPacked8*
                # process.selectedPatJetsPrunedPFCHS8*
                # process.selectedPatJetsPrunedSubjetsPFCHS8*
                # process.selectedPatJetsPrunedPFCHSPacked8*
                # process.packedPatJetsPFCHS8*
                # process.selectedPatJetsPFCHS15*
                # process.selectedPatJetsSoftDropPFCHS15*
                # process.selectedPatJetsSoftDropSubjetsPFCHS15*
                # process.selectedPatJetsSoftDropPFCHSPacked15*
                # process.selectedPatJetsPrunedPFCHS15*
                # process.selectedPatJetsPrunedSubjetsPFCHS15*
                # process.selectedPatJetsPrunedPFCHSPacked15*
                # process.packedPatJetsPFCHS15 *
                process.QGTagger *
                process.egmGsfElectronIDSequence *
                process.egmPhotonIDSequence *
                process.photonIDValueMapProducer * ## ISO MAP FOR PHOTONS
                process.electronIDValueMapProducer * ## ISO MAP FOR PHOTONS
                process.nero
                )

process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(False),
        allowUnscheduled = cms.untracked.bool(True)
)
