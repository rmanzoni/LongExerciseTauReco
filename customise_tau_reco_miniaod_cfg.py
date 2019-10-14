import FWCore.ParameterSet.Config as cms

# trick to make python know where to look for the imports
import sys
sys.path.append('..')

# for example: here
from rerunTauRecoOnMiniAOD import process

updatedTauName = 'slimmedTausNewID'

runSignal = False # Set to False to read in QCD file instead of ZTT

maxEvents = 1000

# readFiles = cms.untracked.vstring()
readFiles = cms.untracked.vstring('file:/eos/cms/store/group/phys_tau/CMSSW_10_6_4_patch1_RelValTenTau_15_500_PU25ns_106X_upgrade2018_realistic_v9_HS_v1/25A52E71-2B61-6747-B2BA-1F51E581E6AC.root')
secFiles  = cms.untracked.vstring()
process.source = cms.Source( 'PoolSource', fileNames=readFiles, secondaryFileNames=secFiles)

print('\t Max events:', process.maxEvents.input.value())

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32( maxEvents ))

# print every event
process.MessageLogger.cerr.FwkReport.reportEvery = 1

##########################################################################################
# single thread for debugging
# process.options.numberOfThreads = cms.untracked.uint32(1)

##########################################################################################
# process.combinatoricRecoTaus.builders[0].signalConeSize = cms.string('max(min(0.1, 3.0/pt()), 0.05)')
process.combinatoricRecoTaus.builders[0].signalConeSize = cms.string('1.e-7')
# process.combinatoricRecoTaus.builders[0].verbosity = cms.int32(3)

##########################################################################################
process.output.outputCommands.append('keep *_'+updatedTauName+'_*_*')
process.output.fileName = cms.untracked.string('miniAOD_rerunTauRECO_newIDs.root')

# Add new TauIDs
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(
    process, 
    cms, 
    debug = True,
    updatedTauName = updatedTauName,
    toKeep = [ 
        'newDM2017v2', 
        'deepTau2017v2', 
    ],
    inputTaus = 'selectedPatTaus',
    )
tauIdEmbedder.runTauID()

# Run DeepTau 
process.TauReco.insert(-1, process.rerunMvaIsolationSequence)
process.TauReco.insert(-1, getattr(process,updatedTauName))


# Schedule definition
# process.schedule = cms.Schedule(deeptau_process.p,process.endjob,process.outpath)





















































## E.G.: change the isolation cone
#process.combinatoricRecoTaus.builders[0].isolationConeSize = cms.double(0.8) # originally 0.5

## E.G.: Remove all decay modes with piZeros
#process.combinatoricRecoTaus.builders[0].decayModes = [dm for dm in process.combinatoricRecoTaus.builders[0].decayModes if dm.nPiZeros==0]

## Some more settings (default values given):
'''
process.combinatoricRecoTaus.builders[0].decayModes = cms.VPSet(
        cms.PSet(
            maxPiZeros = cms.uint32(0),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(0)
        ),
        cms.PSet(
            maxPiZeros = cms.uint32(6),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(1)
        ),
        cms.PSet(
            maxPiZeros = cms.uint32(5),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(1),
            nPiZeros = cms.uint32(2)
        ),
        cms.PSet(
            maxPiZeros = cms.uint32(0),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(2),
            nPiZeros = cms.uint32(0)
        ),
        cms.PSet(
            maxPiZeros = cms.uint32(3),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(2),
            nPiZeros = cms.uint32(1)
        ),
        cms.PSet(
            maxPiZeros = cms.uint32(0),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(3),
            nPiZeros = cms.uint32(0)
        ),
        cms.PSet(
            maxPiZeros = cms.uint32(3),
            maxTracks = cms.uint32(6),
            nCharged = cms.uint32(3),
            nPiZeros = cms.uint32(1)
        )
    )

process.combinatoricRecoTaus.builders[0].qualityCuts = cms.PSet(
            isolationQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.2),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.03),
                minGammaEt = cms.double(1.5),
                minTrackHits = cms.uint32(8),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(1.0),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            leadingTrkOrPFCandOption = cms.string('leadPFCand'),
            primaryVertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices'),
            pvFindingAlgo = cms.string('closestInDeltaZ'),
            recoverLeadingTrk = cms.bool(False),
            signalQualityCuts = cms.PSet(
                maxDeltaZ = cms.double(0.4),
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(1.0),
                minNeutralHadronEt = cms.double(30.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            ),
            vertexTrackFiltering = cms.bool(False),
            vxAssocQualityCuts = cms.PSet(
                maxTrackChi2 = cms.double(100.0),
                maxTransverseImpactParameter = cms.double(0.1),
                minGammaEt = cms.double(1.0),
                minTrackHits = cms.uint32(3),
                minTrackPixelHits = cms.uint32(0),
                minTrackPt = cms.double(0.5),
                minTrackVertexWeight = cms.double(-1.0)
            )
        )

process.hpsPFTauDiscriminationByDecayModeFindingNewDMs.decayModes(
        cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(True),
            mass = cms.bool(True),
            phi = cms.bool(True)
        ),
        assumeStripMass = cms.double(-1.0),
        maxMass = cms.string('1.'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(-1000.0),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(1),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(0),
        nTracksMin = cms.uint32(1)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(True),
            mass = cms.bool(True),
            phi = cms.bool(True)
        ),
        assumeStripMass = cms.double(0.1349),
        maxMass = cms.string('max(1.3, min(1.3*sqrt(pt/100.), 4.2))'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(0.3),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(1),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(1),
        nTracksMin = cms.uint32(1)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(True),
            mass = cms.bool(True),
            phi = cms.bool(True)
        ),
        assumeStripMass = cms.double(0.0),
        maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
        maxPi0Mass = cms.double(0.2),
        minMass = cms.double(0.4),
        minPi0Mass = cms.double(0.05),
        nCharged = cms.uint32(1),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(2),
        nTracksMin = cms.uint32(1)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(False),
            mass = cms.bool(False),
            phi = cms.bool(False)
        ),
        assumeStripMass = cms.double(-1.0),
        maxMass = cms.string('1.2'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(0.0),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(2),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(0),
        nTracksMin = cms.uint32(2)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(False),
            mass = cms.bool(False),
            phi = cms.bool(False)
        ),
        assumeStripMass = cms.double(-1.0),
        maxMass = cms.string('max(1.2, min(1.2*sqrt(pt/100.), 4.0))'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(0.0),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(2),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(1),
        nTracksMin = cms.uint32(2)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(False),
            mass = cms.bool(False),
            phi = cms.bool(False)
        ),
        assumeStripMass = cms.double(-1.0),
        maxMass = cms.string('1.5'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(0.8),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(3),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(0),
        nTracksMin = cms.uint32(2)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(False),
            mass = cms.bool(False),
            phi = cms.bool(False)
        ),
        assumeStripMass = cms.double(-1.0),
        maxMass = cms.string('1.5'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(0.8),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(3),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(0),
        nTracksMin = cms.uint32(2)
    ),
    cms.PSet(
        applyBendCorrection = cms.PSet(
            eta = cms.bool(False),
            mass = cms.bool(False),
            phi = cms.bool(False)
        ),
        assumeStripMass = cms.double(-1.0),
        maxMass = cms.string('1.6'),
        maxPi0Mass = cms.double(1000000000.0),
        minMass = cms.double(0.9),
        minPi0Mass = cms.double(-1000.0),
        nCharged = cms.uint32(3),
        nChargedPFCandsMin = cms.uint32(1),
        nPiZeros = cms.uint32(1),
        nTracksMin = cms.uint32(2)
    )
)

)
'''


# change the output file name, don't overwrite the original file!
# process.output.fileName = cms.untracked.string('{}_miniAOD_rerunTauRECO.root'.format('ZTT' if runSignal else 'QCD'))
