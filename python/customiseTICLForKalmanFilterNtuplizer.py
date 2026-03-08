from Analyzers.Ntuplizer.Ntuplizer_cfi import *

def customiseTICLForKalmanFilterNtuplizer(process):
    process.kfNtuplizer = ntuplizer.clone()
    process.propNtuplizer = ntuplizer.clone(
           KFHits = cms.InputTag("ticlTrackstersStandalonePropagator","KFHits","RECO"),
    )
    process.kfNtuplizerG4e = ntuplizer.clone(
           KFHits = cms.InputTag("ticlTrackstersKalmanFilterG4e","KFHits","RECO")
    )
    process.propNtuplizerG4e = ntuplizer.clone(
           KFHits = cms.InputTag("ticlTrackstersStandalonePropagatorG4e","KFHits","RECO"),
    )
    return process
