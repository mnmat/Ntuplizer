import FWCore.ParameterSet.Config as cms

ntuplizer = cms.EDAnalyzer('Ntuplizer',
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   Tracksters = cms.InputTag("ticlTrackstersMerge","","RECO"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   KFHits = cms.InputTag("ticlTrackstersKalmanFilter","KFHits","RECO"),
   #PropHits = cms.InputTag("ticlTrackstersKalmanFilter","KFHits","RECO"),
   PropHits = cms.InputTag("ticlTrackstersStandalonePropagator","KFHits","RECO"),
   hgcalLayerClusters = cms.InputTag("hgcalLayerClusters", "", "RECO"),
   tracks    = cms.untracked.InputTag('generalTracks'),
   #tracks = cms.untracked.InputTag('standAloneMuons'),
   lcMask = cms.InputTag("ticlTrackstersCLUE3DHigh",""),
   associators = cms.untracked.VInputTag("trackingParticleRecoTrackAsssociation"),
   #associators = cms.untracked.VInputTag("tpToStaMuonAssociation"),
   simVertices = cms.InputTag("g4SimHits"),
   #trackPtMin = cms.double(0.3)
)
