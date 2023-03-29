import FWCore.ParameterSet.Config as cms
import argparse

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D99Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

import sys

eta = sys.argv[2]
energy = sys.argv[3]
#eta = "16"
#energy = "10"
mb = "mb_ngun"
nevents = "500"
cap = "zpos"

fname = '/eos/home-m/mmatthew/Data/KF/CMSSW_13_1_0_pre1/Analyzer/UpdatorStudies/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+fname + 'step3_singlemuon_e'+energy+'GeV_eta'+eta+'_'+cap+'_events'+nevents+'_nopu.root'))


#outfile_ = 'file:/eos/home-m/mmatthew/Data/deleteme.root'
#fname = '/eos/home-m/mmatthew/Data/Analyzer/UpdatorStudies/'+propagator+'/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
#fname = '/eos/home-m/mmatthew/Data/KF/MaterialBudget/Radlen/0_25/'+ cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
fname = '/eos/home-m/mmatthew/Data/KF/CMSSW_13_1_0_pre1/Ntuplizer/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'
outfile_  = 'file:'+fname + 'ttree_singlemuon_e'+energy+'GeV_eta'+eta+'_'+cap+'_events'+nevents+'_nopu.root'

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(outfile_),
                                   closeFileFast = cms.untracked.bool(True)
                               )


process.demo = cms.EDAnalyzer('Ntuplizer',
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   Tracksters = cms.InputTag("ticlTrackstersMerge","","RECO"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   #propagatorEM = cms.InputTag("ticlTrackstersEM","TEST","RECO"),
   #propagatorHAD = cms.InputTag("ticlTrackstersHAD","TEST","RECO"),
   propagatorKF = cms.InputTag("ticlTrackstersKF","Points KF","RECO"),
   xxKF = cms.InputTag("ticlTrackstersKF","xx KF","RECO"),
   xyKF = cms.InputTag("ticlTrackstersKF","xy KF","RECO"),
   yyKF = cms.InputTag("ticlTrackstersKF","yy KF","RECO"),
   xxProp = cms.InputTag("ticlTrackstersKF","xx Prop","RECO"),
   xyProp = cms.InputTag("ticlTrackstersKF","xy Prop","RECO"),
   yyProp = cms.InputTag("ticlTrackstersKF","yy Prop","RECO"),
   abs_fail = cms.InputTag("ticlTrackstersKF","Abs Fail","RECO"),
   kfcharge = cms.InputTag("ticlTrackstersKF","charge KF","RECO"),
   propcharge = cms.InputTag("ticlTrackstersKF","charge Prop","RECO"),
   kfdetid = cms.InputTag("ticlTrackstersKF","detID KF","RECO"),
   propdetid = cms.InputTag("ticlTrackstersKF","detID Prop","RECO"),
   propagator = cms.InputTag("ticlTrackstersKF","Points Prop","RECO"),
   #propagatorTrk = cms.InputTag("ticlTrackstersTrk","TEST","RECO"),
   #propagatorTrkEM = cms.InputTag("ticlTrackstersTrkEM","TEST","RECO"),
   hgcalLayerClusters = cms.InputTag("hgcalLayerClusters", "", "RECO"),
   eta = cms.string(eta),
   energy = cms.string(energy), 
   outdir = cms.string(fname), 
   #trackPtMin = cms.double(0.3)
                              )
process.p = cms.Path(process.demo)

