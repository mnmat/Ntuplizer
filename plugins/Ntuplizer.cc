// -*- C++ -*-
//
// Package:    Analyzer/Ntuplizer
// Class:      Ntuplizer
//
/**\class Ntuplizer Ntuplizer.cc Analyzer/Ntuplizer/plugins/Ntuplizer.cc

 Description: [one line class summary]

 Implementation:f
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Thu, 30 Jun 2022 10:49:35 GMT
//
//

// system include files
#include <memory>
#include <numeric>
#include <sstream>
#include <any>
#include <iomanip>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <fstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TProfile.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/KFHit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/Vector3D.h"


//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include <TVector.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include <algorithm>
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;
using namespace ticl;

class Ntuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Ntuplizer(const edm::ParameterSet&);
  ~Ntuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void clear_arrays();
  virtual void fillHitMap(std::map<DetId, std::pair<const HGCRecHit*,float>>& hitMap, 
      const HGCRecHitCollection& rechitsEE, 
      const HGCRecHitCollection& rechitsFH,
      const HGCRecHitCollection& rechitsBH,
      const std::vector<float>& lcmask) const;
  std::vector<int> matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_);
  LocalError calculateLocalError(DetId detid_, const HGCalDDDConstants* ddd);


  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<KFHit>> KFHitsToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<float>> lcMaskToken_;


  TTree *tree = new TTree("tree","tree");

  std::vector<std::string> detectors, objects, positions, hittypes;

  // variables
  
  int eventnr =0;
  std::string eta_;
  std::string energy_;
  std::string outdir_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;

  // Track

  std::vector<int> track_id;
  std::vector<int> track_charge;
  std::vector<float> track_momentum;
  std::vector<float> track_quality;
  std::vector<float> track_chi2;
  std::vector<float> track_validfraction;
  std::vector<float> track_qoverp;
  std::vector<float> track_algo;


  // KF

  std::vector<float> kf_x;
  std::vector<float> kf_y;
  std::vector<float> kf_z;
  std::vector<float> kf_e;
  std::vector<float> kf_cov_xx;
  std::vector<float> kf_cov_xy;
  std::vector<float> kf_cov_yy;
  std::vector<int> kf_detid;
  std::vector<int> kf_charge;
  std::vector<int> kf_layer;
  std::vector<std::string> kf_dtype;
  std::vector<int> kf_evt;
  std::vector<float> kf_eta;
  std::vector<float> kf_theta;
  std::vector<float> kf_trackid;
  std::vector<float> kf_trackcharge;
  std::vector<float> kf_trackmomentum;
  std::vector<float> kf_trackquality;
  std::vector<float> kf_trackchi2;
  std::vector<float> kf_track_validfraction;
  std::vector<float> kf_track_qoverp;
  std::vector<float> kf_track_algo;

    // RecHits

  std::vector<float> rec_x;
  std::vector<float> rec_y;
  std::vector<float> rec_z;
  std::vector<float> rec_e;
  std::vector<float> rec_cov_xx;
  std::vector<float> rec_cov_xy;
  std::vector<float> rec_cov_yy;
  std::vector<int> rec_detid;
  std::vector<int> rec_layer;
  std::vector<std::string> rec_dtype;
  std::vector<int> rec_evt;
  std::vector<int> rec_kf_compatible;
  std::vector<int> rec_obj_id;
  std::vector<int> rec_pid;
  std::vector<float> rec_mask;


  // SimHits

  std::vector<float> sim_x;
  std::vector<float> sim_y;
  std::vector<float> sim_z;
  std::vector<float> sim_e;
  std::vector<float> sim_cov_xx;
  std::vector<float> sim_cov_xy;
  std::vector<float> sim_cov_yy;
  std::vector<int> sim_detid;
  std::vector<int> sim_layer;
  std::vector<std::string> sim_dtype;
  std::vector<int> sim_evt;
  std::vector<int> sim_kf_compatible;
  std::vector<int> sim_obj_id;
  std::vector<int> sim_pid;
  std::vector<float> sim_mask;


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Ntuplizer::Ntuplizer(const edm::ParameterSet& iConfig) :
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      //abs_failToken_(consumes<float>(iConfig.getParameter<edm::InputTag>("abs_fail"))),
      KFHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("KFHits"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      tracksToken_(consumes<edm::View<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      lcMaskToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("lcMask"))),
      eta_(iConfig.getParameter<std::string>("eta")),
      energy_(iConfig.getParameter<std::string>("energy")),
      outdir_(iConfig.getParameter<std::string>("outdir")){



  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  hittypes = {"Simhits","Rechits","KF"};
  objects = {"Simhits", "Rechits"};        
  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed
  
  usesResource("TFileService");
  edm::Service<TFileService> file;

  tree = file->make<TTree>("tree","tree");

  // SimHits

  tree->Branch("sim_x", &sim_x);
  tree->Branch("sim_y", &sim_y);
  tree->Branch("sim_z", &sim_z);
  tree->Branch("sim_e", &sim_e);
  tree->Branch("sim_layer", &sim_layer);
  tree->Branch("sim_detid", &sim_detid);
  tree->Branch("sim_dtype", &sim_dtype);
  tree->Branch("sim_cov_xx", &sim_cov_xx);
  tree->Branch("sim_cov_xy", &sim_cov_xy);
  tree->Branch("sim_cov_yy", &sim_cov_yy);
  tree->Branch("sim_evt", &sim_evt);
  tree->Branch("sim_kf_compatible", &sim_kf_compatible);
  tree->Branch("sim_obj_id", &sim_obj_id);
  tree->Branch("sim_pid", &sim_pid);
  tree->Branch("sim_mask", &sim_mask);

  // RecHits

  tree->Branch("rec_x", &rec_x);
  tree->Branch("rec_y", &rec_y);
  tree->Branch("rec_z", &rec_z);
  tree->Branch("rec_e", &rec_e);
  tree->Branch("rec_layer", &rec_layer);
  tree->Branch("rec_detid", &rec_detid);
  tree->Branch("rec_dtype", &rec_dtype);
  tree->Branch("rec_cov_xx", &rec_cov_xx);
  tree->Branch("rec_cov_xy", &rec_cov_xy);
  tree->Branch("rec_cov_yy", &rec_cov_yy);
  tree->Branch("rec_evt", &rec_evt);
  tree->Branch("rec_kf_compatible", &rec_kf_compatible);
  tree->Branch("rec_obj_id", &rec_obj_id);
  tree->Branch("rec_pid", &rec_pid);
  tree->Branch("rec_mask", &rec_mask);

  // KF

  tree->Branch("kf_x", &kf_x);
  tree->Branch("kf_y", &kf_y);
  tree->Branch("kf_z", &kf_z);
  tree->Branch("kf_e", &kf_e);
  tree->Branch("kf_layer", &kf_layer);
  tree->Branch("kf_detid", &kf_detid);
  tree->Branch("kf_dtype", &kf_dtype);
  tree->Branch("kf_cov_xx", &kf_cov_xx);
  tree->Branch("kf_cov_xy", &kf_cov_xy);
  tree->Branch("kf_cov_yy", &kf_cov_yy);
  tree->Branch("kf_evt", &kf_evt);
  tree->Branch("kf_eta", &kf_eta);
  tree->Branch("kf_theta", &kf_theta);
  tree->Branch("kf_track_id", &track_id);
  tree->Branch("kf_track_charge", &track_charge);
  tree->Branch("kf_track_momentum", &track_momentum);
  tree->Branch("kf_track_quality", &track_quality);
  tree->Branch("kf_track_chi2", &track_chi2);
  tree->Branch("kf_track_validfraction", &track_validfraction);
  tree->Branch("kf_track_qoverp", &track_qoverp);
  tree->Branch("kf_track_algo", &track_algo);


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

Ntuplizer::~Ntuplizer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty

  tree->Write();

}

//
// member functions
//

// ------------ method called for each event  ------------

LocalError Ntuplizer::calculateLocalError(DetId id, const HGCalDDDConstants* ddd){
  if(recHitTools_.isSilicon(id)){
    float A;
    if(recHitTools_.getSiThickness(id) < 200) A = 1.18; // TODO: replace with non-hardcoded value; hardcoded value from TDR
    else  A = 0.52; // TODO: replace with non-hardcoded value; hardcoded value from TDR
    float a = sqrt(2*A/(3*sqrt(3)));
    double varx = pow(a,4)*5*sqrt(3)/(16*A); // x
    double vary = pow(a,4)*5*sqrt(3)/(16*A); // y 
    return LocalError(varx, 0, vary);
  }
  else{
    const GlobalPoint &pos = recHitTools_.getPosition(id);
    double r = sqrt(pos.x()*pos.x() + pos.y()*pos.y());
    auto radiusLayer = ddd->getRadiusLayer(recHitTools_.getLayer(id));
    int idx = static_cast<int>(std::lower_bound(radiusLayer.begin(), radiusLayer.end(),r)-radiusLayer.begin());
    float rmax = radiusLayer[idx];
    float rmin = radiusLayer[idx-1];

    double phi = recHitTools_.getPhi(id) + M_PI; // radians [0, 2pi]
    double dphi = recHitTools_.getScintDEtaDPhi(id).second; // radians
    double phimin = phi - 0.5*dphi;
    double phimax = phi + 0.5*dphi;

    double A = (rmax*rmax - rmin*rmin)*M_PI*dphi/(2*M_PI);

    double ex2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin - sin(phimin)*cos(phimin) + phimax + sin(phimax)*cos(phimax));
    double ex = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (sin(phimax) - sin(phimin));
    double varx = ex2 - ex*ex;

    double ey2 = 1/(8*A) * (pow(rmax,4) - pow(rmin,4)) * (-phimin + sin(phimin)*cos(phimin) + phimax - sin(phimax)*cos(phimax));
    double ey = 1/(3*A) * (pow(rmax,3) - pow(rmin,3)) * (cos(phimin) - cos(phimax));
    double vary = ey2 - ey*ey;

    double varxy = 1/(16*A)*(pow(rmax,4)-pow(rmin,4))*(cos(2*phimin)-cos(2*phimax)) - ex*ey;
    return LocalError(varx, varxy, vary);
  }
} 

void Ntuplizer::clear_arrays(){

  // SimHit

  sim_x.clear();
  sim_y.clear();
  sim_z.clear();
  sim_e.clear();
  sim_detid.clear();
  sim_layer.clear();
  sim_dtype.clear();
  sim_cov_xx.clear();
  sim_cov_xy.clear();
  sim_cov_yy.clear();
  sim_evt.clear();
  sim_kf_compatible.clear();
  sim_obj_id.clear();
  sim_pid.clear();
  sim_mask.clear();


  // RecHit

  rec_x.clear();
  rec_y.clear();
  rec_z.clear();
  rec_e.clear();
  rec_detid.clear();
  rec_layer.clear();
  rec_dtype.clear();
  rec_cov_xx.clear();
  rec_cov_xy.clear();
  rec_cov_yy.clear();
  rec_evt.clear();
  rec_kf_compatible.clear();
  rec_obj_id.clear();
  rec_pid.clear();
  rec_mask.clear();


  // KF

  kf_x.clear();
  kf_y.clear();
  kf_z.clear();
  kf_e.clear();
  kf_detid.clear();
  kf_layer.clear();
  kf_dtype.clear();
  kf_cov_xx.clear();
  kf_cov_xy.clear();
  kf_cov_yy.clear();
  kf_evt.clear();
  kf_eta.clear();
  kf_theta.clear();
  track_id.clear();
  track_charge.clear();
  track_momentum.clear();
  track_quality.clear();
  track_chi2.clear();
  track_validfraction.clear();
  track_qoverp.clear();
  track_algo.clear();
}

void Ntuplizer::fillHitMap(std::map<DetId, std::pair<const HGCRecHit*, float>>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH,
                                const std::vector<float>& lcmask) const {
  hitMap.clear();
  int counter = 0;
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), std::make_pair(&hit, lcmask[counter]));
    counter++;
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), std::make_pair(&hit, lcmask[counter]));
    counter++;
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), std::make_pair(&hit, lcmask[counter]));
    counter++;
  }
} // end of EfficiencyStudies::fillHitMap


std::vector<int> Ntuplizer::matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_) {
  std::vector<int> matchedIdxs; matchedIdxs.clear();
  for (unsigned int i0=0; i0<rechitdetid_.size(); ++i0) {
    if (detid_ == rechitdetid_[i0]) { matchedIdxs.push_back(i0); }
  }
  return matchedIdxs;
} // end of matchRecHit2CPRecHits


void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //std::cout << "Analyze" << std::endl;

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;
  
  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  edm::Handle<std::vector<float>> lcMaskHandle;
  iEvent.getByToken(lcMaskToken_, lcMaskHandle);

  std::map<DetId, std::pair<const HGCRecHit*, float>> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH, *lcMaskHandle);

  edm::Handle<std::vector<KFHit>> KFHitsHandle;
  iEvent.getByToken(KFHitsToken_, KFHitsHandle);
  const std::vector<KFHit> &kfhits = *KFHitsHandle;

  edm::Handle<edm::View<reco::Track>> tracks_h;
  iEvent.getByToken(tracksToken_,tracks_h);
  const edm::View<reco::Track> & tkx = *(tracks_h.product()); 




  /*
  // Match tracks to simtrack
  edm::Handle<edm::View<reco::Track>> tracks_h;
  iEvent.getByToken(tracksToken_,tracks_h);
  const edm::View<reco::Track> & tkx = *(tracks_h.product());
  
  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  iEvent.getByToken(associatormapRtSs[0], recotosimCollectionH);
  recSimCollP = recotosimCollectionH.product();
  reco::RecoToSimCollection const& recSimColl = *recSimCollP;
  */

  //std::map<float, GlobalPoint> map_gps_kf;
  std::map<int, std::vector<int>> map_detid_kf;
  std::map<float, float> map_xx_kf, map_xy_kf, map_yy_kf;

    // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

  const CaloSubdetectorGeometry *subGeom = geom.getSubdetectorGeometry(DetId::Detector(10), ForwardSubdetector::ForwardEmpty);
  auto geomEE = static_cast<const HGCalGeometry*>(subGeom);
  const HGCalDDDConstants* ddd = &(geomEE->topology().dddConstants());
  auto radiusLayer = ddd->rangeRLayer(8, true);
  std::cout << radiusLayer.first << ", " << radiusLayer.second << std::endl;

  // Test KFHits

/*
  std::cout << "Test KFHits" << std::endl;
  for (const auto& hit: kfhits){
    std::cout << hit.center.x() << "\t" << hit.center.y() << "\t" << hit.center.z() << "\t" << hit.xx << "\t" << hit.xy << "\t" << hit.yy << "\t" << hit.charge << "\t" << hit.detid << "\t" <<std::endl;
    std::cout << hit.test[0][0] << "," << hit.test[1][1] << hit.test[2][2] << hit.test[3][3] << hit.test[4][4] << hit.test[5][5] << std::endl;
  }
  std::cout << "------------------------------------------" << std::endl;
*/
  // KF Hits

  for(int i = 0;i<int(kfhits.size());i++){

    std::map<DetId,std::pair<const HGCRecHit *,float>>::const_iterator itcheck = hitMap.find(kfhits[i].detid);
    float e = 0;
    if (itcheck != hitMap.end()){
      e = hitMap[kfhits[i].detid].first->energy();
    }
    //unsigned int layer_ = recHitTools_.getLayerWithOffset(hits[i].detid);

    std::string detector;
    std::string thickness;
    std::string tmp;
    
    // TODO: if conditions necessary to deal with missing rechits. Not yet possible to determine area of detector properly
    //       Implement it in such a way that we can find the closest detid and from this determine silicon thickness.
    if (kfhits[i].detid==0){
      //auto closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(0))->getClosestCell(hits[i].center);
      //layer_ = recHitTools_.getLayerWithOffset(closest_detid); // Can get rid of this in future iterations as I read it out in KFHits anyway
      detector = "Sc";
      thickness = "None";
      tmp = "Sc";
    }
    else if((kfhits[i].detid==8)||(kfhits[i].detid==9)){
      //auto closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(hits[i].detid))->getClosestCellHex(hits[i].center,true);
      //layer_ = recHitTools_.getLayerWithOffset(closest_detid); // Can get rid of this in future iterations as I read it out in KFHits anyway
      detector = "Si";
      //thickness = std::to_string(int(recHitTools_.getSiThickness(closest_detid))); 
      tmp = detector;
    }
    else if (recHitTools_.isSilicon(kfhits[i].detid)){
      detector = "Si";
      thickness = std::to_string(int(recHitTools_.getSiThickness(kfhits[i].detid))); 
      tmp = detector+" "+thickness;
    }
    else{
      detector = "Sc";
      thickness = "None";
      tmp = "Sc";
    } 

    int tkId = kfhits[i].trackId;
    std::cout << "Track Reco algo: " << tkx[tkId].algoName()<< ", "<< tkx[tkId].algo() << std::endl;
    std::cout << "Track Quality: " <<  tkx[tkId].qualityMask() << std::endl;

    map_detid_kf[kfhits[i].layer].push_back(kfhits[i].detid);
    //map_gps[hits[i].center.z()]=hits[i].center;
    kf_x.push_back(kfhits[i].center.x());
    kf_y.push_back(kfhits[i].center.y());
    kf_z.push_back(kfhits[i].center.z());
    kf_e.push_back(e);
    kf_detid.push_back(kfhits[i].detid);
    kf_layer.push_back(kfhits[i].layer); // 
    kf_dtype.push_back(tmp);
    kf_cov_xx.push_back(kfhits[i].xx);
    kf_cov_xy.push_back(kfhits[i].xy);
    kf_cov_yy.push_back(kfhits[i].yy);
    kf_evt.push_back(eventnr);
    kf_eta.push_back(kfhits[i].eta);
    kf_theta.push_back(kfhits[i].theta);
    track_id.push_back(kfhits[i].trackId);
    track_charge.push_back(kfhits[i].trackCharge);
    track_momentum.push_back(kfhits[i].trackMomentum);
    track_quality.push_back(kfhits[i].trackQuality);
    track_chi2.push_back(kfhits[i].trackChi2);
    track_validfraction.push_back(kfhits[i].trackValidFraction);
    track_qoverp.push_back(kfhits[i].trackQOverP);
    track_algo.push_back(tkx[tkId].algo());
  }
  

  // Loop over Caloparticles 

  std::cout << "Hello" << std::endl;

  std::vector<DetId> tmprechits_; tmprechits_.clear();
  int obj_id = 0;
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    auto pid = cp.particleId();
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){

        DetId detid_ = (it_sc_hae.first);
        std::map<DetId,std::pair<const HGCRecHit *, float>>::const_iterator itcheck = hitMap.find(detid_);
        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);

        std::string detector;
        std::string thickness;
        std::string tmp;
        if(recHitTools_.isSilicon(detid_)){
          detector = "Si";
          thickness = std::to_string(int(recHitTools_.getSiThickness(detid_))); 
          tmp = detector+" "+thickness;
        }
        else{
          detector = "Sc";
          thickness = "None";
          tmp = "Sc";
        } 

        // LocalError

        const CaloSubdetectorGeometry *subGeom = geom.getSubdetectorGeometry(detid_);
        auto geomEE = static_cast<const HGCalGeometry*>(subGeom);
        const HGCalDDDConstants* ddd = &(geomEE->topology().dddConstants());
        auto lerr = calculateLocalError(detid_, ddd); 
        // std::cout << "Layer: " << layer_ << "\t XX: " << lerr.xx() << std::endl;

        // KF 
        int kf_compatible=0;

        int layer = recHitTools_.getLayerWithOffset(detid_);

        /*
        DetId closest_detid;
        if (detector == "Sc") closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCell(gp);
        else closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCellHex(gp, true);
        if(detid_==closest_detid){
          compatible=1;
        }
        */
        if (std::find(map_detid_kf[layer].begin(), map_detid_kf[layer].end(), static_cast<int32_t>(detid_())) != map_detid_kf[layer].end()) {
          kf_compatible=1;
        }
        /*
        if (detid_()==static_cast<uint32_t>(map_detid_kf[layer])){
          kf_compatible=1;
        }
        */
        
        std::cout << (itcheck->second) << std::endl;

        sim_x.push_back(recHitTools_.getPosition(detid_).x());
        sim_y.push_back(recHitTools_.getPosition(detid_).y());
        sim_z.push_back(recHitTools_.getPosition(detid_).z());
        sim_e.push_back(it_sc_hae.second);
        sim_detid.push_back(detid_);
        sim_layer.push_back(layer_);
        sim_dtype.push_back(tmp);
        sim_cov_xx.push_back(lerr.xx());
        sim_cov_xy.push_back(lerr.xy());
        sim_cov_yy.push_back(lerr.yy());
        sim_kf_compatible.push_back(kf_compatible);
        sim_evt.push_back(eventnr);
        sim_obj_id.push_back(obj_id);
        sim_pid.push_back(pid);
        //sim_mask.push_back((itcheck->second));


        if (itcheck != hitMap.end()){
          rec_x.push_back(recHitTools_.getPosition(detid_).x());
          rec_y.push_back(recHitTools_.getPosition(detid_).y());
          rec_z.push_back(recHitTools_.getPosition(detid_).z());
          rec_e.push_back(it_sc_hae.second);
          rec_detid.push_back(detid_);
          rec_layer.push_back(layer_);
          rec_dtype.push_back(tmp);
          rec_cov_xx.push_back(lerr.xx());
          rec_cov_xy.push_back(lerr.xy());
          rec_cov_yy.push_back(lerr.yy()); 
          rec_kf_compatible.push_back(kf_compatible);
          rec_evt.push_back(eventnr);
          rec_obj_id.push_back(obj_id);
          rec_pid.push_back(pid);
          //rec_mask.push_back((itcheck->second));

        }
      }
    }
    obj_id++;
  }
  tree->Fill();

  sim_x.clear();
  sim_y.clear();
  sim_z.clear();
  sim_e.clear();
  sim_detid.clear();
  sim_layer.clear();
  sim_dtype.clear();
  sim_cov_xx.clear();
  sim_cov_xy.clear();
  sim_cov_yy.clear();
  sim_evt.clear();
  sim_kf_compatible.clear();
  sim_obj_id.clear();
  sim_pid.clear();
  sim_mask.clear();

  // RecHit

  rec_x.clear();
  rec_y.clear();
  rec_z.clear();
  rec_e.clear();
  rec_detid.clear();
  rec_layer.clear();
  rec_dtype.clear();
  rec_cov_xx.clear();
  rec_cov_xy.clear();
  rec_cov_yy.clear();
  rec_evt.clear();
  rec_kf_compatible.clear();
  rec_obj_id.clear();
  rec_pid.clear();
  rec_mask.clear();

  // KF

  kf_x.clear();
  kf_y.clear();
  kf_z.clear();
  kf_e.clear();
  kf_detid.clear();
  kf_layer.clear();
  kf_dtype.clear();
  kf_cov_xx.clear();
  kf_cov_xy.clear();
  kf_cov_yy.clear();
  kf_evt.clear();
  kf_eta.clear();
  kf_theta.clear();
  track_id.clear();
  track_charge.clear();
  track_momentum.clear();
  track_quality.clear();
  track_chi2.clear();
  track_validfraction.clear();
  track_qoverp.clear();
  track_algo.clear();
  
  //clear_arrays();
  eventnr=eventnr+1;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void Ntuplizer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void Ntuplizer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Ntuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);