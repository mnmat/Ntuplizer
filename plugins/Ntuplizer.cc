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
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
      const HGCRecHitCollection& rechitsEE, 
      const HGCRecHitCollection& rechitsFH,
      const HGCRecHitCollection& rechitsBH) const;
  std::vector<int> matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_);
  LocalError calculateLocalError(DetId detid_, const HGCalDDDConstants* ddd);


  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<std::vector<ticl::Trackster>> ticlTrackstersMergeToken;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorEMToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorHADToken_;
  /*
  edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorKFToken_;
  edm::EDGetTokenT<std::vector<float>> xxKFToken_;
  edm::EDGetTokenT<std::vector<float>> xyKFToken_;
  edm::EDGetTokenT<std::vector<float>> yyKFToken_;
  edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorToken_;
  edm::EDGetTokenT<std::vector<float>> xxPropToken_;
  edm::EDGetTokenT<std::vector<float>> xyPropToken_;
  edm::EDGetTokenT<std::vector<float>> yyPropToken_;
  edm::EDGetTokenT<std::vector<int>> chargeKFToken_;
  edm::EDGetTokenT<std::vector<int>> chargePropToken_;
  edm::EDGetTokenT<std::vector<int>> detidKFToken_;
  edm::EDGetTokenT<std::vector<int>> detidPropToken_;
  */
  edm::EDGetTokenT<float> abs_failToken_;

  edm::EDGetTokenT<std::vector<KFHit>> KFHitsToken_;
  edm::EDGetTokenT<std::vector<KFHit>> PropHitsToken_;
 //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkToken_;
  //edm::EDGetTokenT<std::vector<Point3DBase<float,GlobalTag>>> propagatorTrkEMToken_;
  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  TTree *tree = new TTree("tree","tree");

  std::vector<std::string> detectors, objects, positions, hittypes;

  // variables
  
  int eventnr =0;
  std::string eta_;
  std::string energy_;
  std::string outdir_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;

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

  // Propagator

  std::vector<float> prop_x;
  std::vector<float> prop_y;
  std::vector<float> prop_z;
  std::vector<float> prop_e;
  std::vector<float> prop_cov_xx;
  std::vector<float> prop_cov_xy;
  std::vector<float> prop_cov_yy;
  std::vector<int> prop_detid;
  std::vector<int> prop_charge;
  std::vector<int> prop_layer;
  std::vector<std::string> prop_dtype;
  std::vector<int> prop_evt;

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
  std::vector<int> rec_prop_compatible;

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
  std::vector<int> sim_prop_compatible;

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
      ticlTrackstersMergeToken(consumes<std::vector<ticl::Trackster> >(iConfig.getParameter<edm::InputTag>("Tracksters"))),
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      abs_failToken_(consumes<float>(iConfig.getParameter<edm::InputTag>("abs_fail"))),
      KFHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("KFHits"))),
      PropHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("PropHits"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      eta_(iConfig.getParameter<std::string>("eta")),
      energy_(iConfig.getParameter<std::string>("energy")),
      outdir_(iConfig.getParameter<std::string>("outdir")){



  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  hittypes = {"Simhits","Rechits","Propagator","KF"};
  objects = {"Simhits", "Rechits"};
  positions = {"Propagator", "KF"};
        
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
  tree->Branch("sim_prop_compatible", &sim_prop_compatible);
  rec_kf_compatible.clear();
  rec_prop_compatible.clear();

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
  tree->Branch("rec_prop_compatible", &rec_prop_compatible);

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
  // Propagator

  tree->Branch("prop_x", &prop_x);
  tree->Branch("prop_y", &prop_y);
  tree->Branch("prop_z", &prop_z);
  tree->Branch("prop_e", &prop_e);
  tree->Branch("prop_layer", &prop_layer);
  tree->Branch("prop_detid", &prop_detid);
  tree->Branch("prop_dtype", &prop_dtype);
  tree->Branch("prop_cov_xx", &prop_cov_xx);
  tree->Branch("prop_cov_xy", &prop_cov_xy);
  tree->Branch("prop_cov_yy", &prop_cov_yy);
  tree->Branch("prop_evt", &prop_evt);

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
  sim_prop_compatible.clear();

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
  rec_prop_compatible.clear();

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

  // Prop

  prop_x.clear();
  prop_y.clear();
  prop_z.clear();
  prop_e.clear();
  prop_detid.clear();
  prop_layer.clear();
  prop_dtype.clear();
  prop_cov_xx.clear();
  prop_cov_xy.clear();
  prop_cov_yy.clear();
  prop_evt.clear();  
}

void Ntuplizer::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), &hit);
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

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

  edm::Handle<std::vector<ticl::Trackster>> ticlTrackstersMerge;
  iEvent.getByToken(ticlTrackstersMergeToken, ticlTrackstersMerge);
  const std::vector<ticl::Trackster>& tracksters = *ticlTrackstersMerge; 

  edm::Handle<float> abs_failHandle_;
  iEvent.getByToken(abs_failToken_,abs_failHandle_);
  const float &abs_fail = *abs_failHandle_;

  edm::Handle<std::vector<KFHit>> KFHitsHandle;
  iEvent.getByToken(KFHitsToken_, KFHitsHandle);
  const std::vector<KFHit> &kfhits = *KFHitsHandle;

  edm::Handle<std::vector<KFHit>> PropHitsHandle;
  iEvent.getByToken(PropHitsToken_, PropHitsHandle);
  const std::vector<KFHit> &prophits = *PropHitsHandle;


  std::map<float, GlobalPoint> map_gps_prop, map_gps_kf;
  std::map<float, float> map_xx_kf, map_xy_kf, map_yy_kf, map_xx_prop, map_xy_prop, map_yy_prop;

  // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

  // Test KFHits

  std::cout << "Test KFHits" << std::endl;
  for (const auto& hit: kfhits){
    std::cout << hit.center.x() << "\t" << hit.center.y() << "\t" << hit.center.z() << "\t" << hit.xx << "\t" << hit.xy << "\t" << hit.yy << "\t" << hit.charge << "\t" << hit.detid << "\t" <<std::endl;
  }
  std::cout << "------------------------------------------" << std::endl;

  // KF & Propagator 
  for (const auto& pos: positions){

    auto &hits = (pos=="KF")? kfhits:prophits;
    auto &map_gps = (pos=="KF")? map_gps_kf:map_gps_prop;

    auto &vec_x = (pos=="KF")? kf_x:prop_x;
    auto &vec_y = (pos=="KF")? kf_y:prop_y;
    auto &vec_z = (pos=="KF")? kf_z:prop_z;
    auto &vec_detid = (pos=="KF")? kf_detid:prop_detid;
    auto &vec_layer = (pos=="KF")? kf_layer:prop_layer;
    auto &vec_dtype = (pos=="KF")? kf_dtype:prop_dtype;
    auto &vec_cov_xx = (pos=="KF")? kf_cov_xx:prop_cov_xx;
    auto &vec_cov_xy = (pos=="KF")? kf_cov_xy:prop_cov_xy;
    auto &vec_cov_yy = (pos=="KF")? kf_cov_yy:prop_cov_yy;
    auto &vec_evt = (pos=="KF")? kf_evt:prop_evt;
    auto &vec_e = (pos=="KF")? kf_e:prop_e;

    for(int i = 0;i<int(hits.size());i++){

      std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(hits[i].detid);
      float e = 0;
      if (itcheck != hitMap.end()){
        e = hitMap[hits[i].detid]->energy();
      }
      unsigned int layer_ = recHitTools_.getLayerWithOffset(hits[i].detid);

      std::string detector;
      std::string thickness;
      std::string tmp;
      if(recHitTools_.isSilicon(hits[i].detid)){
        detector = "Si";
        thickness = std::to_string(int(recHitTools_.getSiThickness(hits[i].detid))); 
        tmp = detector+" "+thickness;
      }
      else{
        detector = "Sc";
        thickness = "None";
        tmp = "Sc";
      } 

      map_gps[hits[i].center.z()]=hits[i].center;
      vec_x.push_back(hits[i].center.x());
      vec_y.push_back(hits[i].center.y());
      vec_z.push_back(hits[i].center.z());
      vec_e.push_back(e);
      vec_detid.push_back(hits[i].detid);
      vec_layer.push_back(layer_);
      vec_dtype.push_back(tmp);
      vec_cov_xx.push_back(hits[i].xx);
      vec_cov_xy.push_back(hits[i].xy);
      vec_cov_yy.push_back(hits[i].yy);
      vec_evt.push_back(eventnr);
    }
  }

  // Loop over Caloparticles 

  std::vector<DetId> tmprechits_; tmprechits_.clear();

  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){

        DetId detid_ = (it_sc_hae.first);
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
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

        std::cout << "Layer: " << layer_ << "\t XX: " << lerr.xx() << std::endl;

        // KF & Propagator 
        int kf_compatible=0;
        int prop_compatible=0;
        for (const auto& pos: positions){
          auto &map_gps = (pos=="KF")? map_gps_kf:map_gps_prop;
          auto &compatible = (pos=="KF")? kf_compatible:prop_compatible;

          auto gp = map_gps[recHitTools_.getPosition(detid_).z()];
          // KF/Propagator compatible?

          DetId closest_detid;
          if (detector == "Sc") closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCell(gp);
          else closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCellHex(gp, true);
          if(detid_==closest_detid){
            compatible=1;
          }
        }
        
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
        sim_prop_compatible.push_back(prop_compatible);
        sim_evt.push_back(eventnr);

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
          rec_prop_compatible.push_back(prop_compatible);
          rec_evt.push_back(eventnr);
        }
      }
    }
  }

  std::cout << "Hits:" << cps.size() << std::endl;
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
  sim_prop_compatible.clear();

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
  rec_prop_compatible.clear();

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

  // Prop

  prop_x.clear();
  prop_y.clear();
  prop_z.clear();
  prop_e.clear();
  prop_detid.clear();
  prop_layer.clear();
  prop_dtype.clear();
  prop_cov_xx.clear();
  prop_cov_xy.clear();
  prop_cov_yy.clear();
  prop_evt.clear();
  
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