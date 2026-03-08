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

#include "Validation/HGCalValidation/interface/CaloParticleSelector.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"

#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/Common/interface/TrackingParticleSelector.h"
#include "CommonTools/RecoAlgos/interface/RecoTrackSelectorBase.h"

#include <Eigen/Dense>

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
  bool inSensorCell(DetId detid_, GlobalPoint point);
  bool isPureHit(DetId detid_, GlobalPoint point,AlgebraicSymMatrix22 cov);
  bool isPureHit_approx(DetId detid_, GlobalPoint point,AlgebraicSymMatrix22 cov);

  void clear_arrays();
  virtual void fillHitMap(std::map<DetId, std::pair<const HGCRecHit*,float>>& hitMap, 
      const HGCRecHitCollection& rechitsEE, 
      const HGCRecHitCollection& rechitsFH,
      const HGCRecHitCollection& rechitsBH,
      const std::vector<float>& lcmask) const;
  std::vector<int> matchRecHit2CPRecHits(DetId detid_, std::vector<DetId> rechitdetid_);
  float calculateAreaHex(DetId id,  const HGCalDDDConstants* ddd) const;
  float calculateAreaScint(DetId id,  const HGCalDDDConstants* ddd) const;
  LocalError calculateLocalError(DetId detid_, const HGCalDDDConstants* ddd);
  
  std::pair<float, float> get_eigenvalues(AlgebraicSymMatrix22 cov);
  float get_theta(AlgebraicSymMatrix22 cov);
  AlgebraicVector2 transform_coordinates(AlgebraicVector2 points, float theta, AlgebraicVector2 center);
  std::pair<float, float> parametrize_line(AlgebraicVector2 p1, AlgebraicVector2 p2);
  bool check_intersection_ellipse_line(std::pair<float,float> eigenvalues, AlgebraicVector2 p1, AlgebraicVector2 p2);
  bool check_intersection_ellipse_point(float axis, float p);
  bool check_intersection(std::pair<float,float> eigenvalues, AlgebraicVector2 p1, AlgebraicVector2 p2);
  std::vector<DetId> getNeighbors(DetId detid_,HGCalTopology& topoEE, HGCalTopology& topoFH, HGCalTopology& topoBH) const;
  std::pair<int,double> getNeighborHitsAndEnergies(std::vector<DetId> neighbors, std::map<DetId, std::pair<const HGCRecHit*, float>>& hitMap) const;
  bool inLC(DetId detid_, std::vector<DetId>& lcdetids) const;
  
  hgcal::RecHitTools recHitTools_;

  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;
  edm::EDGetTokenT<std::vector<KFHit>> KFHitsToken_;
  edm::EDGetTokenT<std::vector<KFHit>> PropHitsToken_;

  edm::EDGetTokenT<reco::CaloClusterCollection> hgcalLayerClustersToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;
  edm::EDGetTokenT<edm::View<reco::Track>> tracksToken_;
  edm::EDGetTokenT<std::vector<float>> lcMaskToken_;
  std::vector<edm::InputTag> associators;
  std::vector<edm::EDGetTokenT<reco::SimToRecoCollection>> associatormapStRs;
  std::vector<edm::EDGetTokenT<reco::RecoToSimCollection>> associatormapRtSs;
  edm::EDGetTokenT<std::vector<SimVertex>> simVerticesToken_;

  std::vector<std::string> detectors, objects, positions, hittypes;


  // variables
  
  // int eventnr =0;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;

  TTree* tree;

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
  std::vector<float> kf_track_id;
  std::vector<float> kf_track_charge;
  std::vector<float> kf_track_momentum;
  std::vector<float> kf_track_quality;
  std::vector<float> kf_track_chi2;
  std::vector<float> kf_track_validfraction;
  std::vector<float> kf_track_qoverp;
  std::vector<float> kf_track_algo;
  std::vector<float> kf_meas_x;
  std::vector<float> kf_meas_y;
  std::vector<float> kf_meas_z;
  std::vector<int> kf_rec_compatible;
  std::vector<int> kf_rec_isclosest;
  std::vector<int> kf_rec_inSensor;
  std::vector<int> kf_rec_pureHit;
  std::vector<int> kf_rec_pureHit_approx;
  std::vector<int> kf_neighbors;
  std::vector<int> kf_inLC;
  std::vector<int> kf_total_neighbors;
  std::vector<float> kf_neighbors_energy;

  // Prop

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
  std::vector<float> prop_eta;
  std::vector<float> prop_theta;
  std::vector<float> prop_track_id;
  std::vector<float> prop_track_charge;
  std::vector<float> prop_track_momentum;
  std::vector<float> prop_track_quality;
  std::vector<float> prop_track_chi2;
  std::vector<float> prop_track_validfraction;
  std::vector<float> prop_track_qoverp;
  std::vector<float> prop_track_algo;
  std::vector<float> prop_meas_x;
  std::vector<float> prop_meas_y;
  std::vector<float> prop_meas_z;
  std::vector<int> prop_rec_compatible;
  std::vector<int> prop_rec_isclosest;
  std::vector<int> prop_rec_inSensor;
  std::vector<int> prop_rec_pureHit;
  std::vector<int> prop_rec_pureHit_approx;
  std::vector<int> prop_neighbors;
  std::vector<int> prop_inLC;
  std::vector<int> prop_total_neighbors;
  std::vector<float> prop_neighbors_energy;

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
  std::vector<int> rec_kf_contained;
  std::vector<int> rec_kf_inSensor;
  std::vector<int> rec_kf_pureHit;
  std::vector<int> rec_kf_pureHit_approx;
  std::vector<int> rec_obj_id;
  std::vector<int> rec_pid;
  std::vector<int> rec_simcluster_id;
  std::vector<int> rec_neighbors;
  std::vector<int> rec_inLC;
  std::vector<int> rec_total_neighbors;
  std::vector<int> rec_neighbors_energy;
  //std::vector<float> rec_mask;


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
  std::vector<int> sim_kf_contained;
  std::vector<int> sim_kf_inSensor;
  std::vector<int> sim_kf_pureHit;
  std::vector<int> sim_kf_pureHit_approx;
  std::vector<int> sim_obj_id;
  std::vector<int> sim_pid;
  std::vector<int> sim_simcluster_id;
  //std::vector<float> sim_mask;


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
      PropHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("PropHits"))),
      hgcalLayerClustersToken_(consumes<reco::CaloClusterCollection>(iConfig.getParameter<edm::InputTag>("hgcalLayerClusters"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      tracksToken_(consumes<edm::View<reco::Track>>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      lcMaskToken_(consumes<std::vector<float>>(iConfig.getParameter<edm::InputTag>("lcMask"))),
      associators(iConfig.getUntrackedParameter<std::vector<edm::InputTag>>("associators")),
      simVerticesToken_(consumes<std::vector<SimVertex>>(iConfig.getParameter<edm::InputTag>("simVertices"))){


  associatormapRtSs.push_back(consumes<reco::RecoToSimCollection>(associators[0]));

  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  hittypes = {"Simhits","Rechits","KF"};
  objects = {"Simhits", "Rechits"};        
  positions = {"KF","Prop"};
  //recHitTools_.reset(new hgcal::RecHitTools());
  //now do what ever initialization is needed

  usesResource("TFileService");
  edm::Service<TFileService> file;
  tree = file->make<TTree>("tree","tree"); // Replace this with an input tree value?

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
  tree->Branch("sim_kf_contained", &sim_kf_contained);
  tree->Branch("sim_kf_inSensor", &sim_kf_inSensor);
  tree->Branch("sim_kf_pureHit", &sim_kf_pureHit);
  tree->Branch("sim_kf_pureHit_approx", &sim_kf_pureHit_approx);
  tree->Branch("sim_obj_id", &sim_obj_id);
  tree->Branch("sim_pid", &sim_pid);
  tree->Branch("sim_simcluster_id", &sim_simcluster_id);
  //tree->Branch("sim_mask", &sim_mask);

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
  tree->Branch("rec_kf_contained", &rec_kf_contained);
  tree->Branch("rec_kf_inSensor", &rec_kf_inSensor);
  tree->Branch("rec_kf_pureHit", &rec_kf_pureHit);
  tree->Branch("rec_kf_pureHit_approx", &rec_kf_pureHit_approx);
  tree->Branch("rec_obj_id", &rec_obj_id);
  tree->Branch("rec_pid", &rec_pid);
  tree->Branch("rec_simcluster_id", &rec_simcluster_id);
  tree->Branch("rec_neighbors", &rec_neighbors);
  tree->Branch("rec_inLC", &rec_inLC);
  tree->Branch("rec_total_neighbors", &rec_total_neighbors);
  tree->Branch("rec_neighbors_energy", &rec_neighbors_energy);
  //tree->Branch("rec_mask", &rec_mask);

  // KF

  tree->Branch("kf_x", &kf_x);
  tree->Branch("kf_y", &kf_y);
  tree->Branch("kf_z", &kf_z);
  tree->Branch("kf_meas_x", &kf_meas_x);
  tree->Branch("kf_meas_y", &kf_meas_y);
  tree->Branch("kf_meas_z", &kf_meas_z);
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
  tree->Branch("kf_track_id", &kf_track_id);
  tree->Branch("kf_track_charge", &kf_track_charge);
  tree->Branch("kf_track_momentum", &kf_track_momentum);
  tree->Branch("kf_track_quality", &kf_track_quality);
  tree->Branch("kf_track_chi2", &kf_track_chi2);
  tree->Branch("kf_track_validfraction", &kf_track_validfraction);
  tree->Branch("kf_track_qoverp", &kf_track_qoverp);
  tree->Branch("kf_track_algo", &kf_track_algo);
  tree->Branch("kf_rec_compatible", &kf_rec_compatible);
  tree->Branch("kf_rec_isclosest", &kf_rec_isclosest);
  tree->Branch("kf_rec_inSensor", &kf_rec_inSensor);
  tree->Branch("kf_rec_pureHit", &kf_rec_pureHit);
  tree->Branch("kf_rec_pureHit_approx", &kf_rec_pureHit_approx);
  tree->Branch("kf_neighbors", &kf_neighbors);
  tree->Branch("kf_inLC", &kf_inLC);
  tree->Branch("kf_total_neighbors", &kf_total_neighbors);
  tree->Branch("kf_neighbors_energy", &kf_neighbors_energy);

    // Prop

  tree->Branch("prop_x", &prop_x);
  tree->Branch("prop_y", &prop_y);
  tree->Branch("prop_z", &prop_z);
  tree->Branch("prop_meas_x", &prop_meas_x);
  tree->Branch("prop_meas_y", &prop_meas_y);
  tree->Branch("prop_meas_z", &prop_meas_z);
  tree->Branch("prop_e", &prop_e);
  tree->Branch("prop_layer", &prop_layer);
  tree->Branch("prop_detid", &prop_detid);
  tree->Branch("prop_dtype", &prop_dtype);
  tree->Branch("prop_cov_xx", &prop_cov_xx);
  tree->Branch("prop_cov_xy", &prop_cov_xy);
  tree->Branch("prop_cov_yy", &prop_cov_yy);
  tree->Branch("prop_evt", &prop_evt);
  tree->Branch("prop_eta", &prop_eta);
  tree->Branch("prop_theta", &prop_theta);
  tree->Branch("prop_track_id", &prop_track_id);
  tree->Branch("prop_track_charge", &prop_track_charge);
  tree->Branch("prop_track_momentum", &prop_track_momentum);
  tree->Branch("prop_track_quality", &prop_track_quality);
  tree->Branch("prop_track_chi2", &prop_track_chi2);
  tree->Branch("prop_track_validfraction", &prop_track_validfraction);
  tree->Branch("prop_track_qoverp", &prop_track_qoverp);
  tree->Branch("prop_track_algo", &prop_track_algo);
  tree->Branch("prop_rec_compatible", &prop_rec_compatible);
  tree->Branch("prop_rec_isclosest", &prop_rec_isclosest);
  tree->Branch("prop_rec_inSensor", &prop_rec_inSensor);
  tree->Branch("prop_rec_pureHit", &prop_rec_pureHit);
  tree->Branch("prop_rec_pureHit_approx", &prop_rec_pureHit_approx);
  tree->Branch("prop_neighbors", &prop_neighbors);
  tree->Branch("prop_inLC", &prop_inLC);
  tree->Branch("prop_total_neighbors", &prop_total_neighbors);
  tree->Branch("prop_neighbors_energy", &prop_neighbors_energy);

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
}

//
// member functions
//

// ------------ method called for each event  ------------

float Ntuplizer::calculateAreaHex(DetId id,  const HGCalDDDConstants* ddd) const {
  if(recHitTools_.getSiThickness(id) < 200) return 1.18; // TODO: replace with non-hardcoded value; hardcoded value from TDR
  else  return 0.52; // TODO: replace with non-hardcoded value; hardcoded value from TDR
}

float Ntuplizer::calculateAreaScint(DetId id, const HGCalDDDConstants* ddd) const{
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
    return A;
}

LocalError Ntuplizer::calculateLocalError(DetId id, const HGCalDDDConstants* ddd){
  if(recHitTools_.isSilicon(id)){
    float A = calculateAreaHex(id,ddd);
    float a = sqrt(2*A/(3*sqrt(3)));
    double varx = pow(a,4)*5*sqrt(3)/(16*A); // x
    double vary = pow(a,4)*5*sqrt(3)/(16*A); // y 
    return LocalError(varx, 0, vary);
  }
  else{
    float A = calculateAreaScint(id, ddd);

    // Exact solution but error prone due to numerical limitations
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

std::pair<float, float> Ntuplizer::get_eigenvalues(AlgebraicSymMatrix22 cov){
  Eigen::Matrix2f ell;
  std::pair<float, float> res;
  ell << cov[0][0], cov[0][1], cov[1][0], cov[1][1];
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2f> eigensolver(ell);
  if (eigensolver.info() != Eigen::Success){
    res.first = cov[0][0];
    res.second =  cov[1][1];
  } else {
    auto eigenvalues = eigensolver.eigenvalues();
    res.first = eigenvalues[0];
    res.second = eigenvalues[1];
  }
  return res;
}

float Ntuplizer::get_theta(AlgebraicSymMatrix22 cov){
  float a = cov[0][0];
  float b = cov[1][0];
  float c = cov[1][1];
  return atan(2*b/(c-a))/2;
}

AlgebraicVector2 Ntuplizer::transform_coordinates(AlgebraicVector2 point, float theta, AlgebraicVector2 center){
  AlgebraicMatrix22 r;
  r[0][0] = cos(-theta);
  r[0][1] = -sin(-theta);
  r[1][0] = sin(theta);
  r[1][1] = cos(theta);
  auto shifted = point - center;
  auto rotated = r*point;
  return shifted;
}

std::pair<float, float> Ntuplizer::parametrize_line(AlgebraicVector2 p1, AlgebraicVector2 p2){
  std::pair<float, float> ret;
  float x1 = p1[0];
  float y1 = p1[1];
  float x2 = p2[0];
  float y2 = p2[1];

  float k = (y2 - y1)/(x2-x1);
  float d = (x2*y1-x1*y2)/(x2-x1);
  ret.first = k;
  ret.second = d;
  return ret;
}

bool Ntuplizer::check_intersection_ellipse_line(std::pair<float,float> eigenvalues, AlgebraicVector2 p1, AlgebraicVector2 p2){
  std::pair<float, float> line = parametrize_line(p1,p2);
  float mbd = std::sqrt(5.991);
  float major = eigenvalues.first*mbd;
  float minor = eigenvalues.second*mbd;
  float major2 = std::pow(major,2);
  float minor2 = std::pow(minor,2);
  float k = line.first;
  float d = line.second;
  float A = major2*std::pow(k,2)+minor2;
  float B = 2*k*d*major2;
  float C = std::pow(d,2)*major2 - major2*minor2;
  if (std::pow(B,2) < 4*A*C){
    return false;
  } else {
    return true;
  }

}

bool Ntuplizer::check_intersection_ellipse_point(float axis, float p){
  if (axis < abs(p)){
    return false;
  } else {
    return true;
  }
}

bool Ntuplizer::check_intersection(std::pair<float,float> eigenvalues, AlgebraicVector2 p1, AlgebraicVector2 p2){
  if (p1[0] == p2[0]){
    return check_intersection_ellipse_point(eigenvalues.first,p1[0]);
  } else {
    return check_intersection_ellipse_line(eigenvalues,p1,p2);
  }
}

bool Ntuplizer::isPureHit(DetId detid_, GlobalPoint point,AlgebraicSymMatrix22 cov){
  // Get Coordinates
  auto x = point.x();
  auto y = point.y();
  AlgebraicVector2 center = AlgebraicVector2(x,y);
  
  // GetCorners
  // std::cout << "Get Corners" << std::endl;
  auto hgcalgeometry = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_));
  auto corners = hgcalgeometry->getCorners(detid_);

  // Calculate eigenvalues
  // std::cout << "Get Eigenvalues" << std::endl;
  std::pair<float, float> eigenvalues = get_eigenvalues(cov);

  // Calculate Rotation Matrix
  
  // std::cout << "Get Theta" << std::endl;
  float theta = get_theta(cov);

  // Coordinate transform of corners
  // std::cout << "Coordinate Transform" << std::endl;
  std::vector<AlgebraicVector2> transformed_corners;
  for (auto& point: corners){
    AlgebraicVector2 p = AlgebraicVector2(point.x(),point.y());
    transformed_corners.push_back(transform_coordinates(p,theta,center));
  }

  // Check for boundary crossing
  // std::cout << "check boundary crossing" << std::endl;
  int cornerSize = 8;
  bool check_intersect = false;
  for (int i = 1; i < cornerSize; i++){
    check_intersect = check_intersection(eigenvalues,transformed_corners[i], transformed_corners[i+1]);
    if (check_intersect) return false;
  }
  return true;
}

bool Ntuplizer::isPureHit_approx(DetId detid_, GlobalPoint point,AlgebraicSymMatrix22 cov){
  
  // Calculate eigenvalues
  // std::cout << "isPureHit_approx: Get Eigenvalues" << std::endl;
  std::pair<float, float> eigenvalues = get_eigenvalues(cov);

  // Calculate Area of Confidence ellipse
  float mbd = std::sqrt(5.991);
  float area_tsos = M_PI*eigenvalues.first*eigenvalues.second*mbd*mbd;

  // Calculate fraction
  // std::cout << "isPureHit_approx: Calculate Fraction" << std::endl;
  auto hgcalgeometry = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_));
  // std::cout << "Get Area" << std::endl;
  float area_sensor = hgcalgeometry->getArea(detid_);
  // std::cout << "Successfully gotten area" << std::endl;
  float frac = area_tsos/area_sensor;

  if (frac>1){
    return false;
  }

  // Check if tsos falls within scaled circle
  // std::cout << "isPureHit_approx: Perform Check" << std::endl;

  auto center = recHitTools_.getPosition(detid_);
  auto diff = point - center;
  float dist = std::sqrt(diff.x()*diff.x() - diff.y()*diff.y());
  float radius = std::sqrt(area_sensor/M_PI*(1-frac));
  if (radius > dist) {
    return true;
  } else {
    return false;
  }
}

bool Ntuplizer::inSensorCell(DetId detid_, GlobalPoint point){
  // inSensorCell checks both hexagonal cells and the trapezoid based on the corners stored in the hgcalGeometry.
  // The corners are stored in a counter clockwise fashion.
  // The hexagonal boundaries are checked using the winding number algorithm.
  // The trapezoid boundaries are checked using simple polar coordinates.

  // Get Coordinates
  auto x = point.x();
  auto y = point.y();
  
  // GetCorners
  auto hgcalgeometry = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_));
  // std::cout << "DetId present in hgcalgeometry? " << hgcalgeometry->present(detid_) << std::endl;
  auto corners = hgcalgeometry->getCorners(detid_);

  /*
  // std::cout << "v1=["<<corners[0].x()<<","<<corners[0].y()<<","<<corners[0].z()<<"]" << std::endl;
  // std::cout << "v2=["<<corners[1].x()<<","<<corners[1].y()<<","<<corners[1].z()<<"]" << std::endl;
  // std::cout << "v3=["<<corners[2].x()<<","<<corners[2].y()<<","<<corners[2].z()<<"]" << std::endl;
  // std::cout << "v4=["<<corners[3].x()<<","<<corners[3].y()<<","<<corners[3].z()<<"]" << std::endl;
  // std::cout << "v5=["<<corners[4].x()<<","<<corners[4].y()<<","<<corners[4].z()<<"]" << std::endl;
  // std::cout << "v6=["<<corners[5].x()<<","<<corners[5].y()<<","<<corners[5].z()<<"]" << std::endl;
  // std::cout << "v7=["<<corners[6].x()<<","<<corners[6].y()<<","<<corners[6].z()<<"]" << std::endl;
  // std::cout << "p=["<<point.x()<<","<<point.y()<<","<<point.z()<<"]" << std::endl;
  */
  bool inside = false;
  if (recHitTools_.isSilicon(detid_)){
    // Check Hexagon
    int cornersSize = 7;
    int j=0;
    for (int i = 1; i < cornersSize; i++) {
      if (((corners[i].y() >+ point.y()) != (corners[j].y() >= point.y())) &&
          (x <= (corners[j].x() - corners[i].x()) * (y - corners[i].y()) / (corners[j].y() - corners[i].y()) + corners[i].x())) {
          inside = !inside;
      }
      j++;
    }
  }
  else {
    // Check Scintillator
    float phi_min = std::atan2(corners[0].y(),corners[0].x()); //radians
    float phi_max = std::atan2(corners[1].y(),corners[1].x()); //radians
    float r_min = std::sqrt(corners[3].x()*corners[3].x() + corners[3].y()*corners[3].y());
    float r_max = std::sqrt(corners[0].x()*corners[0].x() + corners[0].y()*corners[0].y());;

    auto r = std::sqrt(point.x()*point.x()+point.y()*point.y());
    if (point.phi()>=phi_min && point.phi()<=phi_max){
      if (r>=r_min && r<=r_max){
        inside = !inside;
      }
    }
  }
  return inside;
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
  sim_kf_contained.clear();
  sim_kf_inSensor.clear();
  sim_kf_pureHit.clear();
  sim_obj_id.clear();
  sim_pid.clear();
  sim_simcluster_id.clear();
  //sim_mask.clear();

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
  rec_kf_contained.clear();
  rec_kf_inSensor.clear();
  rec_kf_pureHit.clear();
  rec_obj_id.clear();
  rec_pid.clear();
  rec_simcluster_id.clear();
  rec_neighbors.clear();
  rec_inLC.clear();
  rec_total_neighbors.clear();
  rec_neighbors_energy.clear();
  //rec_mask.clear();

  // KF

  kf_x.clear();
  kf_y.clear();
  kf_z.clear();
  kf_meas_x.clear();
  kf_meas_y.clear();
  kf_meas_z.clear();
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
  kf_track_id.clear();
  kf_track_charge.clear();
  kf_track_momentum.clear();
  kf_track_quality.clear();
  kf_track_chi2.clear();
  kf_track_validfraction.clear();
  kf_track_qoverp.clear();
  kf_track_algo.clear();
  kf_rec_compatible.clear();
  kf_rec_isclosest.clear();
  kf_rec_inSensor.clear();
  kf_rec_pureHit.clear();
  kf_rec_pureHit_approx.clear();
  kf_neighbors.clear();
  kf_inLC.clear();
  kf_total_neighbors.clear();
  kf_neighbors_energy.clear();
  // Prop

  prop_x.clear();
  prop_y.clear();
  prop_z.clear();
  prop_meas_x.clear();
  prop_meas_y.clear();
  prop_meas_z.clear();
  prop_e.clear();
  prop_detid.clear();
  prop_layer.clear();
  prop_dtype.clear();
  prop_cov_xx.clear();
  prop_cov_xy.clear();
  prop_cov_yy.clear();
  prop_evt.clear();
  prop_eta.clear();
  prop_theta.clear();
  prop_track_id.clear();
  prop_track_charge.clear();
  prop_track_momentum.clear();
  prop_track_quality.clear();
  prop_track_chi2.clear();
  prop_track_validfraction.clear();
  prop_track_qoverp.clear();
  prop_track_algo.clear();
  prop_rec_compatible.clear();
  prop_rec_isclosest.clear();
  prop_rec_inSensor.clear();
  prop_rec_pureHit.clear();
  prop_rec_pureHit_approx.clear();
  prop_neighbors.clear();
  prop_inLC.clear();
  prop_total_neighbors.clear();
  prop_neighbors_energy.clear();
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


std::vector<DetId> Ntuplizer::getNeighbors(DetId detid_,HGCalTopology& topoEE,HGCalTopology& topoHF,HGCalTopology& topoHB) const {
  // Get neighbors of detid
  std::vector<DetId> neighbors;
  if (detid_.det() == 8){
    neighbors = topoEE.neighbors(detid_);
  } else if (detid_.det() == 9){
    neighbors =  topoHF.neighbors(detid_);
  } else if (detid_.det() == 10){
    neighbors =  topoHB.neighbors(detid_);
  }
  return neighbors;
}

std::pair<int,double> Ntuplizer::getNeighborHitsAndEnergies(std::vector<DetId> neighbors, std::map<DetId, std::pair<const HGCRecHit*, float>>& hitMap) const {
  std::pair<int,float> result;

  // Check if neighboring detids are associated with a RecHit
  int hitNeighbors = 0;
  float totalEnergy = 0;

  for (auto& neighbor: neighbors){
    auto candidate = hitMap.find(neighbor);
    if (candidate != hitMap.end()){
      hitNeighbors++;
      totalEnergy += candidate->second.first->energy();
    }
  }
  result.first = hitNeighbors;
  result.second = totalEnergy;
  return result;
}


bool Ntuplizer::inLC(DetId detid_, std::vector<DetId>& lcs_detids) const {
  auto hit_check = std::find(lcs_detids.begin(),lcs_detids.end(),detid_);
  if (hit_check != lcs_detids.end()){
    return true;
  } else {
    return false;
  }
}

void Ntuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //// std::cout << "Analyze" << std::endl;

  auto const& eventnr = iEvent.id().event();

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

  edm::Handle<std::vector<KFHit>> PropHitsHandle;
  iEvent.getByToken(PropHitsToken_, PropHitsHandle);
  const std::vector<KFHit> &prophits = *PropHitsHandle;

  edm::Handle<edm::View<reco::Track>> tracks_h;
  iEvent.getByToken(tracksToken_,tracks_h);
  const edm::View<reco::Track> & tkx = *(tracks_h.product()); 

  edm::Handle<reco::CaloClusterCollection> layerClusterHandle;
  iEvent.getByToken(hgcalLayerClustersToken_, layerClusterHandle);
  const reco::CaloClusterCollection &lcs = *layerClusterHandle;

  // Create LC map of RecHits
  std::vector<DetId> lcs_detids;
  for (const auto& it_lc : lcs) {
    const std::vector<std::pair<DetId, float>> &hf = it_lc.hitsAndFractions();
    // loop over the rechits of this specific layer cluster
    for (unsigned int j = 0; j < hf.size(); j++) { 
      DetId detid_ = hf[j].first;  
      lcs_detids.push_back(detid_);
    }
  }

  // Match tracks to simtrack
  
  reco::RecoToSimCollection const* recSimCollP = nullptr;
  edm::Handle<reco::RecoToSimCollection> recotosimCollectionH;
  iEvent.getByToken(associatormapRtSs[0], recotosimCollectionH);
  recSimCollP = recotosimCollectionH.product();
  reco::RecoToSimCollection const& recSimColl = *recSimCollP;

  edm::Handle<std::vector<SimVertex>> simVerticesHandle;
  iEvent.getByToken(simVerticesToken_, simVerticesHandle);
  std::vector<SimVertex> const& simVertices = *simVerticesHandle;

  TrackingParticleSelector tpSelector_ = TrackingParticleSelector(0, 100, 1.5, 3,120,280,0,true,false,false,false,{13});
  std::vector<int> signalIdx;




  

  for (edm::View<reco::Track>::size_type i = 0; i < tkx.size(); ++i) {
    RefToBase<reco::Track> track(tracks_h, i);  
        
    if(recSimColl.find(track) == recSimColl.end()) continue;

    std::vector<std::pair<TrackingParticleRef, double> > tps;
    tps = recSimColl[track];
    for (auto tp:tps){
      if(!tpSelector_(*(tp.first))) continue;
      signalIdx.push_back(i); 
    }
  }

  if (signalIdx.size()==0){
    // std::cout << "No Signal found!!!!!!" << std::endl;
    return;
  }

  std::map<int, GlobalPoint> map_gps_kf, map_gps_prop;
  std::map<int, AlgebraicSymMatrix22> map_cov_kf, map_cov_prop;
  std::map<int, std::vector<int>> map_detid_kf,map_detid_prop;
  std::map<float, float> map_xx_kf, map_xy_kf, map_yy_kf, map_xx_prop, map_xy_prop, map_yy_prop;

    // init vars
  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

  const CaloSubdetectorGeometry *geomHB = geom.getSubdetectorGeometry(DetId::Detector(10), ForwardSubdetector::ForwardEmpty);
  auto subgeomHB = static_cast<const HGCalGeometry*>(geomHB);
  auto topoHB = subgeomHB->topology();
  const HGCalDDDConstants* ddd = &(subgeomHB->topology().dddConstants());

  
  const CaloSubdetectorGeometry *geomHF = geom.getSubdetectorGeometry(DetId::Detector(9), ForwardSubdetector::ForwardEmpty);
  auto subgeomHF = static_cast<const HGCalGeometry*>(geomHF);
  auto topoHF = subgeomHF->topology();

  const CaloSubdetectorGeometry *geomEE = geom.getSubdetectorGeometry(DetId::Detector(8), ForwardSubdetector::ForwardEmpty);
  auto subgeomEE = static_cast<const HGCalGeometry*>(geomEE);
  auto topoEE = subgeomEE->topology();

  //auto radiusLayer = ddd->rangeRLayer(8, true);
  //// std::cout << radiusLayer.first << ", " << radiusLayer.second << std::endl;

  // Test KFHits

/*
  // std::cout << "Test KFHits" << std::endl;
  for (const auto& hit: kfhits){
    // std::cout << hit.center.x() << "\t" << hit.center.y() << "\t" << hit.center.z() << "\t" << hit.xx << "\t" << hit.xy << "\t" << hit.yy << "\t" << hit.charge << "\t" << hit.detid << "\t" <<std::endl;
    // std::cout << hit.test[0][0] << "," << hit.test[1][1] << hit.test[2][2] << hit.test[3][3] << hit.test[4][4] << hit.test[5][5] << std::endl;
  }
  // std::cout << "------------------------------------------" << std::endl;
*/
  // KF Hits

  // std::cout << "Fill KFHits" << std::endl;
  for (const auto& pos:positions){
    auto &hits = (pos=="KF")? kfhits:prophits;
    auto &map_gps = (pos=="KF")? map_gps_kf:map_gps_prop;
    auto &map_cov = (pos=="KF")? map_cov_kf:map_cov_prop;

    auto &map_detid = (pos=="KF")? map_detid_kf:map_detid_prop;


    auto &vec_x = (pos=="KF")? kf_x:prop_x;
    auto &vec_y = (pos=="KF")? kf_y:prop_y;
    auto &vec_z = (pos=="KF")? kf_z:prop_z;
    auto &vec_meas_x = (pos=="KF")? kf_meas_x:prop_meas_x;
    auto &vec_meas_y = (pos=="KF")? kf_meas_y:prop_meas_y;
    auto &vec_meas_z = (pos=="KF")? kf_meas_z:prop_meas_z;
    auto &vec_detid = (pos=="KF")? kf_detid:prop_detid;
    auto &vec_layer = (pos=="KF")? kf_layer:prop_layer;
    auto &vec_dtype = (pos=="KF")? kf_dtype:prop_dtype;
    auto &vec_cov_xx = (pos=="KF")? kf_cov_xx:prop_cov_xx;
    auto &vec_cov_xy = (pos=="KF")? kf_cov_xy:prop_cov_xy;
    auto &vec_cov_yy = (pos=="KF")? kf_cov_yy:prop_cov_yy;
    auto &vec_evt = (pos=="KF")? kf_evt:prop_evt;
    auto &vec_e = (pos=="KF")? kf_e:prop_e;
    auto &vec_eta = (pos=="KF")? kf_eta:prop_eta;
    auto &vec_theta = (pos=="KF")? kf_theta:prop_theta;
    auto &vec_track_charge = (pos=="KF")? kf_track_charge:prop_track_charge;
    auto &vec_track_momentum = (pos=="KF")? kf_track_momentum:prop_track_momentum;
    auto &vec_track_quality = (pos=="KF")? kf_track_quality:prop_track_quality;
    auto &vec_track_chi2 = (pos=="KF")? kf_track_chi2:prop_track_chi2;
    auto &vec_track_id = (pos=="KF")? kf_track_id:prop_track_id;
    auto &vec_track_validfraction = (pos=="KF")? kf_track_validfraction:prop_track_validfraction;
    auto &vec_track_qoverp = (pos=="KF")? kf_track_qoverp:prop_track_qoverp;
    auto &vec_track_algo = (pos=="KF")? kf_track_algo:prop_track_algo;
    auto &vec_rec_compatible = (pos=="KF")? kf_rec_compatible:prop_rec_compatible;
    auto &vec_rec_isclosest = (pos=="KF")? kf_rec_isclosest:prop_rec_isclosest;
    auto &vec_rec_inSensor = (pos=="KF")? kf_rec_inSensor:prop_rec_inSensor;
    auto &vec_rec_pureHit = (pos=="KF")? kf_rec_pureHit:prop_rec_pureHit;
    auto &vec_rec_pureHit_approx = (pos=="KF")? kf_rec_pureHit_approx:prop_rec_pureHit_approx;
    auto &vec_neighbors = (pos=="KF")? kf_neighbors:prop_neighbors;
    auto &vec_inLC = (pos=="KF")? kf_inLC:prop_inLC;
    auto &vec_total_neighbors = (pos=="KF")? kf_total_neighbors:prop_total_neighbors;
    auto &vec_neighbors_energy = (pos=="KF")? kf_neighbors_energy:prop_neighbors_energy;


    for(int i = 0;i<int(hits.size());i++){
      if (signalIdx[0] != hits[i].trackId) continue;

      std::map<DetId,std::pair<const HGCRecHit *,float>>::const_iterator itcheck = hitMap.find(hits[i].detid);
      float e = 0;
      if (itcheck != hitMap.end()){
        e = hitMap[hits[i].detid].first->energy();
      }
      //unsigned int layer_ = recHitTools_.getLayerWithOffset(hits[i].detid);

      std::string detector;
      std::string thickness;
      std::string tmp;
      
      // TODO: if conditions necessary to deal with missing rechits. Not yet possible to determine area of detector properly
      //       Implement it in such a way that we can find the closest detid and from this determine silicon thickness.
      auto detid_ = hits[i].detid;
      if (detid_==0){
        //auto closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(0))->getClosestCell(hits[i].center);
        //layer_ = recHitTools_.getLayerWithOffset(closest_detid); // Can get rid of this in future iterations as I read it out in KFHits anyway
        detector = "Sc";
        thickness = "None";
        tmp = "Sc";
      }
      else if((detid_==8)||(detid_==9)){
        //auto closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(hits[i].detid))->getClosestCellHex(hits[i].center,true);
        //layer_ = recHitTools_.getLayerWithOffset(closest_detid); // Can get rid of this in future iterations as I read it out in KFHits anyway
        detector = "Si";
        //thickness = std::to_string(int(recHitTools_.getSiThickness(closest_detid))); 
        tmp = detector;
      }
      else if (recHitTools_.isSilicon(detid_)){
        detector = "Si";
        thickness = std::to_string(int(recHitTools_.getSiThickness(hits[i].detid))); 
        tmp = detector+" "+thickness;
      }
      else{
        detector = "Sc";
        thickness = "None";
        tmp = "Sc";
      } 


      // Neighbors
      
      int num_neighbors = -99;
      float energy = -99;
      int totalNeighbors = -99;
      if (detid_ > 100){
        auto neighbors = getNeighbors(detid_,topoEE,topoHF,topoHB);
        auto neighbors_pair = getNeighborHitsAndEnergies(neighbors, hitMap);
        num_neighbors = neighbors_pair.first;
        energy = neighbors_pair.second;
        energy = e/energy;
        totalNeighbors = neighbors.size();
      }

      // Check if in LC
      bool inLC_ = inLC(detid_,lcs_detids);

      // KF 
      int rec_compatible=0;
      int rec_isclosest=0;
      int rec_inSensor=0;
      int rec_pureHit=0;
      int rec_pureHit_approx=0;
      if (detid_ > 10){
        rec_compatible = 1;
        int layer = recHitTools_.getLayerWithOffset(detid_);

        //
        DetId closest_detid;
        auto gp = hits[i].center;
        AlgebraicSymMatrix22 covariance_matrix;
        covariance_matrix[0][0] = hits[i].xx;
        covariance_matrix[0][1] = hits[i].xy;
        covariance_matrix[1][1] = hits[i].yy;
        [[maybe_unused]] auto cov = AlgebraicSymMatrix22(covariance_matrix);
        auto testgeom = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_));
        if (detector == "Sc"){
          closest_detid = testgeom->getClosestCell(gp);
          // std::cout << "In geometry? " << testgeom->present(closest_detid) << std::endl;
        }
        else {
          closest_detid = testgeom->getClosestCellHex(gp, false);
          // std::cout << "In geometry? " << testgeom->present(closest_detid) << std::endl;
        }
        if(detid_==closest_detid){
          rec_isclosest=1;
        }

        // std::cout << "TSOS conditions" << std::endl;
        // std::cout << "Closest DetId: " << closest_detid.rawId() << ", Det: " << closest_detid.det() << ", Subdet: " << closest_detid.subdetId() << std::endl;
        // std::cout << "DetId: " << DetId(detid_).rawId() << ", Det: " << DetId(detid_).det() << ", Subdet: " << DetId(detid_).subdetId() << std::endl;

        if (inSensorCell(detid_, gp)){
          rec_inSensor=1;
          if (isPureHit(detid_,gp,cov)){
            rec_pureHit=1;
          }
          if (isPureHit_approx(detid_,gp,cov)){
            rec_pureHit_approx=1;
          }
        }
      }

      if (hits[i].detid > 10){
        vec_meas_x.push_back(recHitTools_.getPosition(hits[i].detid).x());
        vec_meas_y.push_back(recHitTools_.getPosition(hits[i].detid).y());
        vec_meas_z.push_back(recHitTools_.getPosition(hits[i].detid).z());
      } else {
        vec_meas_x.push_back(0);
        vec_meas_y.push_back(0);
        vec_meas_z.push_back(0);
      }

      int tkId = hits[i].trackId;

      map_detid[hits[i].layer].push_back(hits[i].detid);
      map_gps[hits[i].layer]=hits[i].center;
      AlgebraicSymMatrix22 covariance_matrix;
      covariance_matrix[0][0] = hits[i].xx;
      covariance_matrix[0][1] = hits[i].xy;
      covariance_matrix[1][1] = hits[i].yy;
      map_cov[hits[i].layer]=AlgebraicSymMatrix22(covariance_matrix);
      vec_x.push_back(hits[i].center.x());
      vec_y.push_back(hits[i].center.y());
      vec_z.push_back(hits[i].center.z());
      vec_e.push_back(e);
      vec_detid.push_back(hits[i].detid);
      vec_layer.push_back(hits[i].layer); // 
      vec_dtype.push_back(tmp);
      vec_cov_xx.push_back(hits[i].xx);
      vec_cov_xy.push_back(hits[i].xy);
      vec_cov_yy.push_back(hits[i].yy);
      vec_evt.push_back(eventnr);
      vec_eta.push_back(hits[i].eta);
      vec_theta.push_back(hits[i].theta);
      vec_track_id.push_back(hits[i].trackId);
      vec_track_charge.push_back(hits[i].trackCharge);
      vec_track_momentum.push_back(hits[i].trackMomentum);
      vec_track_quality.push_back(hits[i].trackQuality);
      vec_track_chi2.push_back(hits[i].trackChi2);
      vec_track_validfraction.push_back(hits[i].trackValidFraction);
      vec_track_qoverp.push_back(hits[i].trackQOverP);
      vec_track_algo.push_back(tkx[tkId].algo());
      vec_rec_compatible.push_back(rec_compatible);
      vec_rec_isclosest.push_back(rec_isclosest);
      vec_rec_inSensor.push_back(rec_inSensor);
      vec_rec_pureHit.push_back(rec_pureHit);
      vec_rec_pureHit_approx.push_back(rec_pureHit_approx);
      vec_neighbors.push_back(num_neighbors);
      vec_inLC.push_back(inLC_);
      vec_total_neighbors.push_back(totalNeighbors);
      vec_neighbors_energy.push_back(energy);
    }
  }
  
  // Define CP selector
  CaloParticleSelector cpSelector_ = CaloParticleSelector(0, 100, 1.5, 3,120,280,0,1000000,true,false,false,false,false,{13},-3.2,3.2);

  // Loop over Caloparticles 

  std::vector<DetId> tmprechits_; tmprechits_.clear();
  int obj_id = 0;
  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 

    //Select only Signal
    if (!cpSelector_(cp,simVertices)) continue;
    auto pid = cp.particleId();
    const SimClusterRefVector& simclusters = cp.simClusters();

    int simClusterId = 0;
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


        // Neighbors
        int num_neighbors = -99;
        float energy = -99;
        int totalNeighbors = -99;
        if (detid_ > 100){
          auto neighbors = getNeighbors(detid_,topoEE,topoHF,topoHB);
          auto neighbors_pair = getNeighborHitsAndEnergies(neighbors, hitMap);
          num_neighbors = neighbors_pair.first;
          energy = neighbors_pair.second;
          energy = it_sc_hae.second/energy;
          totalNeighbors = neighbors.size();
        }
        // Check if in LC
        bool inLC_ = inLC(detid_,lcs_detids);


        // LocalError
        const CaloSubdetectorGeometry *subGeom = geom.getSubdetectorGeometry(detid_);
        auto geomEE = static_cast<const HGCalGeometry*>(subGeom);
        const HGCalDDDConstants* ddd = &(geomEE->topology().dddConstants());
        auto lerr = calculateLocalError(detid_, ddd); 
        // // std::cout << "Layer: " << layer_ << "\t XX: " << lerr.xx() << std::endl;

        // KF 
        int kf_compatible=0;
        int kf_contained=0;
        int kf_inSensor=0;
        int kf_pureHit=0;
        int kf_pureHit_approx=0;

        int layer = recHitTools_.getLayerWithOffset(detid_);

        if (std::find(map_detid_kf[layer].begin(), map_detid_kf[layer].end(), static_cast<int32_t>(detid_())) != map_detid_kf[layer].end()) {
          kf_compatible=1;

          //
          DetId closest_detid;
          auto gp = map_gps_kf[layer];
          auto cov = map_cov_kf[layer];
          if (detector == "Sc") closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCell(gp);
          else closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCellHex(gp, true);
          if(detid_==closest_detid){
            kf_contained=1;
          }

          // std::cout << "Simhit conditions" << std::endl;
          if (inSensorCell(detid_, gp)){
            kf_inSensor=1;

            if (isPureHit(detid_,gp,cov)){
              kf_pureHit=1;
            }
            if (isPureHit_approx(detid_,gp,cov)){
              kf_pureHit_approx=1;
            }
          }
        }
        /*
        if (detid_()==static_cast<uint32_t>(map_detid_kf[layer])){
          kf_compatible=1;
        }
        */
        
        //// std::cout << (itcheck->second) << std::endl;

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
        sim_kf_contained.push_back(kf_contained);
        sim_kf_inSensor.push_back(kf_inSensor);
        sim_kf_pureHit.push_back(kf_pureHit);
        sim_kf_pureHit_approx.push_back(kf_pureHit_approx);
        sim_evt.push_back(eventnr);
        sim_obj_id.push_back(obj_id);
        sim_pid.push_back(pid);
        sim_simcluster_id.push_back(simClusterId);
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
          rec_kf_contained.push_back(kf_contained);
          rec_kf_inSensor.push_back(kf_inSensor);
          rec_kf_pureHit.push_back(kf_pureHit);
          rec_kf_pureHit_approx.push_back(kf_pureHit_approx);
          rec_evt.push_back(eventnr);
          rec_obj_id.push_back(obj_id);
          rec_pid.push_back(pid);
          rec_simcluster_id.push_back(simClusterId);
          rec_neighbors.push_back(num_neighbors);
          rec_inLC.push_back(inLC_);
          rec_total_neighbors.push_back(totalNeighbors);
          rec_neighbors_energy.push_back(energy);
          //rec_mask.push_back((itcheck->second));

        }
      }
      simClusterId++;
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
  sim_kf_contained.clear();
  sim_kf_inSensor.clear();
  sim_kf_pureHit.clear();
  sim_kf_pureHit_approx.clear();
  sim_obj_id.clear();
  sim_pid.clear();
  sim_simcluster_id.clear();
  //sim_mask.clear();

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
  rec_kf_contained.clear();
  rec_kf_inSensor.clear();
  rec_kf_pureHit.clear();
  rec_kf_pureHit_approx.clear();
  rec_obj_id.clear();
  rec_pid.clear();
  rec_simcluster_id.clear();
  rec_neighbors.clear();
  rec_inLC.clear();
  rec_total_neighbors.clear();
  rec_neighbors_energy.clear();
  //rec_mask.clear();

  // KF

  kf_x.clear();
  kf_y.clear();
  kf_z.clear();
  kf_meas_x.clear();
  kf_meas_y.clear();
  kf_meas_z.clear();
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
  kf_track_id.clear();
  kf_track_charge.clear();
  kf_track_momentum.clear();
  kf_track_quality.clear();
  kf_track_chi2.clear();
  kf_track_validfraction.clear();
  kf_track_qoverp.clear();
  kf_track_algo.clear();
  kf_rec_compatible.clear();
  kf_rec_isclosest.clear();
  kf_rec_inSensor.clear();
  kf_rec_pureHit.clear();
  kf_rec_pureHit_approx.clear();
  kf_neighbors.clear();
  kf_inLC.clear();
  kf_total_neighbors.clear();
  kf_neighbors_energy.clear();
  
  // Prop

  prop_x.clear();
  prop_y.clear();
  prop_z.clear();
  prop_meas_x.clear();
  prop_meas_y.clear();
  prop_meas_z.clear();
  prop_e.clear();
  prop_detid.clear();
  prop_layer.clear();
  prop_dtype.clear();
  prop_cov_xx.clear();
  prop_cov_xy.clear();
  prop_cov_yy.clear();
  prop_evt.clear();
  prop_eta.clear();
  prop_theta.clear();
  prop_track_id.clear();
  prop_track_charge.clear();
  prop_track_momentum.clear();
  prop_track_quality.clear();
  prop_track_chi2.clear();
  prop_track_validfraction.clear();
  prop_track_qoverp.clear();
  prop_track_algo.clear();
  prop_rec_compatible.clear();
  prop_rec_isclosest.clear();
  prop_rec_inSensor.clear();
  prop_rec_pureHit.clear();
  prop_rec_pureHit_approx.clear();
  prop_neighbors.clear();
  prop_inLC.clear();
  prop_total_neighbors.clear();
  prop_neighbors_energy.clear();

  //clear_arrays();
  // eventnr=eventnr+1;

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
}

//define this as a plug-in
DEFINE_FWK_MODULE(Ntuplizer);
