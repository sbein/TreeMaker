//
// authors:  Viktor Kutzner, Sam Bein
//
// Produces isolated tracks including those that disappear
//
//         Created:  Wed, 17 May 2017 13:05:08 GMT
//

#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <TLorentzVector.h>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "./DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/DetId/interface/DetId.h"


#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "RecoTracker/DeDx/interface/DeDxTools.h"

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"

#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

//
// class declaration
//

class IsoTrackProducer : public edm::global::EDProducer<> {
public:
  explicit IsoTrackProducer(const edm::ParameterSet&);
   ~IsoTrackProducer() override;  
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
   
private:
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;
    
  float getDeDx(const reco::DeDxHitInfo *hitInfo, bool, bool) const;
  edm::EDGetTokenT<std::vector<reco::Track>> selectedTracksToken;
  edm::EDGetTokenT<std::vector<reco::GsfElectron>> electronsToken;
  edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken;
  edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken;

  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEBToken;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEEToken;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsESToken;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> reducedHcalRecHitsHBToken;
  edm::EDGetTokenT<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >> reducedHcalRecHitsHFToken;
  edm::EDGetTokenT<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >> reducedHcalRecHitsHOToken;
  edm::EDGetTokenT<std::vector<reco::CaloJet>> caloJetsToken;
  edm::EDGetTokenT<reco::VertexCollection> PrimVtxToken;
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandsToken;
  
  
  //from pat de/dx
    //Pat iso track
          edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > gt2dedxStrip_;
          edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > gt2dedxPixel_;
          edm::EDGetTokenT<reco::DeDxHitInfoAss> gt2dedxHitInfo_;
          edm::EDGetTokenT<edm::ValueMap<int> >  gt2dedxHitInfoPrescale_; 
  //from pat de/dx

};

//
// constants, enums and typedefs
//
  double minTrackPt;
  double maxTrackEta;
  double coneRelIsoDR;
  double conePtSumMaxPtFraction;
  double minTrackJetDR;
  double minTrackLeptonDR;
  double maxChargedPFCandSumDR;
  double maxNeutralPFCandSumDR;
  int RequireNumberOfValidPixelHits;
  int RequireNumberOfValidTrackerHits;
  double maxDxy;
  double maxDz;
  int minMissingOuterHits;
  int pixelLayersWithMeasurement;
  int trackerLayersWithMeasurement;
  bool trackQualityUndef;
  bool trackQualityLoose;
  bool trackQualityTight;
  bool trackQualityHighPurity;
  bool trackQualityConfirmed;
  bool trackQualityGoodIterative;
  bool trackQualityLooseSetWithPV;
  bool trackQualityHighPuritySetWithPV;
  bool trackQualityDiscarded;
  bool trackQualitySize;
  double caloEnergyDepositionMaxDR;
  double caloEnergyDepositionMaxE;
  double deadNoisyDR;     
  bool useCaloJetsInsteadOfHits;
  bool doDeDx = false;
  std::string dEdxEstimator;

  bool passExo16044GapsVeto;
  bool passExo16044DeadNoisyECALVeto;
  bool passPFCandVeto;
  
  //Pat iso track          
          bool addPrescaledDeDxTracks_;
          bool usePrecomputedDeDxStrip_;
          bool usePrecomputedDeDxPixel_;
  //Pat iso track
//
// static data member definitions
//

//
// constructors and destructor
//
IsoTrackProducer::IsoTrackProducer(const edm::ParameterSet& iConfig) 
{   
  //register your product    
  produces<std::vector<TLorentzVector> > ("tracks");
  produces<std::vector<double> >         ("tracks@dxyVtx");
  produces<std::vector<double> >         ("tracks@dzVtx");
  produces<std::vector<int> >            ("tracks@nMissingOuterHits");
  produces<std::vector<int> >            ("tracks@pixelLayersWithMeasurement");
  produces<std::vector<int> >            ("tracks@trackerLayersWithMeasurement");
  produces<std::vector<int> >            ("tracks@nMissingInnerHits");
  produces<std::vector<int> >            ("tracks@nMissingMiddleHits");
  produces<std::vector<int> >            ("tracks@nValidPixelHits");
  produces<std::vector<int> >            ("tracks@nValidTrackerHits");
  produces<std::vector<double> >         ("tracks@ptError");
  produces<std::vector<double> >         ("tracks@chi2perNdof");
  produces<std::vector<double> >         ("tracks@trkRelIso");
  produces<std::vector<double> >         ("tracks@trkMiniRelIso");
  produces<std::vector<double> >         ("tracks@minDrLepton");
  produces<std::vector<double> >         ("tracks@trackJetIso");
  produces<std::vector<double> >         ("tracks@matchedCaloEnergy");
  produces<std::vector<double> >         ("tracks@matchedCaloEnergyJets");
  produces<std::vector<double> >         ("tracks@deDxHarmonic2strips");
  produces<std::vector<double> >         ("tracks@deDxHarmonic2pixel");  
  produces<std::vector<double> >         ("tracks@chargedPtSum");
  produces<std::vector<double> >         ("tracks@neutralPtSum");
  produces<std::vector<double> >         ("tracks@neutralWithoutGammaPtSum");
  produces<std::vector<bool> >           ("tracks@passExo16044GapsVeto");
  produces<std::vector<bool> >           ("tracks@passExo16044DeadNoisyECALVeto");
  produces<std::vector<bool> >           ("tracks@passExo16044Tag");
  produces<std::vector<bool> >           ("tracks@passExo16044JetIso");
  produces<std::vector<bool> >           ("tracks@passExo16044LepIso");
  produces<std::vector<bool> >           ("tracks@passPFCandVeto");
  produces<std::vector<bool> >           ("tracks@trackQualityUndef");
  produces<std::vector<bool> >           ("tracks@trackQualityLoose");
  produces<std::vector<bool> >           ("tracks@trackQualityTight");
  produces<std::vector<bool> >           ("tracks@trackQualityHighPurity");
  produces<std::vector<bool> >           ("tracks@trackQualityConfirmed");
  produces<std::vector<bool> >           ("tracks@trackQualityGoodIterative");
  produces<std::vector<bool> >           ("tracks@trackQualityLooseSetWithPV");
  produces<std::vector<bool> >           ("tracks@trackQualityHighPuritySetWithPV");
  produces<std::vector<bool> >           ("tracks@trackQualityDiscarded");
  produces<std::vector<bool> >           ("tracks@trackQualitySize");
  produces<std::vector<int> >         ("tracks@charge");


  PrimVtxToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("PrimaryVertex"));

  selectedTracksToken = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("selectedTracks"));
  electronsToken = consumes<std::vector<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("selectedElectrons"));
  muonsToken = consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("selectedMuons"));
  jetsToken = consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("selectedPFJets"));
  selectedTracksToken = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("selectedTracks"));
  reducedEcalRecHitsEBToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedEcalRecHitsEB"));
  reducedEcalRecHitsEEToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedEcalRecHitsEE"));
  reducedEcalRecHitsESToken = consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >>(iConfig.getParameter<edm::InputTag>("selectedEcalRecHitsES"));
  reducedHcalRecHitsHBToken = consumes<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >>(edm::InputTag("reducedHcalRecHits:hbhereco"));
  reducedHcalRecHitsHFToken = consumes<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >>(edm::InputTag("reducedHcalRecHits:hfreco"));
  reducedHcalRecHitsHOToken = consumes<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >>(edm::InputTag("reducedHcalRecHits:horeco"));
  caloJetsToken = consumes<std::vector<reco::CaloJet>>(iConfig.getParameter<edm::InputTag>("selectedCaloJets"));   
  pfCandsToken = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("selectedPFCand"));

  
  //edm::InputTag estimatorTag("dedxHarmonic2");
  //consumes<edm::ValueMap<reco::DeDxData>>(estimatorTag);

  
  
  edm::InputTag dEdxTag("dedxHitInfo");
  consumes<reco::DeDxHitInfoAss>(dEdxTag);

  minTrackPt = iConfig.getParameter<double>("minTrackPt");
  maxTrackEta = iConfig.getParameter<double>("maxTrackEta");
  coneRelIsoDR = iConfig.getParameter<double>("coneRelIsoDR");
  conePtSumMaxPtFraction = iConfig.getParameter<double>("conePtSumMaxPtFraction");
  minTrackJetDR = iConfig.getParameter<double>("minTrackJetDR");
  minTrackLeptonDR = iConfig.getParameter<double>("minTrackLeptonDR");
  RequireNumberOfValidPixelHits = iConfig.getParameter<int>("RequireNumberOfValidPixelHits");
  RequireNumberOfValidTrackerHits = iConfig.getParameter<int>("RequireNumberOfValidTrackerHits");
  maxDxy = iConfig.getParameter<double>("maxDxy");
  maxDz = iConfig.getParameter<double>("maxDz");

  maxChargedPFCandSumDR = iConfig.getParameter<double>("maxChargedPFCandSumDR");
  maxNeutralPFCandSumDR = iConfig.getParameter<double>("maxNeutralPFCandSumDR");

  minMissingOuterHits = iConfig.getParameter<int>("minMissingOuterHits");
  caloEnergyDepositionMaxDR = iConfig.getParameter<double>("caloEnergyDepositionMaxDR");
  caloEnergyDepositionMaxE = iConfig.getParameter<double>("caloEnergyDepositionMaxE");
  deadNoisyDR = iConfig.getParameter<double>("deadNoisyDR");

  useCaloJetsInsteadOfHits = iConfig.getParameter<bool>("useCaloJetsInsteadOfHits");
  doDeDx = iConfig.getParameter<bool>("doDeDx");
  //dEdxEstimator = iConfig.getParameter<std::string>("dEdxEstimator");//old original, pre-pat
  
  
  //from PAT iso track
  gt2dedxStrip_ = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("dEdxDataStrip"));
  gt2dedxPixel_ = consumes<edm::ValueMap<reco::DeDxData> >(iConfig.getParameter<edm::InputTag>("dEdxDataPixel"));
  gt2dedxHitInfo_= consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dEdxHitInfo"));
  addPrescaledDeDxTracks_ = iConfig.getParameter<bool>("addPrescaledDeDxTracks");
  gt2dedxHitInfoPrescale_ = addPrescaledDeDxTracks_ ? consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("dEdxHitInfoPrescale")) : edm::EDGetTokenT<edm::ValueMap<int>>();
  usePrecomputedDeDxStrip_ = iConfig.getParameter<bool>("usePrecomputedDeDxStrip");
  usePrecomputedDeDxPixel_ = iConfig.getParameter<bool>("usePrecomputedDeDxPixel");
  //from PAT iso track
}


IsoTrackProducer::~IsoTrackProducer()
{}

template <typename T1, typename T2>
double getParticleIsolation(T1 particle, T2 pcollection){

  // get DR of a selected particle with respect to all other particles in a given collection

  double minDR = 1e6;
  for( const auto& item : *pcollection){
    if ((deltaR(particle, item) < minDR) && (deltaR(particle, item) != 0)) minDR = deltaR(particle, item);
  }
  return minDR;
}

template <typename T1, typename T2>
bool checkParticleIsolation(T1 particle, T2 pcollection, double RequiredMinDR){

  // checks particle isolation with respect to another particle collection,
  // returns true if two particles are not inside a cone with minimal DR.

  bool isolated = false;
  double minDR = 1e6;
  for( const auto& item : *pcollection){
    if ((deltaR(particle, item) < minDR) && (deltaR(particle, item) != 0)) minDR = deltaR(particle, item);
  }
  if (minDR > RequiredMinDR) isolated = true;

  return isolated;
}


template <typename T>
double LoopOverRecHits(T hitcollection, reco::Track track, edm::ESHandle<CaloGeometry> CaloGeomHandle, double caloEnergyDepositionMaxDR){

  // get deposited energy for a single recHit collection inside a cone of max DR around a selected track
   
  double energyDepositedPerCalorimeter = 0;      
  for( const auto& hit : *hitcollection){ 
          
    double hitEnergy = hit.energy();
    double cellEta = CaloGeomHandle.product()->getPosition(hit.detid()).eta();
    double cellPhi = CaloGeomHandle.product()->getPosition(hit.detid()).phi();
        
    if (deltaR(cellEta, cellPhi, track.eta(), track.phi()) < caloEnergyDepositionMaxDR) {
      energyDepositedPerCalorimeter += hitEnergy;
    }
  }
       
  return energyDepositedPerCalorimeter;
}


template <typename T>
bool checkNoDeadNoisyECALInTrackCone(T hitcollection, reco::Track track, edm::ESHandle<CaloGeometry> CaloGeomHandle, double DeadNoisyDR){

  // check if selected track points to dead or noisy ECAL channels within a cone of max DR:
   
  double noDeadNoisyCellsInTrackCone = true;

  for( const auto& hit : *hitcollection){
    double cellEta = CaloGeomHandle.product()->getPosition(hit.detid()).eta();
    double cellPhi = CaloGeomHandle.product()->getPosition(hit.detid()).phi();

    std::vector<int> badFlags(2);
    badFlags[0] = EcalRecHit::kDead;
    badFlags[1] = EcalRecHit::kNoisy;

    if (hit.checkFlags(badFlags)) {
      if (deltaR(cellEta, cellPhi, track.eta(), track.phi()) < DeadNoisyDR) {
           noDeadNoisyCellsInTrackCone = false;
      }
    }
  }

  return noDeadNoisyCellsInTrackCone;
}

//
// member functions
//


float IsoTrackProducer::getDeDx(const reco::DeDxHitInfo *hitInfo, bool doPixel, bool doStrip) const
{
    if(hitInfo == nullptr){
        return -1;
    }

    std::vector<float> charge_vec;
    for(unsigned int ih=0; ih<hitInfo->size(); ih++){

        bool isPixel = (hitInfo->pixelCluster(ih) != nullptr);
        bool isStrip = (hitInfo->stripCluster(ih) != nullptr);

        if(isPixel && !doPixel) continue;
        if(isStrip && !doStrip) continue;

        // probably shouldn't happen
        if(!isPixel && !isStrip) continue;

        // shape selection for strips
        if(isStrip && !DeDxTools::shapeSelection(*(hitInfo->stripCluster(ih))))
            continue;
        
        float Norm=0;       
        if(isPixel)
            Norm = 3.61e-06; //compute the normalization factor to get the energy in MeV/mm
        if(isStrip)
            Norm = 3.61e-06 * 265;
            
        //save the dE/dx in MeV/mm to a vector.  
        charge_vec.push_back(Norm*hitInfo->charge(ih)/hitInfo->pathlength(ih));           
    }

    int size = charge_vec.size();
    float result = 0.0;

    //build the harmonic 2 dE/dx estimator
    float expo = -2;
    for(int i=0; i<size; i++){
        result += pow(charge_vec[i], expo);
    }
    result = (size>0) ? pow(result/size, 1./expo) : 0.0;

    return result;
}



// ------------ method called to produce the data  ------------
void IsoTrackProducer::produce(edm::StreamID, edm::Event& iEvent, const edm::EventSetup& iSetup) const
{
  using namespace edm;
  using namespace std;  
  
  // register input collections:
  Handle<std::vector<reco::Track>> selectedTracks;
  iEvent.getByToken( selectedTracksToken, selectedTracks);
  Handle<std::vector<reco::GsfElectron>> electrons;
  iEvent.getByToken( electronsToken, electrons);
  Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken( muonsToken, muons);
  Handle<std::vector<reco::PFJet>> jets;
  iEvent.getByToken( jetsToken, jets);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEB;
  iEvent.getByToken( reducedEcalRecHitsEBToken, reducedEcalRecHitsEB);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsEE;
  iEvent.getByToken( reducedEcalRecHitsEEToken, reducedEcalRecHitsEE);
  Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >> reducedEcalRecHitsES;
  iEvent.getByToken( reducedEcalRecHitsESToken, reducedEcalRecHitsES);
  Handle<edm::SortedCollection<HBHERecHit,edm::StrictWeakOrdering<HBHERecHit> >> reducedHcalRecHitsHB;
  iEvent.getByToken( reducedHcalRecHitsHBToken, reducedHcalRecHitsHB);
  Handle<edm::SortedCollection<HFRecHit,edm::StrictWeakOrdering<HFRecHit> >> reducedHcalRecHitsHF;
  iEvent.getByToken( reducedHcalRecHitsHFToken, reducedHcalRecHitsHF);
  Handle<edm::SortedCollection<HORecHit,edm::StrictWeakOrdering<HORecHit> >> reducedHcalRecHitsHO;
  iEvent.getByToken( reducedHcalRecHitsHOToken, reducedHcalRecHitsHO);
  Handle<std::vector<reco::CaloJet>> caloJets;
  iEvent.getByToken( caloJetsToken, caloJets);
  Handle<std::vector<reco::PFCandidate>> pfCands;
  iEvent.getByToken( pfCandsToken, pfCands);

 
  edm::Handle<edm::ValueMap<reco::DeDxData>> dEdxTrackHandle;
  edm::ValueMap<reco::DeDxData> dEdxTrack;
  edm::Handle<reco::DeDxHitInfoAss> dedxCollH;//Toby  
  

    // associate generalTracks with their DeDx data (estimator for strip dE/dx)
    edm::Handle<edm::ValueMap<reco::DeDxData> > gt2dedxStrip;
    iEvent.getByToken(gt2dedxStrip_, gt2dedxStrip);

    // associate generalTracks with their DeDx data (estimator for pixel dE/dx)
    edm::Handle<edm::ValueMap<reco::DeDxData> > gt2dedxPixel;
    iEvent.getByToken(gt2dedxPixel_, gt2dedxPixel);

    // associate generalTracks with their DeDx hit info (used to estimate pixel dE/dx)
    edm::Handle<reco::DeDxHitInfoAss> gt2dedxHitInfo;
    iEvent.getByToken(gt2dedxHitInfo_, gt2dedxHitInfo);
    edm::Handle<edm::ValueMap<int>> gt2dedxHitInfoPrescale;
    if (addPrescaledDeDxTracks_) {
        iEvent.getByToken(gt2dedxHitInfoPrescale_, gt2dedxHitInfoPrescale);
}

  //edm::Handle<std::vector<reco::DeDxHitInfo>> dedxCollH;//Toby   modified
  //edm::Handle<vector<reco::DeDxHitInfo>> dedxCollH;//Toby  
  //edm::Handle<edm::Association<reco::DeDxHitInfoCollection>> dedxCollH;//Toby  modified
  //edm::Handle<edm::Association<vector<reco::DeDxHitInfo> >> dedxCollH;//Toby modified
  
/*  if (doDeDx)
{
  iEvent.getByLabel(dEdxEstimator, dEdxTrackHandle);
  dEdxTrack = *dEdxTrackHandle.product();
  iEvent.getByLabel("dedxHitInfo", dedxCollH);//Toby
}//Original de/dx pre-pat*/ 

  edm::ESHandle<CaloGeometry> CaloGeomHandle;
  iSetup.get<CaloGeometryRecord>().get(CaloGeomHandle);

  edm::Handle<reco::VertexCollection> vtx_h;
  iEvent.getByToken(PrimVtxToken, vtx_h);
  const reco::Vertex vtx = vtx_h->at(0);

  // register output:
  std::unique_ptr<std::vector<TLorentzVector> > tracks(new std::vector<TLorentzVector>);

  std::unique_ptr<std::vector<double> > selectedTracks_dxyVtx(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_dzVtx(new std::vector<double>);
  std::unique_ptr<std::vector<int> > selectedTracks_trackQuality(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_nMissingOuterHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_pixelLayersWithMeasurement(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_trackerLayersWithMeasurement(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_nMissingInnerHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_nMissingMiddleHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_nValidPixelHits(new std::vector<int>);
  std::unique_ptr<std::vector<int> > selectedTracks_nValidTrackerHits(new std::vector<int>);
  std::unique_ptr<std::vector<double> > selectedTracks_ptError(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_chi2perNdof(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_trkRelIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_trkMiniRelIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_trackJetIso(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_minDrLepton(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_matchedCaloEnergy(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_matchedCaloEnergyJets(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_deDxHarmonic2strips(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_deDxHarmonic2pixel(new std::vector<double>);  
  std::unique_ptr<std::vector<double> > selectedTracks_chargedPtSum(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_neutralPtSum(new std::vector<double>);
  std::unique_ptr<std::vector<double> > selectedTracks_neutralWithoutGammaPtSum(new std::vector<double>);
  std::unique_ptr<std::vector<bool> > selectedTracks_passExo16044GapsVeto(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_passExo16044DeadNoisyECALVeto(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_passPFCandVeto(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_passExo16044JetIso(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_passExo16044LepIso(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_passExo16044Tag(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityUndef(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityLoose(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityTight(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityHighPurity(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityConfirmed(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityGoodIterative(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityLooseSetWithPV(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityHighPuritySetWithPV(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualityDiscarded(new std::vector<bool>);
  std::unique_ptr<std::vector<bool> > selectedTracks_trackQualitySize(new std::vector<bool>);
  std::unique_ptr<std::vector<int> > selectedTracks_charge(new std::vector<int>);

  int itrack = -1;
  for( const auto& track : *selectedTracks){
    itrack+=1;

    // track kinematics:
    bool trackKinematics = ((track.pt() > minTrackPt) && (std::abs(track.eta())) < maxTrackEta);
    if (!(trackKinematics)) continue;

    bool passExo16044Kinematics = track.pt() > 55 && std::abs(track.eta()) < 2.1;

    // vertex consistency:
    bool vertex_match = std::abs(track.dxy(vtx.position())) < maxDxy && std::abs(track.dz(vtx.position())) < maxDz;
    if (!(vertex_match)) continue;

    //This completes the loose candidate selection.

    // check pt sum of tracker activity within a cone around current track:
    bool passedTrackTrackerIso = false;
    double relIso = -1;
    double miniIso = -1;
    double conePtSum_rel  = 0;
    double conePtSum_mini = 0;
    double miniConeDr = 0.2;
    if (track.pt()>50) miniConeDr = 10.0/track.pt();
    else if (track.pt()>200) miniConeDr = 0.05;
    for (const auto& othertrack : *selectedTracks){
        if (!(std::abs(othertrack.dxy(vtx.position()))<0.03 && std::abs(othertrack.dz(vtx.position()))<0.05)) continue;
        if (deltaR(track, othertrack) < 0.00001) continue;
        if (deltaR(track, othertrack) < coneRelIsoDR) conePtSum_rel += othertrack.pt();
        if (deltaR(track, othertrack) < miniConeDr) conePtSum_mini += othertrack.pt();
    }
    if (conePtSum_rel / track.pt() <= conePtSumMaxPtFraction) 
    {	
    	passedTrackTrackerIso = true;
    }
    relIso = conePtSum_rel / track.pt();
    miniIso = conePtSum_mini / track.pt();
    if (!(relIso<0.2)) continue;

    selectedTracks_trkRelIso->push_back(relIso);
    selectedTracks_trkMiniRelIso->push_back(miniIso);


    // check track/jet isolation:
    double trackIsoDR = getParticleIsolation(track, jets);
    selectedTracks_trackJetIso->push_back(trackIsoDR);

    bool passedTrackJetIso = false;
    if (trackIsoDR > minTrackJetDR) passedTrackJetIso = true;
    selectedTracks_passExo16044JetIso->push_back(passedTrackJetIso);


    // check track / lepton isolation:
    double minDrLepton = std::min(getParticleIsolation(track, electrons), getParticleIsolation(track, muons));
    selectedTracks_minDrLepton->push_back(minDrLepton);

    bool passedminDrLepton = false;
    if (minDrLepton > minTrackLeptonDR) passedminDrLepton = true;
    selectedTracks_passExo16044LepIso->push_back(passedminDrLepton);
    
    
    
    // check gaps veto:
    passExo16044GapsVeto = true;
    if ((std::abs(track.eta()) > 0.15) && (std::abs(track.eta()) < 0.35) &&
        (std::abs(track.eta()) > 1.55) && (std::abs(track.eta()) < 1.85) &&
        (std::abs(track.eta()) > 1.42) && (std::abs(track.eta()) < 1.65)) passExo16044GapsVeto = false;
    selectedTracks_passExo16044GapsVeto->push_back(passExo16044GapsVeto);

    // check dead or noisy calo cells:
    passExo16044DeadNoisyECALVeto = checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEB, track, CaloGeomHandle, deadNoisyDR) &&
                                    checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsEE, track, CaloGeomHandle, deadNoisyDR) &&
                                    checkNoDeadNoisyECALInTrackCone(reducedEcalRecHitsES, track, CaloGeomHandle, deadNoisyDR);
    selectedTracks_passExo16044DeadNoisyECALVeto->push_back(passExo16044DeadNoisyECALVeto);

    //now fill the products depending only on track object
    TLorentzVector tlv;
    tlv.SetPxPyPzE(track.px(),track.py(),track.pz(),track.pt()*cosh(track.eta()));
    tracks->push_back(tlv);
    selectedTracks_dxyVtx->push_back(std::abs(track.dxy(vtx.position())));
    selectedTracks_dzVtx->push_back(std::abs(track.dz(vtx.position())));
    reco::HitPattern hitpattern = track.hitPattern();
    selectedTracks_pixelLayersWithMeasurement->push_back(hitpattern.pixelLayersWithMeasurement());
    selectedTracks_trackerLayersWithMeasurement->push_back(hitpattern.trackerLayersWithMeasurement());
    selectedTracks_nMissingOuterHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS));
    selectedTracks_nMissingInnerHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS));
    selectedTracks_nMissingMiddleHits->push_back(hitpattern.trackerLayersWithoutMeasurement(hitpattern.TRACK_HITS));
    selectedTracks_nValidPixelHits->push_back(hitpattern.numberOfValidPixelHits());
    selectedTracks_nValidTrackerHits->push_back(hitpattern.numberOfValidTrackerHits());
    selectedTracks_ptError->push_back(track.ptError());
    selectedTracks_chi2perNdof->push_back(1.0*track.chi2()/track.ndof());
    selectedTracks_trackQualityUndef->push_back(track.quality(track.undefQuality));
    selectedTracks_trackQualityLoose->push_back(track.quality(track.loose));
    selectedTracks_trackQualityTight->push_back(track.quality(track.tight));
    selectedTracks_trackQualityHighPurity->push_back(track.quality(track.highPurity));
    selectedTracks_trackQualityConfirmed->push_back(track.quality(track.confirmed));
    selectedTracks_trackQualityGoodIterative->push_back(track.quality(track.goodIterative));
    selectedTracks_trackQualityLooseSetWithPV->push_back(track.quality(track.looseSetWithPV));
    selectedTracks_trackQualityHighPuritySetWithPV->push_back(track.quality(track.highPuritySetWithPV));
    selectedTracks_trackQualityDiscarded->push_back(track.quality(track.discarded));
    selectedTracks_trackQualitySize->push_back(track.quality(track.qualitySize));
    selectedTracks_charge->push_back(track.charge());
    // save charged/neutral pt sum within DR cone (without selected track):
    double chargedPtSum = 0;
    double neutralPtSum = 0;
    double neutralWithoutGammaPtSum = 0;
    passPFCandVeto = true;
    for( const auto& pfCand : *pfCands) {
        double dR = deltaR(track, pfCand);
        if ( (dR < maxChargedPFCandSumDR) && (pfCand.charge() != 0) )
                chargedPtSum += pfCand.pt();

        if ( (dR < maxNeutralPFCandSumDR) && (pfCand.charge() == 0) ) {
                neutralPtSum += pfCand.pt();
                if (pfCand.pdgId() != 22) neutralWithoutGammaPtSum += pfCand.pt();
        }
        if (dR < 0.01) passPFCandVeto = false;
    }
    selectedTracks_chargedPtSum->push_back(chargedPtSum);
    selectedTracks_neutralPtSum->push_back(neutralPtSum);
    selectedTracks_neutralWithoutGammaPtSum->push_back(neutralWithoutGammaPtSum);
    selectedTracks_passPFCandVeto->push_back(passPFCandVeto);

    // now perform the tight selection:

    // check if fake track:
    bool passExo16044NoFakes = false;
    if ( (hitpattern.numberOfValidPixelHits() >= RequireNumberOfValidPixelHits) &&
         (hitpattern.numberOfValidTrackerHits() >= RequireNumberOfValidTrackerHits) ) {
          // checked for minimal number of hits, now require no missing inner or middle hits:
          if ( hitpattern.trackerLayersWithoutMeasurement(hitpattern.MISSING_INNER_HITS) == 0 &&
               hitpattern.trackerLayersWithoutMeasurement(hitpattern.TRACK_HITS) == 0) {
                // check track displacement:
                if ( vertex_match ) passExo16044NoFakes = true;
          }
    }
        
    // check calo energy deposition:
    double energyDeposited = 0;
    energyDeposited += LoopOverRecHits(reducedEcalRecHitsEB, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
    energyDeposited += LoopOverRecHits(reducedEcalRecHitsEE, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
    energyDeposited += LoopOverRecHits(reducedEcalRecHitsES, track, CaloGeomHandle, caloEnergyDepositionMaxDR);
    energyDeposited += LoopOverRecHits(reducedHcalRecHitsHB, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    energyDeposited += LoopOverRecHits(reducedHcalRecHitsHF, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    energyDeposited += LoopOverRecHits(reducedHcalRecHitsHO, track, CaloGeomHandle, caloEnergyDepositionMaxDR); //check input
    bool passedCaloEnergy = false;
    if (energyDeposited < caloEnergyDepositionMaxE) passedCaloEnergy = true;
    selectedTracks_matchedCaloEnergy->push_back(energyDeposited);

    // optional cross-check, loop over calojets:
    double energyDepositedJets = 0;
    for( const auto& jet : *caloJets){ 
         if (deltaR(jet.detectorP4().eta(), jet.detectorP4().phi(), track.eta(), track.phi()) < caloEnergyDepositionMaxDR) {
             energyDepositedJets += jet.maxEInEmTowers() + jet.maxEInHadTowers();
         }
    }
    selectedTracks_matchedCaloEnergyJets->push_back(energyDepositedJets);

    // check minimum of missing outer hits:
    bool passedMissingOuterHits = false;
    int missingOuterHits = hitpattern.stripLayersWithoutMeasurement(hitpattern.MISSING_OUTER_HITS);
    if (missingOuterHits > minMissingOuterHits) passedMissingOuterHits = true;

  
    // save dE/dx:
    
    float dEdxPixel=-1, dEdxStrip=-1;//from pat
            

    if (doDeDx) 
    {	
        reco::TrackRef trackref = reco::TrackRef(selectedTracks, itrack);
    	    	
    	//Toby, Paul
    	/*
		const reco::DeDxHitInfo *dedxHits = NULL;
		reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(trackref.key());
		if (!dedxHitsRef.isNull())
        dedxHits = &(*dedxHitsRef);
		std::cout << itrack <<  " de/dx info, track with pt= " <<  track.pt() << std::endl;        
		if (dedxHits) {
		for (unsigned int h = 0; h < dedxHits->size(); h++) {
		DetId detid(dedxHits->detId(h)); //detid.subdetId() < 3 Pixel, detid.subdetId() >=3 Strip		
    	std::cout << itrack <<  "detid.subdetId()=" << detid.subdetId() << " charge(h)=" << dedxHits->charge(h) << std::endl;
    	}}
    	*/
    	//Toby, Paul 	

	//from pat
	        // if no dEdx info exists, just store -1

        if(usePrecomputedDeDxStrip_ && gt2dedxStrip.isValid() && gt2dedxStrip->contains(trackref.id())){
            dEdxStrip = (*gt2dedxStrip)[trackref].dEdx();
        }else if(gt2dedxHitInfo.isValid() && gt2dedxHitInfo->contains(trackref.id())){
            const reco::DeDxHitInfo* hitInfo = (*gt2dedxHitInfo)[trackref].get();
            //dEdxStrip = pat::PATIsolatedTrackProducer::getDeDx(hitInfo, false, true);
            dEdxStrip = getDeDx(hitInfo, false, true);            
        }
        if(usePrecomputedDeDxPixel_ && gt2dedxPixel.isValid() && gt2dedxPixel->contains(trackref.id())){
            dEdxPixel = (*gt2dedxPixel)[trackref].dEdx();
        }else if(gt2dedxHitInfo.isValid() && gt2dedxHitInfo->contains(trackref.id())){
            const reco::DeDxHitInfo* hitInfo = (*gt2dedxHitInfo)[trackref].get();
            //dEdxPixel = pat::PATIsolatedTrackProducer::getDeDx(hitInfo, true, false);
            dEdxPixel = getDeDx(hitInfo, true, false);            
        }
	//end: from pat        
        
        
    }
    	
	//std::cout <<  itrack << " summary of dedx; dEdxPixel=" << dEdxPixel << " dEdxStrip="<<dEdxStrip<< " and original dedx=" << std::endl;
    selectedTracks_deDxHarmonic2strips->push_back(dEdxStrip);        
    selectedTracks_deDxHarmonic2pixel->push_back(dEdxPixel);            

    // check if EXO-16-044 selection criteria are satisfied:
    if (passExo16044Kinematics && passExo16044GapsVeto && passExo16044DeadNoisyECALVeto && passedTrackTrackerIso &&
        passExo16044NoFakes && passedminDrLepton && passedTrackJetIso && 
        passedCaloEnergy && passedMissingOuterHits) selectedTracks_passExo16044Tag->push_back(true);
    else selectedTracks_passExo16044Tag->push_back(false);

  }
  


  

  // save track collection in event:
  iEvent.put(std::move(tracks),                                 "tracks");
  iEvent.put(std::move(selectedTracks_dxyVtx),                          "tracks@dxyVtx");
  iEvent.put(std::move(selectedTracks_dzVtx),                           "tracks@dzVtx");
  iEvent.put(std::move(selectedTracks_pixelLayersWithMeasurement),      "tracks@pixelLayersWithMeasurement");
  iEvent.put(std::move(selectedTracks_trackerLayersWithMeasurement),    "tracks@trackerLayersWithMeasurement");
  iEvent.put(std::move(selectedTracks_nMissingOuterHits),               "tracks@nMissingOuterHits");
  iEvent.put(std::move(selectedTracks_nMissingInnerHits),               "tracks@nMissingInnerHits");
  iEvent.put(std::move(selectedTracks_nMissingMiddleHits),              "tracks@nMissingMiddleHits");
  iEvent.put(std::move(selectedTracks_nValidPixelHits),                 "tracks@nValidPixelHits");
  iEvent.put(std::move(selectedTracks_nValidTrackerHits),               "tracks@nValidTrackerHits");
  iEvent.put(std::move(selectedTracks_ptError),                         "tracks@ptError");
  iEvent.put(std::move(selectedTracks_chi2perNdof),                     "tracks@chi2perNdof");
  iEvent.put(std::move(selectedTracks_trkRelIso),                       "tracks@trkRelIso");
  iEvent.put(std::move(selectedTracks_trkMiniRelIso),                   "tracks@trkMiniRelIso");
  iEvent.put(std::move(selectedTracks_trackJetIso),                     "tracks@trackJetIso");
  iEvent.put(std::move(selectedTracks_minDrLepton),                  "tracks@minDrLepton");
  iEvent.put(std::move(selectedTracks_matchedCaloEnergy),               "tracks@matchedCaloEnergy");
  iEvent.put(std::move(selectedTracks_matchedCaloEnergyJets),           "tracks@matchedCaloEnergyJets");
  iEvent.put(std::move(selectedTracks_deDxHarmonic2strips),                   "tracks@deDxHarmonic2strips");
  iEvent.put(std::move(selectedTracks_deDxHarmonic2pixel),                   "tracks@deDxHarmonic2pixel");  
  iEvent.put(std::move(selectedTracks_chargedPtSum),                    "tracks@chargedPtSum");
  iEvent.put(std::move(selectedTracks_neutralPtSum),                    "tracks@neutralPtSum");
  iEvent.put(std::move(selectedTracks_neutralWithoutGammaPtSum),        "tracks@neutralWithoutGammaPtSum");
  iEvent.put(std::move(selectedTracks_passExo16044GapsVeto),            "tracks@passExo16044GapsVeto");
  iEvent.put(std::move(selectedTracks_passExo16044DeadNoisyECALVeto),   "tracks@passExo16044DeadNoisyECALVeto");
  iEvent.put(std::move(selectedTracks_passPFCandVeto),                  "tracks@passPFCandVeto");
  iEvent.put(std::move(selectedTracks_passExo16044JetIso),              "tracks@passExo16044JetIso");
  iEvent.put(std::move(selectedTracks_passExo16044LepIso),              "tracks@passExo16044LepIso");
  iEvent.put(std::move(selectedTracks_passExo16044Tag),                 "tracks@passExo16044Tag");
  iEvent.put(std::move(selectedTracks_trackQualityUndef),               "tracks@trackQualityUndef");
  iEvent.put(std::move(selectedTracks_trackQualityLoose),               "tracks@trackQualityLoose");
  iEvent.put(std::move(selectedTracks_trackQualityTight),               "tracks@trackQualityTight");
  iEvent.put(std::move(selectedTracks_trackQualityHighPurity),          "tracks@trackQualityHighPurity");
  iEvent.put(std::move(selectedTracks_trackQualityConfirmed),           "tracks@trackQualityConfirmed");
  iEvent.put(std::move(selectedTracks_trackQualityGoodIterative),       "tracks@trackQualityGoodIterative");
  iEvent.put(std::move(selectedTracks_trackQualityLooseSetWithPV),      "tracks@trackQualityLooseSetWithPV");
  iEvent.put(std::move(selectedTracks_trackQualityHighPuritySetWithPV), "tracks@trackQualityHighPuritySetWithPV");
  iEvent.put(std::move(selectedTracks_trackQualityDiscarded),           "tracks@trackQualityDiscarded");
  iEvent.put(std::move(selectedTracks_trackQualitySize),                "tracks@trackQualitySize");
  iEvent.put(std::move(selectedTracks_charge),                          "tracks@charge");
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void IsoTrackProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(IsoTrackProducer);
