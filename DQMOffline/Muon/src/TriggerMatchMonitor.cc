/*
 *  See header file for a description of this class.
 *
 *  \author Bibhuprasad Mahakud (Purdue University, West Lafayette, USA)
 */
#include "DQMOffline/Muon/interface/TriggerMatchMonitor.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "TLorentzVector.h"

#include <string>
#include <TMath.h>
using namespace std;
using namespace edm;

//#define DEBUG

TriggerMatchMonitor::TriggerMatchMonitor(const edm::ParameterSet& pSet) {
  LogTrace(metname)<<"[TriggerMatchMonitor] Parameters initialization";
  
  parameters = pSet;

  // the services
  theService = new MuonServiceProxy(parameters.getParameter<ParameterSet>("ServiceParameters"));
  theDbe = edm::Service<DQMStore>().operator->();
  
  beamSpotToken_ = consumes<reco::BeamSpot >(parameters.getUntrackedParameter<edm::InputTag>("offlineBeamSpot")),
  primaryVerticesToken_ = consumes<std::vector<reco::Vertex> >(parameters.getUntrackedParameter<edm::InputTag>("offlinePrimaryVertices")),
  theMuonCollectionLabel_ = consumes<edm::View<reco::Muon> >  (parameters.getParameter<edm::InputTag>("MuonCollection"));
  thePATMuonCollectionLabel_ = consumes<edm::View<pat::Muon> >  (parameters.getParameter<edm::InputTag>("patMuonCollection"));
  theVertexLabel_          = consumes<reco::VertexCollection>(parameters.getParameter<edm::InputTag>("VertexLabel"));
  theBeamSpotLabel_        = mayConsume<reco::BeamSpot>      (parameters.getParameter<edm::InputTag>("BeamSpotLabel"));
  triggerResultsToken_ = consumes<edm::TriggerResults>(parameters.getUntrackedParameter<edm::InputTag>("triggerResults"));
  triggerObjects_ = consumes<std::vector<pat::TriggerObjectStandAlone>>(parameters.getParameter<edm::InputTag>("triggerObjects"));

  // Parameters

  theFolder = parameters.getParameter<string>("folder");
}
TriggerMatchMonitor::~TriggerMatchMonitor() { 
  delete theService;
}

void TriggerMatchMonitor::bookHistograms(DQMStore::IBooker & ibooker,
					  edm::Run const & /*iRun*/,
					  edm::EventSetup const& /*iSetup*/){
    ibooker.cd();
    ibooker.setCurrentFolder(theFolder);

    // monitoring of trigger match parameter
   
    matchHists.push_back(ibooker.book1D("DelR_HLT_IsoMu24_v1", "DeltaR_(offline,HLT)_triggerPass(IsoMu24)", 500, 0.0, 0.5));
    matchHists.push_back(ibooker.book1D("DelR_HLT_IsoMu24_v2", "DeltaR_(offline,HLT)_triggerPass(IsoMu24)", 1000, 0.0, 1.0));
    matchHists.push_back(ibooker.book1D("DelR_HLT_IsoMu24_v3", "DeltaR_(offline,HLT)_triggerPass(IsoMu24)", 5000, 0.0, 5.0));
    matchHists.push_back(ibooker.book1D("PtRatio_HLT_IsoMu24", "PtRatio_(HLTPt/OfflinePt)_IsoMu24", 200, -5., 5.0));

    matchHists.push_back(ibooker.book1D("DelR_L1_IsoMu24_v1", "DeltaR_(offline, L1)_triggerPass(IsoMu24)", 500, 0.0, 1.0));
    matchHists.push_back(ibooker.book1D("DelR_L1_IsoMu24_v2", "DeltaR_(offline, L1)_triggerPass(IsoMu24)", 1000, 0.0, 2.0));
    matchHists.push_back(ibooker.book1D("PtRatio_L1_IsoMu24", "PtRatio_(HLTPt/OfflinePt)_IsoMu24", 200, -5., 5.0));
   
    matchHists.push_back(ibooker.book1D("DelR_HLT_Mu50_v1", "DeltaR_(offline,HLT)_triggerPass(Mu50)", 500, 0.0, 0.5));
    matchHists.push_back(ibooker.book1D("DelR_HLT_Mu50_v2", "DeltaR_(offline,HLT)_triggerPass(Mu50)", 1000, 0.0, 1.0));
    matchHists.push_back(ibooker.book1D("DelR_HLT_Mu50_v3", "DeltaR_(offline,HLT)_triggerPass(Mu50)", 5000, 0.0, 5.0));
    matchHists.push_back(ibooker.book1D("PtRatio_HLT_Mu50", "PtRatio_(HLTPt/OfflinePt)_Mu50", 200, -5., 5.0));

    matchHists.push_back(ibooker.book1D("DelR_L1_Mu50_v1", "DeltaR_(offline, L1)_triggerPass(Mu50)", 500, 0.0, 1.0));
    matchHists.push_back(ibooker.book1D("DelR_L1_Mu50_v2", "DeltaR_(offline, L1)_triggerPass(Mu50)", 1000, 0.0, 2.0));
    matchHists.push_back(ibooker.book1D("PtRatio_L1_Mu50", "PtRatio_(HLTPt/OfflinePt)_Mu50", 200, -5., 5.0));

 




 


}
void TriggerMatchMonitor::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
   LogTrace(metname)<<"[TriggerMatchMonitor] Analyze the mu in different eta regions";
   theService->update(iSetup);
  
   edm::Handle<edm::View<reco::Muon> > muons; 
   iEvent.getByToken(theMuonCollectionLabel_,muons);

   edm::Handle<edm::View<pat::Muon> > PATmuons;
   iEvent.getByToken(thePATMuonCollectionLabel_,PATmuons);

   edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjects;
   iEvent.getByToken(triggerObjects_, triggerObjects);

   Handle<edm::TriggerResults> triggerResults;
   iEvent.getByToken(triggerResultsToken_, triggerResults);


   reco::Vertex::Point posVtx;
   reco::Vertex::Error errVtx;
   Handle<std::vector<reco::Vertex> > recVtxs;
   iEvent.getByToken(primaryVerticesToken_,recVtxs);
   unsigned int theIndexOfThePrimaryVertex = 999.;
   for (unsigned int ind = 0; ind < recVtxs->size(); ++ind) {
     if ( (*recVtxs)[ind].isValid() && !((*recVtxs)[ind].isFake()) ) {
       theIndexOfThePrimaryVertex = ind;
       break;
     }
   }
   if (theIndexOfThePrimaryVertex<100) {
     posVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).position();
     errVtx = ((*recVtxs)[theIndexOfThePrimaryVertex]).error();
   }
   else {
     LogInfo("RecoMuonValidator") << "reco::PrimaryVertex not found, use BeamSpot position instead\n";
     Handle<reco::BeamSpot> recoBeamSpotHandle;
     iEvent.getByToken(beamSpotToken_,recoBeamSpotHandle);
     reco::BeamSpot bs = *recoBeamSpotHandle;
     posVtx = bs.position();
     errVtx(0,0) = bs.BeamWidthX();
     errVtx(1,1) = bs.BeamWidthY();
     errVtx(2,2) = bs.sigmaZ();
   }
   const reco::Vertex thePrimaryVertex(posVtx,errVtx);





   const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);

   if(PATmuons.isValid()){//valid pat Muon

      for(const auto & patMuon : *PATmuons){//pat muon loop
      bool Isolated=patMuon.pfIsolationR04().sumChargedHadronPt + TMath::Max(0., patMuon.pfIsolationR04().sumNeutralHadronEt + patMuon.pfIsolationR04().sumPhotonEt - 0.5*patMuon.pfIsolationR04().sumPUPt)  / patMuon.pt() < 0.25;
      

      if(patMuon.isGlobalMuon() && Isolated && patMuon.isTightMuon(thePrimaryVertex)){//isolated tight muon

      TLorentzVector offlineMuon;
      offlineMuon.SetPtEtaPhiM(patMuon.pt(), patMuon.eta(), patMuon.phi(),0.0);

      char array24[] = "HLT_IsoMu24_v*";  // 
      char *ptrmu24;        // 
      ptrmu24 = array24;      //

      char array50[] = "HLT_Mu50_v*";  // 
      char *ptrmu50;        // 
      ptrmu50 = array50;      //





      if(patMuon.triggered(ptrmu24)){

      try{


      TLorentzVector hltMuon;
      hltMuon.SetPtEtaPhiM(patMuon.hltObject()->pt(),patMuon.hltObject()->eta(),patMuon.hltObject()->phi(),0.0);
   
      double DelRrecoHLT=offlineMuon.DeltaR(hltMuon);
        

      matchHists[0]->Fill(DelRrecoHLT);   
      matchHists[1]->Fill(DelRrecoHLT);
      matchHists[2]->Fill(DelRrecoHLT);
      matchHists[3]->Fill(patMuon.hltObject()->pt()/patMuon.pt());

      if(patMuon.l1Object() !=nullptr){
      TLorentzVector L1Muon;
      L1Muon.SetPtEtaPhiM(patMuon.l1Object()->pt(),patMuon.l1Object()->eta(),patMuon.l1Object()->phi(),0.0);
      double DelRrecoL1=offlineMuon.DeltaR(L1Muon);
      matchHists[4]->Fill(DelRrecoL1);
      matchHists[5]->Fill(DelRrecoL1);
      matchHists[6]->Fill(patMuon.l1Object()->pt()/patMuon.pt());     


      }

       }catch(...){}

      }


      /// Mu50 test
      if(patMuon.triggered(ptrmu50)){

      try{


      TLorentzVector hltMuon_;
      hltMuon_.SetPtEtaPhiM(patMuon.hltObject()->pt(),patMuon.hltObject()->eta(),patMuon.hltObject()->phi(),0.0);

      double DelRrecoHLT_=offlineMuon.DeltaR(hltMuon_);


      matchHists[7]->Fill(DelRrecoHLT_);
      matchHists[8]->Fill(DelRrecoHLT_);
      matchHists[9]->Fill(DelRrecoHLT_);
      matchHists[10]->Fill(patMuon.hltObject()->pt()/patMuon.pt());

      if(patMuon.l1Object() !=nullptr){
      TLorentzVector L1Muon_;
      L1Muon_.SetPtEtaPhiM(patMuon.l1Object()->pt(),patMuon.l1Object()->eta(),patMuon.l1Object()->phi(),0.0);
      double DelRrecoL1_=offlineMuon.DeltaR(L1Muon_);
      matchHists[11]->Fill(DelRrecoL1_);
      matchHists[12]->Fill(DelRrecoL1_);
      matchHists[13]->Fill(patMuon.l1Object()->pt()/patMuon.pt());


      }

       }catch(...){}

      }











      }//isolated tight muon
      } //pat muon loop
      } //valid pat muon







 
}
