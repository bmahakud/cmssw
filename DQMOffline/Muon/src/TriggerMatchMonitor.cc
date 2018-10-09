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
  etaBin = parameters.getParameter<int>("etaBin");
  etaMin = parameters.getParameter<double>("etaMin");
  etaMax = parameters.getParameter<double>("etaMax");

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

    // monitoring of eta parameter
    matchHists.push_back(ibooker.book1D("TriggerMatchDR_IsoMu27", "TriggerMatchDR_IsoMu27", 600, 0.0, 6.0));
    matchHists.push_back(ibooker.book1D("TriggerMatchDR_Mu50", "TriggerMatchDR_Mu50", 600, 0.0, 6.0));


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
    //for(unsigned int trigIndex = 0; trigIndex < trigNames.size(); trigIndex++){
    //std::string testTriggerName = trigNames.triggerName(trigIndex);
    //std::cout<<"Trigger Names: "<<testTriggerName<<std::endl;
    //   }

    if(PATmuons.isValid()){//valid pat Muon

      for(const auto & patMuon : *PATmuons){//pat muon loop


      bool Isolated=patMuon.pfIsolationR04().sumChargedHadronPt + TMath::Max(0., patMuon.pfIsolationR04().sumNeutralHadronEt + patMuon.pfIsolationR04().sumPhotonEt - 0.5*patMuon.pfIsolationR04().sumPUPt)  / patMuon.pt() < 0.25;
      
      std::cout<<"Offline Muon:............................. "<<std::endl;

      if(patMuon.isGlobalMuon() && Isolated && patMuon.isTightMuon(thePrimaryVertex)){//isolated tight muon

      TLorentzVector offlineMuon;
      offlineMuon.SetPtEtaPhiM(patMuon.pt(), patMuon.eta(), patMuon.phi(),0.0);


      for (pat::TriggerObjectStandAlone object : *triggerObjects) {
           object.unpackPathNames(trigNames);        

       if(object.pt() < 10.0) continue;//look for HLT objects with Pt > 10 GeV only.

      //std::cout << "\t   Collection: " << object.collection() << std::endl;
      if (object.collection() !="hltIterL3MuonCandidates::HLT")continue;

      std::cout << "\tTrigger object:  pt " << object.pt() << ", eta " << object.eta() << ", phi " << object.phi() << std::endl;


      const std::vector<std::string>& pathNamesLast = object.pathNames(true); 

      const std::vector<std::string>& pathNamesAll = object.pathNames(false);



      for (const auto& pathName : pathNamesLast){
            //bool isBoth = object.hasPathName( pathName, true, true );//object is associated with l3 filter and associated to the last filter of a successful path. this object caused the trigger to fire.
            const std::string& path_i = pathName;


            std::cout<<"pathName: "<<path_i<<std::endl;

            bool isBoth = object.hasPathName( pathName, true, true );
            bool isL3   = object.hasPathName( pathName, false, true );
            bool isLF   = object.hasPathName( pathName, true, false );
            bool isNone = object.hasPathName( pathName, false, false );

            std::cout<<"isBoth: "<<isBoth<<std::endl;
            std::cout<<"isL3: "<<isL3<<std::endl;
            std::cout<<"isLF: "<<isLF<<std::endl;
            std::cout<<"isNone: "<<isNone<<std::endl;



        if(isBoth && path_i.find("HLT_IsoMu27_v")!= std::string::npos){//iso mu 27


        TLorentzVector TriObjStandaloneMuon;
        TriObjStandaloneMuon.SetPtEtaPhiM(object.pt(), object.eta(), object.phi(),0.0);

        double DeltaR=offlineMuon.DeltaR(TriObjStandaloneMuon);

        matchHists[0]->Fill(DeltaR);

        }//iso mu 27


       if(isBoth && path_i.find("HLT_Mu50_v13") != std::string::npos){//mu 50
         
       TLorentzVector TriObjStandaloneMuon;
       TriObjStandaloneMuon.SetPtEtaPhiM(object.pt(), object.eta(), object.phi(),0.0);
       double DeltaR=offlineMuon.DeltaR(TriObjStandaloneMuon);
       matchHists[1]->Fill(DeltaR);


       }// mu 50






      }



      }//trigger object standalone 


      }//isolated tight muon
     } //pat muon loop
    } //valid pat muon







 
}
