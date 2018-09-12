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
    matchHists.push_back(ibooker.book1D("hltMatchedRecoEff", "hltMatchedRecoEff", etaBin, etaMin, etaMax));
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



  const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResults);



  for(unsigned int trigIndex = 0; trigIndex < trigNames.size(); trigIndex++){
    std::string testTriggerName = trigNames.triggerName(trigIndex);
    std::cout<<"Trigger Names: "<<testTriggerName<<std::endl;



   }



  std::cout<<"unpacked path names ..."<<std::endl;
  for (pat::TriggerObjectStandAlone object : *triggerObjects) {
      object.unpackPathNames(trigNames);
      const std::vector<std::string>& pathNamesAll = object.pathNames(true);
      for (const auto& pathName : pathNamesAll){
      std::cout<<"Unpacked PathName: "<<pathName<<std::endl;
  
     }
   }

  if(PATmuons.isValid()){
  int numMuons=0;

  for(const auto & patMuon : *PATmuons){
   cout<<"pat Muon Pt: "<<patMuon.pt()<<",   Eta:  "<<patMuon.eta()<<endl; 
   numMuons++;


   



   }
  cout<<"Number of pat muons: "<<numMuons<<endl;

   }//valid pat muon







 
  if(muons.isValid()){//check if the muon handle is valid

         // get the track using only the tracker data
	  matchHists[0]->Fill(0.2);

    }//check if muon handle is valid
}
