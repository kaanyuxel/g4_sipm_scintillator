#include "EventAction.hh"

#include <G4SDManager.hh>
#include <G4THitsMap.hh>
#include <G4SystemOfUnits.hh>
#include <G4Event.hh>

#include <G4DigiManager.hh>
#include "digi/G4SipmVoltageTraceDigi.hh"

#include "Analysis.hh"
#include <algorithm>

using namespace std;

void EventAction::EndOfEventAction(const G4Event* event)
{
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();

    G4DigiManager* digiManager = G4DigiManager::GetDMpointer();
    G4DCtable* dcTable = digiManager->GetDCtable();
    
    for(int i = 0; i < dcTable->entries(); i++){
      G4String dmName = dcTable->GetDMname(i);
      G4VDigitizerModule* dm = digiManager->FindDigitizerModule(dmName);
      if(dm){
        dm->Digitize();
      }
    }

    G4DCofThisEvent* dCof = event->GetDCofThisEvent();
    G4int plane;
    
    if(dCof != NULL){
      for(int d = 0; d < dCof->GetNumberOfCollections()/2; d++){
        G4VDigiCollection* dc = dCof->GetDC((d*2)+1);
        if(dc != NULL) {
          if(d == 0)
            plane = 2;
          else
            plane = 1;    
          for(size_t i = 0; i < dc->GetSize(); i++) {
            G4SipmVoltageTraceDigi* digi = (G4SipmVoltageTraceDigi*) dc->GetDigi(i);
            const double tMin = std::max(digi->getTMin(), 0.0);
            const size_t jMin = digi->index(tMin);
            const double tMax = std::min(digi->getTMax(), 1000.0);
            const size_t jMax = std::min(digi->size(), digi->index(tMax));
            for(size_t j = jMin; j < jMax; j++) {
              double time = digi->time(j);
              double voltage = digi->atMeasured(j);
              analysis->FillNtupleDColumn(0,time); 
              analysis->FillNtupleDColumn(1,voltage); 
              analysis->FillNtupleDColumn(2,plane); 
              analysis->FillNtupleDColumn(3,event->GetEventID());
              analysis->AddNtupleRow();
            }
          }    
        } 
      }         
    }   
}
