#include "RunAction.hh"
#include <G4Gamma.hh>
#include <G4Electron.hh>
#include <G4AccumulableManager.hh>
#include <G4SystemOfUnits.hh>
#include "Analysis.hh"

RunAction::RunAction() :
  G4UserRunAction()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->CreateNtuple("test", "test");
  analysisManager->CreateNtupleDColumn("time");
  analysisManager->CreateNtupleDColumn("voltage");
  analysisManager->CreateNtupleDColumn("plane");
  analysisManager->CreateNtupleDColumn("eventID");
  analysisManager->FinishNtuple();
}


void RunAction::BeginOfRunAction(const G4Run*)
{
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  auto analysisManager = G4AnalysisManager::Instance();

  G4String fileName = "time";
  analysisManager->OpenFile(fileName);  
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();

  if (nofEvents == 0) return;

  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

}

RunAction::~RunAction()
{
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}


