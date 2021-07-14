#include "PhysicsList.hh"

#include <G4EmStandardPhysics.hh>
#include <G4DecayPhysics.hh> 
#include <G4ProductionCutsTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4EmExtraPhysics.hh>
#include <G4HadronPhysicsFTFP_BERT.hh>
#include <G4HadronElasticPhysics.hh>

PhysicsList::PhysicsList()
{
  RegisterPhysics(new G4EmStandardPhysics());
  RegisterPhysics(new G4DecayPhysics());
    
  RegisterPhysics(new G4EmExtraPhysics());

  RegisterPhysics(new G4HadronElasticPhysics());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
}


void PhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts(); 
  DumpCutValuesTable();
}
