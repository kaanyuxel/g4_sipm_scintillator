/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4ProcessTable.hh>

#include <PhysicsList.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to register the particles to Geant4:
 *  - register the particles of all physical processes defined in the constructor to Geant4
 *  - displays the registered particles
 */
void PhysicsList::ConstructParticle()   // will be called in the initialisation process
{
	// register the particles of all physical processes defined in the constructor:
	G4VModularPhysicsList::ConstructParticle();

	G4VModularPhysicsList::DumpList();
}

/**
 *  Function to register the physical processes to Geant4:
 *  - register the physical processes defined in the constructor to Geant4
 *  - displays the registered processes
 */
void PhysicsList::ConstructProcess()   // will be called in the initialisation process
{
	// register the physical processes defined in the constructor:
	G4VModularPhysicsList::ConstructProcess();

	// dump list of the registered processes:
	G4cout << G4endl << "processNameList:" << G4endl;
	G4ProcessTable::G4ProcNameVector* processNameList = G4ProcessTable::GetProcessTable()->GetNameList();
	for (G4int i=0; i<G4int(processNameList->size()); i++){
		G4cout << "   " << (*processNameList)[i] << G4endl;
	}
	G4cout << G4endl;
}

/**
 *  Function to set the cuts in Geant4:
 *  - sets the default cuts of Geant4
 *  - displays the cut values
 */
void PhysicsList::SetCuts()   // will be called in the initialisation process
{
	SetCutsWithDefault();
	if (this->verboseLevel > 0) {
		DumpCutValuesTable();
	}
}
