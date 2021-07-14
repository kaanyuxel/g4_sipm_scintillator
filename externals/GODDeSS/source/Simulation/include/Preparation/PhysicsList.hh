/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef OPTICALPHYSICSLIST_HH_
#define OPTICALPHYSICSLIST_HH_

#include <G4VModularPhysicsList.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Registers particles and physics processes, using functions inherited from G4VModularPhysicsList.
class PhysicsList: public G4VModularPhysicsList
{
public:

	/**
	 *  Constructor:
	 *  - passes the verbosity level to Geant4
	 */
	PhysicsList( int verbose = 0   /**< verbosity level for the physics list (default: 0) */
		   )
	{
		// specify the particles and physics processes to be created in the initialisation process:
		//
		// --- EMPTY ---
		//


		SetVerboseLevel(verbose);
	}

	/**
	 *  Destructor (empty)
	 */
	virtual ~PhysicsList()
	{
	}



	// the following function definitions are only needed to overwrite the corresponding G4VModularPhysicsList functions (if this is to be done, e.g. for dumping informations):
	virtual void ConstructParticle();
	virtual void ConstructProcess();
	virtual void SetCuts();
};

#endif
