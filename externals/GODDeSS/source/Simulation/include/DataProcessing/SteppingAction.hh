/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef SteppingAction_h
#define SteppingAction_h 1

#include <G4OpBoundaryProcess.hh>
#include <G4UserSteppingAction.hh>
#include <globals.hh>

#include <SimulationMessenger.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Defines actions to be taken for every step, using functions inherited from G4UserSteppingAction. (A Step is a part of the full trajectory of a particle, confined by two interactions. The sum of all steps of one particle gives its full trajectory.)
class SteppingAction : public G4UserSteppingAction
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	SteppingAction( SimulationMessenger * simulationMessenger   /**< class storing and providing variables which are needed in different parts of the simulation. */
		      )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: Messenger(simulationMessenger)
	{
	}

	/**
	 *  Destructor (empty)
	 */
	~SteppingAction()
	{
	}



	void UserSteppingAction(const G4Step*);

private:
	SimulationMessenger* Messenger;
};

#endif
