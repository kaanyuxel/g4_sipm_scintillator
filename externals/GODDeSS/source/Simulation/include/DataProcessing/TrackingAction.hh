/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef TrackingAction_h
#define TrackingAction_h 1

#include <G4UserTrackingAction.hh>
#include <globals.hh>

#include <SimulationMessenger.hh>
#include <UserEventInformation.hh>
#include <UserRunInformation.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Defines actions to be taken for every track, using functions inherited from G4UserTrackingAction. (A Track is the full trajectory of a particle, it can be divided into steps (which is the part of the trajectory between two interactions).)
class TrackingAction : public G4UserTrackingAction
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	TrackingAction( SimulationMessenger * simulationMessenger   /**< class storing and providing variables which are needed in different parts of the simulation. */
		      )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: Messenger(simulationMessenger)
	, RunInformation(simulationMessenger->GetRunInformation())
	{
	}

	/**
	 *  Destructor (empty)
	 */
	virtual ~TrackingAction()
	{
	}



	virtual void PreUserTrackingAction(const G4Track*);
	virtual void PostUserTrackingAction(const G4Track*);

private:
	SimulationMessenger * Messenger;
	UserRunInformation* RunInformation;
	UserEventInformation* EventInformation;
};

#endif
