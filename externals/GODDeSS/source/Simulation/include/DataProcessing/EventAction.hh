/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef EventAction_h
#define EventAction_h 1

#include <G4UserEventAction.hh>
#include <G4Event.hh>
#include <globals.hh>

#include <fstream>
#include <iostream>

#include <UserEventInformation.hh>
#include <SimulationMessenger.hh>

#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Defines actions to be taken for every event, using functions inherited from G4UserEventAction. (An Event covers everything happening when the particle gun is fired once.) It uses the GODDeSS_DataStorage to obtain the data from the GODDeSS objects.
class EventAction : public G4UserEventAction
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	EventAction( SimulationMessenger * simulationMessenger   /**< class storing and providing variables which are needed in different parts of the simulation. */
		   )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: Messenger(simulationMessenger)
	, RunInformation(Messenger->GetRunInformation())
	, GoddessDataStorage(Messenger->GetGoddessMessenger()->GetDataStorage())
	{
	}

	/**
	 *  Destructor (empty)
	 */
	~EventAction()
	{
	}



	void BeginOfEventAction(const G4Event*);
	void EndOfEventAction(const G4Event*);

private:
	void saveDataToDataFile();
	void saveDataToControlFile();
	void printEventSummary();

	G4int EventID;
	SimulationMessenger * Messenger;
	UserRunInformation * RunInformation;

	UserEventInformation * EventInformation;
	GODDeSS_DataStorage * GoddessDataStorage;
};

#endif
