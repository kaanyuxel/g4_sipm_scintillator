/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef RunAction_h
#define RunAction_h 1

#include <G4UserRunAction.hh>
#include <G4Timer.hh>
#include <G4Run.hh>
#include <globals.hh>

#include <UserRunInformation.hh>
#include <SimulationMessenger.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Defines actions to be taken for every run, using functions inherited from G4UserRunAction. (A Run is started by /run/beamOn and may consist of one or more events.)
class RunAction : public G4UserRunAction
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 *  - creates an object of the UserRunInformation and the G4Timer classes
	 */
	RunAction( SimulationMessenger * simulationMessenger   /**< class storing and providing variables which are needed in different parts of the simulation. */
		 )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: Messenger(simulationMessenger)
	, RunInformation(new UserRunInformation())
	, Timer(new G4Timer)
	{
		Messenger->SetRunInformation(RunInformation);
	}

	/**
	 *  Destructor:
	 *  - deletes the objects which have been created in the constructor (RunAction())
	 */
	~RunAction()
	{
		delete RunInformation;
		delete Timer;
	}



	void BeginOfRunAction(const G4Run* aRun);
	void EndOfRunAction(const G4Run* aRun);

private:
	template <class T> void FindMinMaxValues(T newValue, T & variable_min, T & variable_max) {
		if(!variable_min && !variable_max)
		{
			variable_min = newValue;
			variable_max = newValue;
		}
		else if(variable_min > newValue) variable_min = newValue;
		else if(variable_max < newValue) variable_max = newValue;
	}

	void saveGPSSummary(G4int runID);
	void saveDataSummary(G4int runID);
	void saveControlSummary(G4int runID);

	SimulationMessenger* Messenger;
	UserRunInformation* RunInformation;
	G4Timer* Timer;
};

#endif
