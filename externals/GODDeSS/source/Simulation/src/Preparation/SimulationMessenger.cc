/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <boost/regex.hpp>

#include <SimulationMessenger.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to define the action of the commands.\n
 *  <b> It will be executed by GEANT4, when a command (created by this class) is called in the interactive mode. </b>
 */
void SimulationMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == TrackPhotonsCmd) TrackPhotons = TrackPhotonsCmd->GetNewBoolValue(newValue);
	else if(command == OutputFileCmd)
	{
		DataFile = newValue;
		ControlFile = newValue;
	}
}
