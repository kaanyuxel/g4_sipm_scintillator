/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Event.hh>
#include <G4RunManager.hh>

#include <SteppingAction.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function that will be called by Geant4 for each step.
 */
void SteppingAction::UserSteppingAction(const G4Step* theStep)
{
	G4Track* theTrack = theStep->GetTrack();

	// if the particle is no optical photon, stop here
	if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return;


	G4bool trackOpticalPhotons = Messenger->GetTrackPhotons();
	if(!trackOpticalPhotons) theTrack->SetTrackStatus(fStopAndKill);

	Messenger->GetGoddessMessenger()->GetDataStorage()->SavePreviousStepOpticalPhotonData(theStep);
}
