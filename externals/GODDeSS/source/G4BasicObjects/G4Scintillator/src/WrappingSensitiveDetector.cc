/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Event.hh>
#include <G4EventManager.hh>
// #include <G4OpticalPhoton.hh>

#include "WrappingSensitiveDetector.hh"



// class variables begin with capital letters, local variables with small letters



/**
 * The data processing depends on the hitting particle:
 * - <b> if the particle is a primary particle and the following data has not been saved before: </b> save the time and position of the sensitive detector's hit as well as the particle's momentum when hitting the sensitive detector (ScintillatorTileDataStorage)
 */
G4bool WrappingSensitiveDetector::ProcessHits(G4Step* theStep, G4TouchableHistory*)
{
	G4Track* theTrack = theStep->GetTrack();


	// if the particle is no primary particle, stop here
	if(theTrack->GetParentID() != 0) return true;

	DataStorage->FillEmptyEntriesUpToCurrentTrack(theTrack);


	G4int trackID = theTrack->GetTrackID();

	// if the hit properties have already been set, stop here
	if(! isnan(DataStorage->GetWrappingHitTime(trackID))) return true;



	G4StepPoint * thePreStepPoint  = theStep->GetPreStepPoint();
	G4VPhysicalVolume * volumeThatWasHit = thePreStepPoint->GetPhysicalVolume();

	G4ThreeVector primaryParticleHitPoint = thePreStepPoint->GetPosition();
	primaryParticleHitPoint -= volumeThatWasHit->GetTranslation();
	primaryParticleHitPoint.transform(volumeThatWasHit->GetObjectRotationValue().inverse());
	DataStorage->SetWrappingHitPoint(trackID, primaryParticleHitPoint);

	G4ThreeVector primaryParticleHitMomentum = thePreStepPoint->GetMomentum();
	primaryParticleHitMomentum.transform(volumeThatWasHit->GetObjectRotationValue().inverse());
	DataStorage->SetWrappingHitMomentum(trackID, primaryParticleHitMomentum);

	G4double wrappingHitTime = thePreStepPoint->GetGlobalTime() / CLHEP::ns;
	DataStorage->SetWrappingHitTime(trackID, wrappingHitTime);

	return true;
}
