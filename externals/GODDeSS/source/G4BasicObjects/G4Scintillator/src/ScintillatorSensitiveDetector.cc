/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>

#include "ScintillatorSensitiveDetector.hh"



// class variables begin with capital letters, local variables with small letters



/**
 * The data processing depends on the hitting particle:
 * - <b> if the particle is an optical photons: </b> do nothing
 * - <b> else: </b>
 *   - add up the energy deposit (ScintillatorTileDataStorage)
 *   - <b> if the particle is a primary particle and the following data has not been saved before: </b> save the time and position of the sensitive detector's hit as well as the particle's momentum when hitting the sensitive detector (ScintillatorTileDataStorage)
 */
G4bool ScintillatorSensitiveDetector::ProcessHits(G4Step* theStep, G4TouchableHistory*)
{
	G4Track* theTrack = theStep->GetTrack();


	// if the particle is an optical photon and no primary particle, do nothing
	if(theTrack->GetParentID() != 0 && theTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) return true;


	DataStorage->FillEmptyEntriesUpToCurrentTrack(theTrack);


	G4int trackID = theTrack->GetTrackID();

	if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
	{
		G4double energyDeposit = theStep->GetTotalEnergyDeposit();
		// GetDeltaEnergy()						(NEGATIVE!!!) total energy lost by a particle (to the medium and to secondaries)
		// GetTotalEnergyDeposit()					energy lost to the medium (e.g. by scintillation or delta ray production, bremsstrahlung below the energy cut-off), NOT to secondaries
		// GetTotalEnergyDeposit() - GetNonIonizingEnergyDeposit()	energy available for ionization
		// http://hypernews.slac.stanford.edu/HyperNews/geant4/get/eventtrackmanage/1043.html

		DataStorage->AddUpEnergyDepositionInScintillator(trackID, energyDeposit);
	}


	// if the particle is no primary particle, stop here
	if(theTrack->GetParentID() != 0) return true;


	G4double stepLength = theStep->GetStepLength();
	DataStorage->AddUpPathLengthInScintillator(trackID, stepLength);


	// if the hit properties have already been set, stop here
	if(! isnan(DataStorage->GetScintillatorHitTime(trackID))) return true;


	G4StepPoint * thePreStepPoint = theStep->GetPreStepPoint();
	G4VPhysicalVolume * volumeThatWasHit = thePreStepPoint->GetPhysicalVolume();

	G4ThreeVector primaryParticleHitPoint = thePreStepPoint->GetPosition();
	primaryParticleHitPoint -= volumeThatWasHit->GetTranslation();
	primaryParticleHitPoint.transform(volumeThatWasHit->GetObjectRotationValue().inverse());
	DataStorage->SetScintillatorHitPoint(trackID, primaryParticleHitPoint);

	G4ThreeVector primaryParticleHitMomentum = thePreStepPoint->GetMomentum();
	primaryParticleHitMomentum.transform(volumeThatWasHit->GetObjectRotationValue().inverse());
	DataStorage->SetScintillatorHitMomentum(trackID, primaryParticleHitMomentum);

	G4double scintiHitTime = thePreStepPoint->GetGlobalTime() / CLHEP::ns;
	DataStorage->SetScintillatorHitTime(trackID, scintiHitTime);

	return true;
}
