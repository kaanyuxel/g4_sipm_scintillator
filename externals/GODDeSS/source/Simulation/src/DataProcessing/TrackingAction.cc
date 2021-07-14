/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4Track.hh>
#include <G4OpticalPhoton.hh>
#include <G4ParticleTable.hh>

#include <TrackingAction.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function that will be called by Geant4 at the beginning of each track.
 */
void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	// Visualisation: create trajectory only for primaries and some secundaries
	G4int trackID = aTrack->GetTrackID();
	G4int parentTrackID = aTrack->GetParentID();
	G4String partName = aTrack->GetParticleDefinition()->GetParticleName(); // one cannot use the particle ID, as it is NOT biunique ("0" is used for several Geant4 specific particles (optical photon, geantino,...)

	if(   parentTrackID == 0
	   || partName != G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName()
	   || ((trackID % 2500) >= 0 && (trackID % 2500) < 25)
	  )
	{
		fpTrackingManager->SetStoreTrajectory(true);
	}
	else
	{
		fpTrackingManager->SetStoreTrajectory(false);
	}

	// collecting data:
	// save "global" properties of the track
	EventInformation = (UserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();
// 	TrackIDExistsInMap = EventInformation->TrackIDHasEntries(trackID);   // do this before anything of the current G4Track has been saved to eventInformation (for suppressing doublecounting of particles with more than one G4Track)
	EventInformation->FillEmptyEntriesUpToTrackID(trackID);

	EventInformation->SetParticleName(trackID, partName);

	G4int partID = G4ParticleTable::GetParticleTable()->FindParticle(partName)->GetPDGEncoding();
	EventInformation->SetParticleID(trackID, partID);

	if(aTrack->GetCurrentStepNumber() == 0)
	{
		G4double creationTime = aTrack->GetGlobalTime();
		EventInformation->SetGlobalCreationTime(trackID, creationTime);
	}

	G4double partMass = aTrack->GetParticleDefinition()->GetPDGMass();
	RunInformation->SetParticleMass(partName, partMass);

	G4ThreeVector initialPosition = aTrack->GetVertexPosition();
	EventInformation->SetInitialPosition(trackID, initialPosition);

	G4ThreeVector initialMomentumDirection = aTrack->GetVertexMomentumDirection();
	G4double initEnergy = aTrack->GetVertexKineticEnergy();
	G4ThreeVector initialMomentum = initialMomentumDirection * sqrt( pow(initEnergy, 2) + 2 * initEnergy * partMass );
	EventInformation->SetInitialMomentum(trackID, initialMomentum);

	if(parentTrackID == 0)
	{
		EventInformation->SetParticleIsPrimary(trackID, true);
	}
	else
	{
		EventInformation->SetParticleIsPrimary(trackID, false);

		EventInformation->SetParent(trackID, parentTrackID);
	}
}



/**
 *  Function that will be called by Geant4 at the end of each track.
 */
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
	G4int trackID = aTrack->GetTrackID();
	G4int parentTrackID = aTrack->GetParentID();

	if(parentTrackID != 0)
	{
		G4String productionMechanism = aTrack->GetCreatorProcess()->GetProcessName();
		EventInformation->SetProductionMechanism(trackID, productionMechanism);
	}

	// look only at optical photons
	if(aTrack->GetParticleDefinition()->GetParticleName() == G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName())
	{
// 		if(!TrackIDExistsInMap)   // suppress doublecounting of particles with more than one G4Track
// 		{
			G4String productionMechanism = EventInformation->GetProductionMechanism(trackID);
			G4bool IsScintiPhoton = (productionMechanism == "Scintillation");
			G4bool IsCerenkovPhoton = (productionMechanism == "Cerenkov");
			G4bool IsWLSPhoton = (productionMechanism == "OpWLS");


			G4int numberOfPrimaryParticles = EventInformation->GetNumberOfPrimaryParticles();

			if(IsScintiPhoton)
			{
				EventInformation->IncScintiPhotonCount();

				if(parentTrackID <= numberOfPrimaryParticles) EventInformation->IncScintiPhotonByPrimaryCount();
				else EventInformation->IncScintiPhotonBySecondaryCount();
			}
			else if(IsCerenkovPhoton)
			{
				EventInformation->IncCerenkovPhotonCount();

				if(parentTrackID <= numberOfPrimaryParticles) EventInformation->IncCerenkovPhotonByPrimaryCount();
				else EventInformation->IncCerenkovPhotonBySecondaryCount();
			}
			else if(IsWLSPhoton)
			{
				EventInformation->IncWLSPhotonCount();
			}

			G4StepPoint * thePostStepPoint = aTrack->GetStep()->GetPostStepPoint();

			// find absorption process
			if(thePostStepPoint->GetPhysicalVolume())   //the particle did not leave the world volume
			{
				if(isnan(EventInformation->GetGlobalAbsorptionTime(trackID)))
				{
					G4double globalPhotonAbsorptionTime = thePostStepPoint->GetGlobalTime();
					EventInformation->SetGlobalAbsorptionTime(trackID, globalPhotonAbsorptionTime);
				}

				if(thePostStepPoint->GetProcessDefinedStep()->GetProcessName() == "OpWLS")
				{
					EventInformation->SetAbsorbedIn(trackID, "WLS");

					// count the optical photons absorbed in fibre
					EventInformation->IncAbsorbedInFibreCount();

					if(IsScintiPhoton) EventInformation->IncAbsorbedScintiPhotonInFibreCount();
					else if(IsCerenkovPhoton) EventInformation->IncAbsorbedCerenkovPhotonInFibreCount();
					else if(IsWLSPhoton) EventInformation->IncAbsorbedWLSPhotonInFibreCount();
				}
				else if(Messenger->GetGoddessMessenger()->GetDataStorage()->PhotonDetectorWasHit(trackID))
				{
					G4String nameOfVolumeThatWasHit = Messenger->GetGoddessMessenger()->GetDataStorage()->GetNameOfPhotonDetectorThatWasHit(trackID);
					EventInformation->SetAbsorbedIn(trackID, nameOfVolumeThatWasHit);

					// count the optical photons absorbed in SiPMs
					EventInformation->IncAbsorbedInPhotonDetectorCount();

					if(IsScintiPhoton) EventInformation->IncAbsorbedScintiPhotonInPhotonDetectorCount();
					else if(IsCerenkovPhoton) EventInformation->IncAbsorbedCerenkovPhotonInPhotonDetectorCount();
					else if(IsWLSPhoton) EventInformation->IncAbsorbedWLSPhotonInPhotonDetectorCount();
				}
				else
				{
					EventInformation->SetAbsorbedIn(trackID, "scintiTile");

					// count the optical photons absorbed in scintillator, wrapping, optical cement,...
					EventInformation->IncAbsorbedInScintiCount();

					if(IsScintiPhoton) EventInformation->IncAbsorbedScintiPhotonInScintiCount();
					else if(IsCerenkovPhoton) EventInformation->IncAbsorbedCerenkovPhotonInScintiCount();
					else if(IsWLSPhoton) EventInformation->IncAbsorbedWLSPhotonInScintiCount();
				}
			}
// 		}
	}
}
