/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Event.hh>
#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>
#include <G4PhysicalVolumeStore.hh>

#include "PhotonSensitiveDetector.hh"



// class variables begin with capital letters, local variables with small letters



/**
 * The data processing depends on the hitting particle:
 * - mark that the sensitive detector has been hit
 * - <b> if the particle is an optical photons AND has been created outside the sensitive detector: </b>
 *   - save the time and position of the sensitive detector's hit as well as the particle's momentum when hitting the sensitive detector and the name of the sensitive detector's volume that was hit (PhotonDetectorDataStorage)
 *   - write the data needed to resimulate a photon that hit the sensitive detector (start time, start position, start momentum, polarisation) into a text file
 * - stop the particle
 */
G4bool PhotonSensitiveDetector::ProcessHits(G4Step* theStep, G4TouchableHistory*)
{
	G4Track* theTrack = theStep->GetTrack();
	DataStorage->FillEmptyEntriesUpToCurrentTrack(theTrack);

	// if the particle is no optical photon, stop here
	if(theTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
	{
		G4int trackID = theTrack->GetTrackID();
		DataStorage->SetPhotonDetectorWasHit(trackID);

		// Stop track.
		theTrack->SetTrackStatus(fStopAndKill);

		return true;
	}

	// if optical photon was created inside the sensitive detector, stop here
	G4StepPoint * thePreStepPoint = theStep->GetPreStepPoint();
	if(theTrack->GetLogicalVolumeAtVertex() == thePreStepPoint->GetPhysicalVolume()->GetLogicalVolume())
	{
		// Stop track.
		theTrack->SetTrackStatus(fStopAndKill);

		return true;
	}


	G4int trackID = theTrack->GetTrackID();
	DataStorage->SetPhotonDetectorWasHit(trackID);


	G4VPhysicalVolume * volumeThatWasHit = thePreStepPoint->GetPhysicalVolume();
	G4VPhysicalVolume * photonDetectorThatWasHit = G4PhysicalVolumeStore::GetInstance()->GetVolume(volumeThatWasHit->GetName() + "_coating", false);


	// if optical photon was reflected at the surface of the sensitive detector (and if reflected optical photons are not to be counted), stop here (but do not kill it...)
	// but then also the material properties (especially the refractive index) of the sensitive detector volume have to be defined correctly!!!
// 	if(theStep->GetPostStepPoint()->GetPhysicalVolume() != photonDetectorThatWasHit) return true;


	G4String hitFileName = DataStorage->GetPhotonDetectorHitFile();
	hitFileName = hitFileName.substr(0, hitFileName.rfind(".")) + "_" + volumeThatWasHit->GetName() + hitFileName.substr(hitFileName.rfind("."));
	std::ofstream hitFile;
	hitFile.open(hitFileName.c_str(), std::ios_base::out | std::ios_base::app);

	G4double previousStepStartTime = DataStorage->GetPreviousPreStepPointTime();
	hitFile << "t/ns:\t\t" << previousStepStartTime / CLHEP::ns << "\n";

	G4ThreeVector previousStepStartPos = DataStorage->GetPreviousPreStepPointPosition();
	hitFile << "pos/mm:\t\t" << previousStepStartPos / CLHEP::mm << "\n";

	G4ThreeVector hitMomentum = DataStorage->GetPreviousPreStepPointMomentum();
	hitFile << "momentum/MeV:\t" << hitMomentum / CLHEP::MeV << "\n";

	G4ThreeVector previousStepStartPolarisation = DataStorage->GetPreviousPreStepPointPolarisation();
	hitFile << "pol:\t\t" << previousStepStartPolarisation << "\n\n";

	hitFile.close();


	DataStorage->SetNameOfPhotonDetectorThatWasHit(trackID, volumeThatWasHit->GetName());

	G4ThreeVector hitPoint = thePreStepPoint->GetPosition();
	hitPoint -= photonDetectorThatWasHit->GetTranslation();
	hitPoint.transform(photonDetectorThatWasHit->GetObjectRotationValue().inverse());
	DataStorage->SetPhotonDetectorHitPoint(trackID, hitPoint);

	hitMomentum.transform(photonDetectorThatWasHit->GetObjectRotationValue().inverse());
	DataStorage->SetPhotonDetectorHitMomentum(trackID, hitMomentum);

	// Stop track.
	theTrack->SetTrackStatus(fStopAndKill);

	return true;
}
