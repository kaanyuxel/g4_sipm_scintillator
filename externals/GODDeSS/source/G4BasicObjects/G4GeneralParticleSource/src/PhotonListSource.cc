/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include "PhotonListSource.hh"

#include <G4ParticleTable.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to generate the primary particles.\n
 *  It will be executed by GEANT4, when the simulation is started with the "/run/beamOn" command.\n
 *  If a list with <b> start time </b>, <b> start position </b>, <b> start momentum </b>, and <b> start polarisation </b> of optical photons is provided by PhotonListSourceMessenger, the corresponding optical photons will be generated.
 */
void PhotonListSource::GeneratePrimaries(G4Event* event)
{
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	G4String particleName = "opticalphoton";
	ParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName));

	ParticleGun->SetNumberOfParticles(1);

	// the read in data have to be processed reversely to get the vertices processed in the correct order
	while(Messenger->IsPhotonDataAvailable())
	{
		// get properties of the next photon
		G4double start_time_ns = Messenger->GetNextPhotonStartTime_ns();
		G4ThreeVector start_position_mm = Messenger->GetNextPhotonStartPosition_mm();
		G4ThreeVector start_momentum_MeV = Messenger->GetNextPhotonStartMomentum_MeV();
		G4ThreeVector start_polarisation = Messenger->GetNextPhotonStartPolarisation();

		G4double momentum = sqrt(start_momentum_MeV[0] * start_momentum_MeV[0] + start_momentum_MeV[1] * start_momentum_MeV[1] + start_momentum_MeV[2] * start_momentum_MeV[2]);
		G4ThreeVector start_momentumDirection = start_momentum_MeV;
		start_momentumDirection[0] /= momentum;
		start_momentumDirection[1] /= momentum;
		start_momentumDirection[2] /= momentum;

		// remove the current photon from the messenger
		Messenger->RemoveNextPhotonData();

		// set properties
		ParticleGun->SetParticleTime(start_time_ns);
		ParticleGun->SetParticlePosition(start_position_mm);
		ParticleGun->SetParticleEnergy(momentum);
		ParticleGun->SetParticleMomentumDirection(start_momentumDirection);
		ParticleGun->SetParticlePolarization(start_polarisation);

		// Fire!
		ParticleGun->GeneratePrimaryVertex(event);
	}

	//remove the processed event
	if(Messenger->IsEventAvailable()) Messenger->RemoveEvent();
}
