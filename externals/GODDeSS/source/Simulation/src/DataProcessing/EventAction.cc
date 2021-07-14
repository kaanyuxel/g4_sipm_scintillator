/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4EventManager.hh>
#include <G4OpticalPhoton.hh>

#include <UserRunInformation.hh>

#include <EventAction.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function that will be called by Geant4 at the beginning of each event.
 */
void EventAction::BeginOfEventAction(const G4Event* theEvent)
{
	EventID = theEvent->GetEventID();
	G4cout << G4endl << "-----------------------------------------------------------------------------" << G4endl;
	G4cout << "EventID = " << EventID << G4endl << G4endl;

	Messenger->GetGoddessMessenger()->GetPhotonDetectorConstructor()->WriteEventIDToHitFile(EventID);

	G4EventManager::GetEventManager()->SetUserInformation(new UserEventInformation);

	EventInformation = (UserEventInformation*)theEvent->GetUserInformation();
	EventInformation->clean();
	EventInformation->SetNumberOfPrimaryParticles(theEvent->GetNumberOfPrimaryVertex());

	GoddessDataStorage->clean();
}

/**
 *  Function that will be called by Geant4 at the end of each event.
 */
void EventAction::EndOfEventAction(const G4Event* theEvent)
{
	EventInformation = (UserEventInformation*)theEvent->GetUserInformation();

	saveDataToDataFile();

	saveDataToControlFile();

	printEventSummary();
}



void EventAction::saveDataToDataFile()
{
	G4String dataFileName = Messenger->GetDataFileName();
	std::ofstream dataFile;
	dataFile.open(dataFileName.c_str(), std::ios_base::out | std::ios_base::app);
	dataFile << "EventID:\t\t" << EventID << "\n\n";

	G4int CerenkovPhotonsFromPrimary = 0;
	G4int scintiPhotonsFromPrimary = 0;
	G4int WLSPhotonsFromPrimary = 0;
	G4int CerenkovPhotonsFromSecondary = 0;
	G4int scintiPhotonsFromSecondary = 0;
	G4int WLSPhotonsFromSecondary = 0;
	G4double deltaE = NAN;
	G4double globalHitTime = NAN;

	G4int numParticles = EventInformation->GetNumberOfParticles();
	for (int iter = 1; iter < numParticles + 1; iter++)
	{
		G4String particleName = EventInformation->GetParticleName(iter);
		G4bool isPrimary = EventInformation->GetParticleIsPrimary(iter);

		if(isPrimary)
		{
			dataFile << "primaryParticleName:\t\t\t" << particleName << "\n";

			globalHitTime = GoddessDataStorage->GetScintillatorHitTime(iter);
			dataFile << "t_hit/ns:\t\t\t\t" << globalHitTime / CLHEP::ns << "\n";

			G4ThreeVector scintiHitPoint = GoddessDataStorage->GetScintillatorHitPoint(iter);
			dataFile << "pos_hit/mm:\t\t\t\t" << scintiHitPoint / CLHEP::mm << "\n";
			RunInformation->SetMinMaxValue_ScintiHitPointX_mm(scintiHitPoint.x() / CLHEP::mm);
			RunInformation->SetMinMaxValue_ScintiHitPointY_mm(scintiHitPoint.y() / CLHEP::mm);
			RunInformation->SetMinMaxValue_ScintiHitPointZ_mm(scintiHitPoint.z() / CLHEP::mm);

			G4ThreeVector scintiHitMomentum = GoddessDataStorage->GetScintillatorHitMomentum(iter);
			dataFile << "momentum_hit/MeV:\t\t\t" << scintiHitMomentum / CLHEP::MeV << "\n";

			deltaE = GoddessDataStorage->GetEnergyDepositionInScintillator(iter);
			dataFile << "E_depos/MeV:\t\t\t\t" << deltaE / CLHEP::MeV << "\n\n";
		}
		else if(particleName == G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName())
		{
			G4bool isParentPrimary = EventInformation->GetParentIsPrimary(iter);
			G4String productionMechanism = EventInformation->GetProductionMechanism(iter);

			if(isParentPrimary)
			{
				if(productionMechanism == "Cerenkov") CerenkovPhotonsFromPrimary++;
				if(productionMechanism == "Scintillation") scintiPhotonsFromPrimary++;
				if(productionMechanism == "OpWLS") WLSPhotonsFromPrimary++;
			}
			else
			{
				if(productionMechanism == "Cerenkov") CerenkovPhotonsFromSecondary++;
				if(productionMechanism == "Scintillation") scintiPhotonsFromSecondary++;
				if(productionMechanism == "OpWLS") WLSPhotonsFromSecondary++;
			}
		}
	}

	G4int num_opticalPhotons = CerenkovPhotonsFromPrimary + scintiPhotonsFromPrimary + WLSPhotonsFromPrimary + CerenkovPhotonsFromSecondary + scintiPhotonsFromSecondary + WLSPhotonsFromSecondary;
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(num_opticalPhotons);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition( ((G4double) num_opticalPhotons) / (deltaE / CLHEP::MeV)  );

	dataFile << "Cerenkov_photons_from_primary:\t\t" << CerenkovPhotonsFromPrimary << "\n";
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(CerenkovPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition(((G4double) CerenkovPhotonsFromPrimary) / (deltaE / CLHEP::MeV));

	dataFile << "scinti_photons_from_primary:\t\t" << scintiPhotonsFromPrimary << "\n";
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(scintiPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition(((G4double) scintiPhotonsFromPrimary) / (deltaE / CLHEP::MeV));

	dataFile << "WLS_photons_from_primary:\t\t" << WLSPhotonsFromPrimary << "\n";
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(WLSPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition(((G4double) WLSPhotonsFromPrimary) / (deltaE / CLHEP::MeV));

	dataFile << "Cerenkov_photons_from_secondary:\t" << CerenkovPhotonsFromSecondary << "\n";
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(CerenkovPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition(((G4double) CerenkovPhotonsFromSecondary) / (deltaE / CLHEP::MeV));

	dataFile << "scinti_photons_from_secondary:\t\t" << scintiPhotonsFromSecondary << "\n";
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(scintiPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition(((G4double) scintiPhotonsFromSecondary) / (deltaE / CLHEP::MeV));

	dataFile << "WLS_photons_from_secondary:\t\t" << WLSPhotonsFromSecondary << "\n";
	RunInformation->SetMinMaxValue_opticalPhotonsPerPrimary(WLSPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_opticalPhotonsPerEnergyDeposition(((G4double) WLSPhotonsFromSecondary) / (deltaE / CLHEP::MeV));

	G4int num_AbsorbedScinti_inPhotonDetector = EventInformation->GetAbsorbedScintiPhotonInPhotonDetectorCount();
	G4int num_AbsorbedCer_inPhotonDetector    = EventInformation->GetAbsorbedCerenkovPhotonInPhotonDetectorCount();
	G4int num_AbsorbedWLS_inPhotonDetector    = EventInformation->GetAbsorbedWLSPhotonInPhotonDetectorCount();

	RunInformation->SetMinMaxValue_opticalPhotonsAbsorbedInPhotonDetector(num_AbsorbedScinti_inPhotonDetector + num_AbsorbedCer_inPhotonDetector + num_AbsorbedWLS_inPhotonDetector);

// 	dataFile << "Cerenkov_photons_absorbed_in_SiPM:\t" << num_AbsorbedScinti_inPhotonDetector << "\n";
// 	dataFile << "scinti_photons_absorbed_in_SiPM:\t" << num_AbsorbedCer_inPhotonDetector << "\n";
// 	dataFile << "WLS_photons_absorbed_in_SiPM:\t\t" << num_AbsorbedWLS_inPhotonDetector << "\n";

	for (int iter = 1; iter < numParticles + 1; iter++)
	{
		if(EventInformation->GetParticleName(iter) == G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName() && GoddessDataStorage->PhotonDetectorWasHit(iter))
		{
			dataFile << "\n";

			G4String absorbedIn = GoddessDataStorage->GetNameOfPhotonDetectorThatWasHit(iter);
			dataFile << "abs_vol:\t\t" << absorbedIn << "\n";

			G4double globalAbsorptionTime = EventInformation->GetGlobalAbsorptionTime(iter);
			dataFile << "t_abs/ns:\t\t" << globalAbsorptionTime / CLHEP::ns << "\n";
			if(!isnan(globalHitTime)) RunInformation->SetMinMaxValue_DeltaTimeHitAbsorption_ns(globalAbsorptionTime / CLHEP::ns - globalHitTime / CLHEP::ns);

			G4ThreeVector PhotonDetectorHitPoint = GoddessDataStorage->GetPhotonDetectorHitPoint(iter);
			dataFile << "pos_SiPM_hit/mm:\t" << PhotonDetectorHitPoint / CLHEP::mm << "\n";
			RunInformation->SetMinMaxValue_PhotonDetectorHitPointX_mm(PhotonDetectorHitPoint.x() / CLHEP::mm);
			RunInformation->SetMinMaxValue_PhotonDetectorHitPointY_mm(PhotonDetectorHitPoint.y() / CLHEP::mm);
// 			RunInformation->SetMinMaxValue_PhotonDetectorHitPointZ_mm(PhotonDetectorHitPoint.z() / CLHEP::mm);

			G4ThreeVector PhotonDetectorHitMomentum = GoddessDataStorage->GetPhotonDetectorHitMomentum(iter);
			dataFile << "SiPM_hit_momentum/MeV:\t" << PhotonDetectorHitMomentum / CLHEP::MeV << "\n";
			RunInformation->SetMinMaxValue_OpticalPhotonHitEnergy_MeV(PhotonDetectorHitMomentum.mag() / CLHEP::MeV);
		}
	}

	dataFile << "\n";
	dataFile.close();
}



void EventAction::saveDataToControlFile()
{
	G4int mumberOfEventsToBeSaved = (G4int) ceil(RunInformation->GetNumberOfEvents() / 10.);
	if(mumberOfEventsToBeSaved % 2) mumberOfEventsToBeSaved++;

	G4String controlFileName = Messenger->GetControlFileName();
	std::ofstream controlFile;
	controlFile.open(controlFileName.c_str(), std::ios_base::out | std::ios_base::app);
	if(EventID < mumberOfEventsToBeSaved)       controlFile << "\nEventID:\t\t" << EventID << "\n\n";
	else if(EventID == mumberOfEventsToBeSaved) controlFile << "\n##### from this event on, only the event summary is saved (only 10% of the events are fully saved) #####\n\n";
	else                                        controlFile << "\n";

	G4double control_deltaE_primary_MeV = 0;
	G4double scintiPathLength_primary = 0;
	G4int control_CerenkovPhotonsFromPrimary = 0;
	G4int control_scintiPhotonsFromPrimary = 0;
	G4int control_WLSPhotonsFromPrimary = 0;
	G4int control_CerenkovPhotonsFromSecondary = 0;
	G4int control_scintiPhotonsFromSecondary = 0;
	G4int control_WLSPhotonsFromSecondary = 0;

	G4String name_primary = "";
	G4int ID_primary = 0;
	G4ThreeVector initialPosition_primary;
	G4ThreeVector initialMomentum_primary;
	G4ThreeVector scintiHitPoint_primary;
	G4ThreeVector scintiHitMomentum_primary;
	G4double globalHitTime = -1.;
	std::map<G4String, G4int> names_secondary;
	std::map<G4String, G4int> names_processes;

	G4int numParticles = EventInformation->GetNumberOfParticles();
	for (int iter = 1; iter < numParticles + 1; iter++)
	{
		G4String particleName = EventInformation->GetParticleName(iter);
		G4bool isPrimary = EventInformation->GetParticleIsPrimary(iter);

		if(EventID < mumberOfEventsToBeSaved)
		{
			controlFile << "primary:\t\t" << isPrimary << "\n";   /// This has to be the first entry!!!
			controlFile << "name:\t\t\t" << particleName << "\n";   /// This has to come before productionMechanism and initialEnergy!!!
		}

		if(isPrimary)
		{
			name_primary = particleName;
		}
		else
		{
			try
			{
				G4int number = names_secondary.at(particleName);
				number++;
				names_secondary[particleName] = number;
			}
			catch (...)
			{
				names_secondary[particleName] = 1;
			}
		}

		G4int particleID = EventInformation->GetParticleID(iter);
		if(isPrimary)
		{
			if(!isnan(particleID))
			{
				RunInformation->SetMinMaxValue_Control_PrimaryParticleID(particleID);
				ID_primary = particleID;
			}

			if(EventID < mumberOfEventsToBeSaved) controlFile << "ID:\t\t\t" << particleID << "\n";
		}
		else if(EventID < mumberOfEventsToBeSaved)
		{
			controlFile << "ID:\t\t\t" << particleID << "\n";
			if(!isnan(particleID)) RunInformation->SetMinMaxValue_Control_SecondaryParticleID(particleID);
		}

		if(!isPrimary)
		{
			G4bool isParentPrimary = EventInformation->GetParentIsPrimary(iter);
			if(EventID < mumberOfEventsToBeSaved) controlFile << "primary_p:\t\t" << isParentPrimary << "\n";

			G4String productionMechanism = EventInformation->GetProductionMechanism(iter);
			if(EventID < mumberOfEventsToBeSaved) controlFile << "prod:\t\t\t" << productionMechanism << "\n";

			try
			{
				G4int number = names_processes.at(productionMechanism);
				number++;
				names_processes[productionMechanism] = number;
			}
			catch (...)
			{
				names_processes[productionMechanism] = 1;
			}

			if(particleName == G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName())
			{
				if(isParentPrimary)
				{
					if(productionMechanism == "Cerenkov") control_CerenkovPhotonsFromPrimary++;
					if(productionMechanism == "Scintillation") control_scintiPhotonsFromPrimary++;
					if(productionMechanism == "OpWLS") control_WLSPhotonsFromPrimary++;
				}
				else
				{
					if(productionMechanism == "Cerenkov") control_CerenkovPhotonsFromSecondary++;
					if(productionMechanism == "Scintillation") control_scintiPhotonsFromSecondary++;
					if(productionMechanism == "OpWLS") control_WLSPhotonsFromSecondary++;
				}
			}
		}

		if(isPrimary)
		{
			G4ThreeVector initialPosition = EventInformation->GetInitialPosition(iter);
			initialPosition_primary = initialPosition;
			RunInformation->SetMinMaxValue_Control_PrimaryParticleInitialPositionX_mm(initialPosition.x() / CLHEP::mm);
			RunInformation->SetMinMaxValue_Control_PrimaryParticleInitialPositionY_mm(initialPosition.y() / CLHEP::mm);
			RunInformation->SetMinMaxValue_Control_PrimaryParticleInitialPositionZ_mm(initialPosition.z() / CLHEP::mm);
			if(EventID < mumberOfEventsToBeSaved) controlFile << "pos/mm:\t\t\t" << initialPosition / CLHEP::mm << "\n";

			G4ThreeVector initialMomentum = EventInformation->GetInitialMomentum(iter);
			initialMomentum_primary = initialMomentum;
			if(EventID < mumberOfEventsToBeSaved) controlFile << "momentum/MeV:\t\t" << initialMomentum / CLHEP::MeV << "\n";

			G4double particleMass = RunInformation->GetParticleMass(particleName);
			G4double initialEnergy = sqrt( pow(particleMass, 2) + pow(initialMomentum.mag(), 2) ) - particleMass;
			RunInformation->SetMinMaxValue_Control_PrimaryParticleInitialEnergy_MeV(initialEnergy / CLHEP::MeV);
			if(particleMass)
			{
				G4double initialBetaGamma = initialMomentum.mag() / particleMass;
				RunInformation->SetMinMaxValue_Control_PrimaryParticleInitialBetaGamma(initialBetaGamma);
			}

			globalHitTime = GoddessDataStorage->GetScintillatorHitTime(iter);
			if(!isnan(globalHitTime)) RunInformation->SetMinMaxValue_Control_GlobalScintiHitTime_ns(globalHitTime / CLHEP::ns);
			if(EventID < mumberOfEventsToBeSaved) controlFile << "t_hit/ns:\t\t" << globalHitTime / CLHEP::ns << "\n";

			G4ThreeVector scintiHitPoint = GoddessDataStorage->GetScintillatorHitPoint(iter);
			scintiHitPoint_primary = scintiHitPoint;
			RunInformation->SetMinMaxValue_Control_ScintiHitPointX_mm(scintiHitPoint.x() / CLHEP::mm);
			RunInformation->SetMinMaxValue_Control_ScintiHitPointY_mm(scintiHitPoint.y() / CLHEP::mm);
			RunInformation->SetMinMaxValue_Control_ScintiHitPointZ_mm(scintiHitPoint.z() / CLHEP::mm);
			if(EventID < mumberOfEventsToBeSaved) controlFile << "hit_pos/mm:\t\t" << scintiHitPoint / CLHEP::mm << "\n";

			G4ThreeVector scintiHitMomentum = GoddessDataStorage->GetScintillatorHitMomentum(iter);
			scintiHitMomentum_primary = scintiHitMomentum;
			if(EventID < mumberOfEventsToBeSaved) controlFile << "hit_momentum/MeV:\t" << scintiHitMomentum / CLHEP::MeV << "\n";

			if(particleName != G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName())
			{
				G4double deltaE = GoddessDataStorage->GetEnergyDepositionInScintillator(iter);
				RunInformation->SetMinMaxValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV(deltaE / CLHEP::MeV);
				control_deltaE_primary_MeV += deltaE;
				if(EventID < mumberOfEventsToBeSaved) controlFile << "E_depos/MeV:\t\t" << deltaE / CLHEP::MeV << "\n";
			}

			G4double deltaX = GoddessDataStorage->GetPathLengthInScintillator(iter);
			RunInformation->SetMinMaxValue_Control_PrimaryParticlePathLengthInScintillator_MeV(deltaX / CLHEP::mm);
			scintiPathLength_primary = deltaX;
			if(EventID < mumberOfEventsToBeSaved) controlFile << "DeltaX/mm:\t\t" << deltaX / CLHEP::mm << "\n";

			if(EventID < mumberOfEventsToBeSaved) controlFile << "\n";
		}
		else if(EventID < mumberOfEventsToBeSaved)
		{
			G4double globalCreationTime = -1;
			globalCreationTime = EventInformation->GetGlobalCreationTime(iter);
			controlFile << "t_crea/ns:\t\t" << globalCreationTime / CLHEP::ns << "\n";

			G4ThreeVector initialMomentum = EventInformation->GetInitialMomentum(iter);
			controlFile << "momentum/MeV:\t\t" << initialMomentum / CLHEP::MeV << "\n";

			G4double particleMass = RunInformation->GetParticleMass(particleName);
			G4double initialEnergy = sqrt( pow(particleMass, 2) + pow(initialMomentum.mag(), 2) ) - particleMass;

			if(particleName == G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName())
			{
				RunInformation->SetMinMaxValue_Control_OpticalPhotonInitialEnergy_MeV(initialEnergy / CLHEP::MeV);
			}
			else
			{
				RunInformation->SetMinMaxValue_Control_SecondaryParticleInitialEnergy_MeV(initialEnergy / CLHEP::MeV);

				if(particleMass)
				{
					G4double initialBetaGamma = initialMomentum.mag() / particleMass;
					RunInformation->SetMinMaxValue_Control_SecondaryParticleInitialBetaGamma(initialBetaGamma);
				}
			}

			if(particleName != G4OpticalPhoton::OpticalPhotonDefinition()->GetParticleName())
			{
				G4double deltaE = GoddessDataStorage->GetEnergyDepositionInScintillator(iter);
				controlFile << "E_depos/MeV:\t\t" << deltaE / CLHEP::MeV << "\n";
				RunInformation->SetMinMaxValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV(deltaE / CLHEP::MeV);

				if(GoddessDataStorage->PhotonDetectorWasHit(iter))
				{
					controlFile << "SiPM hit:\t\ttrue\n";
				}

				RunInformation->SetMinMaxValue_Control_DeltaTimeHitCreation_secondary_ns(globalCreationTime / CLHEP::ns - globalHitTime / CLHEP::ns);
			}
			else
			{
				G4double globalAbsorptionTime = EventInformation->GetGlobalAbsorptionTime(iter);
				controlFile << "t_abs/ns:\t\t" << globalAbsorptionTime / CLHEP::ns << "\n";

				RunInformation->SetMinMaxValue_Control_DeltaTimeHitCreation_ns(globalCreationTime / CLHEP::ns - globalHitTime / CLHEP::ns);
				RunInformation->SetMinMaxValue_Control_DeltaTimeCreationAbsorption_ns(globalAbsorptionTime / CLHEP::ns - globalCreationTime / CLHEP::ns);
				RunInformation->SetMinMaxValue_Control_DeltaTimeHitAbsorption_ns(globalAbsorptionTime / CLHEP::ns - globalHitTime / CLHEP::ns);

				G4bool isParentPrimary = EventInformation->GetParentIsPrimary(iter);
				if(isParentPrimary)
				{
					RunInformation->SetMinMaxValue_Control_DeltaTimeHitCreation_parentPrimary_ns(globalCreationTime / CLHEP::ns - globalHitTime / CLHEP::ns);
					RunInformation->SetMinMaxValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns(globalAbsorptionTime / CLHEP::ns - globalCreationTime / CLHEP::ns);
					RunInformation->SetMinMaxValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns(globalAbsorptionTime / CLHEP::ns - globalHitTime / CLHEP::ns);
				}
				else
				{
					RunInformation->SetMinMaxValue_Control_DeltaTimeHitCreation_parentSecondary_ns(globalCreationTime / CLHEP::ns - globalHitTime / CLHEP::ns);
					RunInformation->SetMinMaxValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns(globalAbsorptionTime / CLHEP::ns - globalCreationTime / CLHEP::ns);
					RunInformation->SetMinMaxValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns(globalAbsorptionTime / CLHEP::ns - globalHitTime / CLHEP::ns);
				}

				if( ! isnan(globalAbsorptionTime) )
				{
					G4String absorbedIn = EventInformation->GetAbsorbedIn(iter);
					controlFile << "abs_vol:\t\t" << absorbedIn << "\n";

					if(GoddessDataStorage->PhotonDetectorWasHit(iter))
					{
						G4ThreeVector PhotonDetectorHitPoint = GoddessDataStorage->GetPhotonDetectorHitPoint(iter);
						controlFile << "hit_pos/mm:\t\t" << PhotonDetectorHitPoint / CLHEP::mm << "\n";

						G4ThreeVector PhotonDetectorHitMomentum = GoddessDataStorage->GetPhotonDetectorHitMomentum(iter);
						controlFile << "hit_momentum/MeV:\t" << PhotonDetectorHitMomentum / CLHEP::MeV << "\n";
					}
				}
			}

			controlFile << "\n";
		}
	}

	// create event summaries
	G4int num_primary = EventInformation->GetNumberOfPrimaryParticles();
	RunInformation->IncNumberOfPrimaryParticles(num_primary);

	controlFile << "##### summary of event " << EventID << " #####" << "\n";
	controlFile << "number of primary particles:\t\t" << num_primary << "\n";
	controlFile << "primaryParticleName:\t\t\t" << name_primary << "\n";
	G4cout << "primary particles:\t\t\t" << name_primary << " (" << num_primary << ")" << G4endl;

	controlFile << "primaryParticleID:\t\t\t" << ID_primary << "\n";
	G4cout << "primaryParticleID:\t\t\t" << ID_primary << G4endl;

	controlFile << "primaryParticle_pos/mm:\t\t\t" << initialPosition_primary / CLHEP::mm << "\n";
	G4cout << "primaryParticle_pos/mm:\t\t\t" << initialPosition_primary / CLHEP::mm << G4endl;

	controlFile << "primaryParticle_momentum/MeV:\t\t" << initialMomentum_primary / CLHEP::MeV << "\n";
	G4cout << "primaryParticle_momentum/MeV:\t\t" << initialMomentum_primary / CLHEP::MeV << G4endl;

	controlFile << "primaryParticle_PathLength/mm:\t\t" << scintiPathLength_primary / CLHEP::mm << "\n";
	G4cout << "primaryParticle_PathLength/mm:\t\t" << scintiPathLength_primary / CLHEP::mm << G4endl;

	controlFile << "primaryParticle_E_depos/MeV:\t\t" << control_deltaE_primary_MeV / CLHEP::MeV << "\n";
	G4cout << "primaryParticle_E_depos/MeV:\t\t" << control_deltaE_primary_MeV / CLHEP::MeV << G4endl;

	controlFile << "primaryParticle_hit_time/ns:\t\t" << globalHitTime / CLHEP::ns << "\n";
	G4cout << "primaryParticle_hit_time/ns:\t\t" << globalHitTime / CLHEP::ns << G4endl;

	controlFile << "primaryParticle_hit_pos/mm:\t\t" << scintiHitPoint_primary / CLHEP::mm << "\n";
	G4cout << "primaryParticle_hit_pos/mm:\t\t" << scintiHitPoint_primary / CLHEP::mm << G4endl;

	controlFile << "primaryParticle_hit_momentum/MeV:\t" << scintiHitMomentum_primary / CLHEP::MeV << "\n";
	G4cout << "primaryParticle_hit_momentum/MeV:\t\t" << scintiHitMomentum_primary / CLHEP::MeV << G4endl;

	controlFile << "secondary particles:\t\t\t";
	G4cout << "secondary particles:\t\t\t";
	std::map<G4String, G4int>::iterator iter = names_secondary.begin();
	while(iter != names_secondary.end())
	{
		G4String name = iter->first;
		G4int number = iter->second;
		controlFile << name << " (" << number << ")";
		G4cout << name << " (" << number << ")";
		iter++;
		if(iter != names_secondary.end())
		{
			controlFile << "; ";
			G4cout << "; ";
		}
	}
	controlFile << "\n";
	G4cout << G4endl;

	controlFile << "production processes:\t\t\t";
	G4cout << "production processes:\t\t\t";
	iter = names_processes.begin();
	while(iter != names_processes.end())
	{
		G4String name = iter->first;
		G4int number = iter->second;
		controlFile << name << " (" << number << ")";
		G4cout << name << " (" << number << ")";
		iter++;
		if(iter != names_processes.end())
		{
			controlFile << "; ";
			G4cout << "; ";
		}
	}
	controlFile << "\n";
	G4cout << G4endl;

	int opticalPhotonsFromPrimary = control_CerenkovPhotonsFromPrimary + control_scintiPhotonsFromPrimary + control_WLSPhotonsFromPrimary;
	int opticalPhotonsFromSecondary = control_CerenkovPhotonsFromSecondary + control_scintiPhotonsFromSecondary + control_WLSPhotonsFromSecondary;

	RunInformation->SetMinMaxValue_Control_opticalPhotons(opticalPhotonsFromPrimary + opticalPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromPrimary(opticalPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromSecondary(opticalPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_CerenkovPhotons(control_CerenkovPhotonsFromPrimary + control_CerenkovPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_ScintiPhotons(control_scintiPhotonsFromPrimary + control_scintiPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_WLSPhotons(control_WLSPhotonsFromPrimary + control_WLSPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_CerenkovPhotonsFromPrimary(control_CerenkovPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_Control_ScintiPhotonsFromPrimary(control_scintiPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_Control_WLSPhotonsFromPrimary(control_WLSPhotonsFromPrimary);
	RunInformation->SetMinMaxValue_Control_CerenkovPhotonsFromSecondary(control_CerenkovPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_ScintiPhotonsFromSecondary(control_scintiPhotonsFromSecondary);
	RunInformation->SetMinMaxValue_Control_WLSPhotonsFromSecondary(control_WLSPhotonsFromSecondary);

	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition(((G4double) opticalPhotonsFromPrimary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition(((G4double) opticalPhotonsFromSecondary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition(((G4double) control_CerenkovPhotonsFromPrimary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition(((G4double) control_scintiPhotonsFromPrimary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition(((G4double) control_WLSPhotonsFromPrimary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition(((G4double) control_CerenkovPhotonsFromSecondary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition(((G4double) control_scintiPhotonsFromSecondary) / (control_deltaE_primary_MeV / CLHEP::MeV));
	RunInformation->SetMinMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition(((G4double) control_WLSPhotonsFromSecondary) / (control_deltaE_primary_MeV / CLHEP::MeV));

	controlFile << "CerenkovPhotons_from_primaryParticle:\t" << control_CerenkovPhotonsFromPrimary << "\n";
	controlFile << "scintiPhotons_from_primaryParticle:\t" << control_scintiPhotonsFromPrimary << "\n";
	controlFile << "WLSPhotons_from_primaryParticle:\t" << control_WLSPhotonsFromPrimary << "\n";
	controlFile << "CerenkovPhotons_from_secondaryParticle:\t" << control_CerenkovPhotonsFromSecondary << "\n";
	controlFile << "scintiPhotons_from_secondaryParticle:\t" << control_scintiPhotonsFromSecondary << "\n";
	controlFile << "WLSPhotons_from_secondaryParticle:\t" << control_WLSPhotonsFromSecondary << "\n";

	G4int num_AbsorbedScinti_inFibre            = EventInformation->GetAbsorbedScintiPhotonInFibreCount();
	G4int num_AbsorbedCer_inFibre               = EventInformation->GetAbsorbedCerenkovPhotonInFibreCount();
	G4int num_AbsorbedWLS_inFibre               = EventInformation->GetAbsorbedWLSPhotonInFibreCount();
	G4int num_AbsorbedScinti_inPhotonDetector = EventInformation->GetAbsorbedScintiPhotonInPhotonDetectorCount();
	G4int num_AbsorbedCer_inPhotonDetector    = EventInformation->GetAbsorbedCerenkovPhotonInPhotonDetectorCount();
	G4int num_AbsorbedWLS_inPhotonDetector    = EventInformation->GetAbsorbedWLSPhotonInPhotonDetectorCount();
	G4int num_AbsorbedScinti_inScinti         = EventInformation->GetAbsorbedScintiPhotonInScintiCount();
	G4int num_AbsorbedCer_inScinti            = EventInformation->GetAbsorbedCerenkovPhotonInScintiCount();
	G4int num_AbsorbedWLS_inScinti            = EventInformation->GetAbsorbedWLSPhotonInScintiCount();

	G4int num_Absorbed_inFibre = num_AbsorbedScinti_inFibre + num_AbsorbedCer_inFibre + num_AbsorbedWLS_inFibre;
	G4int num_Absorbed_inPhotonDetector = num_AbsorbedScinti_inPhotonDetector + num_AbsorbedCer_inPhotonDetector + num_AbsorbedWLS_inPhotonDetector;
	G4int num_Absorbed_inScinti = num_AbsorbedScinti_inScinti + num_AbsorbedCer_inScinti + num_AbsorbedWLS_inScinti;
	G4int num_ScintiAbsorbed = num_AbsorbedScinti_inFibre + num_AbsorbedScinti_inPhotonDetector + num_AbsorbedScinti_inScinti;
	G4int num_CerAbsorbed = num_AbsorbedCer_inFibre + num_AbsorbedCer_inPhotonDetector + num_AbsorbedCer_inScinti;
	G4int num_WLSAbsorbed = num_AbsorbedWLS_inFibre + num_AbsorbedWLS_inPhotonDetector + num_AbsorbedWLS_inScinti;

	RunInformation->SetMinMaxValue_Control_opticalPhotonsAbsorbed_volume(num_Absorbed_inScinti, num_Absorbed_inFibre, num_Absorbed_inPhotonDetector, opticalPhotonsFromPrimary + opticalPhotonsFromSecondary);

	RunInformation->SetMinMaxValue_Control_opticalPhotonsAbsorbed_process(num_ScintiAbsorbed, control_scintiPhotonsFromPrimary + control_scintiPhotonsFromSecondary, num_CerAbsorbed, control_CerenkovPhotonsFromPrimary + control_CerenkovPhotonsFromSecondary, num_WLSAbsorbed, control_WLSPhotonsFromPrimary + control_WLSPhotonsFromSecondary);

	controlFile << "Cerenkov_photons_absorbed_in_WLS:\t" << num_AbsorbedCer_inFibre << "\n";
	controlFile << "scinti_photons_absorbed_in_WLS:\t\t" << num_AbsorbedScinti_inFibre << "\n";
	controlFile << "WLS_photons_absorbed_in_WLS:\t\t" << num_AbsorbedWLS_inFibre << "\n";
	controlFile << "Cerenkov_photons_absorbed_in_SiPM:\t" << num_AbsorbedCer_inPhotonDetector << "\n";
	controlFile << "scinti_photons_absorbed_in_SiPM:\t" << num_AbsorbedScinti_inPhotonDetector << "\n";
	controlFile << "WLS_photons_absorbed_in_SiPM:\t\t" << num_AbsorbedWLS_inPhotonDetector << "\n";
	controlFile << "Cerenkov_photons_absorbed_in_Scinti:\t" << num_AbsorbedCer_inScinti << "\n";
	controlFile << "scinti_photons_absorbed_in_Scinti:\t" << num_AbsorbedScinti_inScinti << "\n";
	controlFile << "WLS_photons_absorbed_in_Scinti:\t\t" << num_AbsorbedWLS_inScinti << "\n";

	controlFile << "#####\n";

	controlFile.close();
}



void EventAction::printEventSummary()
{
	G4int num_ScintiPhotons                       = EventInformation->GetScintiPhotonCount();
	G4int num_ScintiPhotons_byPrimary             = EventInformation->GetScintiPhotonByPrimaryCount();
	G4int num_ScintiPhotons_bySecondary           = EventInformation->GetScintiPhotonBySecondaryCount();
	G4int num_CerPhotons                          = EventInformation->GetCerenkovPhotonCount();
	G4int num_CerPhotons_byPrimary                = EventInformation->GetCerenkovPhotonByPrimaryCount();
	G4int num_CerPhotons_bySecondary              = EventInformation->GetCerenkovPhotonBySecondaryCount();
	G4int num_WLSPhotons                          = EventInformation->GetWLSPhotonCount();
	G4int num_Absorbed                            = EventInformation->GetAbsorbedCount();
	G4int num_AbsorbedScinti                      = EventInformation->GetAbsorbedScintiPhotonCount();
	G4int num_AbsorbedCer                         = EventInformation->GetAbsorbedCerenkovPhotonCount();
	G4int num_AbsorbedWLS                         = EventInformation->GetAbsorbedWLSPhotonCount();
	G4int num_Absorbed_inFibre                    = EventInformation->GetAbsorbedInFibreCount();
	G4int num_AbsorbedScinti_inFibre              = EventInformation->GetAbsorbedScintiPhotonInFibreCount();
	G4int num_AbsorbedCer_inFibre                 = EventInformation->GetAbsorbedCerenkovPhotonInFibreCount();
	G4int num_AbsorbedWLS_inFibre                 = EventInformation->GetAbsorbedWLSPhotonInFibreCount();
	G4int num_Absorbed_inPhotonDetector           = EventInformation->GetAbsorbedInPhotonDetectorCount();
	G4int num_AbsorbedScinti_inPhotonDetector     = EventInformation->GetAbsorbedScintiPhotonInPhotonDetectorCount();
	G4int num_AbsorbedCer_inPhotonDetector        = EventInformation->GetAbsorbedCerenkovPhotonInPhotonDetectorCount();
	G4int num_AbsorbedWLS_inPhotonDetector        = EventInformation->GetAbsorbedWLSPhotonInPhotonDetectorCount();
	G4int num_Absorbed_inScinti                   = EventInformation->GetAbsorbedInScintiCount();
	G4int num_AbsorbedScinti_inScinti             = EventInformation->GetAbsorbedScintiPhotonInScintiCount();
	G4int num_AbsorbedCer_inScinti                = EventInformation->GetAbsorbedCerenkovPhotonInScintiCount();
	G4int num_AbsorbedWLS_inScinti                = EventInformation->GetAbsorbedWLSPhotonInScintiCount();

	G4cout << "Number of optical photons:" << G4endl;
	G4cout << "total:                 " << num_ScintiPhotons + num_CerPhotons + num_WLSPhotons << G4endl;
	G4cout << "by scintillation:      " << num_ScintiPhotons << "\t(by primary: " << num_ScintiPhotons_byPrimary << "\tby secondary: " << num_ScintiPhotons_bySecondary << ")" << G4endl;
	G4cout << "by Cerenkov radiation: " << num_CerPhotons << "\t(by primary: " << num_CerPhotons_byPrimary << "\tby secondary: " << num_CerPhotons_bySecondary << ")" << G4endl;
	G4cout << "by WLS:                " << num_WLSPhotons << G4endl << G4endl;

	G4cout << "Number of optical photons absorbed:" << G4endl;
	G4cout << "total:           " << num_Absorbed << "\t(scintillation photon: " << num_AbsorbedScinti << "\tCerenkov photon: " << num_AbsorbedCer << "\tWLS photon: " << num_AbsorbedWLS << ")" << G4endl;
	G4cout << "in fibre:        " << num_Absorbed_inFibre << "\t(scintillation photon: " << num_AbsorbedScinti_inFibre << "\tCerenkov photon: " << num_AbsorbedCer_inFibre << "\tWLS photon: " << num_AbsorbedWLS_inFibre << ")" << G4endl;
	G4cout << "in SiPM:         " << num_Absorbed_inPhotonDetector << "\t(scintillation photon: " << num_AbsorbedScinti_inPhotonDetector << "\tCerenkov photon: " << num_AbsorbedCer_inPhotonDetector << "\tWLS photon: " << num_AbsorbedWLS_inPhotonDetector << ")" << G4endl;
	G4cout << "in scintillator, wrapping, optical cement,...: " << num_Absorbed_inScinti << "\t(scintillation photon: " << num_AbsorbedScinti_inScinti << "\tCerenkov photon: " << num_AbsorbedCer_inScinti << "\tWLS photon: " << num_AbsorbedWLS_inScinti << ")" << G4endl << G4endl;
	G4cout << "-----------------------------------------------------------------------------" << G4endl << G4endl;
}
