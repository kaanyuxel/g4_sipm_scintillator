/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4EventManager.hh>

#include <fstream>
#include <iostream>

#include <RunAction.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function that will be called by Geant4 at the beginning of each run.
 */
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	RunInformation->ResetDefaults();
	RunInformation->SetNumberOfEvents(aRun->GetNumberOfEventToBeProcessed());

	G4int runID = aRun->GetRunID();
	G4cout << "### Run " << runID << " start." << G4endl;

	G4String outFileName = Messenger->GetDataFileName();
	G4String controlOutFileName = Messenger->GetControlFileName();

	std::ofstream outFile;
	std::ofstream controlOutFile;
	std::ofstream hitFile;

	outFile.open(outFileName.c_str(), std::ios_base::out | std::ios_base::app);
	controlOutFile.open(controlOutFileName.c_str(), std::ios_base::out | std::ios_base::app);

	if (runID == 0) outFile << "simulation data\n" << "--------------------\n";
	if (runID == 0) controlOutFile << "simulation data\n" << "--------------------\n";

	outFile << "\nRunID:\t\t" << runID << "\n";
	controlOutFile << "\nRunID:\t\t" << runID << "\n";

	outFile.close();
	controlOutFile.close();

	Messenger->GetGoddessMessenger()->GetPhotonDetectorConstructor()->WriteRunIDToHitFile(runID);

	Timer->Start();
}

/**
 *  Function that will be called by Geant4 at the end of each run.
 */
void RunAction::EndOfRunAction(const G4Run* aRun)
{
	G4int runID = aRun->GetRunID();

	saveGPSSummary(runID);
	saveDataSummary(runID);
	saveControlSummary(runID);


	G4String outFileName = Messenger->GetDataFileName();
	G4String controlOutFileName = Messenger->GetControlFileName();

	std::ofstream outFile;
	std::ofstream controlOutFile;

	outFile.open(outFileName.c_str(), std::ios_base::out | std::ios_base::app);
	controlOutFile.open(controlOutFileName.c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\n--------------------\n" << "all events processed\n" << "--------------------\n";
	controlOutFile << "\n--------------------\n" << "all events processed\n" << "--------------------\n";

	outFile.close();
	controlOutFile.close();


	G4int numberOfEventsProcessed = aRun->GetNumberOfEvent();

	G4cout << G4endl << "number of events processed = " << numberOfEventsProcessed << G4endl;

	Timer->Stop();
	G4cout << "required: " << Timer->GetRealElapsed() << "s" << G4endl << G4endl;
}



void RunAction::saveGPSSummary(G4int runID)
{
	G4int gps_ParticleID_min = 0;
	G4int gps_ParticleID_max = 0;
	G4double gps_initialEnergy_MeV_min = 0.;
	G4double gps_initialEnergy_MeV_max = 0.;
	G4double gps_initialBetaGamma_MeV_min = 0.;
	G4double gps_initialBetaGamma_MeV_max = 0.;
	G4double gps_xPositionOnSourcePlain_mm_min = 0.;
	G4double gps_xPositionOnSourcePlain_mm_max = 0.;
	G4double gps_yPositionOnSourcePlain_mm_min = 0.;
	G4double gps_yPositionOnSourcePlain_mm_max = 0.;
	G4double gps_globalXPosition_mm_min = 0.;
	G4double gps_globalXPosition_mm_max = 0.;
	G4double gps_globalYPosition_mm_min = 0.;
	G4double gps_globalYPosition_mm_max = 0.;
	G4double gps_globalZPosition_mm_min = 0.;
	G4double gps_globalZPosition_mm_max = 0.;

	char line[1000];
	double rfloat;
	double v1;
	double v2;
	double v3;

	// read in the file and get the minimal and maximal values
	G4String inputFileName = Messenger->GetGPSFileName();

	FILE * inFile;
	inFile = fopen(inputFileName.c_str(), "r");

	if(inFile)
	{
		while(!feof(inFile))
		{
			//find key words and read in the data
			fgets(line, 1000, inFile);

			if(strncmp(line, "particleID:", strlen("particleID:")) == 0)
			{
				sscanf(line, "%*s %lf", &rfloat);

				FindMinMaxValues((G4int) rfloat, gps_ParticleID_min, gps_ParticleID_max);
			}
			else if(strncmp(line, "initial_energy_/_MeV:", strlen("initial_energy_/_MeV:")) == 0)
			{
				sscanf(line, "%*s %lf", &rfloat);

				FindMinMaxValues(rfloat, gps_initialEnergy_MeV_min, gps_initialEnergy_MeV_max);
			}
			else if(strncmp(line, "initial_beta*gamma:", strlen("initial_beta*gamma:")) == 0)
			{
				sscanf(line, "%*s %lf", &rfloat);

				FindMinMaxValues(rfloat, gps_initialBetaGamma_MeV_min, gps_initialBetaGamma_MeV_max);
			}
			else if(strncmp(line, "x_position_in_source_plain_/_mm:", strlen("x position in source plain / mm:")) == 0)
			{
				sscanf(line, "%*s %lf", &rfloat);

				FindMinMaxValues(rfloat, gps_xPositionOnSourcePlain_mm_min, gps_xPositionOnSourcePlain_mm_max);
			}
			else if(strncmp(line, "y_position_in_source_plain_/_mm:", strlen("y position in source plain / mm:")) == 0)
			{
				sscanf(line, "%*s %lf", &rfloat);

				FindMinMaxValues(rfloat, gps_yPositionOnSourcePlain_mm_min, gps_yPositionOnSourcePlain_mm_max);
			}
			else if(strncmp(line, "global_initial_particle_position_/_mm:", strlen("global initial particle position / mm:")) == 0)
			{
				sscanf(line, "%*s (%lf ,%lf,%lf)", &v1,&v2,&v3);

				FindMinMaxValues(v1, gps_globalXPosition_mm_min, gps_globalXPosition_mm_max);
				FindMinMaxValues(v2, gps_globalYPosition_mm_min, gps_globalYPosition_mm_max);
				FindMinMaxValues(v3, gps_globalZPosition_mm_min, gps_globalZPosition_mm_max);
			}
		}

		fclose(inFile);
	}

	// save the minimal and maximal values to a summary file
	std::ofstream outFile;
	outFile.open((inputFileName.substr(0, inputFileName.find(".")) + ".sum").c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\n##### summary of Run " << runID << " #####\n\n";

	outFile << "minimal and maximal values of different variables (necessary for creating plots):\n";

	outFile << "particleID_min:\t\t\t\t\t" << gps_ParticleID_min << "\n";
	outFile << "particleID_max:\t\t\t\t\t" << gps_ParticleID_max << "\n";

	outFile << "initial_energy_min_/_MeV:\t\t\t" << gps_initialEnergy_MeV_min << "\n";
	outFile << "initial_energy_max_/_MeV:\t\t\t" << gps_initialEnergy_MeV_max << "\n";

	outFile << "initial_beta*gamma_min:\t\t\t\t" << gps_initialBetaGamma_MeV_min << "\n";
	outFile << "initial_beta*gamma_max:\t\t\t\t" << gps_initialBetaGamma_MeV_max << "\n";

	outFile << "x_position_in_source_plain_min_/_mm:\t\t" << gps_xPositionOnSourcePlain_mm_min << "\n";
	outFile << "x_position_in_source_plain_max_/_mm:\t\t" << gps_xPositionOnSourcePlain_mm_max << "\n";

	outFile << "y_position_in_source_plain_min_/_mm:\t\t" << gps_yPositionOnSourcePlain_mm_min << "\n";
	outFile << "y_position_in_source_plain_max_/_mm:\t\t" << gps_yPositionOnSourcePlain_mm_max << "\n";

	outFile << "global_initial_particle_position_X_min_/_mm:\t" << gps_globalXPosition_mm_min << "\n";
	outFile << "global_initial_particle_position_X_max_/_mm:\t" << gps_globalXPosition_mm_max << "\n";

	outFile << "global_initial_particle_position_Y_min_/_mm:\t" << gps_globalYPosition_mm_min << "\n";
	outFile << "global_initial_particle_position_Y_max_/_mm:\t" << gps_globalYPosition_mm_max << "\n";

	outFile << "global_initial_particle_position_Z_min_/_mm:\t" << gps_globalZPosition_mm_min << "\n";
	outFile << "global_initial_particle_position_Z_max_/_mm:\t" << gps_globalZPosition_mm_max << "\n";

	outFile << "\n#####\n";

	outFile.close();
}



void RunAction::saveDataSummary(G4int runID)
{
	G4String dataFileName = Messenger->GetDataFileName();
	dataFileName = dataFileName.substr(0, dataFileName.find(".")) + ".sum";
	std::ofstream outFile;
	outFile.open(dataFileName.c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\n##### summary of Run " << runID << " #####\n\n";

	outFile << "minimal and maximal values of different variables (necessary for creating plots):\n";

	outFile << "pos_hit_x_min/mm:\t\t\t" << RunInformation->GetMinValue_ScintiHitPointX_mm() << "\n";
	outFile << "pos_hit_x_max/mm:\t\t\t" << RunInformation->GetMaxValue_ScintiHitPointX_mm() << "\n";

	outFile << "pos_hit_y_min/mm:\t\t\t" << RunInformation->GetMinValue_ScintiHitPointY_mm() << "\n";
	outFile << "pos_hit_y_max/mm:\t\t\t" << RunInformation->GetMaxValue_ScintiHitPointY_mm() << "\n";

	outFile << "pos_hit_z_min/mm:\t\t\t" << RunInformation->GetMinValue_ScintiHitPointZ_mm() << "\n";
	outFile << "pos_hit_z_max/mm:\t\t\t" << RunInformation->GetMaxValue_ScintiHitPointZ_mm() << "\n";

	outFile << "optical_photons_per_primary_min:\t" << RunInformation->GetMinValue_opticalPhotonsPerPrimary() << "\n";
	outFile << "optical_photons_per_primary_max:\t" << RunInformation->GetMaxValue_opticalPhotonsPerPrimary() << "\n";

	outFile << "optical_photons_per_E_depos/MeV_min:\t" << RunInformation->GetMinValue_opticalPhotonsPerEnergyDeposition() << "\n";
	outFile << "optical_photons_per_E_depos/MeV_max:\t" << RunInformation->GetMaxValue_opticalPhotonsPerEnergyDeposition() << "\n";

	outFile << "optical_photons_absorbed_in_SiPM_min:\t" << RunInformation->GetMinValue_opticalPhotonsAbsorbedInPhotonDetector() << "\n";
	outFile << "optical_photons_absorbed_in_SiPM_max:\t" << RunInformation->GetMaxValue_opticalPhotonsAbsorbedInPhotonDetector() << "\n";

	outFile << "pos_SiPM_hit_x_min/mm:\t\t\t" << RunInformation->GetMinValue_PhotonDetectorHitPointX_mm() << "\n";
	outFile << "pos_SiPM_hit_x_max/mm:\t\t\t" << RunInformation->GetMaxValue_PhotonDetectorHitPointX_mm() << "\n";

	outFile << "pos_SiPM_hit_y_min/mm:\t\t\t" << RunInformation->GetMinValue_PhotonDetectorHitPointY_mm() << "\n";
	outFile << "pos_SiPM_hit_y_max/mm:\t\t\t" << RunInformation->GetMaxValue_PhotonDetectorHitPointY_mm() << "\n";

// 	outFile << "pos_SiPM_hit_z_min/mm:\t\t\t" << RunInformation->GetMinValue_PhotonDetectorHitPointZ_mm() << "\n";
// 	outFile << "pos_SiPM_hit_z_max/mm:\t\t\t" << RunInformation->GetMaxValue_PhotonDetectorHitPointZ_mm() << "\n";

	outFile << "opticalPhotonHitEnergy_min/MeV:\t\t" << RunInformation->GetMinValue_OpticalPhotonHitEnergy_MeV() << "\n";
	outFile << "opticalPhotonHitEnergy_max/MeV:\t\t" << RunInformation->GetMaxValue_OpticalPhotonHitEnergy_MeV() << "\n";

	outFile << "Delta_t_hit-absorption_min/ns:\t\t" << RunInformation->GetMinValue_DeltaTimeHitAbsorption_ns() << "\n";
	outFile << "Delta_t_hit-absorption_max/ns:\t\t" << RunInformation->GetMaxValue_DeltaTimeHitAbsorption_ns() << "\n";

	outFile << "total_number_of_primary_particles:\t" << RunInformation->GetNumberOfPrimaryParticles() << "\n";

	outFile.close();


	RunInformation->DumpParticleMassesToFile(dataFileName);


	outFile.open(dataFileName.c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\n#####\n";

	outFile.close();
}



void RunAction::saveControlSummary(G4int runID)
{
	G4String controlFileName = Messenger->GetControlFileName();
	controlFileName = controlFileName.substr(0, controlFileName.find(".")) + ".sum";
	std::ofstream outFile;
	outFile.open(controlFileName.c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\n##### summary of Run " << runID << " #####\n\n";

	outFile << "minimal and maximal values of different variables (necessary for creating plots):\n";

	outFile << "primaryParticleID_min:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleID() << "\n";
	outFile << "primaryParticleID_max:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleID() << "\n";

	outFile << "secondaryParticleID_min:\t\t\t\t\t" << RunInformation->GetMinValue_Control_SecondaryParticleID() << "\n";
	outFile << "secondaryParticleID_max:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_SecondaryParticleID() << "\n";

	outFile << "pos_init_x_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleInitialPositionX_mm() << "\n";
	outFile << "pos_init_x_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleInitialPositionX_mm() << "\n";

	outFile << "pos_init_y_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleInitialPositionY_mm() << "\n";
	outFile << "pos_init_y_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleInitialPositionY_mm() << "\n";

	outFile << "pos_init_z_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleInitialPositionZ_mm() << "\n";
	outFile << "pos_init_z_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleInitialPositionZ_mm() << "\n";

	outFile << "pos_hit_x_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_ScintiHitPointX_mm() << "\n";
	outFile << "pos_hit_x_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_ScintiHitPointX_mm() << "\n";

	outFile << "pos_hit_y_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_ScintiHitPointY_mm() << "\n";
	outFile << "pos_hit_y_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_ScintiHitPointY_mm() << "\n";

	outFile << "pos_hit_z_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_ScintiHitPointZ_mm() << "\n";
	outFile << "pos_hit_z_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_ScintiHitPointZ_mm() << "\n";

	outFile << "primaryParticleInitialEnergy_min/MeV:\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleInitialEnergy_MeV() << "\n";
	outFile << "primaryParticleInitialEnergy_max/MeV:\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleInitialEnergy_MeV() << "\n";

	outFile << "secondaryParticleInitialEnergy_min/MeV:\t\t\t\t" << RunInformation->GetMinValue_Control_SecondaryParticleInitialEnergy_MeV() << "\n";
	outFile << "secondaryParticleInitialEnergy_max/MeV:\t\t\t\t" << RunInformation->GetMaxValue_Control_SecondaryParticleInitialEnergy_MeV() << "\n";

	outFile << "opticalPhotonInitialEnergy_min/MeV:\t\t\t\t" << RunInformation->GetMinValue_Control_OpticalPhotonInitialEnergy_MeV() << "\n";
	outFile << "opticalPhotonInitialEnergy_max/MeV:\t\t\t\t" << RunInformation->GetMaxValue_Control_OpticalPhotonInitialEnergy_MeV() << "\n";

	outFile << "primaryParticleInitialBetaGamma_min:\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleInitialBetaGamma() << "\n";
	outFile << "primaryParticleInitialBetaGamma_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleInitialBetaGamma() << "\n";

	outFile << "secondaryParticleInitialBetaGamma_min:\t\t\t\t" << RunInformation->GetMinValue_Control_SecondaryParticleInitialBetaGamma() << "\n";
	outFile << "secondaryParticleInitialBetaGamma_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_SecondaryParticleInitialBetaGamma() << "\n";

	outFile << "E_depos_primary_min/MeV:\t\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV() << "\n";
	outFile << "E_depos_primary_max/MeV:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV() << "\n";

	outFile << "DeltaX_primary_min/mm:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_PrimaryParticlePathLengthInScintillator_MeV() << "\n";
	outFile << "DeltaX_primary_max/mm:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_PrimaryParticlePathLengthInScintillator_MeV() << "\n";

	outFile << "E_depos_secondary_min/MeV:\t\t\t\t\t" << RunInformation->GetMinValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV() << "\n";
	outFile << "E_depos_secondary_max/MeV:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV() << "\n";

	outFile << "t_hit_min/ns:\t\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_GlobalScintiHitTime_ns() << "\n";
	outFile << "t_hit_max/ns:\t\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_GlobalScintiHitTime_ns() << "\n";

	outFile << "Delta_t_hit-creation_secondary_min/ns:\t\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitCreation_secondary_ns() << "\n";
	outFile << "Delta_t_hit-creation_secondary_max/ns:\t\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitCreation_secondary_ns() << "\n";

	outFile << "Delta_t_hit-creation_min/ns:\t\t\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitCreation_ns() << "\n";
	outFile << "Delta_t_hit-creation_max/ns:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitCreation_ns() << "\n";

	outFile << "Delta_t_creation-absorption_min/ns:\t\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeCreationAbsorption_ns() << "\n";
	outFile << "Delta_t_creation-absorption_max/ns:\t\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeCreationAbsorption_ns() << "\n";

	outFile << "Delta_t_hit-absorption_min/ns:\t\t\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitAbsorption_ns() << "\n";
	outFile << "Delta_t_hit-absorption_max/ns:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitAbsorption_ns() << "\n";

	outFile << "Delta_t_hit-creation_parentPrimary_min/ns:\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitCreation_parentPrimary_ns() << "\n";
	outFile << "Delta_t_hit-creation_parentPrimary_max/ns:\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitCreation_parentPrimary_ns() << "\n";

	outFile << "Delta_t_creation-absorption_parentPrimary_min/ns:\t\t" << RunInformation->GetMinValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns() << "\n";
	outFile << "Delta_t_creation-absorption_parentPrimary_max/ns:\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns() << "\n";

	outFile << "Delta_t_hit-absorption_parentPrimary_min/ns:\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns() << "\n";
	outFile << "Delta_t_hit-absorption_parentPrimary_max/ns:\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns() << "\n";

	outFile << "Delta_t_hit-creation_parentSecondary_min/ns:\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitCreation_parentSecondary_ns() << "\n";
	outFile << "Delta_t_hit-creation_parentSecondary_max/ns:\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitCreation_parentSecondary_ns() << "\n";

	outFile << "Delta_t_creation-absorption_parentSecondary_min/ns:\t\t" << RunInformation->GetMinValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns() << "\n";
	outFile << "Delta_t_creation-absorption_parentSecondary_max/ns:\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns() << "\n";

	outFile << "Delta_t_hit-absorption_parentSecondary_min/ns:\t\t\t" << RunInformation->GetMinValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns() << "\n";
	outFile << "Delta_t_hit-absorption_parentSecondary_max/ns:\t\t\t" << RunInformation->GetMaxValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns() << "\n";

	outFile << "optical_photons_min:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_opticalPhotons() << "\n";
	outFile << "optical_photons_max:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_opticalPhotons() << "\n";

	outFile << "optical_photons_from_primary_min:\t\t\t\t" << RunInformation->GetMinValue_Control_opticalPhotonsFromPrimary() << "\n";
	outFile << "optical_photons_from_primary_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_opticalPhotonsFromPrimary() << "\n";

	outFile << "optical_photons_from_secondary_min:\t\t\t\t" << RunInformation->GetMinValue_Control_opticalPhotonsFromSecondary() << "\n";
	outFile << "optical_photons_from_secondary_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_opticalPhotonsFromSecondary() << "\n";

	outFile << "Cerenkov_photons_min:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_CerenkovPhotons() << "\n";
	outFile << "Cerenkov_photons_max:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_CerenkovPhotons() << "\n";

	outFile << "scinti_photons_min:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_ScintiPhotons() << "\n";
	outFile << "scinti_photons_max:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_ScintiPhotons() << "\n";

	outFile << "WLS_photons_min:\t\t\t\t\t\t" << RunInformation->GetMinValue_Control_WLSPhotons() << "\n";
	outFile << "WLS_photons_max:\t\t\t\t\t\t" << RunInformation->GetMaxValue_Control_WLSPhotons() << "\n";

	outFile << "Cerenkov_photons_from_primary_min:\t\t\t\t" << RunInformation->GetMinValue_Control_CerenkovPhotonsFromPrimary() << "\n";
	outFile << "Cerenkov_photons_from_primary_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_CerenkovPhotonsFromPrimary() << "\n";

	outFile << "scinti_photons_from_primary_min:\t\t\t\t" << RunInformation->GetMinValue_Control_ScintiPhotonsFromPrimary() << "\n";
	outFile << "scinti_photons_from_primary_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_ScintiPhotonsFromPrimary() << "\n";

	outFile << "WLS_photons_from_primary_min:\t\t\t\t\t" << RunInformation->GetMinValue_Control_WLSPhotonsFromPrimary() << "\n";
	outFile << "WLS_photons_from_primary_max:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_WLSPhotonsFromPrimary() << "\n";

	outFile << "Cerenkov_photons_from_secondary_min:\t\t\t\t" << RunInformation->GetMinValue_Control_CerenkovPhotonsFromSecondary() << "\n";
	outFile << "Cerenkov_photons_from_secondary_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_CerenkovPhotonsFromSecondary() << "\n";

	outFile << "scinti_photons_from_secondary_min:\t\t\t\t" << RunInformation->GetMinValue_Control_ScintiPhotonsFromSecondary() << "\n";
	outFile << "scinti_photons_from_secondary_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_ScintiPhotonsFromSecondary() << "\n";

	outFile << "WLS_photons_from_secondary_min:\t\t\t\t\t" << RunInformation->GetMinValue_Control_WLSPhotonsFromSecondary() << "\n";
	outFile << "WLS_photons_from_secondary_max:\t\t\t\t\t" << RunInformation->GetMaxValue_Control_WLSPhotonsFromSecondary() << "\n";

	outFile << "optical_photons_from_primary_per_E_depos_primary/MeV_min:\t" << RunInformation->GetMinValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition() << "\n";
	outFile << "optical_photons_from_primary_per_E_depos_primary/MeV_max:\t" << RunInformation->GetMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition() << "\n";

	outFile << "optical_photons_from_secondary_per_E_depos_primary/MeV_min:\t" << RunInformation->GetMinValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition() << "\n";
	outFile << "optical_photons_from_secondary_per_E_depos_primary/MeV_max:\t" << RunInformation->GetMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition() << "\n";

	outFile << "optical_photons_absorbed_volume_min:\t\t\t\t" << RunInformation->GetMinValue_Control_opticalPhotonsAbsorbed_volume() << "\n";
	outFile << "optical_photons_absorbed_volume_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_opticalPhotonsAbsorbed_volume() << "\n";

	outFile << "optical_photons_absorbed_process_min:\t\t\t\t" << RunInformation->GetMinValue_Control_opticalPhotonsAbsorbed_process() << "\n";
	outFile << "optical_photons_absorbed_process_max:\t\t\t\t" << RunInformation->GetMaxValue_Control_opticalPhotonsAbsorbed_process() << "\n";

	outFile.close();


	RunInformation->DumpParticleMassesToFile(controlFileName);


	outFile.open(controlFileName.c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\n#####\n";

	outFile.close();
}
