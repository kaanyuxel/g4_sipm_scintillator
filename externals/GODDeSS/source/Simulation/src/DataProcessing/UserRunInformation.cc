/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <fstream>

#include <UserRunInformation.hh>
#include <CLHEP/Units/SystemOfUnits.h>


// class variables begin with capital letters, local variables with small letters



/**
 *  Function to set the class variables to default values:
 *  G4int = -1 (except number of events (= 0) and particle IDs (= 999999999))
 *  G4double = NAN
 */
void UserRunInformation::ResetDefaults()
{
	NumberOfEvents = 0;
	NumberOfPrimaryParticles = 0;
// for data plots
	ScintiHitPointX_mm_min = NAN;
	ScintiHitPointX_mm_max = NAN;
	ScintiHitPointY_mm_min = NAN;
	ScintiHitPointY_mm_max = NAN;
	ScintiHitPointZ_mm_min = NAN;
	ScintiHitPointZ_mm_max = NAN;
	OpticalPhotonsPerPrimary_min = -1;
	OpticalPhotonsPerPrimary_max = -1;
	OpticalPhotonsPerEnergyDeposition_min = -1;
	OpticalPhotonsPerEnergyDeposition_max = -1;
	OpticalPhotonsAbsorbedInPhotonDetector_min = -1;
	OpticalPhotonsAbsorbedInPhotonDetector_max = -1;
	PhotonDetectorHitPointX_mm_min = NAN;
	PhotonDetectorHitPointX_mm_max = NAN;
	PhotonDetectorHitPointY_mm_min = NAN;
	PhotonDetectorHitPointY_mm_max = NAN;
	PhotonDetectorHitPointZ_mm_min = NAN;
	PhotonDetectorHitPointZ_mm_max = NAN;
	OpticalPhotonHitEnergy_MeV_min = NAN;
	OpticalPhotonHitEnergy_MeV_max = NAN;
	DeltaTimeHitAbsorption_ns_min = NAN;
	DeltaTimeHitAbsorption_ns_max = NAN;

// for control plots
	Control_PrimaryParticleID_min = 999999999;
	Control_PrimaryParticleID_max = 999999999;
	Control_SecondaryParticleID_min = 999999999;
	Control_SecondaryParticleID_max = 999999999;
	Control_PrimaryParticleInitialPositionX_mm_min = NAN;
	Control_PrimaryParticleInitialPositionX_mm_max = NAN;
	Control_PrimaryParticleInitialPositionY_mm_min = NAN;
	Control_PrimaryParticleInitialPositionY_mm_max = NAN;
	Control_PrimaryParticleInitialPositionZ_mm_min = NAN;
	Control_PrimaryParticleInitialPositionZ_mm_max = NAN;
	Control_ScintiHitPointX_mm_min = NAN;
	Control_ScintiHitPointX_mm_max = NAN;
	Control_ScintiHitPointY_mm_min = NAN;
	Control_ScintiHitPointY_mm_max = NAN;
	Control_ScintiHitPointZ_mm_min = NAN;
	Control_ScintiHitPointZ_mm_max = NAN;
	Control_PrimaryParticleInitialEnergy_MeV_min = NAN;
	Control_PrimaryParticleInitialEnergy_MeV_max = NAN;
	Control_SecondaryParticleInitialEnergy_MeV_min = NAN;
	Control_SecondaryParticleInitialEnergy_MeV_max = NAN;
	Control_OpticalPhotonInitialEnergy_MeV_min = NAN;
	Control_OpticalPhotonInitialEnergy_MeV_max = NAN;
	Control_PrimaryParticleInitialBetaGamma_min = NAN;
	Control_PrimaryParticleInitialBetaGamma_max = NAN;
	Control_SecondaryParticleInitialBetaGamma_min = NAN;
	Control_SecondaryParticleInitialBetaGamma_max = NAN;
	Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min = NAN;
	Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max = NAN;
	Control_PrimaryParticlePathLengthInScintillator_MeV_min = NAN;
	Control_PrimaryParticlePathLengthInScintillator_MeV_max = NAN;
	Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min = NAN;
	Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max = NAN;
	Control_GlobalScintiHitTime_ns_min = NAN;
	Control_GlobalScintiHitTime_ns_max = NAN;
	Control_DeltaTimeHitCreation_secondary_ns_min = NAN;
	Control_DeltaTimeHitCreation_secondary_ns_max = NAN;
	Control_DeltaTimeHitCreation_ns_min = NAN;
	Control_DeltaTimeHitCreation_ns_max = NAN;
	Control_DeltaTimeCreationAbsorption_ns_min = NAN;
	Control_DeltaTimeCreationAbsorption_ns_max = NAN;
	Control_DeltaTimeHitAbsorption_ns_min = NAN;
	Control_DeltaTimeHitAbsorption_ns_max = NAN;
	Control_DeltaTimeHitCreation_parentPrimary_ns_min = NAN;
	Control_DeltaTimeHitCreation_parentPrimary_ns_max = NAN;
	Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min = NAN;
	Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max = NAN;
	Control_DeltaTimeHitAbsorption_parentPrimary_ns_min = NAN;
	Control_DeltaTimeHitAbsorption_parentPrimary_ns_max = NAN;
	Control_DeltaTimeHitCreation_parentSecondary_ns_min = NAN;
	Control_DeltaTimeHitCreation_parentSecondary_ns_max = NAN;
	Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min = NAN;
	Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max = NAN;
	Control_DeltaTimeHitAbsorption_parentSecondary_ns_min = NAN;
	Control_DeltaTimeHitAbsorption_parentSecondary_ns_max = NAN;
	Control_OpticalPhotons_min = -1;
	Control_OpticalPhotons_max = -1;
	Control_OpticalPhotonsFromPrimary_min = -1;
	Control_OpticalPhotonsFromPrimary_max = -1;
	Control_OpticalPhotonsFromSecondary_min = -1;
	Control_OpticalPhotonsFromSecondary_max = -1;
	Control_CerenkovPhotons_min = -1;
	Control_CerenkovPhotons_max = -1;
	Control_ScintiPhotons_min = -1;
	Control_ScintiPhotons_max = -1;
	Control_WLSPhotons_min = -1;
	Control_WLSPhotons_max = -1;
	Control_CerenkovPhotonsFromPrimary_min = -1;
	Control_CerenkovPhotonsFromPrimary_max = -1;
	Control_ScintiPhotonsFromPrimary_min = -1;
	Control_ScintiPhotonsFromPrimary_max = -1;
	Control_WLSPhotonsFromPrimary_min = -1;
	Control_WLSPhotonsFromPrimary_max = -1;
	Control_CerenkovPhotonsFromSecondary_min = -1;
	Control_CerenkovPhotonsFromSecondary_max = -1;
	Control_ScintiPhotonsFromSecondary_min = -1;
	Control_ScintiPhotonsFromSecondary_max = -1;
	Control_WLSPhotonsFromSecondary_min = -1;
	Control_WLSPhotonsFromSecondary_max = -1;
	Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min = -1;
	Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max = -1;
	Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min = -1;
	Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max = -1;
	Control_OpticalPhotonsAbsorbed_volume_min = -1;
	Control_OpticalPhotonsAbsorbed_volume_max = -1;
	Control_OpticalPhotonsAbsorbed_process_min = -1;
	Control_OpticalPhotonsAbsorbed_process_max = -1;
}



/**
 *  @return particle mass
 */
G4double UserRunInformation::GetParticleMass(G4String particleName) const
{
	if(MassMap.count(particleName)) return MassMap.at(particleName);
	else return -1.;
}

/**
 *  Function to all saved particle masses into a file.
 */
void UserRunInformation::DumpParticleMassesToFile(G4String pathToFile)
{
	std::ofstream outFile;
	outFile.open(pathToFile.c_str(), std::ios_base::out | std::ios_base::app);

	outFile << "\noccuring particles and their masses/MeV (necessary for getting beta*gamma and energies from the momenta):\n";

	std::map<G4String, G4double>::iterator iter = MassMap.begin();
	while(iter != MassMap.end())
	{
		G4String name = iter->first;
		G4double mass = iter->second;

		outFile << "  " << name << "\t" << mass / CLHEP::MeV << "\n";
		iter++;
	}

	outFile.close();
}


// for data plots

/**
 *  Function to set the minimal and maximal X value of the scintillator hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_ScintiHitPointX_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(ScintiHitPointX_mm_min) && isnan(ScintiHitPointX_mm_max))
	{
		ScintiHitPointX_mm_min = newValue;
		ScintiHitPointX_mm_max = newValue;
	}
	else if(ScintiHitPointX_mm_min > newValue) ScintiHitPointX_mm_min = newValue;
	else if(ScintiHitPointX_mm_max < newValue) ScintiHitPointX_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Y value of the scintillator hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_ScintiHitPointY_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(ScintiHitPointY_mm_min) && isnan(ScintiHitPointY_mm_max))
	{
		ScintiHitPointY_mm_min = newValue;
		ScintiHitPointY_mm_max = newValue;
	}
	else if(ScintiHitPointY_mm_min > newValue) ScintiHitPointY_mm_min = newValue;
	else if(ScintiHitPointY_mm_max < newValue) ScintiHitPointY_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Z value of the scintillator hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_ScintiHitPointZ_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(ScintiHitPointZ_mm_min) && isnan(ScintiHitPointZ_mm_max))
	{
		ScintiHitPointZ_mm_min = newValue;
		ScintiHitPointZ_mm_max = newValue;
	}
	else if(ScintiHitPointZ_mm_min > newValue) ScintiHitPointZ_mm_min = newValue;
	else if(ScintiHitPointZ_mm_max < newValue) ScintiHitPointZ_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons per primary particles (= per event).
 */
void UserRunInformation::SetMinMaxValue_opticalPhotonsPerPrimary(G4int newValue)
{
	if(OpticalPhotonsPerPrimary_min == -1 && OpticalPhotonsPerPrimary_max == -1)
	{
		OpticalPhotonsPerPrimary_min = newValue;
		OpticalPhotonsPerPrimary_max = newValue;
	}
	else if(OpticalPhotonsPerPrimary_min > newValue) OpticalPhotonsPerPrimary_min = newValue;
	else if(OpticalPhotonsPerPrimary_max < newValue) OpticalPhotonsPerPrimary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons per primary particle's energy deposition.
 */
void UserRunInformation::SetMinMaxValue_opticalPhotonsPerEnergyDeposition(G4double newValue)
{
	if(isnan(newValue)) return;

	if(OpticalPhotonsPerEnergyDeposition_min == -1 && OpticalPhotonsPerEnergyDeposition_max == -1)
	{
		OpticalPhotonsPerEnergyDeposition_min = newValue;
		OpticalPhotonsPerEnergyDeposition_max = newValue;
	}
	else if(OpticalPhotonsPerEnergyDeposition_min > newValue) OpticalPhotonsPerEnergyDeposition_min = newValue;
	else if(OpticalPhotonsPerEnergyDeposition_max < newValue) OpticalPhotonsPerEnergyDeposition_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons that have been absorbed in a photon detector (per event).
 */
void UserRunInformation::SetMinMaxValue_opticalPhotonsAbsorbedInPhotonDetector(G4int newValue)
{
	if(OpticalPhotonsAbsorbedInPhotonDetector_min == -1 && OpticalPhotonsAbsorbedInPhotonDetector_max == -1)
	{
		OpticalPhotonsAbsorbedInPhotonDetector_min = newValue;
		OpticalPhotonsAbsorbedInPhotonDetector_max = newValue;
	}
	else if(OpticalPhotonsAbsorbedInPhotonDetector_min > newValue) OpticalPhotonsAbsorbedInPhotonDetector_min = newValue;
	else if(OpticalPhotonsAbsorbedInPhotonDetector_max < newValue) OpticalPhotonsAbsorbedInPhotonDetector_max = newValue;
}

/**
 *  Function to set the minimal and maximal X value of the photon detector hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_PhotonDetectorHitPointX_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(PhotonDetectorHitPointX_mm_min) && isnan(PhotonDetectorHitPointX_mm_max))
	{
		PhotonDetectorHitPointX_mm_min = newValue;
		PhotonDetectorHitPointX_mm_max = newValue;
	}
	else if(PhotonDetectorHitPointX_mm_min > newValue) PhotonDetectorHitPointX_mm_min = newValue;
	else if(PhotonDetectorHitPointX_mm_max < newValue) PhotonDetectorHitPointX_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Y value of the photon detector hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_PhotonDetectorHitPointY_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(PhotonDetectorHitPointY_mm_min) && isnan(PhotonDetectorHitPointY_mm_max))
	{
		PhotonDetectorHitPointY_mm_min = newValue;
		PhotonDetectorHitPointY_mm_max = newValue;
	}
	else if(PhotonDetectorHitPointY_mm_min > newValue) PhotonDetectorHitPointY_mm_min = newValue;
	else if(PhotonDetectorHitPointY_mm_max < newValue) PhotonDetectorHitPointY_mm_max = newValue;
}

// void UserRunInformation::SetMinMaxValue_PhotonDetectorHitPointZ_mm(G4double newValue)
// {
// 	if(isnan(newValue)) return;
//
// 	if(isnan(PhotonDetectorHitPointZ_mm_min) && isnan(PhotonDetectorHitPointZ_mm_max))
// 	{
// 		PhotonDetectorHitPointZ_mm_min = newValue;
// 		PhotonDetectorHitPointZ_mm_max = newValue;
// 	}
// 	else if(PhotonDetectorHitPointZ_mm_min > newValue) PhotonDetectorHitPointZ_mm_min = newValue;
// 	else if(PhotonDetectorHitPointZ_mm_max < newValue) PhotonDetectorHitPointZ_mm_max = newValue;
// }

/**
 *  Function to set the minimal and maximal energy (MeV) of optical photons that have been absorbed in a photon detector.
 */
void UserRunInformation::SetMinMaxValue_OpticalPhotonHitEnergy_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(OpticalPhotonHitEnergy_MeV_min) && isnan(OpticalPhotonHitEnergy_MeV_max))
	{
		OpticalPhotonHitEnergy_MeV_min = newValue;
		OpticalPhotonHitEnergy_MeV_max = newValue;
	}
	else if(OpticalPhotonHitEnergy_MeV_min > newValue) OpticalPhotonHitEnergy_MeV_min = newValue;
	else if(OpticalPhotonHitEnergy_MeV_max < newValue) OpticalPhotonHitEnergy_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the optical photon's absorption time (ns) relative to the primary particle's scintillator hit time.
 */
void UserRunInformation::SetMinMaxValue_DeltaTimeHitAbsorption_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(DeltaTimeHitAbsorption_ns_min) && isnan(DeltaTimeHitAbsorption_ns_max))
	{
		DeltaTimeHitAbsorption_ns_min = newValue;
		DeltaTimeHitAbsorption_ns_max = newValue;
	}
	else if(DeltaTimeHitAbsorption_ns_min > newValue) DeltaTimeHitAbsorption_ns_min = newValue;
	else if(DeltaTimeHitAbsorption_ns_max < newValue) DeltaTimeHitAbsorption_ns_max = newValue;
}


// for control plots
/**
 *  Function to set the minimal and maximal value of the primary particles' particle ID.
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleID(G4int newValue)
{
	if(Control_PrimaryParticleID_min == 999999999 && Control_PrimaryParticleID_max == 999999999)
	{
		Control_PrimaryParticleID_min = newValue;
		Control_PrimaryParticleID_max = newValue;
	}
	else if(Control_PrimaryParticleID_min > newValue) Control_PrimaryParticleID_min = newValue;
	else if(Control_PrimaryParticleID_max < newValue) Control_PrimaryParticleID_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the secondary particles' particle ID.
 */
void UserRunInformation::SetMinMaxValue_Control_SecondaryParticleID(G4int newValue)
{
	if(Control_SecondaryParticleID_min == 999999999 && Control_SecondaryParticleID_max == 999999999)
	{
		Control_SecondaryParticleID_min = newValue;
		Control_SecondaryParticleID_max = newValue;
	}
	else if(Control_SecondaryParticleID_min > newValue) Control_SecondaryParticleID_min = newValue;
	else if(Control_SecondaryParticleID_max < newValue) Control_SecondaryParticleID_max = newValue;
}

/**
 *  Function to set the minimal and maximal X value of the primary particles' initial position (mm).
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleInitialPositionX_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticleInitialPositionX_mm_min) && isnan(Control_PrimaryParticleInitialPositionX_mm_max))
	{
		Control_PrimaryParticleInitialPositionX_mm_min = newValue;
		Control_PrimaryParticleInitialPositionX_mm_max = newValue;
	}
	else if(Control_PrimaryParticleInitialPositionX_mm_min > newValue) Control_PrimaryParticleInitialPositionX_mm_min = newValue;
	else if(Control_PrimaryParticleInitialPositionX_mm_max < newValue) Control_PrimaryParticleInitialPositionX_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Y value of the primary particles' initial position (mm).
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleInitialPositionY_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticleInitialPositionY_mm_min) && isnan(Control_PrimaryParticleInitialPositionY_mm_max))
	{
		Control_PrimaryParticleInitialPositionY_mm_min = newValue;
		Control_PrimaryParticleInitialPositionY_mm_max = newValue;
	}
	else if(Control_PrimaryParticleInitialPositionY_mm_min > newValue) Control_PrimaryParticleInitialPositionY_mm_min = newValue;
	else if(Control_PrimaryParticleInitialPositionY_mm_max < newValue) Control_PrimaryParticleInitialPositionY_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Z value of the primary particles' initial position (mm).
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleInitialPositionZ_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticleInitialPositionZ_mm_min) && isnan(Control_PrimaryParticleInitialPositionZ_mm_max))
	{
		Control_PrimaryParticleInitialPositionZ_mm_min = newValue;
		Control_PrimaryParticleInitialPositionZ_mm_max = newValue;
	}
	else if(Control_PrimaryParticleInitialPositionZ_mm_min > newValue) Control_PrimaryParticleInitialPositionZ_mm_min = newValue;
	else if(Control_PrimaryParticleInitialPositionZ_mm_max < newValue) Control_PrimaryParticleInitialPositionZ_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal X value of the primary particle's scintillator hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_Control_ScintiHitPointX_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_ScintiHitPointX_mm_min) && isnan(Control_ScintiHitPointX_mm_max))
	{
		Control_ScintiHitPointX_mm_min = newValue;
		Control_ScintiHitPointX_mm_max = newValue;
	}
	else if(Control_ScintiHitPointX_mm_min > newValue) Control_ScintiHitPointX_mm_min = newValue;
	else if(Control_ScintiHitPointX_mm_max < newValue) Control_ScintiHitPointX_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Y value of the primary particle's scintillator hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_Control_ScintiHitPointY_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_ScintiHitPointY_mm_min) && isnan(Control_ScintiHitPointY_mm_max))
	{
		Control_ScintiHitPointY_mm_min = newValue;
		Control_ScintiHitPointY_mm_max = newValue;
	}
	else if(Control_ScintiHitPointY_mm_min > newValue) Control_ScintiHitPointY_mm_min = newValue;
	else if(Control_ScintiHitPointY_mm_max < newValue) Control_ScintiHitPointY_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal Z value of the primary particle's scintillator hit position (mm).
 */
void UserRunInformation::SetMinMaxValue_Control_ScintiHitPointZ_mm(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_ScintiHitPointZ_mm_min) && isnan(Control_ScintiHitPointZ_mm_max))
	{
		Control_ScintiHitPointZ_mm_min = newValue;
		Control_ScintiHitPointZ_mm_max = newValue;
	}
	else if(Control_ScintiHitPointZ_mm_min > newValue) Control_ScintiHitPointZ_mm_min = newValue;
	else if(Control_ScintiHitPointZ_mm_max < newValue) Control_ScintiHitPointZ_mm_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the primary particles' initial energy (MeV).
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleInitialEnergy_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticleInitialEnergy_MeV_min) && isnan(Control_PrimaryParticleInitialEnergy_MeV_max))
	{
		Control_PrimaryParticleInitialEnergy_MeV_min = newValue;
		Control_PrimaryParticleInitialEnergy_MeV_max = newValue;
	}
	else if(Control_PrimaryParticleInitialEnergy_MeV_min > newValue) Control_PrimaryParticleInitialEnergy_MeV_min = newValue;
	else if(Control_PrimaryParticleInitialEnergy_MeV_max < newValue) Control_PrimaryParticleInitialEnergy_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the secondary particles' initial energy (MeV).
 */
void UserRunInformation::SetMinMaxValue_Control_SecondaryParticleInitialEnergy_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_SecondaryParticleInitialEnergy_MeV_min) && isnan(Control_SecondaryParticleInitialEnergy_MeV_max))
	{
		Control_SecondaryParticleInitialEnergy_MeV_min = newValue;
		Control_SecondaryParticleInitialEnergy_MeV_max = newValue;
	}
	else if(Control_SecondaryParticleInitialEnergy_MeV_min > newValue) Control_SecondaryParticleInitialEnergy_MeV_min = newValue;
	else if(Control_SecondaryParticleInitialEnergy_MeV_max < newValue) Control_SecondaryParticleInitialEnergy_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the optical photons' initial energy (MeV).
 */
void UserRunInformation::SetMinMaxValue_Control_OpticalPhotonInitialEnergy_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_OpticalPhotonInitialEnergy_MeV_min) && isnan(Control_OpticalPhotonInitialEnergy_MeV_max))
	{
		Control_OpticalPhotonInitialEnergy_MeV_min = newValue;
		Control_OpticalPhotonInitialEnergy_MeV_max = newValue;
	}
	else if(Control_OpticalPhotonInitialEnergy_MeV_min > newValue) Control_OpticalPhotonInitialEnergy_MeV_min = newValue;
	else if(Control_OpticalPhotonInitialEnergy_MeV_max < newValue) Control_OpticalPhotonInitialEnergy_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the primary particles' initial beta * gamma.
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleInitialBetaGamma(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticleInitialBetaGamma_min) && isnan(Control_PrimaryParticleInitialBetaGamma_max))
	{
		Control_PrimaryParticleInitialBetaGamma_min = newValue;
		Control_PrimaryParticleInitialBetaGamma_max = newValue;
	}
	else if(Control_PrimaryParticleInitialBetaGamma_min > newValue) Control_PrimaryParticleInitialBetaGamma_min = newValue;
	else if(Control_PrimaryParticleInitialBetaGamma_max < newValue) Control_PrimaryParticleInitialBetaGamma_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the secondary particles' initial beta * gamma.
 */
void UserRunInformation::SetMinMaxValue_Control_SecondaryParticleInitialBetaGamma(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_SecondaryParticleInitialBetaGamma_min) && isnan(Control_SecondaryParticleInitialBetaGamma_max))
	{
		Control_SecondaryParticleInitialBetaGamma_min = newValue;
		Control_SecondaryParticleInitialBetaGamma_max = newValue;
	}
	else if(Control_SecondaryParticleInitialBetaGamma_min > newValue) Control_SecondaryParticleInitialBetaGamma_min = newValue;
	else if(Control_SecondaryParticleInitialBetaGamma_max < newValue) Control_SecondaryParticleInitialBetaGamma_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the primary particles' energy deposition (MeV) in scintillator material.
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min) && isnan(Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max))
	{
		Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min = newValue;
		Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max = newValue;
	}
	else if(Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min > newValue) Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min = newValue;
	else if(Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max < newValue) Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the primary particles' path length (mm) in scintillator material.
 */
void UserRunInformation::SetMinMaxValue_Control_PrimaryParticlePathLengthInScintillator_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_PrimaryParticlePathLengthInScintillator_MeV_min) && isnan(Control_PrimaryParticlePathLengthInScintillator_MeV_max))
	{
		Control_PrimaryParticlePathLengthInScintillator_MeV_min = newValue;
		Control_PrimaryParticlePathLengthInScintillator_MeV_max = newValue;
	}
	else if(Control_PrimaryParticlePathLengthInScintillator_MeV_min > newValue) Control_PrimaryParticlePathLengthInScintillator_MeV_min = newValue;
	else if(Control_PrimaryParticlePathLengthInScintillator_MeV_max < newValue) Control_PrimaryParticlePathLengthInScintillator_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the secondary particles' energy deposition (MeV) in scintillator material.
 */
void UserRunInformation::SetMinMaxValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min) && isnan(Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max))
	{
		Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min = newValue;
		Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max = newValue;
	}
	else if(Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min > newValue) Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min = newValue;
	else if(Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max < newValue) Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of primary particle's scintillator hit time (ns).
 */
void UserRunInformation::SetMinMaxValue_Control_GlobalScintiHitTime_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_GlobalScintiHitTime_ns_min) && isnan(Control_GlobalScintiHitTime_ns_max))
	{
		Control_GlobalScintiHitTime_ns_min = newValue;
		Control_GlobalScintiHitTime_ns_max = newValue;
	}
	else if(Control_GlobalScintiHitTime_ns_min > newValue) Control_GlobalScintiHitTime_ns_min = newValue;
	else if(Control_GlobalScintiHitTime_ns_max < newValue) Control_GlobalScintiHitTime_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the creation of secondary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitCreation_secondary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitCreation_secondary_ns_min) && isnan(Control_DeltaTimeHitCreation_secondary_ns_max))
	{
		Control_DeltaTimeHitCreation_secondary_ns_min = newValue;
		Control_DeltaTimeHitCreation_secondary_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitCreation_secondary_ns_min > newValue) Control_DeltaTimeHitCreation_secondary_ns_min = newValue;
	else if(Control_DeltaTimeHitCreation_secondary_ns_max < newValue) Control_DeltaTimeHitCreation_secondary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitCreation_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitCreation_ns_min) && isnan(Control_DeltaTimeHitCreation_ns_max))
	{
		Control_DeltaTimeHitCreation_ns_min = newValue;
		Control_DeltaTimeHitCreation_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitCreation_ns_min > newValue) Control_DeltaTimeHitCreation_ns_min = newValue;
	else if(Control_DeltaTimeHitCreation_ns_max < newValue) Control_DeltaTimeHitCreation_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the creation and the absorption of optical photons.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeCreationAbsorption_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeCreationAbsorption_ns_min) && isnan(Control_DeltaTimeCreationAbsorption_ns_max))
	{
		Control_DeltaTimeCreationAbsorption_ns_min = newValue;
		Control_DeltaTimeCreationAbsorption_ns_max = newValue;
	}
	else if(Control_DeltaTimeCreationAbsorption_ns_min > newValue) Control_DeltaTimeCreationAbsorption_ns_min = newValue;
	else if(Control_DeltaTimeCreationAbsorption_ns_max < newValue) Control_DeltaTimeCreationAbsorption_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitAbsorption_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitAbsorption_ns_min) && isnan(Control_DeltaTimeHitAbsorption_ns_max))
	{
		Control_DeltaTimeHitAbsorption_ns_min = newValue;
		Control_DeltaTimeHitAbsorption_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitAbsorption_ns_min > newValue) Control_DeltaTimeHitAbsorption_ns_min = newValue;
	else if(Control_DeltaTimeHitAbsorption_ns_max < newValue) Control_DeltaTimeHitAbsorption_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons that have been created by primary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitCreation_parentPrimary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitCreation_parentPrimary_ns_min) && isnan(Control_DeltaTimeHitCreation_parentPrimary_ns_max))
	{
		Control_DeltaTimeHitCreation_parentPrimary_ns_min = newValue;
		Control_DeltaTimeHitCreation_parentPrimary_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitCreation_parentPrimary_ns_min > newValue) Control_DeltaTimeHitCreation_parentPrimary_ns_min = newValue;
	else if(Control_DeltaTimeHitCreation_parentPrimary_ns_max < newValue) Control_DeltaTimeHitCreation_parentPrimary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the creation and the absorption of optical photons that have been created by primary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min) && isnan(Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max))
	{
		Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min = newValue;
		Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max = newValue;
	}
	else if(Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min > newValue) Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min = newValue;
	else if(Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max < newValue) Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons that have been created by primary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitAbsorption_parentPrimary_ns_min) && isnan(Control_DeltaTimeHitAbsorption_parentPrimary_ns_max))
	{
		Control_DeltaTimeHitAbsorption_parentPrimary_ns_min = newValue;
		Control_DeltaTimeHitAbsorption_parentPrimary_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitAbsorption_parentPrimary_ns_min > newValue) Control_DeltaTimeHitAbsorption_parentPrimary_ns_min = newValue;
	else if(Control_DeltaTimeHitAbsorption_parentPrimary_ns_max < newValue) Control_DeltaTimeHitAbsorption_parentPrimary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons that have been created by secondary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitCreation_parentSecondary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitCreation_parentSecondary_ns_min) && isnan(Control_DeltaTimeHitCreation_parentSecondary_ns_max))
	{
		Control_DeltaTimeHitCreation_parentSecondary_ns_min = newValue;
		Control_DeltaTimeHitCreation_parentSecondary_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitCreation_parentSecondary_ns_min > newValue) Control_DeltaTimeHitCreation_parentSecondary_ns_min = newValue;
	else if(Control_DeltaTimeHitCreation_parentSecondary_ns_max < newValue) Control_DeltaTimeHitCreation_parentSecondary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the creation and the absorption of optical photons that have been created by secondary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min) && isnan(Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max))
	{
		Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min = newValue;
		Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max = newValue;
	}
	else if(Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min > newValue) Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min = newValue;
	else if(Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max < newValue) Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons that have been created by secondary particles.
 */
void UserRunInformation::SetMinMaxValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns(G4double newValue)
{
	if(isnan(newValue)) return;

	if(isnan(Control_DeltaTimeHitAbsorption_parentSecondary_ns_min) && isnan(Control_DeltaTimeHitAbsorption_parentSecondary_ns_max))
	{
		Control_DeltaTimeHitAbsorption_parentSecondary_ns_min = newValue;
		Control_DeltaTimeHitAbsorption_parentSecondary_ns_max = newValue;
	}
	else if(Control_DeltaTimeHitAbsorption_parentSecondary_ns_min > newValue) Control_DeltaTimeHitAbsorption_parentSecondary_ns_min = newValue;
	else if(Control_DeltaTimeHitAbsorption_parentSecondary_ns_max < newValue) Control_DeltaTimeHitAbsorption_parentSecondary_ns_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons per event.
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotons(G4int newValue)
{
	if(Control_OpticalPhotons_min == -1 && Control_OpticalPhotons_max == -1)
	{
		Control_OpticalPhotons_min = newValue;
		Control_OpticalPhotons_max = newValue;
	}
	else if(Control_OpticalPhotons_min > newValue) Control_OpticalPhotons_min = newValue;
	else if(Control_OpticalPhotons_max < newValue) Control_OpticalPhotons_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons that have been created by primary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotonsFromPrimary(G4int newValue)
{
	if(Control_OpticalPhotonsFromPrimary_min == -1 && Control_OpticalPhotonsFromPrimary_max == -1)
	{
		Control_OpticalPhotonsFromPrimary_min = newValue;
		Control_OpticalPhotonsFromPrimary_max = newValue;
	}
	else if(Control_OpticalPhotonsFromPrimary_min > newValue) Control_OpticalPhotonsFromPrimary_min = newValue;
	else if(Control_OpticalPhotonsFromPrimary_max < newValue) Control_OpticalPhotonsFromPrimary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons that have been created by secondary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotonsFromSecondary(G4int newValue)
{
	if(Control_OpticalPhotonsFromSecondary_min == -1 && Control_OpticalPhotonsFromSecondary_max == -1)
	{
		Control_OpticalPhotonsFromSecondary_min = newValue;
		Control_OpticalPhotonsFromSecondary_max = newValue;
	}
	else if(Control_OpticalPhotonsFromSecondary_min > newValue) Control_OpticalPhotonsFromSecondary_min = newValue;
	else if(Control_OpticalPhotonsFromSecondary_max < newValue) Control_OpticalPhotonsFromSecondary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of Cerenkov photons per event.
 */
void UserRunInformation::SetMinMaxValue_Control_CerenkovPhotons(G4int newValue)
{
	if(Control_CerenkovPhotons_min == -1 && Control_CerenkovPhotons_max == -1)
	{
		Control_CerenkovPhotons_min = newValue;
		Control_CerenkovPhotons_max = newValue;
	}
	else if(Control_CerenkovPhotons_min > newValue) Control_CerenkovPhotons_min = newValue;
	else if(Control_CerenkovPhotons_max < newValue) Control_CerenkovPhotons_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of scintillation photons per event.
 */
void UserRunInformation::SetMinMaxValue_Control_ScintiPhotons(G4int newValue)
{
	if(Control_ScintiPhotons_min == -1 && Control_ScintiPhotons_max == -1)
	{
		Control_ScintiPhotons_min = newValue;
		Control_ScintiPhotons_max = newValue;
	}
	else if(Control_ScintiPhotons_min > newValue) Control_ScintiPhotons_min = newValue;
	else if(Control_ScintiPhotons_max < newValue) Control_ScintiPhotons_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons from a wavelength-shifting process per event.
 */
void UserRunInformation::SetMinMaxValue_Control_WLSPhotons(G4int newValue)
{
	if(Control_WLSPhotons_min == -1 && Control_WLSPhotons_max == -1)
	{
		Control_WLSPhotons_min = newValue;
		Control_WLSPhotons_max = newValue;
	}
	else if(Control_WLSPhotons_min > newValue) Control_WLSPhotons_min = newValue;
	else if(Control_WLSPhotons_max < newValue) Control_WLSPhotons_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of Cerenkov photons that have been created by primary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_CerenkovPhotonsFromPrimary(G4int newValue)
{
	if(Control_CerenkovPhotonsFromPrimary_min == -1 && Control_CerenkovPhotonsFromPrimary_max == -1)
	{
		Control_CerenkovPhotonsFromPrimary_min = newValue;
		Control_CerenkovPhotonsFromPrimary_max = newValue;
	}
	else if(Control_CerenkovPhotonsFromPrimary_min > newValue) Control_CerenkovPhotonsFromPrimary_min = newValue;
	else if(Control_CerenkovPhotonsFromPrimary_max < newValue) Control_CerenkovPhotonsFromPrimary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of scintillation photons that have been created by primary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_ScintiPhotonsFromPrimary(G4int newValue)
{
	if(Control_ScintiPhotonsFromPrimary_min == -1 && Control_ScintiPhotonsFromPrimary_max == -1)
	{
		Control_ScintiPhotonsFromPrimary_min = newValue;
		Control_ScintiPhotonsFromPrimary_max = newValue;
	}
	else if(Control_ScintiPhotonsFromPrimary_min > newValue) Control_ScintiPhotonsFromPrimary_min = newValue;
	else if(Control_ScintiPhotonsFromPrimary_max < newValue) Control_ScintiPhotonsFromPrimary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons from a wavelength-shifting process that have been created by primary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_WLSPhotonsFromPrimary(G4int newValue)
{
	if(Control_WLSPhotonsFromPrimary_min == -1 && Control_WLSPhotonsFromPrimary_max == -1)
	{
		Control_WLSPhotonsFromPrimary_min = newValue;
		Control_WLSPhotonsFromPrimary_max = newValue;
	}
	else if(Control_WLSPhotonsFromPrimary_min > newValue) Control_WLSPhotonsFromPrimary_min = newValue;
	else if(Control_WLSPhotonsFromPrimary_max < newValue) Control_WLSPhotonsFromPrimary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of Cerenkov photons that have been created by secondary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_CerenkovPhotonsFromSecondary(G4int newValue)
{
	if(Control_CerenkovPhotonsFromSecondary_min == -1 && Control_CerenkovPhotonsFromSecondary_max == -1)
	{
		Control_CerenkovPhotonsFromSecondary_min = newValue;
		Control_CerenkovPhotonsFromSecondary_max = newValue;
	}
	else if(Control_CerenkovPhotonsFromSecondary_min > newValue) Control_CerenkovPhotonsFromSecondary_min = newValue;
	else if(Control_CerenkovPhotonsFromSecondary_max < newValue) Control_CerenkovPhotonsFromSecondary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of scintillation photons that have been created by secondary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_ScintiPhotonsFromSecondary(G4int newValue)
{
	if(Control_ScintiPhotonsFromSecondary_min == -1 && Control_ScintiPhotonsFromSecondary_max == -1)
	{
		Control_ScintiPhotonsFromSecondary_min = newValue;
		Control_ScintiPhotonsFromSecondary_max = newValue;
	}
	else if(Control_ScintiPhotonsFromSecondary_min > newValue) Control_ScintiPhotonsFromSecondary_min = newValue;
	else if(Control_ScintiPhotonsFromSecondary_max < newValue) Control_ScintiPhotonsFromSecondary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons from a wavelength-shifting process that have been created by secondary particles (per event).
 */
void UserRunInformation::SetMinMaxValue_Control_WLSPhotonsFromSecondary(G4int newValue)
{
	if(Control_WLSPhotonsFromSecondary_min == -1 && Control_WLSPhotonsFromSecondary_max == -1)
	{
		Control_WLSPhotonsFromSecondary_min = newValue;
		Control_WLSPhotonsFromSecondary_max = newValue;
	}
	else if(Control_WLSPhotonsFromSecondary_min > newValue) Control_WLSPhotonsFromSecondary_min = newValue;
	else if(Control_WLSPhotonsFromSecondary_max < newValue) Control_WLSPhotonsFromSecondary_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons that have been created by primary particles (per primary particle's energy deposition).
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition(G4double newValue)
{
	if(isnan(newValue)) return;

	if(Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min == -1 && Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max == -1)
	{
		Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min = newValue;
		Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max = newValue;
	}
	else if(Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min > newValue) Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min = newValue;
	else if(Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max < newValue) Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max = newValue;
}

/**
 *  Function to set the minimal and maximal number of optical photons that have been created by secondary particles (per primary particle's energy deposition).
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition(G4double newValue)
{
	if(isnan(newValue)) return;

	if(Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min == -1 && Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max == -1)
	{
		Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min = newValue;
		Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max = newValue;
	}
	else if(Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min > newValue) Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min = newValue;
	else if(Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max < newValue) Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max = newValue;
}

/**
 *  Function to set the minimal and maximal fraction of optical photons that have been absorbed per event (sorted by the volumes they were absorbed in).
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotonsAbsorbed_volume(G4int num_absorbed_inScinti, G4int num_absorbed_inFibre, G4int num_absorbed_inPhotonDetector, G4int num_total)
{
	if(isnan(num_absorbed_inScinti) || isnan(num_absorbed_inFibre) || isnan(num_absorbed_inPhotonDetector) || isnan(num_total)) return;


	G4double frac_min = ((G4double) num_absorbed_inScinti) / ((G4double) num_total);
	if(frac_min > ((G4double) num_absorbed_inFibre) / ((G4double) num_total)) frac_min = ((G4double) num_absorbed_inFibre) / ((G4double) num_total);
	if(frac_min > ((G4double) num_absorbed_inPhotonDetector) / ((G4double) num_total)) frac_min = ((G4double) num_absorbed_inPhotonDetector) / ((G4double) num_total);

	if(Control_OpticalPhotonsAbsorbed_volume_min == -1 || Control_OpticalPhotonsAbsorbed_volume_min > frac_min) Control_OpticalPhotonsAbsorbed_volume_min = frac_min;


	G4double frac_max = ((G4double) (num_absorbed_inScinti + num_absorbed_inFibre + num_absorbed_inPhotonDetector)) / ((G4double) num_total);

	if(Control_OpticalPhotonsAbsorbed_volume_max == -1 || Control_OpticalPhotonsAbsorbed_volume_max < frac_max) Control_OpticalPhotonsAbsorbed_volume_max = frac_max;
}

/**
 *  Function to set the minimal and maximal fraction of optical photons that have been absorbed per event (sorted by their production process).
 */
void UserRunInformation::SetMinMaxValue_Control_opticalPhotonsAbsorbed_process(G4int num_Scinti_absorbed, G4int num_Scinti_total, G4int num_Cerenkov_absorbed, G4int num_Cerenkov_total, G4int num_WLS_absorbed, G4int num_WLS_total)
{
	if(isnan(num_Scinti_absorbed) || isnan(num_Scinti_total) || isnan(num_Cerenkov_absorbed) || isnan(num_Cerenkov_total) || isnan(num_WLS_absorbed) || isnan(num_WLS_total)) return;

	G4int num_absorbed = num_Scinti_absorbed + num_Cerenkov_absorbed + num_WLS_absorbed;
	G4int num_total = num_Scinti_total + num_Cerenkov_total + num_WLS_total;

	G4double frac_min = -1.;
	if(num_Scinti_total) frac_min = ((G4double) num_Scinti_absorbed) / ((G4double) num_Scinti_total);
	if(num_Cerenkov_total && frac_min > ((G4double) num_Cerenkov_absorbed) / ((G4double) num_Cerenkov_total)) frac_min = ((G4double) num_Cerenkov_absorbed) / ((G4double) num_Cerenkov_total);
	if(num_WLS_total && frac_min > ((G4double) num_WLS_absorbed) / ((G4double) num_WLS_total)) frac_min = ((G4double) num_WLS_absorbed) / ((G4double) num_WLS_total);
	if(num_total && frac_min > ((G4double) num_absorbed) / ((G4double) num_total)) frac_min = ((G4double) num_absorbed) / ((G4double) num_total);

	if(Control_OpticalPhotonsAbsorbed_process_min == -1 || Control_OpticalPhotonsAbsorbed_process_min > frac_min) Control_OpticalPhotonsAbsorbed_process_min = frac_min;

	G4double frac_max = -1.;
	if(num_Scinti_total) frac_max = ((G4double) num_Scinti_absorbed) / ((G4double) num_Scinti_total);
	if(num_Cerenkov_total && frac_max < ((G4double) num_Cerenkov_absorbed) / ((G4double) num_Cerenkov_total)) frac_max = ((G4double) num_Cerenkov_absorbed) / ((G4double) num_Cerenkov_total);
	if(num_WLS_total && frac_max < ((G4double) num_WLS_absorbed) / ((G4double) num_WLS_total)) frac_max = ((G4double) num_WLS_absorbed) / ((G4double) num_WLS_total);
	if(num_total && frac_max < ((G4double) num_absorbed) / ((G4double) num_total)) frac_max = ((G4double) num_absorbed) / ((G4double) num_total);

	if(Control_OpticalPhotonsAbsorbed_process_max == -1 || Control_OpticalPhotonsAbsorbed_process_max < frac_max) Control_OpticalPhotonsAbsorbed_process_max = frac_max;
}
