/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef UserRunInformation_h
#define UserRunInformation_h 1

#include <globals.hh>

#include <map>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Class to save the data which belonging to the current run.
class UserRunInformation
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	UserRunInformation()
	{
		ResetDefaults();
		MassMap.clear();
	}

	/**
	 *  Destructor (empty)
	 */
	~UserRunInformation()
	{
	}

	void ResetDefaults();


	/**
	 *  Function to set the number of events.
	 */
	void SetNumberOfEvents(G4int newValue)
	{ NumberOfEvents = newValue; }
	/**
	 *  @return number of events
	 */
	G4int GetNumberOfEvents() const
	{ return NumberOfEvents; }

	/**
	 *  Function to increment the total number of primary particles in the run.
	 */
	void IncNumberOfPrimaryParticles(G4int npp)
	{ NumberOfPrimaryParticles += npp; }
	/**
	 *  @return total number of primary particles in the run
	 */
	G4int GetNumberOfPrimaryParticles() const
	{ return NumberOfPrimaryParticles; }

	/**
	 *  Function to save the particle mass.
	 */
	void SetParticleMass(G4String particleName, G4double particleMass)
	{ MassMap[particleName] = particleMass; }
	G4double GetParticleMass(G4String particleName) const;
	void DumpParticleMassesToFile(G4String pathToFile);

// for data plots
	void SetMinMaxValue_ScintiHitPointX_mm(G4double newValue);
	/**
	 *  @return minimal X value of the scintillator hit position (mm)
	 */
	G4double GetMinValue_ScintiHitPointX_mm() const
	{ return ScintiHitPointX_mm_min; }
	/**
	 *  @return maximal X value of the scintillator hit position (mm)
	 */
	G4double GetMaxValue_ScintiHitPointX_mm() const
	{ return ScintiHitPointX_mm_max; }

	void SetMinMaxValue_ScintiHitPointY_mm(G4double newValue);
	/**
	 *  @return minimal Y value of the scintillator hit position (mm)
	 */
	G4double GetMinValue_ScintiHitPointY_mm() const
	{ return ScintiHitPointY_mm_min; }
	/**
	 *  @return maximal Y value of the scintillator hit position (mm)
	 */
	G4double GetMaxValue_ScintiHitPointY_mm() const
	{ return ScintiHitPointY_mm_max; }

	void SetMinMaxValue_ScintiHitPointZ_mm(G4double newValue);
	/**
	 *  @return minimal Z value of the scintillator hit position (mm)
	 */
	G4double GetMinValue_ScintiHitPointZ_mm() const
	{ return ScintiHitPointZ_mm_min; }
	/**
	 *  @return maximal Z value of the scintillator hit position (mm)
	 */
	G4double GetMaxValue_ScintiHitPointZ_mm() const
	{ return ScintiHitPointZ_mm_max; }

	void SetMinMaxValue_opticalPhotonsPerPrimary(G4int newValue);
	/**
	 *  @return minimal number of optical photons per primary particles (= per event)
	 */
	G4int GetMinValue_opticalPhotonsPerPrimary() const
	{ return OpticalPhotonsPerPrimary_min; }
	/**
	 *  @return maximal number of optical photons per primary particles (= per event)
	 */
	G4int GetMaxValue_opticalPhotonsPerPrimary() const
	{ return OpticalPhotonsPerPrimary_max; }

	void SetMinMaxValue_opticalPhotonsPerEnergyDeposition(G4double newValue);
	/**
	 *  @return minimal number of optical photons per primary particle's energy deposition
	 */
	G4double GetMinValue_opticalPhotonsPerEnergyDeposition() const
	{ return OpticalPhotonsPerEnergyDeposition_min; }
	/**
	 *  @return maximal number of optical photons per primary particle's energy deposition
	 */
	G4double GetMaxValue_opticalPhotonsPerEnergyDeposition() const
	{ return OpticalPhotonsPerEnergyDeposition_max; }

	void SetMinMaxValue_opticalPhotonsAbsorbedInPhotonDetector(G4int newValue);
	/**
	 *  @return minimal number of optical photons that have been absorbed in a photon detector (per event)
	 */
	G4int GetMinValue_opticalPhotonsAbsorbedInPhotonDetector() const
	{ return OpticalPhotonsAbsorbedInPhotonDetector_min; }
	/**
	 *  @return maximal number of optical photons that have been absorbed in a photon detector (per event)
	 */
	G4int GetMaxValue_opticalPhotonsAbsorbedInPhotonDetector() const
	{ return OpticalPhotonsAbsorbedInPhotonDetector_max; }

	void SetMinMaxValue_PhotonDetectorHitPointX_mm(G4double newValue);
	/**
	 *  @return minimal X value of the photon detector hit position (mm)
	 */
	G4double GetMinValue_PhotonDetectorHitPointX_mm() const
	{ return PhotonDetectorHitPointX_mm_min; }
	/**
	 *  @return maximal X value of the photon detector hit position (mm)
	 */
	G4double GetMaxValue_PhotonDetectorHitPointX_mm() const
	{ return PhotonDetectorHitPointX_mm_max; }

	void SetMinMaxValue_PhotonDetectorHitPointY_mm(G4double newValue);
	/**
	 *  @return minimal Y value of the photon detector hit position (mm)
	 */
	G4double GetMinValue_PhotonDetectorHitPointY_mm() const
	{ return PhotonDetectorHitPointY_mm_min; }
	/**
	 *  @return maximal Y value of the photon detector hit position (mm)
	 */
	G4double GetMaxValue_PhotonDetectorHitPointY_mm() const
	{ return PhotonDetectorHitPointY_mm_max; }

// 	void SetMinMaxValue_PhotonDetectorHitPointZ_mm(G4double newValue);
// 	G4double GetMinValue_PhotonDetectorHitPointZ_mm() const
// 	{ return PhotonDetectorHitPointZ_mm_min; }
// 	G4double GetMaxValue_PhotonDetectorHitPointZ_mm() const
// 	{ return PhotonDetectorHitPointZ_mm_max; }

	void SetMinMaxValue_OpticalPhotonHitEnergy_MeV(G4double newValue);
	/**
	 *  @return minimal energy (MeV) of optical photons that have been absorbed in a photon detector
	 */
	G4double GetMinValue_OpticalPhotonHitEnergy_MeV() const
	{ return OpticalPhotonHitEnergy_MeV_min; }
	/**
	 *  @return maximal energy (MeV) of optical photons that have been absorbed in a photon detector
	 */
	G4double GetMaxValue_OpticalPhotonHitEnergy_MeV() const
	{ return OpticalPhotonHitEnergy_MeV_max; }

	void SetMinMaxValue_DeltaTimeHitAbsorption_ns(G4double newValue);
	/**
	 *  @return minimal value of the optical photon's absorption time (ns) relative to the primary particle's scintillator hit time
	 */
	G4double GetMinValue_DeltaTimeHitAbsorption_ns() const
	{ return DeltaTimeHitAbsorption_ns_min; }
	/**
	 *  @return maximal value of the optical photon's absorption time (ns) relative to the primary particle's scintillator hit time
	 */
	G4double GetMaxValue_DeltaTimeHitAbsorption_ns() const
	{ return DeltaTimeHitAbsorption_ns_max; }

// for control plots
	void SetMinMaxValue_Control_PrimaryParticleID(G4int newValue);
	/**
	 *  @return minimal value of the primary particles' particle ID
	 */
	G4int GetMinValue_Control_PrimaryParticleID() const
	{ return Control_PrimaryParticleID_min; }
	/**
	 *  @return maximal value of the primary particles' particle ID
	 */
	G4int GetMaxValue_Control_PrimaryParticleID() const
	{ return Control_PrimaryParticleID_max; }

	void SetMinMaxValue_Control_SecondaryParticleID(G4int newValue);
	/**
	 *  @return minimal value of the secondary particles' particle ID
	 */
	G4int GetMinValue_Control_SecondaryParticleID() const
	{ return Control_SecondaryParticleID_min; }
	/**
	 *  @return maximal value of the secondary particles' particle ID
	 */
	G4int GetMaxValue_Control_SecondaryParticleID() const
	{ return Control_SecondaryParticleID_max; }

	void SetMinMaxValue_Control_PrimaryParticleInitialPositionX_mm(G4double newValue);
	/**
	 *  @return minimal X value of the primary particles' initial position (mm)
	 */
	G4double GetMinValue_Control_PrimaryParticleInitialPositionX_mm() const
	{ return Control_PrimaryParticleInitialPositionX_mm_min; }
	/**
	 *  @return maximal X value of the primary particles' initial position (mm)
	 */
	G4double GetMaxValue_Control_PrimaryParticleInitialPositionX_mm() const
	{ return Control_PrimaryParticleInitialPositionX_mm_max; }

	void SetMinMaxValue_Control_PrimaryParticleInitialPositionY_mm(G4double newValue);
	/**
	 *  @return minimal Y value of the primary particles' initial position (mm)
	 */
	G4double GetMinValue_Control_PrimaryParticleInitialPositionY_mm() const
	{ return Control_PrimaryParticleInitialPositionY_mm_min; }
	/**
	 *  @return maximal Y value of the primary particles' initial position (mm)
	 */
	G4double GetMaxValue_Control_PrimaryParticleInitialPositionY_mm() const
	{ return Control_PrimaryParticleInitialPositionY_mm_max; }

	void SetMinMaxValue_Control_PrimaryParticleInitialPositionZ_mm(G4double newValue);
	/**
	 *  @return minimal Z value of the primary particles' initial position (mm)
	 */
	G4double GetMinValue_Control_PrimaryParticleInitialPositionZ_mm() const
	{ return Control_PrimaryParticleInitialPositionZ_mm_min; }
	/**
	 *  @return maximal Z value of the primary particles' initial position (mm)
	 */
	G4double GetMaxValue_Control_PrimaryParticleInitialPositionZ_mm() const
	{ return Control_PrimaryParticleInitialPositionZ_mm_max; }

	void SetMinMaxValue_Control_ScintiHitPointX_mm(G4double newValue);
	/**
	 *  @return minimal X value of the primary particle's scintillator hit position (mm)
	 */
	G4double GetMinValue_Control_ScintiHitPointX_mm() const
	{ return Control_ScintiHitPointX_mm_min; }
	/**
	 *  @return maximal X value of the primary particle's scintillator hit position (mm)
	 */
	G4double GetMaxValue_Control_ScintiHitPointX_mm() const
	{ return Control_ScintiHitPointX_mm_max; }

	void SetMinMaxValue_Control_ScintiHitPointY_mm(G4double newValue);
	/**
	 *  @return minimal Y value of the primary particle's scintillator hit position (mm)
	 */
	G4double GetMinValue_Control_ScintiHitPointY_mm() const
	{ return Control_ScintiHitPointY_mm_min; }
	/**
	 *  @return maximal Y value of the primary particle's scintillator hit position (mm)
	 */
	G4double GetMaxValue_Control_ScintiHitPointY_mm() const
	{ return Control_ScintiHitPointY_mm_max; }

	void SetMinMaxValue_Control_ScintiHitPointZ_mm(G4double newValue);
	/**
	 *  @return minimal Z value of the primary particle's scintillator hit position (mm)
	 */
	G4double GetMinValue_Control_ScintiHitPointZ_mm() const
	{ return Control_ScintiHitPointZ_mm_min; }
	/**
	 *  @return maximal Z value of the primary particle's scintillator hit position (mm)
	 */
	G4double GetMaxValue_Control_ScintiHitPointZ_mm() const
	{ return Control_ScintiHitPointZ_mm_max; }

	void SetMinMaxValue_Control_PrimaryParticleInitialEnergy_MeV(G4double newValue);
	/**
	 *  @return minimal value of the primary particles' initial energy (MeV)
	 */
	G4double GetMinValue_Control_PrimaryParticleInitialEnergy_MeV() const
	{ return Control_PrimaryParticleInitialEnergy_MeV_min; }
	/**
	 *  @return maximal value of the primary particles' initial energy (MeV)
	 */
	G4double GetMaxValue_Control_PrimaryParticleInitialEnergy_MeV() const
	{ return Control_PrimaryParticleInitialEnergy_MeV_max; }

	void SetMinMaxValue_Control_SecondaryParticleInitialEnergy_MeV(G4double newValue);
	/**
	 *  @return minimal value of the secondary particles' initial energy (MeV)
	 */
	G4double GetMinValue_Control_SecondaryParticleInitialEnergy_MeV() const
	{ return Control_SecondaryParticleInitialEnergy_MeV_min; }
	/**
	 *  @return maximal value of the secondary particles' initial energy (MeV)
	 */
	G4double GetMaxValue_Control_SecondaryParticleInitialEnergy_MeV() const
	{ return Control_SecondaryParticleInitialEnergy_MeV_max; }

	void SetMinMaxValue_Control_OpticalPhotonInitialEnergy_MeV(G4double newValue);
	/**
	 *  @return minimal value of the optical photons' initial energy (MeV)
	 */
	G4double GetMinValue_Control_OpticalPhotonInitialEnergy_MeV() const
	{ return Control_OpticalPhotonInitialEnergy_MeV_min; }
	/**
	 *  @return maximal value of the optical photons' initial energy (MeV)
	 */
	G4double GetMaxValue_Control_OpticalPhotonInitialEnergy_MeV() const
	{ return Control_OpticalPhotonInitialEnergy_MeV_max; }

	void SetMinMaxValue_Control_PrimaryParticleInitialBetaGamma(G4double newValue);
	/**
	 *  @return minimal value of the primary particles' initial beta * gamma
	 */
	G4double GetMinValue_Control_PrimaryParticleInitialBetaGamma() const
	{ return Control_PrimaryParticleInitialBetaGamma_min; }
	/**
	 *  @return maximal value of the primary particles' initial beta * gamma
	 */
	G4double GetMaxValue_Control_PrimaryParticleInitialBetaGamma() const
	{ return Control_PrimaryParticleInitialBetaGamma_max; }

	void SetMinMaxValue_Control_SecondaryParticleInitialBetaGamma(G4double newValue);
	/**
	 *  @return minimal value of the secondary particles' initial beta * gamma
	 */
	G4double GetMinValue_Control_SecondaryParticleInitialBetaGamma() const
	{ return Control_SecondaryParticleInitialBetaGamma_min; }
	/**
	 *  @return maximal value of the secondary particles' initial beta * gamma
	 */
	G4double GetMaxValue_Control_SecondaryParticleInitialBetaGamma() const
	{ return Control_SecondaryParticleInitialBetaGamma_max; }

	void SetMinMaxValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV(G4double newValue);
	/**
	 *  @return minimal value of the primary particles' energy deposition (MeV) in scintillator material
	 */
	G4double GetMinValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV() const
	{ return Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min; }
	/**
	 *  @return maximal value of the primary particles' energy deposition (MeV) in scintillator material
	 */
	G4double GetMaxValue_Control_PrimaryParticleEnergyDepositionInScintillator_MeV() const
	{ return Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max; }

	void SetMinMaxValue_Control_PrimaryParticlePathLengthInScintillator_MeV(G4double newValue);
	/**
	 *  @return minimal value of the primary particles' energy deposition (MeV) in scintillator material
	 */
	G4double GetMinValue_Control_PrimaryParticlePathLengthInScintillator_MeV() const
	{ return Control_PrimaryParticlePathLengthInScintillator_MeV_min; }
	/**
	 *  @return maximal value of the primary particles' energy deposition (MeV) in scintillator material
	 */
	G4double GetMaxValue_Control_PrimaryParticlePathLengthInScintillator_MeV() const
	{ return Control_PrimaryParticlePathLengthInScintillator_MeV_max; }

	void SetMinMaxValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV(G4double newValue);
	/**
	 *  @return minimal value of the secondary particles' energy deposition (MeV) in scintillator material
	 */
	G4double GetMinValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV() const
	{ return Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min; }
	/**
	 *  @return maximal value of the secondary particles' energy deposition (MeV) in scintillator material
	 */
	G4double GetMaxValue_Control_SecondaryParticleEnergyDepositionInScintillator_MeV() const
	{ return Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max; }

	void SetMinMaxValue_Control_GlobalScintiHitTime_ns(G4double newValue);
	/**
	 *  @return minimal value of primary particle's scintillator hit time (ns)
	 */
	G4double GetMinValue_Control_GlobalScintiHitTime_ns() const
	{ return Control_GlobalScintiHitTime_ns_min; }
	/**
	 *  @return maximal value of primary particle's scintillator hit time (ns)
	 */
	G4double GetMaxValue_Control_GlobalScintiHitTime_ns() const
	{ return Control_GlobalScintiHitTime_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitCreation_secondary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the creation of secondary particles
	 */
	G4double GetMinValue_Control_DeltaTimeHitCreation_secondary_ns() const
	{ return Control_DeltaTimeHitCreation_secondary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the creation of secondary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeHitCreation_secondary_ns() const
	{ return Control_DeltaTimeHitCreation_secondary_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitCreation_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons
	 */
	G4double GetMinValue_Control_DeltaTimeHitCreation_ns() const
	{ return Control_DeltaTimeHitCreation_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons
	 */
	G4double GetMaxValue_Control_DeltaTimeHitCreation_ns() const
	{ return Control_DeltaTimeHitCreation_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeCreationAbsorption_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the creation and the absorption of optical photons
	 */
	G4double GetMinValue_Control_DeltaTimeCreationAbsorption_ns() const
	{ return Control_DeltaTimeCreationAbsorption_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the creation and the absorption of optical photons
	 */
	G4double GetMaxValue_Control_DeltaTimeCreationAbsorption_ns() const
	{ return Control_DeltaTimeCreationAbsorption_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitAbsorption_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons
	 */
	G4double GetMinValue_Control_DeltaTimeHitAbsorption_ns() const
	{ return Control_DeltaTimeHitAbsorption_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons
	 */
	G4double GetMaxValue_Control_DeltaTimeHitAbsorption_ns() const
	{ return Control_DeltaTimeHitAbsorption_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitCreation_parentPrimary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons that have been created by primary particles
	 */
	G4double GetMinValue_Control_DeltaTimeHitCreation_parentPrimary_ns() const
	{ return Control_DeltaTimeHitCreation_parentPrimary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons that have been created by primary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeHitCreation_parentPrimary_ns() const
	{ return Control_DeltaTimeHitCreation_parentPrimary_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the creation and the absorption of optical photons that have been created by primary particles
	 */
	G4double GetMinValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns() const
	{ return Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the creation and the absorption of optical photons that have been created by primary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeCreationAbsorption_parentPrimary_ns() const
	{ return Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons that have been created by primary particles
	 */
	G4double GetMinValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns() const
	{ return Control_DeltaTimeHitAbsorption_parentPrimary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons that have been created by primary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeHitAbsorption_parentPrimary_ns() const
	{ return Control_DeltaTimeHitAbsorption_parentPrimary_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitCreation_parentSecondary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons that have been created by secondary particles
	 */
	G4double GetMinValue_Control_DeltaTimeHitCreation_parentSecondary_ns() const
	{ return Control_DeltaTimeHitCreation_parentSecondary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the creation of optical photons that have been created by secondary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeHitCreation_parentSecondary_ns() const
	{ return Control_DeltaTimeHitCreation_parentSecondary_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the creation and the absorption of optical photons that have been created by secondary particles
	 */
	G4double GetMinValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns() const
	{ return Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the creation and the absorption of optical photons that have been created by secondary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeCreationAbsorption_parentSecondary_ns() const
	{ return Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max; }

	void SetMinMaxValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns(G4double newValue);
	/**
	 *  @return minimal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons that have been created by secondary particles
	 */
	G4double GetMinValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns() const
	{ return Control_DeltaTimeHitAbsorption_parentSecondary_ns_min; }
	/**
	 *  @return maximal value of the time (ns) between the primary particle's scintillator hit and the absorption of optical photons that have been created by secondary particles
	 */
	G4double GetMaxValue_Control_DeltaTimeHitAbsorption_parentSecondary_ns() const
	{ return Control_DeltaTimeHitAbsorption_parentSecondary_ns_max; }

	void SetMinMaxValue_Control_opticalPhotons(G4int newValue);
	/**
	 *  @return minimal number of optical photons per event
	 */
	G4int GetMinValue_Control_opticalPhotons() const
	{ return Control_OpticalPhotons_min; }
	/**
	 *  @return maximal number of optical photons per event
	 */
	G4int GetMaxValue_Control_opticalPhotons() const
	{ return Control_OpticalPhotons_max; }

	void SetMinMaxValue_Control_opticalPhotonsFromPrimary(G4int newValue);
	/**
	 *  @return minimal number of optical photons that have been created by primary particles (per event)
	 */
	G4int GetMinValue_Control_opticalPhotonsFromPrimary() const
	{ return Control_OpticalPhotonsFromPrimary_min; }
	/**
	 *  @return maximal number of optical photons that have been created by primary particles (per event)
	 */
	G4int GetMaxValue_Control_opticalPhotonsFromPrimary() const
	{ return Control_OpticalPhotonsFromPrimary_max; }

	void SetMinMaxValue_Control_opticalPhotonsFromSecondary(G4int newValue);
	/**
	 *  @return minimal number of optical photons that have been created by secondary particles (per event)
	 */
	G4int GetMinValue_Control_opticalPhotonsFromSecondary() const
	{ return Control_OpticalPhotonsFromSecondary_min; }
	/**
	 *  @return maximal number of optical photons that have been created by secondary particles (per event)
	 */
	G4int GetMaxValue_Control_opticalPhotonsFromSecondary() const
	{ return Control_OpticalPhotonsFromSecondary_max; }

	void SetMinMaxValue_Control_CerenkovPhotons(G4int newValue);
	/**
	 *  @return minimal number of Cerenkov photons per event
	 */
	G4int GetMinValue_Control_CerenkovPhotons() const
	{ return Control_CerenkovPhotons_min; }
	/**
	 *  @return maximal number of Cerenkov photons per event
	 */
	G4int GetMaxValue_Control_CerenkovPhotons() const
	{ return Control_CerenkovPhotons_max; }

	void SetMinMaxValue_Control_ScintiPhotons(G4int newValue);
	/**
	 *  @return minimal number of scintillation photons per event
	 */
	G4int GetMinValue_Control_ScintiPhotons() const
	{ return Control_ScintiPhotons_min; }
	/**
	 *  @return maximal number of scintillation photons per event
	 */
	G4int GetMaxValue_Control_ScintiPhotons() const
	{ return Control_ScintiPhotons_max; }

	void SetMinMaxValue_Control_WLSPhotons(G4int newValue);
	/**
	 *  @return minimal number of optical photons from a wavelength-shifting process per event
	 */
	G4int GetMinValue_Control_WLSPhotons() const
	{ return Control_WLSPhotons_min; }
	/**
	 *  @return maximal number of optical photons from a wavelength-shifting process per event
	 */
	G4int GetMaxValue_Control_WLSPhotons() const
	{ return Control_WLSPhotons_max; }

	void SetMinMaxValue_Control_CerenkovPhotonsFromPrimary(G4int newValue);
	/**
	 *  @return minimal number of Cerenkov photons that have been created by primary particles (per event)
	 */
	G4int GetMinValue_Control_CerenkovPhotonsFromPrimary() const
	{ return Control_CerenkovPhotonsFromPrimary_min; }
	/**
	 *  @return maximal number of Cerenkov photons that have been created by primary particles (per event)
	 */
	G4int GetMaxValue_Control_CerenkovPhotonsFromPrimary() const
	{ return Control_CerenkovPhotonsFromPrimary_max; }

	void SetMinMaxValue_Control_ScintiPhotonsFromPrimary(G4int newValue);
	/**
	 *  @return minimal number of scintillation photons that have been created by primary particles (per event)
	 */
	G4int GetMinValue_Control_ScintiPhotonsFromPrimary() const
	{ return Control_ScintiPhotonsFromPrimary_min; }
	/**
	 *  @return maximal number of scintillation photons that have been created by primary particles (per event)
	 */
	G4int GetMaxValue_Control_ScintiPhotonsFromPrimary() const
	{ return Control_ScintiPhotonsFromPrimary_max; }

	void SetMinMaxValue_Control_WLSPhotonsFromPrimary(G4int newValue);
	/**
	 *  @return minimal number of optical photons from a wavelength-shifting process that have been created by primary particles (per event)
	 */
	G4int GetMinValue_Control_WLSPhotonsFromPrimary() const
	{ return Control_WLSPhotonsFromPrimary_min; }
	/**
	 *  @return minimal number of optical photons from a wavelength-shifting process that have been created by primary particles (per event)
	 */
	G4int GetMaxValue_Control_WLSPhotonsFromPrimary() const
	{ return Control_WLSPhotonsFromPrimary_max; }

	void SetMinMaxValue_Control_CerenkovPhotonsFromSecondary(G4int newValue);
	/**
	 *  @return minimal number of Cerenkov photons that have been created by secondary particles (per event)
	 */
	G4int GetMinValue_Control_CerenkovPhotonsFromSecondary() const
	{ return Control_CerenkovPhotonsFromSecondary_min; }
	/**
	 *  @return maximal number of Cerenkov photons that have been created by secondary particles (per event)
	 */
	G4int GetMaxValue_Control_CerenkovPhotonsFromSecondary() const
 	{ return Control_CerenkovPhotonsFromSecondary_max; }

	void SetMinMaxValue_Control_ScintiPhotonsFromSecondary(G4int newValue);
	/**
	 *  @return minimal number of scintillation photons that have been created by secondary particles (per event)
	 */
	G4int GetMinValue_Control_ScintiPhotonsFromSecondary() const
	{ return Control_ScintiPhotonsFromSecondary_min; }
	/**
	 *  @return maximal number of scintillation photons that have been created by secondary particles (per event)
	 */
	G4int GetMaxValue_Control_ScintiPhotonsFromSecondary() const
	{ return Control_ScintiPhotonsFromSecondary_max; }

	void SetMinMaxValue_Control_WLSPhotonsFromSecondary(G4int newValue);
	/**
	 *  @return minimal number of optical photons from a wavelength-shifting process that have been created by secondary particles (per event)
	 */
	G4int GetMinValue_Control_WLSPhotonsFromSecondary() const
	{ return Control_WLSPhotonsFromSecondary_min; }
	/**
	 *  @return maximal number of optical photons from a wavelength-shifting process that have been created by secondary particles (per event)
	 */
	G4int GetMaxValue_Control_WLSPhotonsFromSecondary() const
	{ return Control_WLSPhotonsFromSecondary_max; }

	void SetMinMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition(G4double newValue);
	/**
	 *  @return minimal number of optical photons that have been created by primary particles (per primary particle's energy deposition)
	 */
	G4double GetMinValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition() const
	{ return Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min; }
	/**
	 *  @return maximal number of optical photons that have been created by primary particles (per primary particle's energy deposition)
	 */
	G4double GetMaxValue_Control_opticalPhotonsFromPrimaryPerPrimaryEnergyDeposition() const
	{ return Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max; }

	void SetMinMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition(G4double newValue);
	/**
	 *  @return minimal number of optical photons that have been created by secondary particles (per primary particle's energy deposition)
	 */
	G4double GetMinValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition() const
	{ return Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min; }
	/**
	 *  @return maximal number of optical photons that have been created by secondary particles (per primary particle's energy deposition)
	 */
	G4double GetMaxValue_Control_opticalPhotonsFromSecondaryPerPrimaryEnergyDeposition() const
	{ return Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max; }

	void SetMinMaxValue_Control_opticalPhotonsAbsorbed_volume( G4int num_absorbed_inScinti,		/**< number of optical photons absorbed in scintillator material */
								   G4int num_absorbed_inFibre,		/**< number of optical photons absorbed in a fibre */
								   G4int num_absorbed_inPhotonDetector,	/**< number of optical photons absorbed in a photon detector */
								   G4int num_total			/**< total number of optical photons */
								 );
	/**
	 *  @return minimal fraction of optical photons that have been absorbed per event (sorted by the volumes they were absorbed in)
	 */
	G4double GetMinValue_Control_opticalPhotonsAbsorbed_volume() const
	{ return Control_OpticalPhotonsAbsorbed_volume_min; }
	/**
	 *  @return maximal fraction of optical photons that have been absorbed per event (sorted by the volumes they were absorbed in)
	 */
	G4double GetMaxValue_Control_opticalPhotonsAbsorbed_volume() const
	{ return Control_OpticalPhotonsAbsorbed_volume_max; }

	void SetMinMaxValue_Control_opticalPhotonsAbsorbed_process( G4int num_Scinti_absorbed,	/**< number of absorbed scintillation photons */
								    G4int num_Scinti_total,	/**< total number of scintillation photons */
								    G4int num_Cerenkov_absorbed,	/**< number of absorbed Cerenkov photons */
								    G4int num_Cerenkov_total,	/**< total number of Cerenkov photons */
								    G4int num_WLS_absorbed,	/**< number of absorbed optical photons from a wavelength-shifting process */
								    G4int num_WLS_total		/**< total number of optical photons from a wavelength-shifting process */
								  );
	/**
	 *  @return minimal fraction of optical photons that have been absorbed per event (sorted by their production process)
	 */
	G4double GetMinValue_Control_opticalPhotonsAbsorbed_process() const
	{ return Control_OpticalPhotonsAbsorbed_process_min; }
	/**
	 *  @return maximal fraction of optical photons that have been absorbed per event (sorted by their production process)
	 */
	G4double GetMaxValue_Control_opticalPhotonsAbsorbed_process() const
	{ return Control_OpticalPhotonsAbsorbed_process_max; }

private:
	G4int NumberOfEvents;
	G4int NumberOfPrimaryParticles;
	std::map<G4String, G4double> MassMap;

// for data plots
	G4double ScintiHitPointX_mm_min;
	G4double ScintiHitPointX_mm_max;
	G4double ScintiHitPointY_mm_min;
	G4double ScintiHitPointY_mm_max;
	G4double ScintiHitPointZ_mm_min;
	G4double ScintiHitPointZ_mm_max;
	G4int OpticalPhotonsPerPrimary_min;
	G4int OpticalPhotonsPerPrimary_max;
	G4double OpticalPhotonsPerEnergyDeposition_min;
	G4double OpticalPhotonsPerEnergyDeposition_max;
	G4int OpticalPhotonsAbsorbedInPhotonDetector_min;
	G4int OpticalPhotonsAbsorbedInPhotonDetector_max;
	G4double PhotonDetectorHitPointX_mm_min;
	G4double PhotonDetectorHitPointX_mm_max;
	G4double PhotonDetectorHitPointY_mm_min;
	G4double PhotonDetectorHitPointY_mm_max;
	G4double PhotonDetectorHitPointZ_mm_min;
	G4double PhotonDetectorHitPointZ_mm_max;
	G4double OpticalPhotonHitEnergy_MeV_min;
	G4double OpticalPhotonHitEnergy_MeV_max;
	G4double DeltaTimeHitAbsorption_ns_min;
	G4double DeltaTimeHitAbsorption_ns_max;

// for control plots
	G4int Control_PrimaryParticleID_min;
	G4int Control_PrimaryParticleID_max;
	G4int Control_SecondaryParticleID_min;
	G4int Control_SecondaryParticleID_max;
	G4double Control_PrimaryParticleInitialPositionX_mm_min;
	G4double Control_PrimaryParticleInitialPositionX_mm_max;
	G4double Control_PrimaryParticleInitialPositionY_mm_min;
	G4double Control_PrimaryParticleInitialPositionY_mm_max;
	G4double Control_PrimaryParticleInitialPositionZ_mm_min;
	G4double Control_PrimaryParticleInitialPositionZ_mm_max;
	G4double Control_ScintiHitPointX_mm_min;
	G4double Control_ScintiHitPointX_mm_max;
	G4double Control_ScintiHitPointY_mm_min;
	G4double Control_ScintiHitPointY_mm_max;
	G4double Control_ScintiHitPointZ_mm_min;
	G4double Control_ScintiHitPointZ_mm_max;
	G4double Control_PrimaryParticleInitialEnergy_MeV_min;
	G4double Control_PrimaryParticleInitialEnergy_MeV_max;
	G4double Control_SecondaryParticleInitialEnergy_MeV_min;
	G4double Control_SecondaryParticleInitialEnergy_MeV_max;
	G4double Control_OpticalPhotonInitialEnergy_MeV_min;
	G4double Control_OpticalPhotonInitialEnergy_MeV_max;
	G4double Control_PrimaryParticleInitialBetaGamma_min;
	G4double Control_PrimaryParticleInitialBetaGamma_max;
	G4double Control_SecondaryParticleInitialBetaGamma_min;
	G4double Control_SecondaryParticleInitialBetaGamma_max;
	G4double Control_PrimaryParticleEnergyDepositionInScintillator_MeV_min;
	G4double Control_PrimaryParticleEnergyDepositionInScintillator_MeV_max;
	G4double Control_PrimaryParticlePathLengthInScintillator_MeV_min;
	G4double Control_PrimaryParticlePathLengthInScintillator_MeV_max;
	G4double Control_SecondaryParticleEnergyDepositionInScintillator_MeV_min;
	G4double Control_SecondaryParticleEnergyDepositionInScintillator_MeV_max;
	G4double Control_GlobalScintiHitTime_ns_min;
	G4double Control_GlobalScintiHitTime_ns_max;
	G4double Control_DeltaTimeHitCreation_secondary_ns_min;
	G4double Control_DeltaTimeHitCreation_secondary_ns_max;
	G4double Control_DeltaTimeHitCreation_ns_min;
	G4double Control_DeltaTimeHitCreation_ns_max;
	G4double Control_DeltaTimeCreationAbsorption_ns_min;
	G4double Control_DeltaTimeCreationAbsorption_ns_max;
	G4double Control_DeltaTimeHitAbsorption_ns_min;
	G4double Control_DeltaTimeHitAbsorption_ns_max;
	G4double Control_DeltaTimeHitCreation_parentPrimary_ns_min;
	G4double Control_DeltaTimeHitCreation_parentPrimary_ns_max;
	G4double Control_DeltaTimeCreationAbsorption_parentPrimary_ns_min;
	G4double Control_DeltaTimeCreationAbsorption_parentPrimary_ns_max;
	G4double Control_DeltaTimeHitAbsorption_parentPrimary_ns_min;
	G4double Control_DeltaTimeHitAbsorption_parentPrimary_ns_max;
	G4double Control_DeltaTimeHitCreation_parentSecondary_ns_min;
	G4double Control_DeltaTimeHitCreation_parentSecondary_ns_max;
	G4double Control_DeltaTimeCreationAbsorption_parentSecondary_ns_min;
	G4double Control_DeltaTimeCreationAbsorption_parentSecondary_ns_max;
	G4double Control_DeltaTimeHitAbsorption_parentSecondary_ns_min;
	G4double Control_DeltaTimeHitAbsorption_parentSecondary_ns_max;
	G4int Control_OpticalPhotons_min;
	G4int Control_OpticalPhotons_max;
	G4int Control_OpticalPhotonsFromPrimary_min;
	G4int Control_OpticalPhotonsFromPrimary_max;
	G4int Control_OpticalPhotonsFromSecondary_min;
	G4int Control_OpticalPhotonsFromSecondary_max;
	G4int Control_CerenkovPhotons_min;
	G4int Control_CerenkovPhotons_max;
	G4int Control_ScintiPhotons_min;
	G4int Control_ScintiPhotons_max;
	G4int Control_WLSPhotons_min;
	G4int Control_WLSPhotons_max;
	G4int Control_CerenkovPhotonsFromPrimary_min;
	G4int Control_CerenkovPhotonsFromPrimary_max;
	G4int Control_ScintiPhotonsFromPrimary_min;
	G4int Control_ScintiPhotonsFromPrimary_max;
	G4int Control_WLSPhotonsFromPrimary_min;
	G4int Control_WLSPhotonsFromPrimary_max;
	G4int Control_CerenkovPhotonsFromSecondary_min;
	G4int Control_CerenkovPhotonsFromSecondary_max;
	G4int Control_ScintiPhotonsFromSecondary_min;
	G4int Control_ScintiPhotonsFromSecondary_max;
	G4int Control_WLSPhotonsFromSecondary_min;
	G4int Control_WLSPhotonsFromSecondary_max;
	G4double Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_min;
	G4double Control_OpticalPhotonsFromPrimaryPerPrimaryEnergyDeposition_max;
	G4double Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_min;
	G4double Control_OpticalPhotonsFromSecondaryPerPrimaryEnergyDeposition_max;
	G4double Control_OpticalPhotonsAbsorbed_volume_min;
	G4double Control_OpticalPhotonsAbsorbed_volume_max;
	G4double Control_OpticalPhotonsAbsorbed_process_min;
	G4double Control_OpticalPhotonsAbsorbed_process_max;
};

#endif
