/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef UserEventInformation_h
#define UserEventInformation_h 1


#include <G4VUserEventInformation.hh>
#include <G4ThreeVector.hh>
#include <globals.hh>

#include <map>
#include <vector>

#include <boost/any.hpp>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Class to save the data which belonging to the current event, using functions inherited from G4VUserEventInformation.
class UserEventInformation : public G4VUserEventInformation
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	UserEventInformation()
	: NumberOfPrimaryParticles(0)
	, ScintiPhoton(0)
	, ScintiPhotonByPrimary(0)
	, ScintiPhotonBySecondary(0)
	, CerenkovPhoton(0)
	, CerenkovPhotonByPrimary(0)
	, CerenkovPhotonBySecondary(0)
	, WLSCount(0)
	, AbsorbedInFibre(0)
	, AbsorbedScintiPhotonInFibre(0)
	, AbsorbedCerenkovPhotonInFibre(0)
	, AbsorbedWLSPhotonInFibre(0)
	, AbsorbedInPhotonDetector(0)
	, AbsorbedScintiPhotonInPhotonDetector(0)
	, AbsorbedCerenkovPhotonInPhotonDetector(0)
	, AbsorbedWLSPhotonInPhotonDetector(0)
	, AbsorbedInScinti(0)
	, AbsorbedScintiPhotonInScinti(0)
	, AbsorbedCerenkovPhotonInScinti(0)
	, AbsorbedWLSPhotonInScinti(0)
	, ParticleIDKey("particleID")
	, ParticleNameKey("particleName")
	, IsPrimaryKey("isPrimary")
	, ProductionMechanismKey("productionMechanism")
	, InitialPositionKey("initialPosition")
	, GlobalCreationTimeKey("globalCreationTime")
	, GlobalAbsorptionTimeKey("globalAbsorptionTime")
	, InitialMomentumKey("InitialMomentum")
	, ParentTrackIDKey("parentTrackID")
	, ParentNameKey("parentName")
	, ParentIDKey("parentID")
	, AbsorbedInKey("absorbedIn")
	{
	}

	/**
	 *  Destructor (empty)
	 */
	~UserEventInformation()
	{
	}



	inline void Print()const{};

	/**
	 *  Function to set the number of primary particles.
	 */
	void SetNumberOfPrimaryParticles(G4int npp)
	{ NumberOfPrimaryParticles = npp; }
	/**
	 *  @return number of primary particles
	 */
	G4int GetNumberOfPrimaryParticles() const
	{ return NumberOfPrimaryParticles; }

	/**
	 *  Function to increment the number of scintillation photons.
	 */
	void IncScintiPhotonCount()
	{ ScintiPhoton++; }
	/**
	 *  @return number of scintillation photons
	 */
	G4int GetScintiPhotonCount() const
	{ return ScintiPhoton; }

	/**
	 *  Function to increment the number of scintillation photons that where created by a primary particle.
	 */
	void IncScintiPhotonByPrimaryCount()
	{ ScintiPhotonByPrimary++; }
	/**
	 *  @return number of scintillation photons that where created by a primary particle
	 */
	G4int GetScintiPhotonByPrimaryCount() const
	{ return ScintiPhotonByPrimary; }

	/**
	 *  Function to increment the number of scintillation photons that where created by a secondary particle.
	 */
	void IncScintiPhotonBySecondaryCount()
	{ ScintiPhotonBySecondary++; }
	/**
	 *  @return number of scintillation photons that where created by a secondary particle
	 */
	G4int GetScintiPhotonBySecondaryCount() const
	{ return ScintiPhotonBySecondary; }

	/**
	 *  Function to increment the number of Cerenkov photons.
	 */
	void IncCerenkovPhotonCount()
	{ CerenkovPhoton++; }
	/**
	 *  @return number of Cerenkov photons
	 */
	G4int GetCerenkovPhotonCount() const
	{ return CerenkovPhoton; }

	/**
	 *  Function to increment the number of Cerenkov photons that where created by a primary particle.
	 */
	void IncCerenkovPhotonByPrimaryCount()
	{ CerenkovPhotonByPrimary++; }
	/**
	 *  @return number of Cerenkov photons that where created by a primary particle
	 */
	G4int GetCerenkovPhotonByPrimaryCount() const
	{ return CerenkovPhotonByPrimary; }

	/**
	 *  Function to increment the number of Cerenkov photons that where created by a secondary particle.
	 */
	void IncCerenkovPhotonBySecondaryCount()
	{ CerenkovPhotonBySecondary++; }
	/**
	 *  @return number of Cerenkov photons that where created by a secondary particle
	 */
	G4int GetCerenkovPhotonBySecondaryCount() const
	{ return CerenkovPhotonBySecondary; }

	/**
	 *  Function to increment the number of optical photons from a wavelength-shifting process.
	 */
	void IncWLSPhotonCount()
	{ WLSCount++; }
	/**
	 *  @return number of optical photons from a wavelength-shifting process.
	 */
	G4int GetWLSPhotonCount() const
	{ return WLSCount; }

	/**
	 *  @return number of optical photons that have been absorbed
	 */
	G4int GetAbsorbedCount() const
	{ return (AbsorbedInScinti + AbsorbedInFibre + AbsorbedInPhotonDetector); }

	/**
	 *  @return number of scintillation photons that have been absorbed
	 */
	G4int GetAbsorbedScintiPhotonCount() const
	{ return (AbsorbedScintiPhotonInScinti + AbsorbedScintiPhotonInFibre + AbsorbedScintiPhotonInPhotonDetector); }

	/**
	 *  @return number of Cerenkov photons that have been absorbed
	 */
	G4int GetAbsorbedCerenkovPhotonCount() const
	{ return (AbsorbedCerenkovPhotonInScinti + AbsorbedCerenkovPhotonInFibre + AbsorbedCerenkovPhotonInPhotonDetector); }

	/**
	 *  @return number of optical photons from a wavelength-shifting process that have been absorbed
	 */
	G4int GetAbsorbedWLSPhotonCount() const
	{ return (AbsorbedWLSPhotonInScinti + AbsorbedWLSPhotonInFibre + AbsorbedWLSPhotonInPhotonDetector); }

	/**
	 *  Function to save the name of the volume, the particle with the corresponding track ID was absorbed in.
	 */
	void SetAbsorbedIn(G4int trackID, G4String absorbedIn)
	{ (ParticleVector[trackID])[AbsorbedInKey] = absorbedIn; }
	G4String GetAbsorbedIn(G4int trackID) const;

	/**
	 *  Function to increment the number of optical photons that have been absorbed in a fibre.
	 */
	void IncAbsorbedInFibreCount()
	{ AbsorbedInFibre++; }
	/**
	 *  @return number of optical photons that have been absorbed in a fibre
	 */
	G4int GetAbsorbedInFibreCount() const
	{ return AbsorbedInFibre; }

	/**
	 *  Function to increment the number of scintillation photons that have been absorbed in a fibre.
	 */
	void IncAbsorbedScintiPhotonInFibreCount()
	{ AbsorbedScintiPhotonInFibre++; }
	/**
	 *  @return number of scintillation photons that have been absorbed in a fibre
	 */
	G4int GetAbsorbedScintiPhotonInFibreCount() const
	{ return AbsorbedScintiPhotonInFibre; }

	/**
	 *  Function to increment the number of Cerenkov photons that have been absorbed in a fibre.
	 */
	void IncAbsorbedCerenkovPhotonInFibreCount()
	{ AbsorbedCerenkovPhotonInFibre++; }
	/**
	 *  @return number of Cerenkov photons that have been absorbed in a fibre
	 */
	G4int GetAbsorbedCerenkovPhotonInFibreCount() const
	{ return AbsorbedCerenkovPhotonInFibre; }

	/**
	 *  Function to increment the number of optical photons from a wavelength-shifting process that have been absorbed in a fibre.
	 */
	void IncAbsorbedWLSPhotonInFibreCount()
	{ AbsorbedWLSPhotonInFibre++; }
	/**
	 *  @return number of optical photons from a wavelength-shifting process that have been absorbed in a fibre
	 */
	G4int GetAbsorbedWLSPhotonInFibreCount() const
	{ return AbsorbedWLSPhotonInFibre; }

	/**
	 *  Function to increment the number of optical photons that have been absorbed in a photon detector.
	 */
	void IncAbsorbedInPhotonDetectorCount()
	{ AbsorbedInPhotonDetector++; }
	/**
	 *  @return number of optical photons that have been absorbed in a photon detector
	 */
	G4int GetAbsorbedInPhotonDetectorCount() const
	{ return AbsorbedInPhotonDetector; }

	/**
	 *  Function to increment the number of scintillation photons that have been absorbed in a photon detector.
	 */
	void IncAbsorbedScintiPhotonInPhotonDetectorCount()
	{ AbsorbedScintiPhotonInPhotonDetector++; }
	/**
	 *  @return number of scintillation photons that have been absorbed in a photon detector
	 */
	G4int GetAbsorbedScintiPhotonInPhotonDetectorCount() const
	{ return AbsorbedScintiPhotonInPhotonDetector; }

	/**
	 *  Function to increment the number of Cerenkov photons that have been absorbed in a photon detector.
	 */
	void IncAbsorbedCerenkovPhotonInPhotonDetectorCount()
	{ AbsorbedCerenkovPhotonInPhotonDetector++; }
	/**
	 *  @return number of Cerenkov photons that have been absorbed in a photon detector
	 */
	G4int GetAbsorbedCerenkovPhotonInPhotonDetectorCount() const
	{ return AbsorbedCerenkovPhotonInPhotonDetector; }

	/**
	 *  Function to increment the number of optical photons from a wavelength-shifting process that have been absorbed in a photon detector.
	 */
	void IncAbsorbedWLSPhotonInPhotonDetectorCount()
	{ AbsorbedWLSPhotonInPhotonDetector++; }
	/**
	 *  @return number of optical photons from a wavelength-shifting process that have been absorbed in a photon detector
	 */
	G4int GetAbsorbedWLSPhotonInPhotonDetectorCount() const
	{ return AbsorbedWLSPhotonInPhotonDetector; }

	/**
	 *  Function to increment the number of optical photons that have been absorbed in a scintillator tile.
	 */
	void IncAbsorbedInScintiCount()
	{ AbsorbedInScinti++; }
	/**
	 *  @return number of optical photons that have been absorbed in a scintillator tile
	 */
	G4int GetAbsorbedInScintiCount() const
	{ return AbsorbedInScinti; }

	/**
	 *  Function to increment the number of scintillation photons that have been absorbed in a scintillator tile.
	 */
	void IncAbsorbedScintiPhotonInScintiCount()
	{ AbsorbedScintiPhotonInScinti++; }
	/**
	 *  @return number of scintillation photons that have been absorbed in a scintillator tile
	 */
	G4int GetAbsorbedScintiPhotonInScintiCount() const
	{ return AbsorbedScintiPhotonInScinti; }

	/**
	 *  Function to increment the number of Cerenkov photons that have been absorbed in a scintillator tile.
	 */
	void IncAbsorbedCerenkovPhotonInScintiCount()
	{ AbsorbedCerenkovPhotonInScinti++; }
	/**
	 *  @return number of Cerenkov photons that have been absorbed in a scintillator tile
	 */
	G4int GetAbsorbedCerenkovPhotonInScintiCount() const
	{ return AbsorbedCerenkovPhotonInScinti; }

	/**
	 *  Function to increment the number of optical photons from a wavelength-shifting process that have been absorbed in a scintillator tile.
	 */
	void IncAbsorbedWLSPhotonInScintiCount()
	{ AbsorbedWLSPhotonInScinti++; }
	/**
	 *  @return number of optical photons from a wavelength-shifting process that have been absorbed in a scintillator tile
	 */
	G4int GetAbsorbedWLSPhotonInScintiCount() const
	{ return AbsorbedWLSPhotonInScinti; }

	/**
	 *  @return number of particles that have been saved
	 */
	G4int GetNumberOfParticles() const
	{ return ParticleVector.size() - 1; }   // ParticleVector[0] is always empty as the trackIDs start with 1
	G4bool TrackIDExistsInVector(G4int trackID) const;
// 	G4bool TrackIDHasEntries(G4int trackID) const;
	void FillEmptyEntriesUpToTrackID(G4int trackID);

	/**
	 *  Function to save particle ID of the particle with the corresponding track ID.
	 */
	void SetParticleID(G4int trackID, G4int partID)
	{ (ParticleVector[trackID])[ParticleIDKey] = partID; }
	G4int GetParticleID(G4int trackID) const;

	/**
	 *  Function to save name of the particle with the corresponding track ID.
	 */
	void SetParticleName(G4int trackID, G4String partName)
	{ (ParticleVector[trackID])[ParticleNameKey] = partName; }
	G4String GetParticleName(G4int trackID) const;

	/**
	 *  Function to save whether the particle with the corresponding track ID is a primary particle.
	 */
	void SetParticleIsPrimary(G4int trackID, G4bool partIsPrimary)
	{ (ParticleVector[trackID])[IsPrimaryKey] = partIsPrimary; }
	G4bool GetParticleIsPrimary(G4int trackID) const;

	/**
	 *  Function to save production mechanism of the particle with the corresponding track ID.
	 */
	void SetProductionMechanism(G4int trackID, G4String prodMech)
	{ (ParticleVector[trackID])[ProductionMechanismKey] = prodMech; }
	G4String GetProductionMechanism(G4int trackID) const;

	/**
	 *  Function to save initial position of the particle with the corresponding track ID.
	 */
	void SetInitialPosition(G4int trackID, G4ThreeVector initPos)
	{ (ParticleVector[trackID])[InitialPositionKey] = initPos; }
	G4ThreeVector GetInitialPosition(G4int trackID) const;

	/**
	 *  Function to save creation time of the particle with the corresponding track ID.
	 */
	void SetGlobalCreationTime(G4int trackID, G4double creationTime)
	{ (ParticleVector[trackID])[GlobalCreationTimeKey] = creationTime; }
	G4double GetGlobalCreationTime(G4int trackID) const;

	/**
	 *  Function to save initial momentum of the particle with the corresponding track ID.
	 */
	void SetInitialMomentum(G4int trackID, G4ThreeVector initMomentum)
	{ (ParticleVector[trackID])[InitialMomentumKey] = initMomentum; }
	G4ThreeVector GetInitialMomentum(G4int trackID) const;

	void SetParent(G4int trackID, G4int parentTrackID);
	G4int GetParentTrackID(G4int trackID) const;
	G4String GetParentName(G4int trackID) const;
	G4int GetParentID(G4int trackID) const;
	G4bool GetParentIsPrimary(G4int trackID) const;

	/**
	 *  Function to save absorption time of the particle with the corresponding track ID.
	 */
	void SetGlobalAbsorptionTime(G4int trackID, G4double absorptionTime)
	{ (ParticleVector[trackID])[GlobalAbsorptionTimeKey] = absorptionTime; }
	G4double GetGlobalAbsorptionTime(G4int trackID) const;


	/**
	 *  Function to removed all saved data.
	 */
	void clean()
	{ ParticleVector.clear(); }



private:
	G4int NumberOfPrimaryParticles;
	G4int ScintiPhoton;
	G4int ScintiPhotonByPrimary;
	G4int ScintiPhotonBySecondary;
	G4int CerenkovPhoton;
	G4int CerenkovPhotonByPrimary;
	G4int CerenkovPhotonBySecondary;
	G4int WLSCount;
	G4int AbsorbedInFibre;
	G4int AbsorbedScintiPhotonInFibre;
	G4int AbsorbedCerenkovPhotonInFibre;
	G4int AbsorbedWLSPhotonInFibre;
	G4int AbsorbedInPhotonDetector;
	G4int AbsorbedScintiPhotonInPhotonDetector;
	G4int AbsorbedCerenkovPhotonInPhotonDetector;
	G4int AbsorbedWLSPhotonInPhotonDetector;
	G4int AbsorbedInScinti;
	G4int AbsorbedScintiPhotonInScinti;
	G4int AbsorbedCerenkovPhotonInScinti;
	G4int AbsorbedWLSPhotonInScinti;
	std::vector< std::map<G4String, boost::any> > ParticleVector;
	G4String ParticleIDKey;
	G4String ParticleNameKey;
	G4String IsPrimaryKey;
	G4String ProductionMechanismKey;
	G4String InitialPositionKey;
	G4String GlobalCreationTimeKey;
	G4String GlobalAbsorptionTimeKey;
	G4String InitialMomentumKey;
	G4String ParentTrackIDKey;
	G4String ParentNameKey;
	G4String ParentIDKey;
	G4String AbsorbedInKey;
};

#endif
