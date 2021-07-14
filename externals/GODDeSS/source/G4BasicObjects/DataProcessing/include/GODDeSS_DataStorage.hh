/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef GODDeSS_DataStorage_h
#define GODDeSS_DataStorage_h 1


#include <G4ThreeVector.hh>
#include <G4String.hh>
#include <G4Step.hh>

#include <map>
#include <vector>

#include <boost/any.hpp>



// class variables begin with capital letters, local variables with small letters



///  a class to save the data which was taken by the sensitive detectors
class GODDeSS_DataStorage
{
public:
	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	GODDeSS_DataStorage()
	: PhotonDetectorHitFile("")
	, ScintillatorHitPositionKey("ScintillatorHitPosition")
	, ScintillatorHitMomentumKey("ScintillatorHitMomentum")
	, ScintillatorHitTimeKey("ScintillatorHitTime")
	, EnergyDepositionKey("deltaE_Scinti")
	, PathLengthKey("deltaX_Scinti")
	, WrappingHitPositionKey("WrappingHitPosition")
	, WrappingHitMomentumKey("WrappingHitMomentum")
	, WrappingHitTimeKey("WrappingHitTime")
	, FibreWasHitKey("FibreWasHit")
	, PhotonDetectorWasHitKey("PhotonDetectorWasHit")
	, PhotonDetectorHitPositionKey("PhotonDetectorHitPosition")
	, PhotonDetectorHitMomentumKey("PhotonDetectorHitMomentum")
	, PhotonDetectorHitNameKey("PhotonDetectorHitName")
	, ParticleIDKey("particleID")
	, ParticleNameKey("particleName")
	, IsPrimaryKey("isPrimary")
	, ProductionMechanismKey("productionMechanism")
	{
		SetDefaults();
	}

	/**
	 *  Destructor (empty)
	 */
	~GODDeSS_DataStorage()
	{
	}



	/**
	 *  Function to set the path of the text file, which the data needed to resimulate a photon that hit the sensitive detector is written into.
	 */
	void SetPhotonDetectorHitFile(G4String photonDetectorHitFile)
	{
		PhotonDetectorHitFile = photonDetectorHitFile;
	}
	/**
	 *  @return path of the text file, which the data needed to resimulate a photon that hit the sensitive detector is written into
	 */
	G4String GetPhotonDetectorHitFile() const
	{
		return PhotonDetectorHitFile;
	}


	void SetScintillatorHitPoint(G4int trackID, G4ThreeVector hitPoint);
	G4ThreeVector GetScintillatorHitPoint(G4int trackID) const;

	void SetScintillatorHitMomentum(G4int trackID, G4ThreeVector hitMomentum);
	G4ThreeVector GetScintillatorHitMomentum(G4int trackID) const;

	void SetScintillatorHitTime(G4int trackID, G4double hitTime);
	G4double GetScintillatorHitTime(G4int trackID) const;

	void AddUpEnergyDepositionInScintillator(G4int trackID, G4double deltaE);
	G4double GetEnergyDepositionInScintillator(G4int trackID) const;

	void AddUpPathLengthInScintillator(G4int trackID, G4double deltaX);
	G4double GetPathLengthInScintillator(G4int trackID) const;


	void SetWrappingHitPoint(G4int trackID, G4ThreeVector hitPoint);
	G4ThreeVector GetWrappingHitPoint(G4int trackID) const;

	void SetWrappingHitMomentum(G4int trackID, G4ThreeVector hitMomentum);
	G4ThreeVector GetWrappingHitMomentum(G4int trackID) const;

	void SetWrappingHitTime(G4int trackID, G4double hitTime);
	G4double GetWrappingHitTime(G4int trackID) const;



	void SetFibreWasHit(G4int trackID);
	G4bool FibreWasHit(G4int trackID) const;



	/**
	 *  Function to save the data needed to resimulate a photon that hit the sensitive detector. It has to be processed for every step (i.e. in the SteppingAction).
	 */
	void SavePreviousStepOpticalPhotonData(const G4Step * step)
	{
		G4StepPoint * preStepPoint = step->GetPreStepPoint();
		PreviousPreStepPointMomentum = preStepPoint->GetMomentum();
		PreviousPreStepPointPosition = preStepPoint->GetPosition();
		PreviousPreStepPointTime = preStepPoint->GetGlobalTime();
		PreviousPreStepPointPolarisation = preStepPoint->GetPolarization();

		G4double stepLength = step->GetStepLength();
		if(stepLength > 0.5 * CLHEP::mm)
		{
			PreviousPreStepPointPosition += step->GetDeltaPosition() * (1. - 0.5 * CLHEP::mm / stepLength);
			PreviousPreStepPointTime += step->GetDeltaTime() * (1. - 0.5 * CLHEP::mm / stepLength);
		}
	}
	/**
	 *  @return start time (for resimulation purpose) of a photon hitting the sensitive detector
	 */
	G4double GetPreviousPreStepPointTime() const
	{
		return PreviousPreStepPointTime;
	}
	/**
	 *  @return start position (for resimulation purpose) of a photon hitting the sensitive detector
	 */
	G4ThreeVector GetPreviousPreStepPointPosition() const
	{
		return PreviousPreStepPointPosition;
	}
	/**
	 *  @return start momentum (for resimulation purpose) of a photon hitting the sensitive detector
	 */
	G4ThreeVector GetPreviousPreStepPointMomentum() const
	{
		return PreviousPreStepPointMomentum;
	}
	/**
	 *  @return polarisation (for resimulation purpose) of a photon hitting the sensitive detector.
	 */
	G4ThreeVector GetPreviousPreStepPointPolarisation() const
	{
		return PreviousPreStepPointPolarisation;
	}

	/**
	*  Function to save the names of the existing sensitive detectors volumes of photon detectors.
	*/
	void AddPhotonSensitiveDetectorVolumeName(G4String photonSensitiveDetectorVolumeName)
	{
		PhotonSensitiveDetectorNameVector.push_back(photonSensitiveDetectorVolumeName);
	}
	/**
	*  @return vector with the names of the existing sensitive detectors volumes of photon detectors
	*/
	std::vector<G4String> GetPhotonSensitiveDetectorVolumeNames() const
	{
		return PhotonSensitiveDetectorNameVector;
	}

	void SetPhotonDetectorWasHit(G4int trackID);
	G4bool PhotonDetectorWasHit(G4int trackID) const;

	void SetPhotonDetectorHitPoint(G4int trackID, G4ThreeVector hitPoint);
	G4ThreeVector GetPhotonDetectorHitPoint(G4int trackID) const;

	void SetPhotonDetectorHitMomentum(G4int trackID, G4ThreeVector hitMomentum);
	G4ThreeVector GetPhotonDetectorHitMomentum(G4int trackID) const;

	/**
	 *  Function to save the name of the sensitive detector's volume that was hit by the particle with the corresponding track ID.
	 */
	void SetNameOfPhotonDetectorThatWasHit(G4int trackID, G4String hitVolume)
	{
		(ParticleVector[trackID])[PhotonDetectorHitNameKey] = hitVolume;
	}
	G4String GetNameOfPhotonDetectorThatWasHit(G4int trackID) const;



	/**
	 *  Function to save particle ID of the particle with the corresponding track ID.
	 */
	void SetParticleID(G4int trackID, G4int partID)
	{
		(ParticleVector[trackID])[ParticleIDKey] = partID;
	}
	G4int GetParticleID(G4int trackID) const;

	/**
	 *  Function to save name of the particle with the corresponding track ID.
	 */
	void SetParticleName(G4int trackID, G4String partName)
	{
		(ParticleVector[trackID])[ParticleNameKey] = partName;
	}
	G4String GetParticleName(G4int trackID) const;

	/**
	 *  Function to save whether the particle with the corresponding track ID is a primary particle.
	 */
	void SetParticleIsPrimary(G4int trackID)
	{
		(ParticleVector[trackID])[IsPrimaryKey] = true;
	}
	G4bool GetParticleIsPrimary(G4int trackID) const;

	/**
	 *  Function to save production mechanism of the particle with the corresponding track ID.
	 */
	void SetProductionMechanism(G4int trackID, G4String prodMech)
	{
		(ParticleVector[trackID])[ProductionMechanismKey] = prodMech;
	}
	G4String GetProductionMechanism(G4int trackID) const;



	G4bool TrackIDExistsInVector(G4int trackID) const;
	G4bool TrackIDHasEntries(G4int trackID) const;
	void FillEmptyEntriesUpToCurrentTrack(G4Track * theTrack);



	void SetDefaults();
	void clean() { ParticleVector.clear(); }



private:
	G4String PhotonDetectorHitFile;

	G4ThreeVector PreviousPreStepPointMomentum;
	G4ThreeVector PreviousPreStepPointPosition;
	G4double PreviousPreStepPointTime;
	G4ThreeVector PreviousPreStepPointPolarisation;

	std::vector< std::map<G4String, boost::any> > ParticleVector;
	G4String ScintillatorHitPositionKey;
	G4String ScintillatorHitMomentumKey;
	G4String ScintillatorHitTimeKey;
	G4String EnergyDepositionKey;
	G4String PathLengthKey;
	G4String WrappingHitPositionKey;
	G4String WrappingHitMomentumKey;
	G4String WrappingHitTimeKey;
	G4String FibreWasHitKey;
	std::vector<G4String> PhotonSensitiveDetectorNameVector;
	G4String PhotonDetectorWasHitKey;
	G4String PhotonDetectorHitPositionKey;
	G4String PhotonDetectorHitMomentumKey;
	G4String PhotonDetectorHitNameKey;
	G4String ParticleIDKey;
	G4String ParticleNameKey;
	G4String IsPrimaryKey;
	G4String ProductionMechanismKey;
};

#endif
