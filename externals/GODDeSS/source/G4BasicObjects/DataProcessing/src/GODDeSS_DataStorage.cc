/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4ParticleTable.hh>
#include <G4EventManager.hh>   //somehow needed for theTrack->GetCreatorProcess()->GetProcessName()

#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to save the position of the particle with the corresponding track ID when first hitting a scintillator tile.
 */
void GODDeSS_DataStorage::SetScintillatorHitPoint(G4int trackID, G4ThreeVector hitPoint)
{
	if( ! TrackIDHasEntries(trackID) )						(ParticleVector[trackID])[ScintillatorHitPositionKey] = hitPoint;
	else if( ! (ParticleVector[trackID]).count(ScintillatorHitPositionKey) )	(ParticleVector[trackID])[ScintillatorHitPositionKey] = hitPoint;
}
/**
 *  @return <b> if the entry exists for this particle: </b> position of the particle with the corresponding track ID when first hitting a scintillator tile
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector GODDeSS_DataStorage::GetScintillatorHitPoint(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(ScintillatorHitPositionKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(ScintillatorHitPositionKey) );
	}

	return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  Function to save the momentum of the particle with the corresponding track ID when first hitting a scintillator tile.
 */
void GODDeSS_DataStorage::SetScintillatorHitMomentum(G4int trackID, G4ThreeVector hitMomentum)
{
	if( ! TrackIDHasEntries(trackID) )						(ParticleVector[trackID])[ScintillatorHitMomentumKey] = hitMomentum;
	else if( ! (ParticleVector[trackID]).count(ScintillatorHitMomentumKey) )	(ParticleVector[trackID])[ScintillatorHitMomentumKey] = hitMomentum;
}
/**
 *  @return <b> if the entry exists for this particle: </b> momentum of the particle with the corresponding track ID when first hitting a scintillator tile
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector GODDeSS_DataStorage::GetScintillatorHitMomentum(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(ScintillatorHitMomentumKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(ScintillatorHitMomentumKey) );
	}

	return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  Function to save the time when the particle with the corresponding track ID first hits a scintillator tile.
 */
void GODDeSS_DataStorage::SetScintillatorHitTime(G4int trackID, G4double hitTime)
{
	if( ! TrackIDHasEntries(trackID) )					(ParticleVector[trackID])[ScintillatorHitTimeKey] = hitTime;
	else if( ! (ParticleVector[trackID]).count(ScintillatorHitTimeKey) )	(ParticleVector[trackID])[ScintillatorHitTimeKey] = hitTime;
}
/**
 *  @return <b> if the entry exists for this particle: </b> time when the particle with the corresponding track ID first hits a scintillator tile
 *  @return <b> else: </b> NAN
 */
G4double GODDeSS_DataStorage::GetScintillatorHitTime(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(ScintillatorHitTimeKey) ) return boost::any_cast<G4double>( (ParticleVector[trackID]).at(ScintillatorHitTimeKey) );
	}

	return NAN;
}

/**
 *  Function to add up the energy deposition of the particle with the corresponding track ID in scintillator material.
 */
void GODDeSS_DataStorage::AddUpEnergyDepositionInScintillator(G4int trackID, G4double deltaE)
{
	if( ! TrackIDHasEntries(trackID) )					(ParticleVector[trackID])[EnergyDepositionKey] = deltaE;
	else if( ! (ParticleVector[trackID]).count(ScintillatorHitTimeKey) )	(ParticleVector[trackID])[EnergyDepositionKey] = deltaE;
	else if( isnan(GetEnergyDepositionInScintillator(trackID)) )		(ParticleVector[trackID])[EnergyDepositionKey] = deltaE;
	else
	{
		G4double deltaE_inMap = GetEnergyDepositionInScintillator(trackID);
		(ParticleVector[trackID])[EnergyDepositionKey] = deltaE_inMap + deltaE;
	}
}
/**
 *  @return <b> if the entry exists for this particle: </b> energy deposition of the particle with the corresponding track ID in scintillator material
 *  @return <b> else: </b> NAN
 */
G4double GODDeSS_DataStorage::GetEnergyDepositionInScintillator(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(EnergyDepositionKey) ) return boost::any_cast<G4double>( (ParticleVector[trackID]).at(EnergyDepositionKey) );
	}

	return NAN;
}

/**
 *  Function to add up the path length of the particle with the corresponding track ID in scintillator material.
 */
void GODDeSS_DataStorage::AddUpPathLengthInScintillator(G4int trackID, G4double deltaX)
{
	if( ! TrackIDHasEntries(trackID) )					(ParticleVector[trackID])[PathLengthKey] = deltaX;
	else if( ! (ParticleVector[trackID]).count(ScintillatorHitTimeKey) )	(ParticleVector[trackID])[PathLengthKey] = deltaX;
	else if( isnan(GetPathLengthInScintillator(trackID)) )		(ParticleVector[trackID])[PathLengthKey] = deltaX;
	else
	{
		G4double deltaX_inMap = GetPathLengthInScintillator(trackID);
		(ParticleVector[trackID])[PathLengthKey] = deltaX_inMap + deltaX;
	}
}
/**
 *  @return <b> if the entry exists for this particle: </b> path length of the particle with the corresponding track ID in scintillator material
 *  @return <b> else: </b> NAN
 */
G4double GODDeSS_DataStorage::GetPathLengthInScintillator(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(PathLengthKey) ) return boost::any_cast<G4double>( (ParticleVector[trackID]).at(PathLengthKey) );
	}

	return NAN;
}



/**
 *  Function to save the position of the particle with the corresponding track ID when first hitting a wrapping.
 */
void GODDeSS_DataStorage::SetWrappingHitPoint(G4int trackID, G4ThreeVector hitPoint)
{
	if( ! TrackIDHasEntries(trackID) )					(ParticleVector[trackID])[WrappingHitPositionKey] = hitPoint;
	else if( ! (ParticleVector[trackID]).count(WrappingHitPositionKey) )	(ParticleVector[trackID])[WrappingHitPositionKey] = hitPoint;
}
/**
 *  @return <b> if the entry exists for this particle: </b> position of the particle with the corresponding track ID when first hitting a wrapping
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector GODDeSS_DataStorage::GetWrappingHitPoint(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(WrappingHitPositionKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(WrappingHitPositionKey) );
	}

	return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  Function to save the momentum of the particle with the corresponding track ID when first hitting a wrapping.
 */
void GODDeSS_DataStorage::SetWrappingHitMomentum(G4int trackID, G4ThreeVector hitMomentum)
{
	if( ! TrackIDHasEntries(trackID) )					(ParticleVector[trackID])[WrappingHitMomentumKey] = hitMomentum;
	else if( ! (ParticleVector[trackID]).count(WrappingHitMomentumKey) )	(ParticleVector[trackID])[WrappingHitMomentumKey] = hitMomentum;
}
/**
 *  @return <b> if the entry exists for this particle: </b> momentum of the particle with the corresponding track ID when first hitting a wrapping
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector GODDeSS_DataStorage::GetWrappingHitMomentum(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(WrappingHitMomentumKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(WrappingHitMomentumKey) );
	}

	return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  Function to save the time when the particle with the corresponding track ID first hits a wrapping.
 */
void GODDeSS_DataStorage::SetWrappingHitTime(G4int trackID, G4double hitTime)
{
	if( ! TrackIDHasEntries(trackID) )					(ParticleVector[trackID])[WrappingHitTimeKey] = hitTime;
	else if( ! (ParticleVector[trackID]).count(WrappingHitTimeKey) )	(ParticleVector[trackID])[WrappingHitTimeKey] = hitTime;
}
/**
 *  @return <b> if the entry exists for this particle: </b> time when the particle with the corresponding track ID first hits a wrapping
 *  @return <b> else: </b> NAN
 */
G4double GODDeSS_DataStorage::GetWrappingHitTime(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(WrappingHitTimeKey) ) return boost::any_cast<G4double>( (ParticleVector[trackID]).at(WrappingHitTimeKey) );
	}

	return NAN;
}



/**
 *  Function to change the G4bool, if a fibre was hit by the particle, to "true".
 */
void GODDeSS_DataStorage::SetFibreWasHit(G4int trackID)
{
	(ParticleVector[trackID])[FibreWasHitKey] = true;
}
/**
 *  @return <b> if the entry exists for this particle: </b> G4bool if a fibre was hit by the particle with the corresponding track ID
 *  @return <b> else: </b> false
 */
G4bool GODDeSS_DataStorage::FibreWasHit(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(FibreWasHitKey) ) return boost::any_cast<G4bool>( (ParticleVector[trackID]).at(FibreWasHitKey) );
	}

	return false;
}



/**
 *  Function to change the G4bool, if a sensitive detector was hit by the particle, to "true".
 */
void GODDeSS_DataStorage::SetPhotonDetectorWasHit(G4int trackID)
{
	(ParticleVector[trackID])[PhotonDetectorWasHitKey] = true;
}
/**
 *  @return <b> if the entry exists for this particle: </b> G4bool if a photon detector that was hit by the particle with the corresponding track ID
 *  @return <b> else: </b> false
 */
G4bool GODDeSS_DataStorage::PhotonDetectorWasHit(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(PhotonDetectorWasHitKey) ) return boost::any_cast<G4bool>( (ParticleVector[trackID]).at(PhotonDetectorWasHitKey) );
	}

	return false;
}

/**
 *  Function to save the position on the photon detector, where the particle with the corresponding track ID hit it.
 */
void GODDeSS_DataStorage::SetPhotonDetectorHitPoint(G4int trackID, G4ThreeVector hitPoint)
{
	if( ! TrackIDHasEntries(trackID) )						(ParticleVector[trackID])[PhotonDetectorHitPositionKey] = hitPoint;
	else if( ! (ParticleVector[trackID]).count(PhotonDetectorHitPositionKey) )	(ParticleVector[trackID])[PhotonDetectorHitPositionKey] = hitPoint;
}
/**
 *  @return <b> if the entry exists for this particle: </b> position on the photon detector, where the particle with the corresponding track ID hit it
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector GODDeSS_DataStorage::GetPhotonDetectorHitPoint(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(PhotonDetectorHitPositionKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(PhotonDetectorHitPositionKey) );
	}

	return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  Function to save the momentum of the particle with the corresponding track ID, when hitting the photon detector.
 */
void GODDeSS_DataStorage::SetPhotonDetectorHitMomentum(G4int trackID, G4ThreeVector hitMomentum)
{
	if( ! TrackIDHasEntries(trackID) )						(ParticleVector[trackID])[PhotonDetectorHitMomentumKey] = hitMomentum;
	else if( ! (ParticleVector[trackID]).count(PhotonDetectorHitMomentumKey) )	(ParticleVector[trackID])[PhotonDetectorHitMomentumKey] = hitMomentum;
}
/**
 *  @return <b> if the entry exists for this particle: </b> momentum of the particle with the corresponding track ID, when hitting the photon detector
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector GODDeSS_DataStorage::GetPhotonDetectorHitMomentum(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(PhotonDetectorHitMomentumKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(PhotonDetectorHitMomentumKey) );
	}

	return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  @return <b> if the entry exists for this particle: </b> name of the sensitive detector's volume that was hit by the particle with the corresponding track ID
 *  @return <b> else: </b> "---"
 */
G4String GODDeSS_DataStorage::GetNameOfPhotonDetectorThatWasHit(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(PhotonDetectorHitNameKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(PhotonDetectorHitNameKey) );
	}

	return "---";
}



/**
 *  @return <b> if the entry exists for this particle: </b> particle ID of the particle with the corresponding track ID
 *  @return <b> else: </b> 999999999 (one cannot use "0", as this is used for Geant4 specific particles (optical photon, geantino,...))
 */
G4int GODDeSS_DataStorage::GetParticleID(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(ParticleIDKey) ) return boost::any_cast<G4int>( (ParticleVector[trackID]).at(ParticleIDKey) );
	}

	return 999999999; // one cannot use "0", this is used for Geant4 specific particles (optical photon, geantino,...)
}

/**
 *  @return <b> if the entry exists for this particle: </b> name of the particle with the corresponding track ID
 *  @return <b> else: </b> "---"
 */
G4String GODDeSS_DataStorage::GetParticleName(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(ParticleNameKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(ParticleNameKey) );
	}

	return "---";
}

/**
 *  @return <b> if the entry exists for this particle: </b> G4bool, if the particle with the corresponding track ID is a primary particle
 *  @return <b> else: </b> true
 */
G4bool GODDeSS_DataStorage::GetParticleIsPrimary(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(IsPrimaryKey) ) return boost::any_cast<G4bool>( (ParticleVector[trackID]).at(IsPrimaryKey) );
	}

	return true;
}

/**
 *  @return <b> if the entry exists for this particle: </b> production mechanism of the particle with the corresponding track ID
 *  @return <b> else: </b> "---"
 */
G4String GODDeSS_DataStorage::GetProductionMechanism(G4int trackID) const
{
	if(TrackIDHasEntries(trackID))
	{
		if( (ParticleVector[trackID]).count(ProductionMechanismKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(ProductionMechanismKey) );
	}

	return "---";
}



/**
 *  Function to set the class variables to default values:
 */
void GODDeSS_DataStorage::SetDefaults()
{
	/** - start momentum (for resimulation purpose) of a photon hitting the sensitive detector \code PreviousPreStepPointMomentum = G4ThreeVector(NAN, NAN, NAN); \endcode */
	PreviousPreStepPointMomentum = G4ThreeVector(NAN, NAN, NAN);

	/** - start position (for resimulation purpose) of a photon hitting the sensitive detector \code PreviousPreStepPointPosition = G4ThreeVector(NAN, NAN, NAN); \endcode */
	PreviousPreStepPointPosition = G4ThreeVector(NAN, NAN, NAN);

	/** - start time (for resimulation purpose) of a photon hitting the sensitive detector \code PreviousPreStepPointTime = NAN; \endcode */
	PreviousPreStepPointTime = NAN;

	/** - polarisation (for resimulation purpose) of a photon hitting the sensitive detector \code PreviousPreStepPointPolarisation = G4ThreeVector(NAN, NAN, NAN); \endcode */
	PreviousPreStepPointPolarisation = G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  @return <b> if the entry exists for this track ID: </b> true
 *  @return <b> else: </b> false
 */
G4bool GODDeSS_DataStorage::TrackIDExistsInVector(G4int trackID) const
{
	if(ParticleVector.size() > (unsigned) trackID) return true;
	else return false;
}

/**
 *  @return <b> if the entry exists for this track ID and has entries: </b> true
 *  @return <b> else: </b> false
 */
G4bool GODDeSS_DataStorage::TrackIDHasEntries(G4int trackID) const
{
	if( TrackIDExistsInVector(trackID) && (ParticleVector[trackID]).size() != 0) return true;
	else return false;
}

/**
 *  Function to tell the UserEventInformation, that a particle with the given track ID as been simulated and its properties will be saved. <b> This function has to be called before any information about the particle can be saved! </b>
 */
void GODDeSS_DataStorage::FillEmptyEntriesUpToCurrentTrack(G4Track * theTrack)
{
	G4int trackID = theTrack->GetTrackID();

	if(!TrackIDExistsInVector(trackID))
	{
		for(int iter = ParticleVector.size(); iter < trackID + 1; iter++)
		{
			ParticleVector.push_back(std::map<G4String, boost::any>());
		}

		G4String partName = theTrack->GetParticleDefinition()->GetParticleName();
		SetParticleName(trackID, partName);
		SetParticleID(trackID, G4ParticleTable::GetParticleTable()->FindParticle(partName)->GetPDGEncoding());
		G4bool isPrimaryParticle = (! theTrack->GetParentID());
		if (isPrimaryParticle) SetParticleIsPrimary(trackID);
		else SetProductionMechanism(trackID, ((G4String) theTrack->GetCreatorProcess()->GetProcessName()));
	}
}
