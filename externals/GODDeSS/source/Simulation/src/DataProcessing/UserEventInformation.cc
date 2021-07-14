/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4OpticalPhoton.hh>

#include <UserEventInformation.hh>
#include <G4OpticalPhoton.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  @return <b> if the entry exists for this track ID: </b> true
 *  @return <b> else: </b> false
 */
G4bool UserEventInformation::TrackIDExistsInVector(G4int trackID) const
{
	if(ParticleVector.size() > (unsigned) trackID) return true;
	else return false;
}

// G4bool UserEventInformation::TrackIDHasEntries(G4int trackID) const
// {
// 	if( TrackIDExistsInVector(trackID) && (ParticleVector[trackID]).size() != 0) return true;
// 	else return false;
// }

/**
 *  Function to tell the UserEventInformation, that a particle with the given track ID as been simulated and its properties will be saved. <b> This function has to be called before any information about the particle can be saved! </b>
 */
void UserEventInformation::FillEmptyEntriesUpToTrackID(G4int trackID)
{
	if(!TrackIDExistsInVector(trackID))
	{
		for(int iter = ParticleVector.size(); iter < trackID + 1; iter++)
		{
			ParticleVector.push_back(std::map<G4String, boost::any>());
		}
	}
}

/**
 *  @return <b> if the entry exists for this particle: </b> particle ID of the particle with the corresponding track ID
 *  @return <b> else: </b> 999999999 (one cannot use "0", as this is used for Geant4 specific particles (optical photon, geantino,...))
 */
G4int UserEventInformation::GetParticleID(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(ParticleIDKey) ) return boost::any_cast<G4int>( (ParticleVector[trackID]).at(ParticleIDKey) );
	else return 999999999; // one cannot use "0", this is used for Geant4 specific particles (optical photon, geantino,...)
}

/**
 *  @return <b> if the entry exists for this particle: </b> name of the particle with the corresponding track ID
 *  @return <b> else: </b> "---"
 */
G4String UserEventInformation::GetParticleName(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(ParticleNameKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(ParticleNameKey) );
	else return "---";
}

/**
 *  @return <b> if the entry exists for this particle: </b> G4bool, if the particle with the corresponding track ID is a primary particle
 *  @return <b> else: </b> true
 */
G4bool UserEventInformation::GetParticleIsPrimary(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(IsPrimaryKey) ) return boost::any_cast<G4bool>( (ParticleVector[trackID]).at(IsPrimaryKey) );
	else return true;
}

/**
 *  @return <b> if the entry exists for this particle: </b> production mechanism of the particle with the corresponding track ID
 *  @return <b> else: </b> "---"
 */
G4String UserEventInformation::GetProductionMechanism(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(ProductionMechanismKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(ProductionMechanismKey) );
	else return "---";
}

/**
 *  @return <b> if the entry exists for this particle: </b> initial position of the particle with the corresponding track ID
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector UserEventInformation::GetInitialPosition(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(InitialPositionKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(InitialPositionKey) );
	else return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  @return <b> if the entry exists for this particle: </b> creation time of the particle with the corresponding track ID
 *  @return <b> else: </b> NAN
 */
G4double UserEventInformation::GetGlobalCreationTime(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(GlobalCreationTimeKey) ) return boost::any_cast<G4double>( (ParticleVector[trackID]).at(GlobalCreationTimeKey) );
	else return NAN;
}

/**
 *  @return <b> if the entry exists for this particle: </b> initial momentum of the particle with the corresponding track ID
 *  @return <b> else: </b> G4ThreeVector(NAN, NAN, NAN)
 */
G4ThreeVector UserEventInformation::GetInitialMomentum(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(InitialMomentumKey) ) return boost::any_cast<G4ThreeVector>( (ParticleVector[trackID]).at(InitialMomentumKey) );
	else return G4ThreeVector(NAN, NAN, NAN);
}

/**
 *  Function to save the parent particle of the particle with the corresponding track ID.
 */
void UserEventInformation::SetParent(G4int trackID, G4int parentTrackID)
{
	(ParticleVector[trackID])[ParentTrackIDKey] = parentTrackID;
	if(parentTrackID != 0)
	{
		(ParticleVector[trackID])[ParentNameKey] = GetParticleName(parentTrackID);
		(ParticleVector[trackID])[ParentIDKey] = GetParticleID(parentTrackID);
	}
}
/**
 *  @return <b> if the entry exists for this particle: </b> parent particle's track ID of the particle with the corresponding track ID
 *  @return <b> else: </b> -1
 */
G4int UserEventInformation::GetParentTrackID(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(ParentTrackIDKey) ) return boost::any_cast<G4int>( (ParticleVector[trackID]).at(ParentTrackIDKey) );
	else return -1;
}
/**
 *  @return <b> if the entry exists for this particle: </b> parent particle's name of the particle with the corresponding track ID
 *  @return <b> else: </b> "---"
 */
G4String UserEventInformation::GetParentName(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(ParentNameKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(ParentNameKey) );
	else return "---";
}
/**
 *  @return <b> if the entry exists for this particle: </b> parent particle's particle ID of the particle with the corresponding track ID
 *  @return <b> else: </b> 999999999
 */
G4int UserEventInformation::GetParentID(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(ParentIDKey) ) return boost::any_cast<G4int>( (ParticleVector[trackID]).at(ParentIDKey) );
	else return 999999999; // one cannot use "0", this is used for Geant4 specific particles (optical photon, geantino,...)
}
/**
 *  @return G4bool, if the parent particle of the particle with the corresponding track ID is a primary particle (if the parent particle is an optical photon, the request will be repeated with its track ID)
 */
G4bool UserEventInformation::GetParentIsPrimary(G4int trackID) const
{
	while(GetParentName(trackID) == G4OpticalPhoton::Definition()->GetParticleName() && ( ! GetParticleIsPrimary(GetParentTrackID(trackID)) )) trackID = GetParentTrackID(trackID);

	return GetParticleIsPrimary(GetParentTrackID(trackID));
}

/**
 *  @return <b> if the entry exists for this particle: </b> absorption time of the particle with the corresponding track ID
 *  @return <b> else: </b> NAN
 */
G4double UserEventInformation::GetGlobalAbsorptionTime(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(GlobalAbsorptionTimeKey) ) return boost::any_cast<G4double>( (ParticleVector[trackID]).at(GlobalAbsorptionTimeKey) );
	else return NAN;
}

/**
 *  @return <b> if the entry exists for this particle: </b> name of the volume, the particle with the corresponding track ID was absorbed in
 *  @return <b> else: </b> "---"
 */
G4String UserEventInformation::GetAbsorbedIn(G4int trackID) const
{
	if( (ParticleVector[trackID]).count(AbsorbedInKey) ) return boost::any_cast<G4String>( (ParticleVector[trackID]).at(AbsorbedInKey) );
	else return "---";
}
