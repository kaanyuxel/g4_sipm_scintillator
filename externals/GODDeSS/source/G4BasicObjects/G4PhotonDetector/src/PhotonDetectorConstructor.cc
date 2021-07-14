/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <boost/lexical_cast.hpp>

#include <PhotonDetectorConstructor.hh>
#include <PhotonDetectorPhysicsList.icc>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to construct a photon detector:
 *  - runs applyVolumeCounter() on \em photon_detector_name
 *  - creates a new G4PhotonDetector object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4PhotonDetector object
 */
G4PhotonDetector * PhotonDetectorConstructor::ConstructPhotonDetector( G4double edge_length,
								       G4VPhysicalVolume* mother_volume,
								       G4String photon_detector_name,
								       G4ThreeVector sensitive_surface_normal_relative_to_reference,
								       G4ThreeVector sensitive_surface_position_relative_to_reference,
								       G4VPhysicalVolume* reference_volume,
								       G4Transform3D reference_transform )
{
	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;

	// modify the name if another volume was named like this before
	applyVolumeCounter(photon_detector_name);

	// create the volumes
	G4PhotonDetector * photonDetector = new G4PhotonDetector( edge_length,
								  mother_volume,
								  photon_detector_name,
								  sensitive_surface_normal_relative_to_reference,
								  sensitive_surface_position_relative_to_reference,
								  reference_volume,
								  reference_transform,
								  SearchOverlaps,
								  PropertyTools,
								  DataStorage );

	// reset to the default values
	SetDefaults();

	// return the volumes
	return photonDetector;
}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
/**
 *  Function to construct construct replications of a photon detector:
 *  - runs applyVolumeCounter() on \em photon_detector_name
 *  - creates a new G4PhotonDetector objects
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return vector with pointers to the created G4Fibre objects
 */
vector<G4PhotonDetector *> PhotonDetectorConstructor::ConstructPhotonDetectors( G4double edge_length,
										 G4VPhysicalVolume* mother_volume,
										 G4int number_of_replications,
										 G4Transform3D replication_transf,
										 G4String photon_detector_name,
										 G4ThreeVector sensitive_surface_normal_relative_to_reference,
										 G4ThreeVector sensitive_surface_position_relative_to_reference,
										 G4VPhysicalVolume* reference_volume,
										 G4Transform3D reference_transform )
{
	std::vector<G4PhotonDetector *> objectsVector;

	if(number_of_replications < 2)
	{
		G4cerr << "##########" << "# WARNING:" << "# Replications of a photon detector is to be build, but the number of replications is to low. No photon detector will be build." << "##########" << G4endl;
		return objectsVector;
	}

	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;


	G4RotationMatrix rotation = G4RotationMatrix();
	G4ThreeVector translation = G4ThreeVector(0, 0, 0);
	G4RotationMatrix rotation_replication = replication_transf.getRotation();
	G4ThreeVector translation_replication = replication_transf.getTranslation();

	for(int iter = 0; iter < number_of_replications; iter++)
	{
		if(iter != 0)
		{
			rotation *= rotation_replication;
			translation += translation_replication;
		}

		replication_transf = G4Transform3D(rotation, translation);

		// modify the name if another volume was named like this before
		G4String temp_photon_detector_name = photon_detector_name;
		applyVolumeCounter(photon_detector_name);

		// create the volumes
		G4PhotonDetector * photonDetector = new G4PhotonDetector( edge_length,
									   mother_volume,
									   temp_photon_detector_name,
									   sensitive_surface_normal_relative_to_reference,
									   sensitive_surface_position_relative_to_reference,
									   reference_volume,
									   reference_transform,
									   SearchOverlaps,
									   PropertyTools,
									   DataStorage,
									   replication_transf );

		objectsVector.push_back(photonDetector);
	}

	// reset to the default values
	SetDefaults();

	// return the volumes
	return objectsVector;
}



/**
 *  Function for avoiding multiple volumes with the same name:
 *  - a string variable containing the name prefix has to be passed to this function
 *  - the name prefix will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector
 *  - the extended name prefix will be passed to the original string variable
 */
void PhotonDetectorConstructor::applyVolumeCounter(G4String & physical_volume_name)
{
	if(NamePrefixVector.size())
	{
		G4int i = 1;

		// if nothing is found, stop looping
		G4bool continueWhileLoop = true;

		// loop as long as something is found
		while(continueWhileLoop)
		{
			continueWhileLoop = false;

			for(unsigned int iter = 0; iter < NamePrefixVector.size(); iter++)
			{
				// if something is found, increment the counter and continue the while loop and continue looping
				if((i == 1 && NamePrefixVector[iter] == physical_volume_name) || NamePrefixVector[iter] == (physical_volume_name + "_" + boost::lexical_cast<std::string>(i)).c_str())
				{
					i++;
					continueWhileLoop = true;
					break;
				}
			}

			// if nothing has found at last, set new physical_volume_name
			if(!continueWhileLoop && i != 1) physical_volume_name += "_" + boost::lexical_cast<std::string>(i);
		}
	}

	NamePrefixVector.push_back(physical_volume_name);
}



/**
 *  Function to set the class variables to default values:
 */
void PhotonDetectorConstructor::SetDefaults()
{
	/** - pointer to the reference volume (physical volume relative to which the G4PhotonDetector will be orientated) \code ReferenceVolume_physical = 0; \endcode */
	ReferenceVolume_physical = 0;

	/** - transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) \code ReferenceVolume_transformation = G4Transform3D(); \endcode */
	ReferenceVolume_transformation = G4Transform3D();

	/** - name of the G4PhotonDetector%'s volumes (it will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector) \code PhotonDetectorName = "photonDetector"; \endcode */
	PhotonDetectorName = "photonDetector";

	/** - surface normal of the G4PhotonDetector%'s sensitive area (relative to the reference volume) \code SurfaceNormal_rel = G4ThreeVector(0., 0., 1.); \endcode */
	SurfaceNormal_rel = G4ThreeVector(0., 0., 1.);

	/** - position of the G4PhotonDetector%'s sensitive area (relative to the reference volume) \code SurfacePos_rel = G4ThreeVector(0., 0., 0.); \endcode */
	SurfacePos_rel = G4ThreeVector(0., 0., 0.);
}



/**
 *  Function to write the event ID into the text file containing the photons that hit a G4PhotonDetector.
 */
void PhotonDetectorConstructor::WriteEventIDToHitFile(G4int eventID)
{
	std::vector<G4String> photonSensitiveDetectorVolumeNames = DataStorage->GetPhotonSensitiveDetectorVolumeNames();

	if(photonSensitiveDetectorVolumeNames.size())
	{
		for(unsigned int iter = 0; iter < photonSensitiveDetectorVolumeNames.size(); iter++)
		{
			G4String hitFileName = DataStorage->GetPhotonDetectorHitFile();
			hitFileName = hitFileName.substr(0, hitFileName.rfind(".")) + "_" + photonSensitiveDetectorVolumeNames[iter] + hitFileName.substr(hitFileName.rfind("."));

			std::ofstream hitFile;
			hitFile.open( hitFileName.c_str(), std::ios_base::out | std::ios_base::app );
			hitFile << "\nEventID:\t" << eventID << "\n\n";
			hitFile.close();
		}
	}
}

/**
 *  Function to write the run ID into the text file containing the photons that hit a G4PhotonDetector.
 */
void PhotonDetectorConstructor::WriteRunIDToHitFile(G4int runID)
{
	std::vector<G4String> photonSensitiveDetectorVolumeNames = DataStorage->GetPhotonSensitiveDetectorVolumeNames();

	if(photonSensitiveDetectorVolumeNames.size())
	{
		for(unsigned int iter = 0; iter < photonSensitiveDetectorVolumeNames.size(); iter++)
		{
			G4String hitFileName = DataStorage->GetPhotonDetectorHitFile();
			hitFileName = hitFileName.substr(0, hitFileName.rfind(".")) + "_" + photonSensitiveDetectorVolumeNames[iter] + hitFileName.substr(hitFileName.rfind("."));

			std::ofstream hitFile;
			hitFile.open( hitFileName.c_str(), std::ios_base::out | std::ios_base::app );
			hitFile << "\nRunID:\t\t" << runID << "\n";
			hitFile.close();
		}
	}
}
