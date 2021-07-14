/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef PHOTONDETECTORCONSTRUCTOR_H
#define PHOTONDETECTORCONSTRUCTOR_H

#include <G4VModularPhysicsList.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VPhysicalVolume.hh>
#include <G4OpticalSurface.hh>

#include <G4PhotonDetector.hh>
#include <PropertyToolsManager.hh>
#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class making it easier and more flexible to create G4PhotonDetector objects. <b> Additionally, it automatically deals with the registration of the needed physics! </b>
class PhotonDetectorConstructor
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 *  - adds the needed physics processes (specified in PhotonDetectorPhysicsList) to the physics list
	 */
	PhotonDetectorConstructor( G4VModularPhysicsList * userPhysicsList,	/**< physics list which is used for the simulation */
				   PropertyToolsManager * propertyTools,	/**< pointer to the PropertyToolsManager that is to be used */
				   GODDeSS_DataStorage * dataStorage,		/**< pointer to the GODDeSS_DataStorage that is to be used */
				   G4bool searchOverlaps,			/**< Geant should search for overlaps when placing the physical volumes of G4PhotonDetector%'s ("true" or "false") */
				   int verbose = 0				/**< verbosity level for adding the needed physics processes (specified in PhotonDetectorPhysicsList) to the physics list (default: 0) */
				 )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	{
		LoadPhysicsList(userPhysicsList, verbose);

		SetDefaults();
	}

	/**
	 *  Destructor (empty)
	 */
	~PhotonDetectorConstructor()
	{
	}



	// These are overloaded functions. They are calling the main function (defined below) using default values for the parameters which are not specified.
	#include <PhotonDetectorConstructor.icc>



	// These are the main function for creating the volumes. They will be overloaded to allow setting the optional parameters via the function call or via set-functions.

	// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	vector<G4PhotonDetector *> ConstructPhotonDetectors( G4double edge_length,						/**< edge length of the G4PhotonDetector%'s sensitive area */
							      G4VPhysicalVolume* mother_volume,					/**< G4PhotonDetector%'s mother volume */
							      G4int number_of_replications,					/**< number of replications that is to be created */
							      G4Transform3D replication_transf,					/**< transformation of each replication, relative to the previous one */
							      G4String photon_detector_name,					/**< name of the G4PhotonDetector%'s volume (it will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector) */
							      G4ThreeVector sensitive_surface_normal_relative_to_reference,	/**< surface normal of the G4PhotonDetector%'s sensitive area (relative to the reference volume) */
							      G4ThreeVector sensitive_surface_position_relative_to_reference,	/**< position of the G4PhotonDetector%'s sensitive area (relative to the reference volume) */
							      G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4PhotonDetector will be orientated (if 0, the mother volume is the reference volume) */
							      G4Transform3D reference_transformation				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
							    );

	// option to only one object
	G4PhotonDetector * ConstructPhotonDetector( G4double edge_length,						/**< edge length of the G4PhotonDetector%'s sensitive area */
						    G4VPhysicalVolume* mother_volume,					/**< G4PhotonDetector%'s mother volume */
						    G4String photon_detector_name,					/**< name of the G4PhotonDetector%'s volume (it will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector) */
						    G4ThreeVector sensitive_surface_normal_relative_to_reference,	/**< surface normal of the G4PhotonDetector%'s sensitive area (relative to the reference volume) */
						    G4ThreeVector sensitive_surface_position_relative_to_reference,	/**< position of the G4PhotonDetector%'s sensitive area (relative to the reference volume) */
						    G4VPhysicalVolume* reference_volume,					/**< reference volume relative to which the G4PhotonDetector will be orientated (if 0, the mother volume is the reference volume) */
						    G4Transform3D reference_transformation				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
						  );



	// These are the set-functions.
	/**
	 *  Function to set the name of the G4PhotonDetector%'s volumes (it will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector).
	 *
	 *  <b> The value set by this function will only apply to the very next G4PhotonDetector object that is created! </b>
	 */
	void SetPhotonDetectorName(G4String photon_detector_name)
	{
		PhotonDetectorName = photon_detector_name;
	}

	/**
	 *  Function to set the surface normal of the G4PhotonDetector%'s sensitive area (relative to the reference volume).
	 *
	 *  <b> The value set by this function will only apply to the very next G4PhotonDetector object that is created! </b>
	 */
	void SetSensitiveSurfaceNormalRelativeToReferenceVolume(G4ThreeVector sensitive_surface_normal_relative_to_reference)
	{
		SurfaceNormal_rel = sensitive_surface_normal_relative_to_reference;
	}

	/**
	 *  Function to set the position of the G4PhotonDetector%'s sensitive area (relative to the reference volume).
	 *
	 *  <b> The value set by this function will only apply to the very next G4PhotonDetector object that is created! </b>
	 */
	void SetSensitiveSurfacePositionRelativeToReferenceVolume(G4ThreeVector sensitive_surface_position_relative_to_reference)
	{
		SurfacePos_rel = sensitive_surface_position_relative_to_reference;
	}

	/**
	 *  Function to set the reference volume relative to which the G4PhotonDetector will be orientated.
	 *
	 *  <b> The value set by this function will only apply to the very next G4PhotonDetector object that is created! </b>
	 */
	void SetPhotonDetectorReferenceVolume(G4VPhysicalVolume* reference_volume)
	{
		ReferenceVolume_physical = reference_volume;
	}

	/**
	 *  Function to set the transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume).
	 *
	 *  <b> The value set by this function will only apply to the very next G4PhotonDetector object that is created! </b>
	 */
	void SetPhotonDetectorReferenceVolumeTransformation(G4Transform3D reference_volume_transformation)
	{
		ReferenceVolume_transformation = reference_volume_transformation;
	}



	void WriteEventIDToHitFile( G4int eventID /**< event ID which is to be written into the text file*/
				  );

	void WriteRunIDToHitFile( G4int runID /**< run ID which is to be written into the text file*/
				);

private:
	void applyVolumeCounter( G4String & physical_volume_name /**< name which is to be used for naming the volumes */
			       );
	void SetDefaults();

	void registerPhysics(G4VModularPhysicsList * physicsList, G4VPhysicsConstructor * physicsConstructor);
	void LoadPhysicsList(G4VModularPhysicsList * physicsList, int verbose = 0);


	G4bool SearchOverlaps;
	PropertyToolsManager * PropertyTools;
	std::vector<G4String> NamePrefixVector;
	GODDeSS_DataStorage * DataStorage;

	G4VPhysicalVolume * ReferenceVolume_physical;
	G4Transform3D ReferenceVolume_transformation;
	G4String PhotonDetectorName;
	G4ThreeVector SurfaceNormal_rel;
	G4ThreeVector SurfacePos_rel;
};

#endif
