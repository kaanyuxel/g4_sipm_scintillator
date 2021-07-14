/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef SCINTILLATORTILECONSTRUCTOR_H
#define SCINTILLATORTILECONSTRUCTOR_H

#include <G4VModularPhysicsList.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VPhysicalVolume.hh>
#include <G4OpticalSurface.hh>

#include <G4ScintillatorTile.hh>
#include <G4Wrapping.hh>
#include <PropertyToolsManager.hh>
#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class making it easier and more flexible to create G4ScintillatorTile and G4Wrapping objects. <b> Additionally, it automatically deals with the registration of the needed physics! </b>
class ScintillatorTileConstructor
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 *  - adds the needed physics processes (specified in ScintillatorTilePhysicsList) to the physics list
	 */
	ScintillatorTileConstructor( G4VModularPhysicsList * userPhysicsList,	/**< physics list which is used for the simulation */
				     PropertyToolsManager * propertyTools,	/**< pointer to the PropertyToolsManager that is to be used */
				     GODDeSS_DataStorage * dataStorage,		/**< pointer to the GODDeSS_DataStorage that is to be used */
				     G4bool searchOverlaps,			/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
				     int verbose = 0				/**< verbosity level for adding the needed physics processes (specified in FibrePhysicsList) to the physics list (default: 0) */
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
	 *  Destructor:
	 *  - deletes the data storage object
	 */
	~ScintillatorTileConstructor()
	{
	}



	// These are overloaded functions. They are calling the main function (defined below) using default values for the parameters which are not specified.
	#include <ScintillatorTileConstructor.icc>



	// These are the main function for creating the volumes. They will be overloaded to allow setting the optional parameters via the function call or via set-functions.

	// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	vector<G4ScintillatorTile *> ConstructScintillators( G4ThreeVector scintillator_dimensions,				/**< dimensions of the G4ScintillatorTile */
							      G4String scintillator_property_file,				/**< path to the file containing the G4ScintillatorTile%'s properties */
							      G4VPhysicalVolume* mother_volume,					/**< G4ScintillatorTile%'s mother volume */
							      G4int number_of_replications,					/**< number of replications that is to be created */
							      G4Transform3D replication_transf,					/**< transformation of each replication, relative to the previous one */
							      G4String scintillator_name,					/**< name of the G4ScintillatorTile%'s volume (it will be extended to distinguish between different G4ScintillatorTile%s and different volumes of one G4ScintillatorTile) */
							      G4Transform3D transf_relative_to_reference,			/**< G4ScintillatorTile%'s transformation relative to the reference volume */
							      G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4ScintillatorTile will be orientated (if 0, the mother volume is the reference volume) */
							      G4Transform3D reference_transformation,				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
							      G4bool cut_dimple,						/**< a dimple is to be cut out out of the G4ScintillatorTile ("true" or "false") */
							      G4double dimple_radius,						/**< radius of the dimple */
							      G4ThreeVector dimple_centre_relative_to_scintillator_centre,	/**< position of the dimple's centre relative to the G4ScintillatorTile%'s centre */
							      G4bool twin_tile,							/**< the G4ScintillatorTile is to be a twin tile ("true" or "false") */
							      G4String twin_tile_property_file,					/**< path to the file containing the properties of the twin tile's reflective layer */
							      G4bool constructSensitiveDetector					/**< a sensitive detector is to be constructed ("true" or "false") */
							    );

	// option to only one object
	G4ScintillatorTile * ConstructScintillator( G4ThreeVector scintillator_dimensions,				/**< dimensions of the G4ScintillatorTile */
						    G4String scintillator_property_file,				/**< path to the file containing the G4ScintillatorTile%'s properties */
						    G4VPhysicalVolume* mother_volume,					/**< G4ScintillatorTile%'s mother volume */
						    G4String scintillator_name,						/**< name of the G4ScintillatorTile%'s volume (it will be extended to distinguish between different G4ScintillatorTile%s and different volumes of one G4ScintillatorTile) */
						    G4Transform3D transf_relative_to_reference,				/**< G4ScintillatorTile%'s transformation relative to the reference volume */
						    G4VPhysicalVolume* reference_volume,					/**< reference volume relative to which the G4ScintillatorTile will be orientated (if 0, the mother volume is the reference volume) */
						    G4Transform3D reference_transformation,				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
						    G4bool cut_dimple,							/**< a dimple is to be cut out out of the G4ScintillatorTile ("true" or "false") */
						    G4double dimple_radius,						/**< radius of the dimple */
						    G4ThreeVector dimple_centre_relative_to_scintillator_centre,	/**< position of the dimple's centre relative to the G4ScintillatorTile%'s centre */
						    G4bool twin_tile,							/**< the G4ScintillatorTile is to be a twin tile ("true" or "false") */
						    G4String twin_tile_property_file,					/**< path to the file containing the properties of the twin tile's reflective layer */
						    G4bool constructSensitiveDetector					/**< a sensitive detector is to be constructed ("true" or "false") */
						  );



	// These are the set-functions.
	/**
	 *  Function to set the name of the G4ScintillatorTile%'s volumes (it will be extended to distinguish between different G4ScintillatorTile%s and different volumes of one G4ScintillatorTile).
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetScintillatorName(G4String scintillator_name)
	{
		ScintillatorName = scintillator_name;
	}

	/**
	 *  Function to set the G4ScintillatorTile%'s transformation relative to the reference volume.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetScintillatorTransformation(G4Transform3D transf_relative_to_reference)
	{
		ScintillatorTransformation = transf_relative_to_reference;
	}

	/**
	 *  Function to set the reference volume relative to which the G4ScintillatorTile will be orientated.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetScintillatorReferenceVolume(G4VPhysicalVolume* reference_volume)
	{
		ScintillatorReferenceVolume_physical = reference_volume;
	}

	/**
	 *  Function to set the transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume).
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetScintillatorReferenceVolumeTransformation(G4Transform3D reference_volume_transformation)
	{
		ReferenceVolume_transformation = reference_volume_transformation;
	}

	/**
	 *  Function to set the dimple's radius and centre position.\n
	 *  A sphere with this specifications will be cut out of the G4ScintillatorTile.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetScintillatorCutDimple(G4double dimple_radius, G4ThreeVector dimple_centre_relative_to_scintillator_centre)
	{
		ScintillatorCutDimple = true;
		ScintillatorDimpleRadius = dimple_radius;
		ScintillatorDimpleCentrePosition = dimple_centre_relative_to_scintillator_centre;
	}

	/**
	 *  Function to set the path to the file containing the properties of the twin tile's reflective layer.\n
	 *  If this function is used, the G4ScintillatorTile will be split into two volumes with a reflective layer between them.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetScintillatorTwinTile(G4String twin_tile_property_file)
	{
		ScintillatorTwinTile = true;
		ScintillatorTwinTilePropertyFile = twin_tile_property_file;
	}

	/**
	 *  Function to make the next G4ScintillatorTile or G4Wrapping being constructed with a sensitive detector.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void ConstructASensitiveDetector()
	{
		ConstructSensitiveDetector = true;
	}



	void ConstructLUTWrapping( G4ScintillatorTile* scinti_to_be_wrapped,	/**< G4ScintillatorTile that is to be wrapped with the G4Wrapping */
				   G4OpticalSurfaceModel surfaceModel,		/**< model for the optical surfaces */
				   G4OpticalSurfaceFinish surfaceFinish,		/**< finish for the optical surfaces */
				   G4SurfaceType surfaceType,			/**< type of the optical surfaces */
				   G4double roughness,				/**< roughness of the wrapping surface (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
				   G4VPhysicalVolume* surrounding_volume = 0	/**< volume sharing the surface with the G4ScintillatorTile that is to be wrapped (if 0, the whole G4ScintillatorTile%'s surface will be wrapped, otherwise only the shared surface will be wrapped) */
				 );



	// option to only one object
	G4Wrapping * ConstructWrapping( G4ScintillatorTile* scinti_to_be_wrapped,	/**< G4ScintillatorTile that is to be wrapped with the G4Wrapping */
					G4String wrapping_property_file,		/**< path to the file containing the G4Wrapping%'s properties */
					G4VPhysicalVolume* mother_volume,		/**< G4Wrapping%'s mother volume */
					G4String wrapping_name,				/**< name of the G4Wrapping%'s volume (it will be extended to distinguish between different G4Wrapping%s and different volumes of one G4Wrapping) */
					std::vector<G4VPhysicalVolume *> cut_volumes,	/**< volumes that is to be cut out of the G4Wrapping */
					G4bool constructSensitiveDetector		/**< a sensitive detector is to be constructed ("true" or "false") */
				      );



	// These are the set-functions.
	/**
	 *  Function to set the G4Wrapping%'s mother volume manually.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetWrappingMotherVolume(G4VPhysicalVolume * wrapping_mother_volume_physical)
	{
		WrappingMotherVolume_physical = wrapping_mother_volume_physical;
	}

	/**
	 *  Function to set the name of the G4Wrapping%'s volumes (it will be extended to distinguish between different G4Wrapping%s and different volumes of one G4Wrapping).
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetWrappingName(G4String wrapping_name)
	{
		WrappingName = wrapping_name;
	}

	/**
	 *  Function to specify which volumes is to be cut out of the G4Wrapping.\n
	 *  If this function is <b> not </b> used, all volumes that have been placed before (exept the world volume) will be cut out.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetWrappingCutVolumes(std::vector<G4VPhysicalVolume *> cut_volumes)
	{
		WrappingCutVolumes = cut_volumes;
		WrappingCutVolumesSet = true;
	}



	/**
	 *  Function to set the optical surface properties of G4ScintillatorTile or G4Wrapping.
	 *
	 *  <b> The value set by this function will only apply to the very next G4ScintillatorTile or G4Wrapping object that is created! </b>
	 */
	void SetOpticalSurfaceProperties( G4OpticalSurfaceModel surface_model,
					   G4OpticalSurfaceFinish surface_finish
					 )
	{
		SurfaceModel = surface_model;
		SurfaceFinish = surface_finish;
	}

private:
	void applyVolumeCounter( G4String & physical_volume_name /**< name which is to be used for naming the volumes */
			       );
	void SetDefaults();

	void registerPhysics(G4VModularPhysicsList * physicsList, G4VPhysicsConstructor * physicsConstructor);
	void LoadPhysicsList(G4VModularPhysicsList * physicsList, int verbose = 0);


	G4bool ConstructSensitiveDetector;
	G4bool SearchOverlaps;
	PropertyToolsManager * PropertyTools;
	std::vector<G4String> NamePrefixVector;
	G4OpticalSurfaceModel SurfaceModel;
	G4OpticalSurfaceFinish SurfaceFinish;
	GODDeSS_DataStorage * DataStorage;

	G4String ScintillatorName;
	G4Transform3D ScintillatorTransformation;
	G4VPhysicalVolume * ScintillatorReferenceVolume_physical;
	G4Transform3D ReferenceVolume_transformation;
	G4bool ScintillatorCutDimple;
	G4double ScintillatorDimpleRadius;
	G4ThreeVector ScintillatorDimpleCentrePosition;
	G4bool ScintillatorTwinTile;
	G4String ScintillatorTwinTilePropertyFile;
	G4VPhysicalVolume * WrappingMotherVolume_physical;
	G4String WrappingName;
	G4bool WrappingCutVolumesSet;
	std::vector<G4VPhysicalVolume *> WrappingCutVolumes;
};

#endif
