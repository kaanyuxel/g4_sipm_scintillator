/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <boost/lexical_cast.hpp>

#include <G4PhysicalVolumeStore.hh>

#include <ScintillatorTileConstructor.hh>
#include <ScintillatorTilePhysicsList.icc>



// class variables begin with capital letters, local variables with small letters



// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
/**
 *  Function to construct construct replications of a scintillator tile:
 *  - runs applyVolumeCounter() on \em scintillator_name
 *  - creates a new G4ScintillatorTile objects
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return vector with pointers to the created G4Fibre objects
 */
vector<G4ScintillatorTile *> ScintillatorTileConstructor::ConstructScintillators( G4ThreeVector scintillator_dimensions,
										   G4String scintillator_properties_file,
										   G4VPhysicalVolume* mother_volume,
										   G4int number_of_replications,
										   G4Transform3D replication_transf,
										   G4String scintillator_name,
										   G4Transform3D transf_relative_to_reference,
										   G4VPhysicalVolume* reference_volume,
										   G4Transform3D reference_transform,
										   G4bool cut_dimple,
										   G4double dimple_radius,
										   G4ThreeVector dimple_centre_relative_to_scintillator_centre,
										   G4bool twin_tile,
										   G4String twin_tile_properties_file,
										   G4bool constructSensitiveDetector )
{
	std::vector<G4ScintillatorTile *> objectsVector;

	if(number_of_replications < 2)
	{
		G4cerr << "##########" << "# WARNING:" << "# Replications of a scintillator tile is to be build, but the number of replications is to low. No scintillator tile will be build." << "##########" << G4endl;
		return objectsVector;
	}

	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;
	if(!cut_dimple)
	{
		ScintillatorDimpleRadius = NAN;
		ScintillatorDimpleCentrePosition = G4ThreeVector(NAN, NAN, NAN);
	}


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
		G4String temp_scintillator_name = scintillator_name;
		applyVolumeCounter(temp_scintillator_name);

		// create the volumes
		G4ScintillatorTile * scintillator = new G4ScintillatorTile( scintillator_dimensions,
									     scintillator_properties_file,
									     mother_volume,
									     temp_scintillator_name,
									     transf_relative_to_reference,
									     reference_volume,
									     reference_transform,
									     cut_dimple,
									     dimple_radius,
									     dimple_centre_relative_to_scintillator_centre,
									     twin_tile,
									     twin_tile_properties_file,
									     constructSensitiveDetector,
									     SearchOverlaps,
									     PropertyTools,
									     SurfaceModel,
									     SurfaceFinish,
									     DataStorage,
									     replication_transf );

		objectsVector.push_back(scintillator);
	}

	// reset to the default values
	SetDefaults();

	// return the volumes
	return objectsVector;
}



/**
 *  Function to construct a scintillator tile:
 *  - runs applyVolumeCounter() on \em scintillator_name
 *  - creates a new G4ScintillatorTile object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4ScintillatorTile object
 */
G4ScintillatorTile * ScintillatorTileConstructor::ConstructScintillator( G4ThreeVector scintillator_dimensions,
									 G4String scintillator_property_file,
									 G4VPhysicalVolume* mother_volume,
									 G4String scintillator_name,
									 G4Transform3D transf_relative_to_reference,
									 G4VPhysicalVolume* reference_volume,
									 G4Transform3D reference_transform,
									 G4bool cut_dimple,
									 G4double dimple_radius,
									 G4ThreeVector dimple_centre_relative_to_scintillator_centre,
									 G4bool twin_tile,
									 G4String twin_tile_property_file,
									 G4bool constructSensitiveDetector )
{
	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;
	if(!cut_dimple)
	{
		ScintillatorDimpleRadius = NAN;
		ScintillatorDimpleCentrePosition = G4ThreeVector(NAN, NAN, NAN);
	}

	// modify the name if another volume was named like this before
	applyVolumeCounter(scintillator_name);

	// create the volumes
	G4ScintillatorTile * scintillator = new G4ScintillatorTile( scintillator_dimensions,
								    scintillator_property_file,
								    mother_volume,
								    scintillator_name,
								    transf_relative_to_reference,
								    reference_volume,
								    reference_transform,
								    cut_dimple,
								    dimple_radius,
								    dimple_centre_relative_to_scintillator_centre,
								    twin_tile,
								    twin_tile_property_file,
								    constructSensitiveDetector,
								    SearchOverlaps,
								    PropertyTools,
								    SurfaceModel,
								    SurfaceFinish,
								    DataStorage );

	// reset to the default values
	SetDefaults();

	// return the volumes
	return scintillator;
}



/**
 *  Function to create a wrapping around a G4ScintillatorTile by using the surface properties and Look Up Tables:
 *  - changes the surface properties (G4Wrapping::ConstructLUTWrapping())
 *  - resets the class variables to default values (using SetDefaults())
 */
void ScintillatorTileConstructor::ConstructLUTWrapping( G4ScintillatorTile* scinti_to_be_wrapped,
							G4OpticalSurfaceModel surface_model,
							G4OpticalSurfaceFinish surface_finish,
							G4SurfaceType surface_type,
							G4double roughness,
							G4VPhysicalVolume* surrounding_volume )
{
	// create the volumes
	G4Wrapping::ConstructLUTWrapping(scinti_to_be_wrapped, surface_model, surface_finish, surface_type, roughness, surrounding_volume);

	// reset to the default values
	SetDefaults();
}



/**
 *  Function to construct a wrapping volume around a G4ScintillatorTile:
 *  - runs applyVolumeCounter() on \em wrapping_name
 *  - creates a new G4Wrapping object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4Wrapping object
 */
G4Wrapping * ScintillatorTileConstructor::ConstructWrapping( G4ScintillatorTile* scinti_to_be_wrapped,
							     G4String wrapping_property_file,
							     G4VPhysicalVolume* mother_volume,
							     G4String wrapping_name,
							     std::vector<G4VPhysicalVolume *> cut_volumes,
							     G4bool constructSensitiveDetector )
{
	// the default cut_volumes are all physical volumes
	if(!WrappingCutVolumesSet) cut_volumes = *G4PhysicalVolumeStore::GetInstance();

	// modify the name if another volume was named like this before
	applyVolumeCounter(wrapping_name);

	// create the volumes
	G4Wrapping * wrapping = new G4Wrapping(scinti_to_be_wrapped,
					       wrapping_property_file,
					       mother_volume,
					       wrapping_name,
					       cut_volumes,
					       constructSensitiveDetector,
					       SearchOverlaps,
					       PropertyTools,
					       SurfaceModel,
					       SurfaceFinish,
					       DataStorage);

	// reset to the default values
	SetDefaults();

	// return the volumes
	return wrapping;
}



/**
 *  Function for avoiding multiple volumes with the same name:
 *  - a string variable containing the name prefix has to be passed to this function
 *  - the name prefix will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre
 *  - the extended name prefix will be passed to the original string variable
 */
void ScintillatorTileConstructor::applyVolumeCounter(G4String & physical_volume_name)
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
void ScintillatorTileConstructor::SetDefaults()
{
	/** - G4bool, whether a sensitive detector is to be constructed \code ConstructSensitiveDetector = false; \endcode */
	ConstructSensitiveDetector = false;

	/** - model for the optical surfaces \code SurfaceModel = unified; \endcode */
	SurfaceModel = unified;

	/** - finish for the optical surfaces \code SurfaceFinish = ground; \endcode */
	SurfaceFinish = ground;

	/** - name of the G4ScintillatorTile%'s volumes (it will be extended to distinguish between different G4ScintillatorTile%s and different volumes of one G4ScintillatorTile) \code ScintillatorName = "scintillator"; \endcode */
	ScintillatorName = "scintillator";

	/** - G4ScintillatorTile%'s transformation (relative to the reference volume) \code ScintillatorTransformation = G4Transform3D(); \endcode */
	ScintillatorTransformation = G4Transform3D();

	/** - pointer to the reference volume (physical volume relative to which the G4ScintillatorTile will be orientated) \code ScintillatorReferenceVolume_physical = 0; \endcode */
	ScintillatorReferenceVolume_physical = 0;

	/** - transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) \code ReferenceVolume_transformation = G4Transform3D(); \endcode */
	ReferenceVolume_transformation = G4Transform3D();

	/** - G4bool, if a dimple is to be cut out out of the G4ScintillatorTile \code ScintillatorCutDimple = false; \endcode */
	ScintillatorCutDimple = false;

	/** - radius of the dimple \code ScintillatorDimpleRadius = NAN; \endcode */
	ScintillatorDimpleRadius = NAN;

	/** - position of the dimple's centre relative to the G4ScintillatorTile%'s centre \code ScintillatorDimpleCentrePosition = G4ThreeVector(NAN, NAN, NAN); \endcode */
	ScintillatorDimpleCentrePosition = G4ThreeVector(NAN, NAN, NAN);

	/** - G4bool, if the G4ScintillatorTile is to be a twin tile \code ScintillatorTwinTile = false; \endcode */
	ScintillatorTwinTile = false;

	/** - path to the file containing the the properties of the twin tile's reflective layer \code ScintillatorTwinTilePropertyFile = ""; \endcode */
	ScintillatorTwinTilePropertyFile = "";

	/** - G4Wrapping%'s mother volume \code WrappingMotherVolume_physical = 0; \endcode */
	WrappingMotherVolume_physical = 0;

	/** - name of the G4Wrapping%'s volumes (it will be extended to distinguish between different G4Wrapping%s and different volumes of one G4Wrapping) \code WrappingName = "wrapping"; \endcode */
	WrappingName = "wrapping";

	/** - G4bool, if the volumes (that is to be cut out of the G4Wrapping) have been set \code WrappingCutVolumesSet = false; \endcode */
	WrappingCutVolumesSet = false;

	/** - volumes that is to be cut out of the G4Wrapping \code if(WrappingCutVolumes.size() == 0) WrappingCutVolumes = std::vector<G4VPhysicalVolume *>();
 else WrappingCutVolumes.clear(); \endcode */
	if(WrappingCutVolumes.size() == 0) WrappingCutVolumes = std::vector<G4VPhysicalVolume *>();
	else WrappingCutVolumes.clear();
}
