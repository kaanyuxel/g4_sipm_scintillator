/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <boost/lexical_cast.hpp>

#include <FibreConstructor.hh>
#include <FibrePhysicsList.icc>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to construct a <b> straight </b> fibre with a given <b> start </b> and <b> end point </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4Fibre object
 */
G4Fibre * FibreConstructor::ConstructFibre( G4String fibre_property_file,
					    G4VPhysicalVolume* mother_volume,
					    G4ThreeVector startPoint_relative_to_reference,
					    G4ThreeVector endPoint_relative_to_reference,
					    G4VPhysicalVolume* reference_volume,
					    G4Transform3D reference_transform,
					    G4String fibre_name_prefix,
					    G4bool onlyInsideMother,
					    G4bool cutAuntVolumesWithDaughters,
					    G4double startPoint_reflectivity,
					    G4double startPoint_sigmaAlpha,
					    G4double endPoint_reflectivity,
					    G4double endPoint_sigmaAlpha,
					    G4bool glued,
					    G4String glue_property_file,
					    G4String embedment_profile,
					    G4String where_is_embedment_to_be_flush_with_fibre,
					    std::vector<G4VPhysicalVolume *> volumes_with_embedment,
					    G4bool constructSensitiveDetector
					  )
{
	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;

	// modify the name if another volume was named like this before
	applyVolumeCounter(fibre_name_prefix);

	// create the volumes
	G4Fibre * fibre = new G4Fibre(fibre_property_file,
				      mother_volume,
				      startPoint_relative_to_reference,
				      endPoint_relative_to_reference,
				      reference_volume,
				      reference_transform,
				      fibre_name_prefix,
				      onlyInsideMother,
				      cutAuntVolumesWithDaughters,
				      startPoint_reflectivity,
				      startPoint_sigmaAlpha,
				      endPoint_reflectivity,
				      endPoint_sigmaAlpha,
				      glued,
				      glue_property_file,
				      embedment_profile,
				      where_is_embedment_to_be_flush_with_fibre,
				      volumes_with_embedment,
				      constructSensitiveDetector,
				      SearchOverlaps,
				      PropertyTools,
				      DataStorage);

	// reset to the default values
	SetDefaults();

	// return the volumes
	return fibre;
}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
/**
 *  Function to construct replications of a <b> straight </b> fibre with a given <b> start </b> and <b> end point </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre objects
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return vector with pointers to the created G4Fibre objects
 */
vector<G4Fibre *> FibreConstructor::ConstructFibres( G4String fibre_property_file,
						     G4VPhysicalVolume* mother_volume,
						     G4ThreeVector startPoint_relative_to_reference,
						     G4ThreeVector endPoint_relative_to_reference,
						     G4int number_of_replications,
						     G4Transform3D replication_transf,
						     G4VPhysicalVolume* reference_volume,
						     G4Transform3D reference_transform,
						     G4String fibre_name_prefix,
						     G4bool onlyInsideMother,
						     G4bool cutAuntVolumesWithDaughters,
						     G4double startPoint_reflectivity,
						     G4double startPoint_sigmaAlpha,
						     G4double endPoint_reflectivity,
						     G4double endPoint_sigmaAlpha,
						     G4bool glued,
						     G4String glue_property_file,
						     G4String embedment_profile,
						     G4String where_is_embedment_to_be_flush_with_fibre,
						     std::vector<G4VPhysicalVolume *> volumes_with_embedment,
						     G4bool constructSensitiveDetector
						   )
{
	std::vector<G4Fibre *> objectsVector;

	if(number_of_replications < 2)
	{
		G4cout << "##########" << "# WARNING:" << "# Replications of a fibre is to be build, but the number of replications is to low. No fibre will be build." << "##########" << G4endl;
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
		G4String temp_fibre_name_prefix = fibre_name_prefix;
		applyVolumeCounter(temp_fibre_name_prefix);

		// create the volumes
		G4Fibre * fibre = new G4Fibre(fibre_property_file,
					      mother_volume,
					      startPoint_relative_to_reference,
					      endPoint_relative_to_reference,
					      reference_volume,
					      reference_transform,
					      temp_fibre_name_prefix,
					      onlyInsideMother,
					      cutAuntVolumesWithDaughters,
					      startPoint_reflectivity,
					      startPoint_sigmaAlpha,
					      endPoint_reflectivity,
					      endPoint_sigmaAlpha,
					      glued,
					      glue_property_file,
					      embedment_profile,
					      where_is_embedment_to_be_flush_with_fibre,
					      volumes_with_embedment,
					      constructSensitiveDetector,
					      SearchOverlaps,
					      PropertyTools,
					      DataStorage,
					      replication_transf);

		objectsVector.push_back(fibre);
	}

	// reset to the default values
	SetDefaults();

	// return the volumes
	return objectsVector;
}

/**
 *  Function to construct a <b> straight </b> fibre with a given <b> length </b> and <b> transformation </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4Fibre object
 */
G4Fibre * FibreConstructor::ConstructFibre( G4String fibre_property_file,
					    G4VPhysicalVolume* mother_volume,
					    G4double length,
					    G4Transform3D transf_relative_to_reference,
					    G4VPhysicalVolume* reference_volume,
					    G4Transform3D reference_transform,
					    G4String fibre_name_prefix,
					    G4bool onlyInsideMother,
					    G4bool cutAuntVolumesWithDaughters,
					    G4double startPoint_reflectivity,
					    G4double startPoint_sigmaAlpha,
					    G4double endPoint_reflectivity,
					    G4double endPoint_sigmaAlpha,
					    G4bool glued,
					    G4String glue_property_file,
					    G4String embedment_profile,
					    G4String where_is_embedment_to_be_flush_with_fibre,
					    std::vector<G4VPhysicalVolume *> volumes_with_embedment,
					    G4bool constructSensitiveDetector
					  )
{
	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;

	// modify the name if another volume was named like this before
	applyVolumeCounter(fibre_name_prefix);

	// create the volumes
	G4Fibre * fibre = new G4Fibre(fibre_property_file,
				      mother_volume,
				      length,
				      transf_relative_to_reference,
				      reference_volume,
				      reference_transform,
				      fibre_name_prefix,
				      onlyInsideMother,
				      cutAuntVolumesWithDaughters,
				      startPoint_reflectivity,
				      startPoint_sigmaAlpha,
				      endPoint_reflectivity,
				      endPoint_sigmaAlpha,
				      glued,
				      glue_property_file,
				      embedment_profile,
				      where_is_embedment_to_be_flush_with_fibre,
				      volumes_with_embedment,
				      constructSensitiveDetector,
				      SearchOverlaps,
				      PropertyTools,
				      DataStorage);

	// reset to the default values
	SetDefaults();

	// return the volumes
	return fibre;
}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
/**
 *  Function to construct replications of a <b> straight </b> fibre with a given <b> length </b> and <b> transformation </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre objects
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return vector with pointers to the created G4Fibre objects
 */
vector<G4Fibre *> FibreConstructor::ConstructFibres( G4String fibre_property_file,
						     G4VPhysicalVolume* mother_volume,
						     G4double length,
						     G4int number_of_replications,
						     G4Transform3D replication_transf,
						     G4Transform3D transf_relative_to_reference,
						     G4VPhysicalVolume* reference_volume,
						     G4Transform3D reference_transform,
						     G4String fibre_name_prefix,
						     G4bool onlyInsideMother,
						     G4bool cutAuntVolumesWithDaughters,
						     G4double startPoint_reflectivity,
						     G4double startPoint_sigmaAlpha,
						     G4double endPoint_reflectivity,
						     G4double endPoint_sigmaAlpha,
						     G4bool glued,
						     G4String glue_property_file,
						     G4String embedment_profile,
						     G4String where_is_embedment_to_be_flush_with_fibre,
						     std::vector<G4VPhysicalVolume *> volumes_with_embedment,
						     G4bool constructSensitiveDetector
						   )
{
	std::vector<G4Fibre *> objectsVector;

	if(number_of_replications < 2)
	{
		G4cout << "##########" << "# WARNING:" << "# Replications of a fibre is to be build, but the number of replications is to low. No fibre will be build." << "##########" << G4endl;
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
		G4String temp_fibre_name_prefix = fibre_name_prefix;
		applyVolumeCounter(temp_fibre_name_prefix);

		// create the volumes
		G4Fibre * fibre = new G4Fibre(fibre_property_file,
					      mother_volume,
					      length,
					      transf_relative_to_reference,
					      reference_volume,
					      reference_transform,
					      temp_fibre_name_prefix,
					      onlyInsideMother,
					      cutAuntVolumesWithDaughters,
					      startPoint_reflectivity,
					      startPoint_sigmaAlpha,
					      endPoint_reflectivity,
					      endPoint_sigmaAlpha,
					      glued,
					      glue_property_file,
					      embedment_profile,
					      where_is_embedment_to_be_flush_with_fibre,
					      volumes_with_embedment,
					      constructSensitiveDetector,
					      SearchOverlaps,
					      PropertyTools,
					      DataStorage,
					      replication_transf);

		objectsVector.push_back(fibre);
	}

	// reset to the default values
	SetDefaults();

	// return the volumes
	return objectsVector;
}

/**
 *  Function to construct a <b> bent </b> fibre with a given <b> bending axis </b>, <b> bending angle </b>, <b> start </b> and <b> end point </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4Fibre object
 */
G4Fibre * FibreConstructor::ConstructFibre( G4String fibre_property_file,
					    G4VPhysicalVolume* mother_volume,
					    G4ThreeVector startPoint_relative_to_reference,
					    G4ThreeVector endPoint_relative_to_reference,
					    G4double bending_delta_angle,
					    G4ThreeVector bending_axis,
					    G4VPhysicalVolume* reference_volume,
					    G4Transform3D reference_transform,
					    G4String fibre_name_prefix,
					    G4bool onlyInsideMother,
					    G4bool cutAuntVolumesWithDaughters,
					    G4double startPoint_reflectivity,
					    G4double startPoint_sigmaAlpha,
					    G4double endPoint_reflectivity,
					    G4double endPoint_sigmaAlpha,
					    G4bool glued,
					    G4String glue_property_file,
					    G4String embedment_profile,
					    G4String where_is_embedment_to_be_flush_with_fibre,
					    std::vector<G4VPhysicalVolume *> volumes_with_embedment,
					    G4bool constructSensitiveDetector
					  )
{
	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;

	// modify the name if another volume was named like this before
	applyVolumeCounter(fibre_name_prefix);

	// create the volumes
	G4Fibre * fibre = new G4Fibre(fibre_property_file,
				      mother_volume,
				      startPoint_relative_to_reference,
				      endPoint_relative_to_reference,
				      bending_delta_angle,
				      bending_axis,
				      reference_volume,
				      reference_transform,
				      fibre_name_prefix,
				      onlyInsideMother,
				      cutAuntVolumesWithDaughters,
				      startPoint_reflectivity,
				      startPoint_sigmaAlpha,
				      endPoint_reflectivity,
				      endPoint_sigmaAlpha,
				      glued,
				      glue_property_file,
				      embedment_profile,
				      where_is_embedment_to_be_flush_with_fibre,
				      volumes_with_embedment,
				      constructSensitiveDetector,
				      SearchOverlaps,
				      PropertyTools,
				      DataStorage);

	// reset to the default values
	SetDefaults();

	// return the volumes
	return fibre;
}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
/**
 *  Function to construct replications of a <b> bent </b> fibre with a given <b> bending axis </b>, <b> bending angle </b>, <b> start </b> and <b> end point </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre objects
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return vector with pointers to the created G4Fibre objects
 */
vector<G4Fibre *> FibreConstructor::ConstructFibres( G4String fibre_property_file,
						     G4VPhysicalVolume* mother_volume,
						     G4ThreeVector startPoint_relative_to_reference,
						     G4ThreeVector endPoint_relative_to_reference,
						     G4double bending_delta_angle,
						     G4ThreeVector bending_axis,
						     G4int number_of_replications,
						     G4Transform3D replication_transf,
						     G4VPhysicalVolume* reference_volume,
						     G4Transform3D reference_transform,
						     G4String fibre_name_prefix,
						     G4bool onlyInsideMother,
						     G4bool cutAuntVolumesWithDaughters,
						     G4double startPoint_reflectivity,
						     G4double startPoint_sigmaAlpha,
						     G4double endPoint_reflectivity,
						     G4double endPoint_sigmaAlpha,
						     G4bool glued,
						     G4String glue_property_file,
						     G4String embedment_profile,
						     G4String where_is_embedment_to_be_flush_with_fibre,
						     std::vector<G4VPhysicalVolume *> volumes_with_embedment,
						     G4bool constructSensitiveDetector
						   )
{
	std::vector<G4Fibre *> objectsVector;

	if(number_of_replications < 2)
	{
		G4cout << "##########" << "# WARNING:" << "# Replications of a fibre is to be build, but the number of replications is to low. No fibre will be build." << "##########" << G4endl;
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
		G4String temp_fibre_name_prefix = fibre_name_prefix;
		applyVolumeCounter(temp_fibre_name_prefix);

		// create the volumes
		G4Fibre * fibre = new G4Fibre(fibre_property_file,
					      mother_volume,
					      startPoint_relative_to_reference,
					      endPoint_relative_to_reference,
					      bending_delta_angle,
					      bending_axis,
					      reference_volume,
					      reference_transform,
					      temp_fibre_name_prefix,
					      onlyInsideMother,
					      cutAuntVolumesWithDaughters,
					      startPoint_reflectivity,
					      startPoint_sigmaAlpha,
					      endPoint_reflectivity,
					      endPoint_sigmaAlpha,
					      glued,
					      glue_property_file,
					      embedment_profile,
					      where_is_embedment_to_be_flush_with_fibre,
					      volumes_with_embedment,
					      constructSensitiveDetector,
					      SearchOverlaps,
					      PropertyTools,
					      DataStorage,
					      replication_transf);

		objectsVector.push_back(fibre);
	}

	// reset to the default values
	SetDefaults();

	// return the volumes
	return objectsVector;
}

/**
 *  Function to construct a <b> bent </b> fibre with a given <b> circular centre </b>, <b> bending axis </b>, <b> bending radius </b>, <b> start </b> and <b> end angle </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre object
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return pointer to the created G4Fibre object
 */
G4Fibre * FibreConstructor::ConstructFibre( G4String fibre_property_file,
					    G4VPhysicalVolume* mother_volume,
					    G4double bending_start_angle,
					    G4double bending_end_angle,
					    G4ThreeVector bending_circular_centre,
					    G4ThreeVector bending_axis,
					    G4double bending_radius,
					    G4VPhysicalVolume* reference_volume,
					    G4Transform3D reference_transform,
					    G4String fibre_name_prefix,
					    G4bool onlyInsideMother,
					    G4bool cutAuntVolumesWithDaughters,
					    G4double startPoint_reflectivity,
					    G4double startPoint_sigmaAlpha,
					    G4double endPoint_reflectivity,
					    G4double endPoint_sigmaAlpha,
					    G4bool glued,
					    G4String glue_property_file,
					    G4String embedment_profile,
					    G4String where_is_embedment_to_be_flush_with_fibre,
					    std::vector<G4VPhysicalVolume *> volumes_with_embedment,
					    G4bool constructSensitiveDetector
					  )
{
	// the "optional" parameters
	if(!reference_volume) reference_volume = mother_volume;

	// modify the name if another volume was named like this before
	applyVolumeCounter(fibre_name_prefix);

	// create the volumes
	G4Fibre * fibre = new G4Fibre(fibre_property_file,
				      mother_volume,
				      bending_start_angle,
				      bending_end_angle,
				      bending_circular_centre,
				      bending_axis,
				      bending_radius,
				      reference_volume,
				      reference_transform,
				      fibre_name_prefix,
				      onlyInsideMother,
				      cutAuntVolumesWithDaughters,
				      startPoint_reflectivity,
				      startPoint_sigmaAlpha,
				      endPoint_reflectivity,
				      endPoint_sigmaAlpha,
				      glued,
				      glue_property_file,
				      embedment_profile,
				      where_is_embedment_to_be_flush_with_fibre,
				      volumes_with_embedment,
				      constructSensitiveDetector,
				      SearchOverlaps,
				      PropertyTools,
				      DataStorage);

	// reset to the default values
	SetDefaults();

	// return the volumes
	return fibre;
}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
/**
 *  Function to construct replications of a <b> bent </b> fibre with a given <b> circular centre </b>, <b> bending axis </b>, <b> bending radius </b>, <b> start </b> and <b> end angle </b>:
 *  - runs applyVolumeCounter() on \em fibre_name_prefix
 *  - creates a new G4Fibre objects
 *  - resets the class variables to default values (using SetDefaults())
 *
 *  @return vector with pointers to the created G4Fibre objects
 */
vector<G4Fibre *> FibreConstructor::ConstructFibres( G4String fibre_property_file,
						     G4VPhysicalVolume* mother_volume,
						     G4double bending_start_angle,
						     G4double bending_end_angle,
						     G4ThreeVector bending_circular_centre,
						     G4ThreeVector bending_axis,
						     G4double bending_radius,
						     G4int number_of_replications,
						     G4Transform3D replication_transf,
						     G4VPhysicalVolume* reference_volume,
						     G4Transform3D reference_transform,
						     G4String fibre_name_prefix,
						     G4bool onlyInsideMother,
						     G4bool cutAuntVolumesWithDaughters,
						     G4double startPoint_reflectivity,
						     G4double startPoint_sigmaAlpha,
						     G4double endPoint_reflectivity,
						     G4double endPoint_sigmaAlpha,
						     G4bool glued,
						     G4String glue_property_file,
						     G4String embedment_profile,
						     G4String where_is_embedment_to_be_flush_with_fibre,
						     std::vector<G4VPhysicalVolume *> volumes_with_embedment,
						     G4bool constructSensitiveDetector
						   )
{
	std::vector<G4Fibre *> objectsVector;

	if(number_of_replications < 2)
	{
		G4cout << "##########" << "# WARNING:" << "# Replications of a fibre is to be build, but the number of replications is to low. No fibre will be build." << "##########" << G4endl;
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
		G4String temp_fibre_name_prefix = fibre_name_prefix;
		applyVolumeCounter(temp_fibre_name_prefix);

		// create the volumes
		G4Fibre * fibre = new G4Fibre(fibre_property_file,
					      mother_volume,
					      bending_start_angle,
					      bending_end_angle,
					      bending_circular_centre,
					      bending_axis,
					      bending_radius,
					      reference_volume,
					      reference_transform,
					      temp_fibre_name_prefix,
					      onlyInsideMother,
					      cutAuntVolumesWithDaughters,
					      startPoint_reflectivity,
					      startPoint_sigmaAlpha,
					      endPoint_reflectivity,
					      endPoint_sigmaAlpha,
					      glued,
					      glue_property_file,
					      embedment_profile,
					      where_is_embedment_to_be_flush_with_fibre,
					      volumes_with_embedment,
					      constructSensitiveDetector,
					      SearchOverlaps,
					      PropertyTools,
					      DataStorage,
					      replication_transf);

		objectsVector.push_back(fibre);
	}

	// reset to the default values
	SetDefaults();

	// return the volumes
	return objectsVector;
}



/**
 *  Function for avoiding multiple volumes with the same name:
 *  - a string variable containing the name prefix has to be passed to this function
 *  - the name prefix will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre
 *  - the extended name prefix will be passed to the original string variable
 */
void FibreConstructor::applyVolumeCounter(G4String & physical_volume_name)
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
void FibreConstructor::SetDefaults()
{
	/** - G4bool, if a sensitive detector is to be created \code ConstructSensitiveDetector = false; \endcode */
	ConstructSensitiveDetector = false;

	/** - pointer to the reference volume (physical volume relative to which the G4Fibre will be orientated) \code ReferenceVolume_physical = 0; \endcode */
	ReferenceVolume_physical = 0;

	/** - transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) \code ReferenceVolume_transformation = G4Transform3D(); \endcode */
	ReferenceVolume_transformation = G4Transform3D();

	/** - G4Fibre%'s transformation (relative to the reference volume) \code FibreTransformation_rel = G4Transform3D(); \endcode */
	FibreTransformation_rel = G4Transform3D();

	/** - prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) \code FibreNamePrefix = "fibre" \endcode */
	FibreNamePrefix = "fibre";

	/** - G4bool, if the G4Fibre is to be placed only into its mother volume, i.e. all parts outside the mother volume will be cut away \code OnlyInsideMother = false; \endcode */
	OnlyInsideMother = false;

	/** - G4bool, if the G4Fibre should NOT be placed into aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves \code CutAuntVolumesWithDaughters = false; \endcode */
	CutAuntVolumesWithDaughters = false;

	/** - reflectivity of the G4Fibre%'s beginning (start point) reflective layer \code Reflectivity_startPoint = NAN; \endcode */
	Reflectivity_startPoint = NAN;

	/** - roughness of the G4Fibre%'s beginning (start point) (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) \code Roughness_startPoint = 0. * deg; \endcode */
	Roughness_startPoint = 0. * CLHEP::deg;

	/** - reflectivity of the G4Fibre%'s end reflective layer \code Reflectivity_endPoint = NAN; \endcode */
	Reflectivity_endPoint = NAN;

	/** - roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) \code Roughness_endPoint = 0. * deg; \endcode */
	Roughness_endPoint = 0. * CLHEP::deg;

	/** - G4bool, if the G4Fibre is to be glued to at least one volume ("true" or "false") \code Glued = false; \endcode */
	Glued = false;

	/** - path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) \code OpticalCementPropertyFile = ""; \endcode */
	OpticalCementPropertyFile = "";

	/** - profile of the glue's volume ("round" or "quadratic") \code EmbedmentProfile = ""; \endcode */
	EmbedmentProfile = "";

	/** - the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) \code WhereShouldEmbedmentBeFlushWithFibre = ""; \endcode */
	WhereIsEmbedmentToBeFlushWithFibre = "";

	VolumesWithEmbedment.clear();
}
