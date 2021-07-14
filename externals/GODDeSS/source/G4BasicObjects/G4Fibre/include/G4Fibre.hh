/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef G4FIBRE_H
#define G4FIBRE_H

#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4Material.hh>
#include <G4VSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4OpticalSurface.hh>
#include <G4Colour.hh>
#include <G4SDManager.hh>

#include <FibreSensitiveDetector.hh>
#include <GoddessProperties.hh>
#include <PropertyToolsManager.hh>

#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class generating the materials, volumes, and optical properties needed for a fibre and allows to access them.
class G4Fibre
{
public:

	/// <b> >>> construct a STRAIGHT fibre with a given START and END POINT </b>
	/**
	 *  Constructor to construct a <b> straight </b> fibre with a given <b> start </b> and <b> end point </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
		 G4ThreeVector endPoint_relative_to_reference,			/**< G4Fibre%'s end point relative to the reference volume */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage				/**< pointer to the GODDeSS_DataStorage that is to be used */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, StartPoint(startPoint_relative_to_reference)
	, EndPoint(endPoint_relative_to_reference)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		Length = NAN;
		FibreTransformation_rel = G4Transform3D();
		ReplicationTransformation = G4Transform3D();

		FibreBent = false;
		BendingCircularCentre = G4ThreeVector(NAN, NAN, NAN);
		BendingAxis = G4ThreeVector(NAN, NAN, NAN);
		BendingRadius = NAN;
		BendingDeltaAngle = NAN;
		BendingStartAngle = NAN;
		BendingEndAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("straight");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("straight");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a STRAIGHT fibre with a given START and END POINT </b>
	/**
	 *  Constructor to construct replications of a <b> straight </b> fibre with a given <b> start </b> and <b> end point </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
		 G4ThreeVector endPoint_relative_to_reference,			/**< G4Fibre%'s end point relative to the reference volume */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage,				/**< pointer to the GODDeSS_DataStorage that is to be used */
		 G4Transform3D replication_transf				/**< transformation of each replication, relative to the previous one */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, StartPoint(startPoint_relative_to_reference)
	, EndPoint(endPoint_relative_to_reference)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	, ReplicationTransformation(replication_transf)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		Length = NAN;
		FibreTransformation_rel = G4Transform3D();

		FibreBent = false;
		BendingCircularCentre = G4ThreeVector(NAN, NAN, NAN);
		BendingAxis = G4ThreeVector(NAN, NAN, NAN);
		BendingRadius = NAN;
		BendingDeltaAngle = NAN;
		BendingStartAngle = NAN;
		BendingEndAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("straight");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("straight");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

	/// <b> >>> construct a STRAIGHT fibre with a given LENGTH and TRANSFORMATION </b>
	/**
	 *  Constructor to construct a <b> straight </b> fibre with a given <b> length </b> and <b> transformation </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4double length,						/**< G4Fibre%'s length */
		 G4Transform3D transf_relative_to_reference,			/**< G4Fibre%'s transformation relative to the reference volume */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage				/**< pointer to the GODDeSS_DataStorage that is to be used */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, Length(length)
	, FibreTransformation_rel(transf_relative_to_reference)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		StartPoint = G4ThreeVector(NAN, NAN, NAN);
		EndPoint = G4ThreeVector(NAN, NAN, NAN);
		ReplicationTransformation = G4Transform3D();

		FibreBent = false;
		BendingCircularCentre = G4ThreeVector(NAN, NAN, NAN);
		BendingAxis = G4ThreeVector(NAN, NAN, NAN);
		BendingRadius = NAN;
		BendingDeltaAngle = NAN;
		BendingStartAngle = NAN;
		BendingEndAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("straight");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("straight");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a STRAIGHT fibre with a given LENGTH and TRANSFORMATION </b>
	/**
	 *  Constructor to construct replications of a <b> straight </b> fibre with a given <b> length </b> and <b> transformation </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4double length,						/**< G4Fibre%'s length */
		 G4Transform3D transf_relative_to_reference,			/**< G4Fibre%'s transformation relative to the reference volume */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage,				/**< pointer to the GODDeSS_DataStorage that is to be used */
		 G4Transform3D replication_transf				/**< transformation of each replication, relative to the previous one */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, Length(length)
	, FibreTransformation_rel(transf_relative_to_reference)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	, ReplicationTransformation(replication_transf)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		StartPoint = G4ThreeVector(NAN, NAN, NAN);
		EndPoint = G4ThreeVector(NAN, NAN, NAN);

		FibreBent = false;
		BendingCircularCentre = G4ThreeVector(NAN, NAN, NAN);
		BendingAxis = G4ThreeVector(NAN, NAN, NAN);
		BendingRadius = NAN;
		BendingDeltaAngle = NAN;
		BendingStartAngle = NAN;
		BendingEndAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("straight");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("straight");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

	/// <b> >>> construct a BENT fibre with a given BENDING AXIS, BENDING ANGLE, START and END POINT </b>
	/**
	 *  Constructor to construct a <b> bent </b> fibre with a given <b> bending axis </b>, <b> bending angle </b>, <b> start </b> and <b> end point </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
		 G4ThreeVector endPoint_relative_to_reference,			/**< G4Fibre%'s end point relative to the reference volume */
		 G4double bending_delta_angle,					/**< G4Fibre%'s bending angle */
		 G4ThreeVector bending_axis,					/**< G4Fibre%'s bending axis */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage				/**< pointer to the GODDeSS_DataStorage that is to be used */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, StartPoint(startPoint_relative_to_reference)
	, EndPoint(endPoint_relative_to_reference)
	, BendingAxis(bending_axis.perpPart(EndPoint - StartPoint) / bending_axis.perpPart(EndPoint - StartPoint).mag())
	, BendingDeltaAngle(bending_delta_angle)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		FibreBent = true;
		BendingStartAngle = NAN;
		BendingEndAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		Length = NAN;
		FibreTransformation_rel = G4Transform3D();
		ReplicationTransformation = G4Transform3D();

		// get bending properties
		G4ThreeVector linkingVector = EndPoint - StartPoint;
		BendingCircularCentre = (EndPoint + StartPoint) / 2. - linkingVector.mag() / 2. / tan(BendingDeltaAngle / 2.) * linkingVector.cross(BendingAxis) / linkingVector.mag();
		BendingRadius = linkingVector.mag() / 2. / sin(BendingDeltaAngle / 2.);
		if(CreateReflectiveStartPointVolumes || CreateReflectiveEndPointVolumes) Reflective_bendingDeltaAngle = 2. * asin(Reflective_thickness / 2. / BendingRadius);
		if(!(CreateReflectiveStartPointVolumes && CreateReflectiveEndPointVolumes)) Rough_bendingDeltaAngle = 2. * asin(Rough_thickness / 2. / BendingRadius);

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("bent");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("bent");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a BENT fibre with a given BENDING AXIS, BENDING ANGLE, START and END POINT </b>
	/**
	 *  Constructor to construct replications of a <b> bent </b> fibre with a given <b> bending axis </b>, <b> bending angle </b>, <b> start </b> and <b> end point </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
		 G4ThreeVector endPoint_relative_to_reference,			/**< G4Fibre%'s end point relative to the reference volume */
		 G4double bending_delta_angle,					/**< G4Fibre%'s bending angle */
		 G4ThreeVector bending_axis,					/**< G4Fibre%'s bending axis */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage,				/**< pointer to the GODDeSS_DataStorage that is to be used */
		 G4Transform3D replication_transf				/**< transformation of each replication, relative to the previous one */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, StartPoint(startPoint_relative_to_reference)
	, EndPoint(endPoint_relative_to_reference)
	, BendingAxis(bending_axis.perpPart(EndPoint - StartPoint) / bending_axis.perpPart(EndPoint - StartPoint).mag())
	, BendingDeltaAngle(bending_delta_angle)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	, ReplicationTransformation(replication_transf)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		FibreBent = true;
		BendingStartAngle = NAN;
		BendingEndAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		Length = NAN;
		FibreTransformation_rel = G4Transform3D();

		// get bending properties
		G4ThreeVector linkingVector = EndPoint - StartPoint;
		BendingCircularCentre = (EndPoint + StartPoint) / 2. - linkingVector.mag() / 2. / tan(BendingDeltaAngle / 2.) * linkingVector.cross(BendingAxis) / linkingVector.mag();
		BendingRadius = linkingVector.mag() / 2. / sin(BendingDeltaAngle / 2.);
		if(CreateReflectiveStartPointVolumes || CreateReflectiveEndPointVolumes) Reflective_bendingDeltaAngle = 2. * asin(Reflective_thickness / 2. / BendingRadius);
		if(!(CreateReflectiveStartPointVolumes && CreateReflectiveEndPointVolumes)) Rough_bendingDeltaAngle = 2. * asin(Rough_thickness / 2. / BendingRadius);

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("bent");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("bent");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

	/// <b> >>> construct a BENT fibre with a given CIRCULAR CENTRE, BENDING AXIS, BENDING RADIUS, START and END ANGLE </b>
	/**
	 *  Constructor to construct a <b> bent </b> fibre with a given <b> circular centre </b>, <b> bending axis </b>, <b> bending radius </b>, <b> start </b> and <b> end angle </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4double bending_start_angle,					/**< G4Fibre%'s start angle (in the bending circle) */
		 G4double bending_end_angle,					/**< G4Fibre%'s end angle (in the bending circle) */
		 G4ThreeVector bending_circular_centre,				/**< circular centre of the G4Fibre%'s bending */
		 G4ThreeVector bending_axis,					/**< G4Fibre%'s bending axis */
		 G4double bending_radius,					/**< G4Fibre%'s bending radius */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage				/**< pointer to the GODDeSS_DataStorage that is to be used */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, BendingCircularCentre(bending_circular_centre)
	, BendingAxis(bending_axis / bending_axis.mag())
	, BendingRadius(bending_radius)
	, BendingStartAngle(bending_start_angle)
	, BendingEndAngle(bending_end_angle)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		if(BendingEndAngle < BendingStartAngle)
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# bending_end_angle (" << bending_end_angle << ") < bending_start_angle (" << bending_start_angle << "). No bent fibre can be build." << "##########" << G4endl;
			return;
		}

		// set default values
		SetDefaults();

		FibreBent = true;
		BendingDeltaAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		StartPoint = G4ThreeVector(NAN, NAN, NAN);
		EndPoint = G4ThreeVector(NAN, NAN, NAN);
		Length = NAN;
		FibreTransformation_rel = G4Transform3D();
		ReplicationTransformation = G4Transform3D();

		// get bending properties
		BendingDeltaAngle = BendingEndAngle - BendingStartAngle;
		if(CreateReflectiveStartPointVolumes || CreateReflectiveEndPointVolumes) Reflective_bendingDeltaAngle = 2. * asin(Reflective_thickness / 2. / BendingRadius);
		if(!(CreateReflectiveStartPointVolumes && CreateReflectiveEndPointVolumes)) Rough_bendingDeltaAngle = 2. * asin(Rough_thickness / 2. / BendingRadius);

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("bent");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("bent");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a BENT fibre with a given CIRCULAR CENTRE, BENDING AXIS, BENDING RADIUS, START and END ANGLE </b>
	/**
	 *  Constructor to construct replications of a <b> bent </b> fibre with a given <b> circular centre </b>, <b> bending axis </b>, <b> bending radius </b>, <b> start </b> and <b> end angle </b>:
	 *  - sets class variables to default values
	 *  - loads the property file and determines the fibres type
	 *  - sets the variables for the fibre properties and its placement (InitialiseVariables() \& GenerateTransformation())
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructVolumes())
	 */
	G4Fibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s property */
		 G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
		 G4double bending_start_angle,					/**< G4Fibre%'s start angle (in the bending circle) */
		 G4double bending_end_angle,					/**< G4Fibre%'s end angle (in the bending circle) */
		 G4ThreeVector bending_circular_centre,				/**< circular centre of the G4Fibre%'s bending */
		 G4ThreeVector bending_axis,					/**< G4Fibre%'s bending axis */
		 G4double bending_radius,					/**< G4Fibre%'s bending radius */
		 G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
		 G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
		 G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
		 G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
		 G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
		 G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
		 G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
		 G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
		 G4String glue_property_file,					/**< path to the file containing the glue's property (the glue may consist of any material composition, also air or vacuum) */
		 G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
		 G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
		 std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
		 G4bool constructSensitiveDetector,				/**< a sensitive detector is to be constructed ("true" or "false") */
		 G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		 PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
		 GODDeSS_DataStorage * dataStorage,				/**< pointer to the GODDeSS_DataStorage that is to be used */
		 G4Transform3D replication_transf				/**< transformation of each replication, relative to the previous one */
	       )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, ReferenceVolume_transformation(reference_transformation)
	, BendingCircularCentre(bending_circular_centre)
	, BendingAxis(bending_axis / bending_axis.mag())
	, BendingRadius(bending_radius)
	, BendingStartAngle(bending_start_angle)
	, BendingEndAngle(bending_end_angle)
	, FibreNamePrefix(fibre_name_prefix)
	, OnlyInsideMother(onlyInsideMother)
	, CutAuntVolumesWithDaughters(cutAuntVolumesWithDaughters)
	, Glued(glued)
	, EmbedmentProfile(embedment_profile)
	, WhereIsEmbedmentToBeFlushWithFibre(where_is_embedment_to_be_flush_with_fibre)
	, VolumesWithEmbedment(volumes_with_embedment)
	, Reflectivity_startPoint(startPoint_reflectivity)
	, Roughness_startPoint(startPoint_sigmaAlpha)
	, Reflectivity_endPoint(endPoint_reflectivity)
	, Roughness_endPoint(endPoint_sigmaAlpha)
	, ReplicationTransformation(replication_transf)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		if(BendingEndAngle < BendingStartAngle)
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# bending_end_angle (" << bending_end_angle << ") < bending_start_angle (" << bending_start_angle << "). No bent fibre can be build." << "##########" << G4endl;
			return;
		}

		// set default values
		SetDefaults();

		FibreBent = true;
		BendingDeltaAngle = NAN;
		Reflective_bendingDeltaAngle = NAN;
		Rough_bendingDeltaAngle = NAN;

		StartPoint = G4ThreeVector(NAN, NAN, NAN);
		EndPoint = G4ThreeVector(NAN, NAN, NAN);
		Length = NAN;
		FibreTransformation_rel = G4Transform3D();

		// get bending properties
		BendingDeltaAngle = BendingEndAngle - BendingStartAngle;
		if(CreateReflectiveStartPointVolumes || CreateReflectiveEndPointVolumes) Reflective_bendingDeltaAngle = 2. * asin(Reflective_thickness / 2. / BendingRadius);
		if(!(CreateReflectiveStartPointVolumes && CreateReflectiveEndPointVolumes)) Rough_bendingDeltaAngle = 2. * asin(Rough_thickness / 2. / BendingRadius);

		// the obligatory parameters
		FibreProperties.load(fibre_property_file);

		if(FibreProperties.containsNumber("lightOutput") || FibreProperties.containsNumber("lightOutput_rel"))
		{
			IsScinti = true;
			FibreNamePrefix = FibreNamePrefix + "_(scinti)";
		}
		else if(PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon"))
		{
			IsWLS = true;
			FibreNamePrefix = FibreNamePrefix + "_(WLS)";
		}

		if(Glued && glue_property_file == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no property file is specified. No volume will be build." << G4endl << "##########" << G4endl;
			Glued = false;
		}
		if(Glued && EmbedmentProfile == "")
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no profile file is specified. A quadratic profile will be used." << G4endl << "##########" << G4endl;
			EmbedmentProfile = "quadratic";
		}
		if( Glued && ! (WhereIsEmbedmentToBeFlushWithFibre == "" || WhereIsEmbedmentToBeFlushWithFibre == "S" || WhereIsEmbedmentToBeFlushWithFibre == "E" || WhereIsEmbedmentToBeFlushWithFibre == "SE") )
		{
			G4cerr << "##########" << G4endl << "# WARNING:" << G4endl << "# An optical cement volume is to be build, but no valid string has been given to specify where the embedment is to be flush with the fibre. \"\" will be used." << G4endl << "##########" << G4endl;
			WhereIsEmbedmentToBeFlushWithFibre = "";
		}

		if(Glued) OpticalCementProperties.load(glue_property_file);

		// calculate variables and create materials:
		InitialiseVariables();
		GenerateTransformation("bent");
		DefineMaterials();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);

		// create the volumes
		ConstructVolumes("bent");

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

	/**
	 *  Destructor (empty)
	 */
	~G4Fibre()
	{
	}



//   #######################   ///
//   #  Getter functions:  #   ///
//   #######################   ///
// basic name of the object:
	/** @return the G4Fibre%'s name prefix */
	G4String GetFibreName()
	{ return FibreNamePrefix; }

// dimensions of the object:
	/**
	 *  @return <b> straight fibres: </b> the fibre's total length
	 *  @return <b> bent fibres: </b> NAN
	 */
	G4double GetFibreLength()
	{ return Length; }

	/**
	 *  @return <b> round fibres: </b> the fibre's total radius
	 *  @return <b> quadratic fibres: </b> NAN
	 */
	G4double GetFibreRadius()
	{ return FibreRadius; }

	/**
	 *  @return <b> round fibres: </b> NAN
	 *  @return <b> quadratic fibres: </b> the fibre's total edge length
	 */
	G4double GetFibreEdgeLength()
	{ return FibreEdgeLength; }

// for bent fibres:
// NOTE: the fibre is originally generated circularly around the z-axis, 0° is parallel to the x-axis and the angle is counted mathematically right-handed around the z-axis
	/**
	 *  @return <b> straight fibres: </b> G4ThreeVector(NAN, NAN, NAN)
	 *  @return <b> bent fibres: </b> the position of the fibre's bending circular centre
	 */
	G4ThreeVector GetCircularCentreOfBentFibre()
	{ return BendingCircularCentre; }

	/**
	 *  @return <b> straight fibres: </b> G4ThreeVector(NAN, NAN, NAN)
	 *  @return <b> bent fibres: </b> the fibre's bending axis
	 */
	G4ThreeVector GetBendingAxisOfBentFibre()
	{ return BendingAxis; }

	/**
	 *  @return <b> straight fibres: </b> NAN
	 *  @return <b> bent fibres: </b> the fibre's bending radius
	 */
	G4double GetBendingRadiusOfBentFibre()
	{ return BendingRadius; }

	/**
	 *  @return <b> straight fibres: </b> NAN
	 *  @return <b> bent fibres: </b> the fibre's bending angle
	 */
	G4double GetBendingDeltaAngleOfBentFibre()
	{ return BendingDeltaAngle; }

	/**
	 *  @return <b> straight fibres: </b> NAN
	 *  @return <b> bent fibres: </b> the fibre's bending start angle
	 */
	G4double GetBendingStartAngleOfBentFibre()
	{ return BendingStartAngle; }

	/**
	 *  @return <b> straight fibres: </b> NAN
	 *  @return <b> bent fibres: </b> the fibre's bending end angle
	 */
	G4double GetBendingEndAngleOfBentFibre()
	{ return BendingEndAngle; }

// relative diameter with respect to the full fibre diameter
// NOTE: relative diameter = relative radius = relative edge length
	/**
	 *  @return <b> if a coating exists: </b> the inner relative diameter of the fibre's coating (compaired to the full diameter)
	 *  @return <b> else: </b> NAN
	 */
	G4double GetInnerRelativeDiameterOfCoating()
	{ return RelativeFibreRadius_CoatingMin; }

	/**
	 *  @return <b> if a coating exists: </b> the outer relative diameter of the fibre's coating (compaired to the full diameter)
	 *  @return <b> else: </b> NAN
	 */
	G4double GetOuterRelativeDiameterOfCoating()
	{ return RelativeFibreRadius_CoatingMax; }

	/**
	 *  @return <b> if a 3. cladding exists: </b> the relative diameter of the fibre's 3. cladding (compaired to the full diameter)
	 *  @return <b> else: </b> NAN
	 */
	G4double GetRelativeDiameterOfThirdCladding()
	{ return RelativeFibreRadius_Cladding3; }

	/**
	 *  @return <b> if a 2. cladding exists: </b> the relative diameter of the fibre's 2. cladding (compaired to the full diameter)
	 *  @return <b> else: </b> NAN
	 */
	G4double GetRelativeDiameterOfSecondCladding()
	{ return RelativeFibreRadius_Cladding2; }

	/**
	 *  @return <b> if a 1. cladding exists: </b> the relative diameter of the fibre's 1. cladding (compaired to the full diameter)
	 *  @return <b> else: </b> NAN
	 */
	G4double GetRelativeDiameterOfFirstCladding()
	{ return RelativeFibreRadius_Cladding1; }

	/** @return the relative diameter of the fibre's core (compaired to the full diameter) */
	G4double GetRelativeDiameterOfFibreCore()
	{ return RelativeFibreRadius_Core; }

// position and orientation of the object:
	/** @return the transformation of the G4Fibre inside the mother volume */
	G4Transform3D GetFibreInsideMotherTransformation()
	{
		return FibreTransformation_insideMother;
	}

	/** @return the transformation of the G4Fibre outside the mother volume in the grandmother volume */
	G4Transform3D GetFibreInsideGrandMotherTransformation()
	{ return FibreTransformation_outsideMother; }

// mother volume of the object:
	/** @return pointer to the G4Fibre%'s mother volume */
	G4VPhysicalVolume * GetMotherVolume_physicalVolume()
	{
		return MotherVolume_physical;
	}

// physical volumes of the object:
	/**
	 *  @return <b> if any G4Fibre%'s physical volume exists inside the mother volume: </b> pointer to the outermost G4Fibre%'s physical volume (in the following order: optical cement, 3. cladding, 2. cladding, 1. cladding, core), existing inside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetOutermostVolumeInsideMother_physicalVolume()
	{
		if(Glued) return FibreEmbedment_physical[0];
		else if(Cladding3Exists) return Cladding3_physical[0];
		else if(Cladding2Exists) return Cladding2_physical[0];
		else if(Cladding1Exists) return Cladding1_physical[0];
		else return FibreCore_physical[0];
	}

	/**
	 *  @return <b> if any G4Fibre%'s physical volume exists outside the mother volume: </b> pointer to the outermost G4Fibre%'s physical volume (in the following order: 3. cladding, 2. cladding, 1. cladding, core), existing outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetOutermostVolumeOutsideMother_physicalVolume()
	{
		if(Cladding3Exists) return Cladding3_physical[Cladding3_physical.size() - 1];
		else if(Cladding2Exists) return Cladding2_physical[Cladding2_physical.size() - 1];
		else if(Cladding1Exists) return Cladding1_physical[Cladding1_physical.size() - 1];
		else return FibreCore_physical[FibreCore_physical.size() - 1];
	}

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s optical cement exists: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetOpticalCement_physicalVolume()
	{ return FibreEmbedment_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s coating exists intside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetCoatingInsideMother_physicalVolume()
	{ return Coating_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s coating exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetCoatingOutsideMother_physicalVolume()
	{ return Coating_physical[Coating_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 3. cladding exists intside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetThirdCladdingInsideMother_physicalVolume()
	{ return Cladding3_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 3. cladding exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetThirdCladdingOutsideMother_physicalVolume()
	{ return Cladding3_physical[Cladding3_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 2. cladding exists intside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetSecondCladdingInsideMother_physicalVolume()
	{ return Cladding2_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 2. cladding exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetSecondCladdingOutsideMother_physicalVolume()
	{ return Cladding2_physical[Cladding2_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 1. cladding exists intside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFirstCladdingInsideMother_physicalVolume()
	{ return Cladding1_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 1. cladding exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFirstCladdingOutsideMother_physicalVolume()
	{ return Cladding1_physical[Cladding1_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s core logical intside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFibreCoreInsideMother_physicalVolume()
	{ return FibreCore_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s core exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFibreCoreOutsideMother_physicalVolume()
	{ return FibreCore_physical[FibreCore_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s reflective start point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetReflectiveStartPointInsideMother_physicalVolume()
	{ return ReflectiveStartPoint_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s reflective start point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetReflectiveStartPointOutsideMother_physicalVolume()
	{ return ReflectiveStartPoint_physical[ReflectiveStartPoint_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s reflective end point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetReflectiveEndPointInsideMother_physicalVolume()
	{ return ReflectiveStartPoint_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s reflective end point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetReflectiveEndPointOutsideMother_physicalVolume()
	{ return ReflectiveStartPoint_physical[ReflectiveStartPoint_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 3. cladding roughened start point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetThirdCladdingRoughenedStartPointInsideMother_physicalVolume()
	{ return RoughenedStartPoint_cladding3_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 3. cladding roughened start point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetThirdCladdingRoughenedStartPointOutsideMother_physicalVolume()
	{ return RoughenedStartPoint_cladding3_physical[RoughenedStartPoint_cladding3_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 2. cladding roughened start point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetSecondCladdingRoughenedStartPointInsideMother_physicalVolume()
	{ return RoughenedStartPoint_cladding2_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 2. cladding roughened start point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetSecondCladdingRoughenedStartPointOutsideMother_physicalVolume()
	{ return RoughenedStartPoint_cladding2_physical[RoughenedStartPoint_cladding2_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 1. cladding roughened start point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFirstCladdingRoughenedStartPointInsideMother_physicalVolume()
	{ return RoughenedStartPoint_cladding1_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 1. cladding roughened start point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFirstCladdingRoughenedStartPointOutsideMother_physicalVolume()
	{ return RoughenedStartPoint_cladding1_physical[RoughenedStartPoint_cladding1_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s core roughened start point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFibreCoreRoughenedStartPointInsideMother_physicalVolume()
	{ return RoughenedStartPoint_fibreCore_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s core roughened start point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFibreCoreRoughenedStartPointOutsideMother_physicalVolume()
	{ return RoughenedStartPoint_fibreCore_physical[RoughenedStartPoint_fibreCore_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 3. cladding roughened end point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetThirdCladdingRoughenedEndPointInsideMother_physicalVolume()
	{ return RoughenedEndPoint_cladding3_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 3. cladding roughened end point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetThirdCladdingRoughenedEndPointOutsideMother_physicalVolume()
	{ return RoughenedEndPoint_cladding3_physical[RoughenedEndPoint_cladding3_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 2. cladding roughened end point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetSecondCladdingRoughenedEndPointInsideMother_physicalVolume()
	{ return RoughenedEndPoint_cladding2_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 2. cladding roughened end point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetSecondCladdingRoughenedEndPointOutsideMother_physicalVolume()
	{ return RoughenedEndPoint_cladding2_physical[RoughenedEndPoint_cladding2_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 1. cladding roughened end point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFirstCladdingRoughenedEndPointInsideMother_physicalVolume()
	{ return RoughenedEndPoint_cladding1_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s 1. cladding roughened end point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFirstCladdingRoughenedEndPointOutsideMother_physicalVolume()
	{ return RoughenedEndPoint_cladding1_physical[RoughenedEndPoint_cladding1_physical.size() - 1]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s core roughened end point exists inside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFibreCoreRoughenedEndPointInsideMother_physicalVolume()
	{ return RoughenedEndPoint_fibreCore_physical[0]; }

	/**
	 *  @return <b> if a physical volume for a G4Fibre%'s core roughened end point exists outside the mother volume: </b> pointer to this physical volume
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetFibreCoreRoughenedEndPointOutsideMother_physicalVolume()
	{ return RoughenedEndPoint_fibreCore_physical[RoughenedEndPoint_fibreCore_physical.size() - 1]; }

// logical volumes of the object:
	/**
	 *  @return <b> if any G4Fibre%'s logical volume exists inside the mother volume: </b> pointer to the outermost G4Fibre%'s logical volume (in the following order: optical cement, 3. cladding, 2. cladding, 1. cladding, core), existing inside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetOutermostVolumeInsideMother_logicalVolume()
	{
		if(GetOutermostVolumeInsideMother_physicalVolume()) return GetOutermostVolumeInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if any G4Fibre%'s logical volume exists outside the mother volume: </b> pointer to the outermost G4Fibre%'s logical volume (in the following order: 3. cladding, 2. cladding, 1. cladding, core), existing outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetOutermostVolumeOutsideMother_logicalVolume()
	{
		if(GetOutermostVolumeOutsideMother_physicalVolume()) return GetOutermostVolumeOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s optical cement exists: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetOpticalCement_logicalVolume()
	{
		if(GetOpticalCement_physicalVolume()) return GetOpticalCement_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 3. cladding exists intside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetThirdCladdingInsideMother_logicalVolume()
	{
		if(GetThirdCladdingInsideMother_physicalVolume()) return GetThirdCladdingInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 3. cladding exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetThirdCladdingOutsideMother_logicalVolume()
	{
		if(GetThirdCladdingOutsideMother_physicalVolume()) return GetThirdCladdingOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 2. cladding exists intside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetSecondCladdingInsideMother_logicalVolume()
	{
		if(GetSecondCladdingInsideMother_physicalVolume()) return GetSecondCladdingInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 2. cladding exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetSecondCladdingOutsideMother_logicalVolume()
	{
		if(GetSecondCladdingOutsideMother_physicalVolume()) return GetSecondCladdingOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 1. cladding exists intside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFirstCladdingInsideMother_logicalVolume()
	{
		if(GetFirstCladdingInsideMother_physicalVolume()) return GetFirstCladdingInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 1. cladding exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFirstCladdingOutsideMother_logicalVolume()
	{
		if(GetFirstCladdingOutsideMother_physicalVolume()) return GetFirstCladdingOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s core logical intside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFibreCoreInsideMother_logicalVolume()
	{
		if(GetFibreCoreInsideMother_physicalVolume()) return GetFibreCoreInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s core exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFibreCoreOutsideMother_logicalVolume()
	{
		if(GetFibreCoreOutsideMother_physicalVolume()) return GetFibreCoreOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s reflective start point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetReflectiveStartPointInsideMother_logicalVolume()
	{
		if(GetReflectiveStartPointInsideMother_physicalVolume()) return GetReflectiveStartPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s reflective start point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetReflectiveStartPointOutsideMother_logicalVolume()
	{
		if(GetReflectiveStartPointOutsideMother_physicalVolume()) return GetReflectiveStartPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s reflective end point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetReflectiveEndPointInsideMother_logicalVolume()
	{
		if(GetReflectiveEndPointInsideMother_physicalVolume()) return GetReflectiveEndPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s reflective end point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetReflectiveEndPointOutsideMother_logicalVolume()
	{
		if(GetReflectiveEndPointOutsideMother_physicalVolume()) return GetReflectiveEndPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 3. cladding roughened start point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetThirdCladdingRoughenedStartPointInsideMother_logicalVolume()
	{
		if(GetThirdCladdingRoughenedStartPointInsideMother_physicalVolume()) return GetThirdCladdingRoughenedStartPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 3. cladding roughened start point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetThirdCladdingRoughenedStartPointOutsideMother_logicalVolume()
	{
		if(GetThirdCladdingRoughenedStartPointOutsideMother_physicalVolume()) return GetThirdCladdingRoughenedStartPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 2. cladding roughened start point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetSecondCladdingRoughenedStartPointInsideMother_logicalVolume()
	{
		if(GetSecondCladdingRoughenedStartPointInsideMother_physicalVolume()) return GetSecondCladdingRoughenedStartPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 2. cladding roughened start point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetSecondCladdingRoughenedStartPointOutsideMother_logicalVolume()
	{
		if(GetSecondCladdingRoughenedStartPointOutsideMother_physicalVolume()) return GetSecondCladdingRoughenedStartPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 1. cladding roughened start point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFirstCladdingRoughenedStartPointInsideMother_logicalVolume()
	{
		if(GetFirstCladdingRoughenedStartPointInsideMother_physicalVolume()) return GetFirstCladdingRoughenedStartPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 1. cladding roughened start point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFirstCladdingRoughenedStartPointOutsideMother_logicalVolume()
	{
		if(GetFirstCladdingRoughenedStartPointOutsideMother_physicalVolume()) return GetFirstCladdingRoughenedStartPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s core roughened start point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFibreCoreRoughenedStartPointInsideMother_logicalVolume()
	{
		if(GetFibreCoreRoughenedStartPointInsideMother_physicalVolume()) return GetFibreCoreRoughenedStartPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s core roughened start point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFibreCoreRoughenedStartPointOutsideMother_logicalVolume()
	{
		if(GetFibreCoreRoughenedStartPointOutsideMother_physicalVolume()) return GetFibreCoreRoughenedStartPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 3. cladding roughened end point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetThirdCladdingRoughenedEndPointInsideMother_logicalVolume()
	{
		if(GetThirdCladdingRoughenedEndPointInsideMother_physicalVolume()) return GetThirdCladdingRoughenedEndPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 3. cladding roughened end point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetThirdCladdingRoughenedEndPointOutsideMother_logicalVolume()
	{
		if(GetThirdCladdingRoughenedEndPointOutsideMother_physicalVolume()) return GetThirdCladdingRoughenedEndPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 2. cladding roughened end point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetSecondCladdingRoughenedEndPointInsideMother_logicalVolume()
	{
		if(GetSecondCladdingRoughenedEndPointInsideMother_physicalVolume()) return GetSecondCladdingRoughenedEndPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 2. cladding roughened end point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetSecondCladdingRoughenedEndPointOutsideMother_logicalVolume()
	{
		if(GetSecondCladdingRoughenedEndPointOutsideMother_physicalVolume()) return GetSecondCladdingRoughenedEndPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 1. cladding roughened end point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFirstCladdingRoughenedEndPointInsideMother_logicalVolume()
	{
		if(GetFirstCladdingRoughenedEndPointInsideMother_physicalVolume()) return GetFirstCladdingRoughenedEndPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s 1. cladding roughened end point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFirstCladdingRoughenedEndPointOutsideMother_logicalVolume()
	{
		if(GetFirstCladdingRoughenedEndPointOutsideMother_physicalVolume()) return GetFirstCladdingRoughenedEndPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s core roughened end point exists inside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFibreCoreRoughenedEndPointInsideMother_logicalVolume()
	{
		if(GetFibreCoreRoughenedEndPointInsideMother_physicalVolume()) return GetFibreCoreRoughenedEndPointInsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/**
	 *  @return <b> if a logical volume for a G4Fibre%'s core roughened end point exists outside the mother volume: </b> pointer to this logical volume
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetFibreCoreRoughenedEndPointOutsideMother_logicalVolume()
	{
		if(GetFibreCoreRoughenedEndPointOutsideMother_physicalVolume()) return GetFibreCoreRoughenedEndPointOutsideMother_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

// solid volumes of the object:
	/**
	 *  @return <b> if any G4Fibre%'s solid volume exists inside the mother volume: </b> pointer to the outermost G4Fibre%'s solid volume (in the following order: optical cement, 3. cladding, 2. cladding, 1. cladding, core), existing inside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetOutermostVolumeInsideMother_solidVolume()
	{
		if(GetOutermostVolumeInsideMother_logicalVolume()) return GetOutermostVolumeInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if any G4Fibre%'s solid volume exists outside the mother volume: </b> pointer to the outermost G4Fibre%'s solid volume (in the following order: 3. cladding, 2. cladding, 1. cladding, core), existing outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetOutermostVolumeOutsideMother_solidVolume()
	{
		if(GetOutermostVolumeOutsideMother_logicalVolume()) return GetOutermostVolumeOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s optical cement exists: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetOpticalCement_solidVolume()
	{
		if(GetOpticalCement_logicalVolume()) return GetOpticalCement_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding exists intside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingInsideMother_solidVolume()
	{
		if(GetThirdCladdingInsideMother_logicalVolume()) return GetThirdCladdingInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingOutsideMother_solidVolume()
	{
		if(GetThirdCladdingOutsideMother_logicalVolume()) return GetThirdCladdingOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladding_basicSolidVolume()
	{
		if(GetThirdCladdingInsideMother_solidVolume()) return GetThirdCladdingInsideMother_solidVolume()->GetConstituentSolid(0);
		else if(GetThirdCladdingOutsideMother_solidVolume()) return GetThirdCladdingOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding exists intside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingInsideMother_solidVolume()
	{
		if(GetSecondCladdingInsideMother_logicalVolume()) return GetSecondCladdingInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingOutsideMother_solidVolume()
	{
		if(GetSecondCladdingOutsideMother_logicalVolume()) return GetSecondCladdingOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladding_basicSolidVolume()
	{
		if(GetSecondCladdingInsideMother_solidVolume()) return GetSecondCladdingInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetSecondCladdingOutsideMother_solidVolume()) return GetSecondCladdingOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding exists intside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingInsideMother_solidVolume()
	{
		if(GetFirstCladdingInsideMother_logicalVolume()) return GetFirstCladdingInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingOutsideMother_solidVolume()
	{
		if(GetFirstCladdingOutsideMother_logicalVolume()) return GetFirstCladdingOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladding_basicSolidVolume()
	{
		if(GetFirstCladdingInsideMother_solidVolume()) return GetFirstCladdingInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetFirstCladdingOutsideMother_solidVolume()) return GetFirstCladdingOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core logical intside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreInsideMother_solidVolume()
	{
		if(GetFibreCoreInsideMother_logicalVolume()) return GetFibreCoreInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreOutsideMother_solidVolume()
	{
		if(GetFibreCoreOutsideMother_logicalVolume()) return GetFibreCoreOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCore_basicSolidVolume()
	{
		if(GetFibreCoreInsideMother_solidVolume()) return GetFibreCoreInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetFibreCoreOutsideMother_solidVolume()) return GetFibreCoreOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s reflective start point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetReflectiveStartPointInsideMother_solidVolume()
	{
		if(GetReflectiveStartPointInsideMother_logicalVolume()) return GetReflectiveStartPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s reflective start point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetReflectiveStartPointOutsideMother_solidVolume()
	{
		if(GetReflectiveStartPointOutsideMother_logicalVolume()) return GetReflectiveStartPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s reflective start point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetReflectiveStartPoint_basicSolidVolume()
	{
		if(GetReflectiveStartPointInsideMother_solidVolume()) return GetReflectiveStartPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetReflectiveStartPointOutsideMother_solidVolume()) return GetReflectiveStartPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s reflective end point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetReflectiveEndPointInsideMother_solidVolume()
	{
		if(GetReflectiveEndPointInsideMother_logicalVolume()) return GetReflectiveEndPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s reflective end point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetReflectiveEndPointOutsideMother_solidVolume()
	{
		if(GetReflectiveEndPointOutsideMother_logicalVolume()) return GetReflectiveEndPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s reflective end point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetReflectiveEndPoint_basicSolidVolume()
	{
		if(GetReflectiveEndPointInsideMother_solidVolume()) return GetReflectiveEndPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetReflectiveEndPointOutsideMother_solidVolume()) return GetReflectiveEndPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding roughened start point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingRoughenedStartPointInsideMother_solidVolume()
	{
		if(GetThirdCladdingRoughenedStartPointInsideMother_logicalVolume()) return GetThirdCladdingRoughenedStartPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding roughened start point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingRoughenedStartPointOutsideMother_solidVolume()
	{
		if(GetThirdCladdingRoughenedStartPointOutsideMother_logicalVolume()) return GetThirdCladdingRoughenedStartPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding roughened start point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingRoughenedStartPoint_basicSolidVolume()
	{
		if(GetThirdCladdingRoughenedStartPointInsideMother_solidVolume()) return GetThirdCladdingRoughenedStartPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetThirdCladdingRoughenedStartPointOutsideMother_solidVolume()) return GetThirdCladdingRoughenedStartPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding roughened start point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingRoughenedStartPointInsideMother_solidVolume()
	{
		if(GetSecondCladdingRoughenedStartPointInsideMother_logicalVolume()) return GetSecondCladdingRoughenedStartPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding roughened start point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingRoughenedStartPointOutsideMother_solidVolume()
	{
		if(GetSecondCladdingRoughenedStartPointOutsideMother_logicalVolume()) return GetSecondCladdingRoughenedStartPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding roughened start point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingRoughenedStartPoint_basicSolidVolume()
	{
		if(GetSecondCladdingRoughenedStartPointInsideMother_solidVolume()) return GetSecondCladdingRoughenedStartPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetSecondCladdingRoughenedStartPointOutsideMother_solidVolume()) return GetSecondCladdingRoughenedStartPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding roughened start point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingRoughenedStartPointInsideMother_solidVolume()
	{
		if(GetFirstCladdingRoughenedStartPointInsideMother_logicalVolume()) return GetFirstCladdingRoughenedStartPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding roughened start point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingRoughenedStartPointOutsideMother_solidVolume()
	{
		if(GetFirstCladdingRoughenedStartPointOutsideMother_logicalVolume()) return GetFirstCladdingRoughenedStartPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding roughened start point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingRoughenedStartPoint_basicSolidVolume()
	{
		if(GetFirstCladdingRoughenedStartPointInsideMother_solidVolume()) return GetFirstCladdingRoughenedStartPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetFirstCladdingRoughenedStartPointOutsideMother_solidVolume()) return GetFirstCladdingRoughenedStartPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core roughened start point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreRoughenedStartPointInsideMother_solidVolume()
	{
		if(GetFibreCoreRoughenedStartPointInsideMother_logicalVolume()) return GetFibreCoreRoughenedStartPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core roughened start point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreRoughenedStartPointOutsideMother_solidVolume()
	{
		if(GetFibreCoreRoughenedStartPointOutsideMother_logicalVolume()) return GetFibreCoreRoughenedStartPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core roughened start point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreRoughenedStartPoint_basicSolidVolume()
	{
		if(GetFibreCoreRoughenedStartPointInsideMother_solidVolume()) return GetFibreCoreRoughenedStartPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetFibreCoreRoughenedStartPointOutsideMother_solidVolume()) return GetFibreCoreRoughenedStartPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding roughened end point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingRoughenedEndPointInsideMother_solidVolume()
	{
		if(GetThirdCladdingRoughenedEndPointInsideMother_logicalVolume()) return GetThirdCladdingRoughenedEndPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding roughened end point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingRoughenedEndPointOutsideMother_solidVolume()
	{
		if(GetThirdCladdingRoughenedEndPointOutsideMother_logicalVolume()) return GetThirdCladdingRoughenedEndPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 3. cladding roughened end point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetThirdCladdingRoughenedEndPoint_basicSolidVolume()
	{
		if(GetThirdCladdingRoughenedEndPointInsideMother_solidVolume()) return GetThirdCladdingRoughenedEndPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetThirdCladdingRoughenedEndPointOutsideMother_solidVolume()) return GetThirdCladdingRoughenedEndPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding roughened end point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingRoughenedEndPointInsideMother_solidVolume()
	{
		if(GetSecondCladdingRoughenedEndPointInsideMother_logicalVolume()) return GetSecondCladdingRoughenedEndPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding roughened end point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingRoughenedEndPointOutsideMother_solidVolume()
	{
		if(GetSecondCladdingRoughenedEndPointOutsideMother_logicalVolume()) return GetSecondCladdingRoughenedEndPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 2. cladding roughened end point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetSecondCladdingRoughenedEndPoint_basicSolidVolume()
	{
		if(GetSecondCladdingRoughenedEndPointInsideMother_solidVolume()) return GetSecondCladdingRoughenedEndPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetSecondCladdingRoughenedEndPointOutsideMother_solidVolume()) return GetSecondCladdingRoughenedEndPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding roughened end point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingRoughenedEndPointInsideMother_solidVolume()
	{
		if(GetFirstCladdingRoughenedEndPointInsideMother_logicalVolume()) return GetFirstCladdingRoughenedEndPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding roughened end point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingRoughenedEndPointOutsideMother_solidVolume()
	{
		if(GetFirstCladdingRoughenedEndPointOutsideMother_logicalVolume()) return GetFirstCladdingRoughenedEndPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s 1. cladding roughened end point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFirstCladdingRoughenedEndPoint_basicSolidVolume()
	{
		if(GetFirstCladdingRoughenedEndPointInsideMother_solidVolume()) return GetFirstCladdingRoughenedEndPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetFirstCladdingRoughenedEndPointOutsideMother_solidVolume()) return GetFirstCladdingRoughenedEndPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core roughened end point exists inside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreRoughenedEndPointInsideMother_solidVolume()
	{
		if(GetFibreCoreRoughenedEndPointInsideMother_logicalVolume()) return GetFibreCoreRoughenedEndPointInsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core roughened end point exists outside the mother volume: </b> pointer to this solid volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreRoughenedEndPointOutsideMother_solidVolume()
	{
		if(GetFibreCoreRoughenedEndPointOutsideMother_logicalVolume()) return GetFibreCoreRoughenedEndPointOutsideMother_logicalVolume()->GetSolid();
		else return 0;
	}

	/**
	 *  @return <b> if a solid volume for a G4Fibre%'s core roughened end point exists: </b> pointer to its basic solid volume, i.e. the solid volume before the G4Fibre was cut into parts inside and outside the mother volume
	 *  @return <b> else: </b> 0
	 */
	G4VSolid * GetFibreCoreRoughenedEndPoint_basicSolidVolume()
	{
		if(GetFibreCoreRoughenedEndPointInsideMother_solidVolume()) return GetFibreCoreRoughenedEndPointInsideMother_solidVolume()->GetConstituentSolid(0);
		if(GetFibreCoreRoughenedEndPointOutsideMother_solidVolume()) return GetFibreCoreRoughenedEndPointOutsideMother_solidVolume()->GetConstituentSolid(0);
		else return 0;
	}

// materials of the object:
	/**
	 *  @return <b> if the G4Fibre is glued: </b> pointer to the material of the optical cement
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetOpticalCementMaterial()
	{ return Material_OpticalCement; }

	/**
	 *  @return <b> if a coating exists: </b> pointer to the material of the G4Fibre%'s coating
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetacketMaterial()
	{ return Material_Coating; }

	/**
	 *  @return <b> if a 3. cladding exists: </b> pointer to the material of the G4Fibre%'s 3. cladding
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetThirdCladdingMaterial()
	{ return Material_Cladding3; }

	/**
	 *  @return <b> if a 2. cladding exists: </b> pointer to the material of the G4Fibre%'s 2. cladding
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetSecondCladdingMaterial()
	{ return Material_Cladding2; }

	/**
	 *  @return <b> if a 1. cladding exists: </b> pointer to the material of the G4Fibre%'s 1. cladding
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetFirstCladdingMaterial()
	{ return Material_Cladding1; }

	/** @return pointer to the material of the G4Fibre%'s core */
	G4Material * GetFibreCoreMaterial()
	{ return Material_FibreCore; }

// optical surfaces of the object:
	/**
	 *  @return <b> if a pointer to the optical surface of the G4Fibre%'s start point exists: </b> pointer to the optical surface of the G4Fibre%'s start point
	 *  @return <b> else: </b> 0
	 */
	G4OpticalSurface * GetReflectiveStartPointOpticalSurface()
	{ return OptSurf_startPoint; }

	/**
	 *  @return <b> if a pointer to the optical surface of the G4Fibre%'s end point exists: </b> pointer to the optical surface of the G4Fibre%'s end point
	 *  @return <b> else: </b> 0
	 */
	G4OpticalSurface * GetReflectiveEndPointOpticalSurface()
	{ return OptSurf_endPoint; }

// sensitive detector:
	/** @return G4bool, if the G4Fibre has a sensitive detector */
	G4bool HasSensitiveDetector()
	{ return ConstructSensitiveDetector; }

	/** @return pointer to the sensitive detector */
	FibreSensitiveDetector * GetSensitiveDetector()
	{ return ((FibreSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("PhotonDetectorSD", false)); }

// other properties of the object:
	/** @return G4bool, if the fibre is bent */
	G4bool IsFibreBent()
	{ return FibreBent; }

	/** @return G4bool, if the fibre is round */
	G4bool IsFibreRound()
	{ return FibreRound; }

	/** @return G4bool, if it is a scintillating fibre */
	G4bool IsScintillatingFibre()
	{ return IsScinti; }

	/** @return G4bool, if it is a wave-length-shifting fibre */
	G4bool IsWLSFibre()
	{ return IsWLS; }

	/** @return G4bool, if it is a light-guiding fibre (neither wave-length-shifting nor scintillating) */
	G4bool IsLightGuidingFibre()
	{
		if(!IsWLS && !IsScinti) return true;
		else return false;
	}

	/** @return G4bool, if the fibre is glued */
	G4bool IsFibreGlued()
	{ return Glued; }

	/** @return G4bool, if a coating exists */
	G4bool DoesCoatingExist()
	{ return CoatingExists; }

	/** @return G4bool, if a 3. cladding exists */
	G4bool DoesThirdCladdingExist()
	{ return Cladding3Exists; }

	/** @return G4bool, if a 2. cladding exists */
	G4bool DoesSecondCladdingExist()
	{ return Cladding2Exists; }

	/** @return G4bool, if a 1. cladding exists */
	G4bool DoesFirstCladdingExist()
	{ return Cladding1Exists; }

	/** @return G4bool, if the fibre has a reflective start point */
	G4bool FibreHasReflectiveStartPoint()
	{ return CreateReflectiveStartPointVolumes; }

	/** @return G4bool, if the fibre has a reflective end point */
	G4bool FibreHasReflectiveEndPoint()
	{ return CreateReflectiveEndPointVolumes; }

	/** @return G4bool, if the fibre has a roughened start point */
	G4bool FibreHasRoughenedStartPoint()
	{
		if(!CreateReflectiveStartPointVolumes && Roughness_startPoint >= 1e-12) return true;
		else return false;
	}

	/** @return G4bool, if the fibre has a roughened end point */
	G4bool FibreHasRoughenedEndPoint()
	{
		if(!CreateReflectiveEndPointVolumes && Roughness_endPoint >= 1e-12) return true;
		else return false;
	}

	/**
	 *  @return <b> if the fibre has a reflective start point: </b> reflectivity of the fibre's reflective start point
	 *  @return <b> else: </b> NAN
	 */
	G4double GetFibreStartPointReflectivity()
	{
		if(CreateReflectiveStartPointVolumes) return Reflectivity_startPoint;
		else return NAN;
	}

	/**
	 *  @return <b> if the fibre has a reflective end point: </b> reflectivity of the fibre's reflective end point
	 *  @return <b> else: </b> NAN
	 */
	G4double GetFibreEndPointReflectivity()
	{
		if(CreateReflectiveEndPointVolumes) return Reflectivity_endPoint;
		else return NAN;
	}

	/**
	 *  @return <b> if the fibre has a no reflective start point: </b> roughness of the fibre's start point
	 *  @return <b> else: </b> NAN
	 */
	G4double	 GetFibreStartPointRoughness()
	{
		if(!CreateReflectiveStartPointVolumes) return Roughness_startPoint;
		else return NAN;
	}

	/**
	 *  @return <b> if the fibre has no reflective end point: </b> roughness of the fibre's end point
	 *  @return <b> else: </b> NAN
	 */
	G4double	 GetFibreEndPointRoughness()
	{
		if(!CreateReflectiveEndPointVolumes) return Roughness_endPoint;
		else return NAN;
	}

private:
	void SetDefaults();

	void InitialiseVariables();
	void GenerateTransformation(G4String fibreType	/**< which type of G4Fibre is to be created (bent or straight) */
				   );

	void DefineMaterials();

	void DefineOpticalCementMaterialProperties();
	void DefineCoreMaterialProperties();
	void DefineCladding1MaterialProperties();
	void DefineCladding2MaterialProperties();
	void DefineCladding3MaterialProperties();
	void DefineCoatingMaterialProperties();

	void ConstructSurface();


	void ConstructVolumes(G4String fibreType	/**< which type of G4Fibre is to be created (bent or straight) */
			     );

	std::vector<G4VPhysicalVolume *>  ConstructFibreLayerPhysical( std::vector<G4VPhysicalVolume *> &mothers,	/**< vector of mother volumes */
								       G4Transform3D cutTransformation,	/**< transformation for cutting the G4Fibre into parts inside and outside the mother volume */
								       G4Transform3D transformation,	/**< transformation for placing the G4Fibre layer into the mother volume */
								       G4String nameBase,		/**< base of the G4Fibre%'s volumes' names (it will be extended to distinguish between different volumes belonging to one G4Fibre) */
								       G4double fibreRadiusMin_rel,	/**< G4Fibre layer's relative inner radius with respect to the G4Fibre%'s full radius */
								       G4double fibreRadiusMax_rel,	/**< G4Fibre layer's relative outer radius with respect to the G4Fibre%'s full radius */
								       G4Material * material,		/**< G4Fibre layer's material */
								       G4Colour colour_in,		/**< G4Fibre layer's colour for parts inside the mother volume (for visualisation) */
								       G4Colour colour_out,		/**< G4Fibre layer's colour for parts outside the mother volume (for visualisation) */
								       G4String fibreType,		/**< which type of G4Fibre is to be created (bent or straight) */
								       G4String fibrePart = ""		/**< which part of the G4Fibre is to be created (roughened end, reflective end or normal part) */
								     );


	std::vector<boost::any> ConstructFibreLayerLogical( G4Transform3D cutTransformation,	/**< transformation for cutting the G4Fibre into parts inside and outside the mother volume */
							    G4Transform3D transformation,	/**< transformation for placing the G4Fibre layer into the mother volume */
							    G4String nameBase,			/**< base of the G4Fibre%'s volumes' names (it will be extended to distinguish between different volumes belonging to one G4Fibre) */
							    G4double fibreRadiusMin_rel,		/**< G4Fibre layer's relative inner radius with respect to the G4Fibre%'s full radius */
							    G4double fibreRadiusMax_rel,		/**< G4Fibre layer's relative outer radius with respect to the G4Fibre%'s full radius */
							    G4Material * material,		/**< G4Fibre layer's material */
							    G4Colour colour_in,			/**< G4Fibre layer's colour for parts inside the mother volume (for visualisation) */
							    G4Colour colour_out,			/**< G4Fibre layer's colour for parts outside the mother volume (for visualisation) */
							    G4String fibreType,			/**< which type of G4Fibre is to be created (bent or straight) */
							    G4String fibrePart			/**< which part of the G4Fibre is to be created (roughened end, reflective end or normal part) */
							  );


	std::vector<G4VPhysicalVolume *>  ConstructEmbedmentPhysical( G4String nameBase,		/**< base of the G4Fibre%'s volumes' names (it will be extended to distinguish between different volumes belonging to one G4Fibre) */
								      G4String fibreType		/**< which type of G4Fibre is to be created (bent or straight) */
								     );


	std::vector<boost::any> ConstructEmbedmentLogical( G4String nameBase,			/**< base of the G4Fibre%'s volumes' names (it will be extended to distinguish between different volumes belonging to one G4Fibre) */
							   G4String fibreType			/**< which type of G4Fibre is to be created (bent or straight) */
							  );



	std::vector<G4VPhysicalVolume *> findGrandMotherAndAuntVolumes(G4bool cutAuntVolumesWithDaughters);



	G4bool CriticalErrorOccured;
	G4bool SearchOverlaps;
	G4bool ConstructSensitiveDetector;
	PropertyToolsManager * PropertyTools;
	GODDeSS_DataStorage * DataStorage;
	G4VPhysicalVolume * MotherVolume_physical;
	G4VSolid * MotherVolume_solid;
	G4VPhysicalVolume * ReferenceVolume_physical;
	G4Transform3D ReferenceVolume_transformation;
	std::vector<G4VPhysicalVolume *> GrandMotherAndAuntVolumes;

// Materials & Elements
	G4Material* Material_OpticalCement;
	G4Material* Material_FibreCore;
	G4Material* Material_Cladding1;
	G4Material* Material_Cladding2;
	G4Material* Material_Cladding3;
	G4Material* Material_Coating;

// Fibre
	GoddessProperties FibreProperties;
	G4ThreeVector StartPoint;
	G4ThreeVector EndPoint;
	G4ThreeVector BendingCircularCentre;
	G4ThreeVector BendingAxis;
	G4double BendingRadius;
	G4double BendingDeltaAngle;
	G4double BendingStartAngle;
	G4double BendingEndAngle;
	G4double Length;
	G4Transform3D FibreTransformation_rel;
	G4String FibreNamePrefix;
	G4bool OnlyInsideMother;
	G4bool CutAuntVolumesWithDaughters;
	G4bool FibreBent;
	G4bool FibreRound;
	G4bool IsScinti;
	G4bool IsWLS;
	G4bool Glued;
	G4bool Cladding1Exists;
	G4bool Cladding2Exists;
	G4bool Cladding3Exists;
	G4bool CoatingExists;
	G4double FibreRadius;
	G4double RelativeFibreRadius_Core;
	G4double RelativeFibreRadius_Cladding1;
	G4double RelativeFibreRadius_Cladding2;
	G4double RelativeFibreRadius_Cladding3;
	G4double RelativeFibreRadius_CoatingMin;
	G4double RelativeFibreRadius_CoatingMax;
	G4double FibreEdgeLength;
	G4Transform3D FibreTransformation_insideMother;
	G4Transform3D FibreTransformation_outsideMother;
	std::vector<G4VPhysicalVolume *> FibreCore_physical;
	std::vector<G4VPhysicalVolume *> Cladding1_physical;
	std::vector<G4VPhysicalVolume *> Cladding2_physical;
	std::vector<G4VPhysicalVolume *> Cladding3_physical;
	std::vector<G4VPhysicalVolume *> Coating_physical;
	std::vector<G4VPhysicalVolume *> FibreCore_mothers_physical;
	std::vector<G4VPhysicalVolume *> Cladding1_mothers_physical;
	std::vector<G4VPhysicalVolume *> Cladding2_mothers_physical;
	std::vector<G4VPhysicalVolume *> Cladding3_mothers_physical;
	std::vector<G4VPhysicalVolume *> Coating_mothers_physical;

// Embedment
	std::vector<G4VPhysicalVolume *> FibreEmbedment_physical;
	GoddessProperties OpticalCementProperties;
	G4String EmbedmentProfile;
	G4String WhereIsEmbedmentToBeFlushWithFibre;
	std::vector<G4VPhysicalVolume *> VolumesWithEmbedment;
	G4double EmbedmentThickness;

// Reflective/roughend/smooth fibre end
	G4OpticalSurface * OptSurf_startPoint;
	G4OpticalSurface * OptSurf_endPoint;
	G4OpticalSurface * OptSurf_perfectlySmooth;

	G4bool CreateReflectiveStartPointVolumes;
	G4double Reflectivity_startPoint;
	G4Transform3D ReflectiveStartPointTransformation_insideMother;
	G4Transform3D ReflectiveStartPointTransformation_outsideMother;
	std::vector<G4VPhysicalVolume *> ReflectiveStartPoint_physical;
	std::vector<G4VPhysicalVolume *> ReflectiveStartPoint_mothers_physical;

	G4double Roughness_startPoint;
	G4Transform3D RoughenedStartPointTransformation_insideMother;
	G4Transform3D RoughenedStartPointTransformation_outsideMother;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_cladding3_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_cladding2_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_cladding1_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_fibreCore_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_cladding3_mothers_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_cladding2_mothers_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_cladding1_mothers_physical;
	std::vector<G4VPhysicalVolume *> RoughenedStartPoint_fibreCore_mothers_physical;

	G4bool CreateReflectiveEndPointVolumes;
	G4double Reflectivity_endPoint;
	G4Transform3D ReflectiveEndPointTransformation_insideMother;
	G4Transform3D ReflectiveEndPointTransformation_outsideMother;
	std::vector<G4VPhysicalVolume *> ReflectiveEndPoint_physical;
	std::vector<G4VPhysicalVolume *> ReflectiveEndPoint_mothers_physical;

	G4double Roughness_endPoint;
	G4Transform3D RoughenedEndPointTransformation_insideMother;
	G4Transform3D RoughenedEndPointTransformation_outsideMother;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_cladding3_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_cladding2_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_cladding1_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_fibreCore_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_cladding3_mothers_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_cladding2_mothers_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_cladding1_mothers_physical;
	std::vector<G4VPhysicalVolume *> RoughenedEndPoint_fibreCore_mothers_physical;

	G4double Reflective_thickness;
	G4double Reflective_bendingDeltaAngle;

	G4double Rough_thickness;
	G4double Rough_bendingDeltaAngle;



	G4Transform3D ReplicationTransformation;
};

#endif
