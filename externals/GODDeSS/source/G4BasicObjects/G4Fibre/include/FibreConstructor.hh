/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef FIBRECONSTRUCTOR_H
#define FIBRECONSTRUCTOR_H

#include <G4VModularPhysicsList.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VPhysicalVolume.hh>

#include <G4Fibre.hh>
#include <PropertyToolsManager.hh>
#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class making it easier and more flexible to create G4Fibre objects. <b> Additionally, it automatically deals with the registration of the needed physics! </b>
class FibreConstructor
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 *  - adds the needed physics processes (specified in FibrePhysicsList) to the physics list
	 */
	FibreConstructor( G4VModularPhysicsList * userPhysicsList,	/**< physics list which is used for the simulation */
			  PropertyToolsManager * propertyTools,		/**< pointer to the PropertyToolsManager that is to be used */
			  GODDeSS_DataStorage * dataStorage,		/**< pointer to the GODDeSS_DataStorage that is to be used */
			  G4bool searchOverlaps,				/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
			  int verbose = 0				/**< verbosity level for adding the needed physics processes (specified in FibrePhysicsList) to the physics list (default: 0) */
			)
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	{
		if(SearchOverlaps)
		{
			SearchOverlaps = false;   // NOTE Due to the previous commit (Automatic distributes of the fibres to the different volumes), the fibres cannot be checked for overlaps anymore...
			G4cerr << "##########" <<
			G4endl << "# WARNING:" <<
			G4endl << "# Checking for overlaps is VERY SLOW for G4Fibres." <<
			G4endl << "# It has therefore been deactivated in the constructor of the FibreConstructor class (source/G4BasicObjects/G4Fibre/include/FibreConstructor.hh)." <<
			G4endl << "# If you desperately want to check G4Fibres for overlaps, comment out the if-statement that contains this warning-message." <<
			G4endl << "##########" << G4endl;
		}

		// adds the needed physics processes to the userPhysicsList
		LoadPhysicsList(userPhysicsList, verbose);

		SetDefaults();
	}

	/**
	 *  Destructor (empty)
	 */
	~FibreConstructor()
	{
	}



	// These are overloaded functions. They are calling the main function (defined below) using default values for the parameters which are not specified.
	#include <FibreConstructor_bent.icc>



	// These are the main function for creating the volumes. They will be overloaded to allow setting the optional parameters via the function call or via set-functions.

	// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a BENT fibre with a given CIRCULAR CENTRE, BENDING AXIS, BENDING RADIUS, START and END ANGLE </b>
	vector<G4Fibre *> ConstructFibres( G4String fibre_property_file,				/**< path to the file containing the G4Fibre%'s properties */
					   G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
					   G4double bending_start_angle,					/**< G4Fibre%'s start angle (in the bending circle) */
					   G4double bending_end_angle,					/**< G4Fibre%'s end angle (in the bending circle) */
					   G4ThreeVector bending_circular_centre,			/**< circular centre of the G4Fibre%'s bending */
					   G4ThreeVector bending_axis,					/**< G4Fibre%'s bending axis */
					   G4double bending_radius,					/**< G4Fibre%'s bending radius */
					   G4int number_of_replications,					/**< number of replications that are to be created */
					   G4Transform3D replication_transf,				/**< transformation of each replication, relative to the previous one */
					   G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
					   G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
					   G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
					   G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
					   G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
					   G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
					   G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
					   G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
					   G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
					   std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
					   G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
					 );

	// option to only one object
	/// <b> >>> construct a BENT fibre with a given CIRCULAR CENTRE, BENDING AXIS, BENDING RADIUS, START and END ANGLE </b>
	G4Fibre * ConstructFibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s properties */
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
				  G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
				  G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
				  G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
				  G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
				  G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
				  std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
				  G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
				);



	// These are overloaded functions. They are calling the main function (defined below) using default values for the parameters which are not specified.
	#include <FibreConstructor_bent_startEnd.icc>



	// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a BENT fibre with a given BENDING AXIS, BENDING ANGLE, START and END POINT </b>
	vector<G4Fibre *> ConstructFibres( G4String fibre_property_file,				/**< path to the file containing the G4Fibre%'s properties */
					   G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
					   G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
					   G4ThreeVector endPoint_relative_to_reference,		/**< G4Fibre%'s end point relative to the reference volume */
					   G4double bending_delta_angle,					/**< G4Fibre%'s bending angle */
					   G4ThreeVector bending_axis,					/**< G4Fibre%'s bending axis */
					   G4int number_of_replications,					/**< number of replications that are to be created */
					   G4Transform3D replication_transf,				/**< transformation of each replication, relative to the previous one */
					   G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
					   G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
					   G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
					   G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
					   G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
					   G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
					   G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
					   G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
					   G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
					   std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
					   G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
					 );

	// option to only one object
	/// <b> >>> construct a BENT fibre with a given BENDING AXIS, BENDING ANGLE, START and END POINT </b>
	G4Fibre * ConstructFibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s properties */
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
				  G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
				  G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
				  G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
				  G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
				  G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
				  std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
				  G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
				);



	// These are overloaded functions. They are calling the main function (defined below) using default values for the parameters which are not specified.
	#include <FibreConstructor_straight.icc>



	// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/// <b> >>> construct replications of a STRAIGHT fibre with a given LENGTH and TRANSFORMATION </b>
	// This is the main function for creating the volumes. It will be overloaded to allow setting the optional parameters via the function call or via set-functions.
	vector<G4Fibre *> ConstructFibres( G4String fibre_property_file,				/**< path to the file containing the G4Fibre%'s properties */
					   G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
					   G4double length,						/**< G4Fibre%'s length */
					   G4int number_of_replications,					/**< number of replications that are to be created */
					   G4Transform3D replication_transf,				/**< transformation of each replication, relative to the previous one */
					   G4Transform3D transf_relative_to_reference,			/**< G4Fibre%'s transformation relative to the reference volume */
					   G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
					   G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
					   G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
					   G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
					   G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
					   G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
					   G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
					   G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
					   G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
					   std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
					   G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
					 );

	// option to only one object
	/// <b> >>> construct a STRAIGHT fibre with a given LENGTH and TRANSFORMATION </b>
	G4Fibre * ConstructFibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s properties */
				  G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
				  G4double length,						/**< G4Fibre%'s length */
				  G4Transform3D transf_relative_to_reference,			/**< G4Fibre%'s transformation relative to the reference volume */
				  G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
				  G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
				  G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
				  G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
				  G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
				  G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
				  G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
				  G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
				  G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
				  G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
				  G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
				  std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
				  G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
				);



	// These are overloaded functions. They are calling the main function (defined below) using default values for the parameters which are not specified.
	#include <FibreConstructor_straight_startEnd.icc>



	// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	// This is the main function for creating the volumes. It will be overloaded to allow setting the optional parameters via the function call or via set-functions.
	/// <b> >>> construct replications of a STRAIGHT fibre with a given START and END POINT </b>
	vector<G4Fibre *> ConstructFibres( G4String fibre_property_file,				/**< path to the file containing the G4Fibre%'s properties */
					   G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
					   G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
					   G4ThreeVector endPoint_relative_to_reference,		/**< G4Fibre%'s end point relative to the reference volume */
					   G4int number_of_replications,					/**< number of replications that are to be created */
					   G4Transform3D replication_transf,				/**< transformation of each replication, relative to the previous one */
					   G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
					   G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
					   G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
					   G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
					   G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
					   G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
					   G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
					   G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
					   G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
					   G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
					   G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
					   std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
					   G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
					 );

	// option to only one object
	/// <b> >>> construct a STRAIGHT fibre with a given START and END POINT </b>
	G4Fibre * ConstructFibre( G4String fibre_property_file,					/**< path to the file containing the G4Fibre%'s properties */
				  G4VPhysicalVolume* mother_volume,				/**< volume that is to be the G4Fibre%'s mother volume */
				  G4ThreeVector startPoint_relative_to_reference,		/**< G4Fibre%'s start point relative to the reference volume */
				  G4ThreeVector endPoint_relative_to_reference,			/**< G4Fibre%'s end point relative to the reference volume */
				  G4VPhysicalVolume* reference_volume,				/**< reference volume relative to which the G4Fibre will be orientated (if 0, the mother volume is the reference volume) */
				  G4Transform3D reference_transformation,			/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
				  G4String fibre_name_prefix,					/**< prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre) */
				  G4bool onlyInsideMother,					/**< if "true", only G4Fibre parts inside the G4Fibre%'s mother volume will be created */
				  G4bool cutAuntVolumesWithDaughters,				/**< if "true", only NO G4Fibre parts will be created inside the G4Fibre%'s aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves */
				  G4double startPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s beginning (start point) reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
				  G4double startPoint_sigmaAlpha,				/**< roughness of the G4Fibre%'s beginning (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4double endPoint_reflectivity,				/**< reflectivity of the G4Fibre%'s end reflective layer (NAN means, that no reflective/absorbing layer is to be added) */
				  G4double endPoint_sigmaAlpha,					/**< roughness of the G4Fibre%'s end (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian - no influence for reflective layers) */
				  G4bool glued,							/**< the G4Fibre is to be glued to at least one volume ("true" or "false") */
				  G4String glue_property_file,					/**< path to the file containing the glue's properties (the glue may consist of any material composition, also air or vacuum) */
				  G4String embedment_profile,					/**< profile that the G4Fibre%'s glue volume should have ("round" or "quadratic") */
				  G4String where_is_embedment_to_be_flush_with_fibre,		/**< the ends at which the G4Fibre%'s glue volume is to be flush with the G4Fibre ("" = none, "S" = start point, "E" = end point or "SE" = start and end point) */
				  std::vector<G4VPhysicalVolume *> volumes_with_embedment,	/**< the volumes to which the G4Fibre is to be glued */
				  G4bool constructSensitiveDetector				/**< a sensitive detector is to be constructed ("true" or "false") */
				);



	// These are the setter-functions.
	/**
	 *  Function to set the G4Fibre%'s transformation relative to the reference volume.
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreTransformationRelativeToReferenceVolume(G4Transform3D transf_relative_to_reference)
	{
		FibreTransformation_rel = transf_relative_to_reference;
	}

	/**
	 *  Function to set the prefix of the G4Fibre%'s volumes' names (it will be extended to distinguish between different G4Fibre%s and different volumes of one G4Fibre).
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreNamePrefix(G4String fibre_name_prefix)
	{
		FibreNamePrefix = fibre_name_prefix;
	}

	/**
	 *  Function to make the G4Fibre to be placed only into its mother volume, i.e. all parts outside the mother volume will be cut away.
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetCreateFibreOnlyInsideMotherVolume()
	{
		OnlyInsideMother = true;
	}

	/**
	 *  Function to make the G4Fibre NOT to be placed into aunt volumes (volumes that have been places into the same volume like the the mother volume) that contain daugther volumes themselves.
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetDoNotCreateFibreInsideAuntVolumesWithDaughters()
	{
		CutAuntVolumesWithDaughters = true;
	}

	/**
	 *  Function to set the reference volume relative to which the G4Fibre will be orientated.
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreReferenceVolume(G4VPhysicalVolume* reference_volume)
	{
		ReferenceVolume_physical = reference_volume;
	}

	/**
	 *  Function to set the transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume).
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreReferenceVolumeTransformation(G4Transform3D reference_volume_transformation)
	{
		ReferenceVolume_transformation = reference_volume_transformation;
	}

	/**
	 *  Function to set the reflectivity of the G4Fibre%'s beginning (start point) reflective layer.\n
	 *  The G4Fibre%'s beginning (start point) will have a reflective layer using this reflectivity.
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreStartPointReflectivity(G4double startPoint_reflectivity)
	{
		Reflectivity_startPoint = startPoint_reflectivity;
	}

	/**
	 *  Function to set roughness of the G4Fibre%'s end.\n
	 *  The parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian.\n
	 *  The G4Fibre%'s end will be roughened using this roughness.
	 *  ! In rare cases (if the thin part of the fibre, which creates the roughened end, is seperated between different volumes AND the fibre is orientated vertically to the dividing surface AND a optical photon is created inside this thin part AND the optical photon is orientated vertically to the dividing surface) this may lead to unphysical behaviour !
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreStartPointRoughness(G4double startPoint_sigmaAlpha)
	{
		Roughness_startPoint = startPoint_sigmaAlpha;
	}

	/**
	 *  Function to set the reflectivity of the G4Fibre%'s beginning (start point) reflective layer.\n
	 *  The G4Fibre%'s beginning (start point) will have a reflective layer using this reflectivity.
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreEndPointReflectivity(G4double endPoint_reflectivity)
	{
		Reflectivity_endPoint = endPoint_reflectivity;
	}

	/**
	 *  Function to set roughness of the G4Fibre%'s end.\n
	 *  The parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian.\n
	 *  The G4Fibre%'s end will be roughened using this roughness.
	 *  ! In rare cases (if the thin part of the fibre, which creates the roughened end, is seperated between different volumes AND the fibre is orientated vertically to the dividing surface AND a optical photon is created inside this thin part AND the optical photon is orientated vertically to the dividing surface) this may lead to unphysical behaviour !
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreEndPointRoughness(G4double endPoint_sigmaAlpha)
	{
		Roughness_endPoint = endPoint_sigmaAlpha;
	}

	/**
	 *  \overload
	 *  \n
	 *
	 *  This functon can be used to define multiple volumes in which the embedment is to be created.
	 */
	void SetFibreGlued(G4String glue_property_file, G4String embedment_profile, std::vector<G4VPhysicalVolume *> volumes_with_embedment, G4String where_is_embedment_to_be_flush_with_fibre = "")
	{
		Glued = true;
		OpticalCementPropertyFile = glue_property_file;
		EmbedmentProfile = embedment_profile;
		VolumesWithEmbedment = volumes_with_embedment;
		WhereIsEmbedmentToBeFlushWithFibre = where_is_embedment_to_be_flush_with_fibre;
	}

	/**
	 *  Function to set the path to the file containing the glue's property file as well as the profile of the glue's volume ("round" or "quadratic") and the volume in which the embedment is to be created (default: Mothervolume of the G4Fibre).\n
	 *  The glue may also be air.\n
	 *  If this function is used, a volume with the specified profile will be placed around the G4Fibre.\n
	 *  Otherwise, the G4Fibre will be "moulded" into the mother volume
	 *
	 *  <b> The value set by this function will only apply to the very next G4Fibre object that is created! </b>
	 */
	void SetFibreGlued(G4String glue_property_file, G4String embedment_profile, G4VPhysicalVolume * volume_with_embedment = 0, G4String where_is_embedment_to_be_flush_with_fibre = "")
	{
		Glued = true;
		OpticalCementPropertyFile = glue_property_file;
		EmbedmentProfile = embedment_profile;

		if(volume_with_embedment)
		{
			VolumesWithEmbedment.clear();
			VolumesWithEmbedment.push_back(volume_with_embedment);
		}

		WhereIsEmbedmentToBeFlushWithFibre = where_is_embedment_to_be_flush_with_fibre;
	}

private:
	void applyVolumeCounter( G4String & physical_volume_name /**< name prefix which is to be used for naming the volumes */
			       );
	void SetDefaults();

	void registerPhysics(G4VModularPhysicsList * physicsList, G4VPhysicsConstructor * physicsConstructor);
	void LoadPhysicsList(G4VModularPhysicsList * physicsList, int verbose = 0);


	G4bool ConstructSensitiveDetector;
	G4bool SearchOverlaps;
	PropertyToolsManager * PropertyTools;
	GODDeSS_DataStorage * DataStorage;
	std::vector<G4String> NamePrefixVector;

	G4VPhysicalVolume * ReferenceVolume_physical;
	G4Transform3D ReferenceVolume_transformation;
	G4Transform3D FibreTransformation_rel;
	G4String FibreNamePrefix;
	G4bool OnlyInsideMother;
	G4bool CutAuntVolumesWithDaughters;
	G4bool Glued;
	G4String OpticalCementPropertyFile;
	G4String EmbedmentProfile;
	G4String WhereIsEmbedmentToBeFlushWithFibre;
	std::vector<G4VPhysicalVolume *> VolumesWithEmbedment;
	G4double Reflectivity_startPoint;
	G4double Roughness_startPoint;
	G4double Reflectivity_endPoint;
	G4double Roughness_endPoint;
};

#endif
