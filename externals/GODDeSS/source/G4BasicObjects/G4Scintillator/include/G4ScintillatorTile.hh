/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef G4SCINTILLATORTile_H
#define G4SCINTILLATORTile_H

#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4VPhysicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4Material.hh>
#include <G4SDManager.hh>

#include <GoddessProperties.hh>
#include <PropertyToolsManager.hh>
#include <ScintillatorSensitiveDetector.hh>

#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class generating the materials, volumes, and optical properties needed for a scintillator tile and allows to access them.
class G4ScintillatorTile
{
public:

	/**
	 *  Constructor to construct a scintillator tile:
	 *  - sets class variables to default values
	 *  - loads the property file
	 *  - sets the variables for the scintillator tile's placement (considering the transformation of the scintillator tile relative to the reference volume and the transformation of the reference volume relative to the mother volume)
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructScintiVolume())
	 */
	G4ScintillatorTile( G4ThreeVector scintillator_dimensions,				/**< dimensions of the G4ScintillatorTile */
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
			    G4bool constructSensitiveDetector,					/**< a sensitive detector is to be constructed ("true" or "false") */
			    G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
			    PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
			    G4OpticalSurfaceModel surfaceModel,					/**< model for the optical surfaces */
			    G4OpticalSurfaceFinish surfaceFinish,				/**< finish for the optical surfaces */
			    GODDeSS_DataStorage * dataStorage					/**< pointer to the GODDeSS_DataStorage that is to be used */
			  )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, SurfaceModel(surfaceModel)
	, SurfaceFinish(surfaceFinish)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, Dimensions(scintillator_dimensions)
	, ScintillatorName(scintillator_name)
	, CutDimple(cut_dimple)
	, DimpleRadius(dimple_radius)
	, DimpleCentrePosition(dimple_centre_relative_to_scintillator_centre)
	, TwinTile(twin_tile)
	{
		if(reference_volume != MotherVolume_physical && reference_volume->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		// the obligatory parameters
		ScintiProperties.load(scintillator_property_file);
		if(TwinTile) TwinTileProperties.load(twin_tile_property_file);


		G4RotationMatrix rotation = MotherVolume_physical->GetObjectRotationValue().inverse() * reference_volume->GetObjectRotationValue() * reference_transformation.getRotation() * transf_relative_to_reference.getRotation();

		// WARNING-NOTE: as transform() changes the vector it is applied to, the following 7 commands are needed instead of a 1 line command
		G4ThreeVector translation = transf_relative_to_reference.getTranslation();
		translation.transform(reference_volume->GetObjectRotationValue());
		translation += reference_volume->GetObjectTranslation();
		translation.transform(reference_transformation.getRotation());
		translation += reference_transformation.getTranslation();
		translation -= MotherVolume_physical->GetObjectTranslation();
		translation.transform(MotherVolume_physical->GetObjectRotationValue().inverse());

		Transformation = G4Transform3D(rotation, translation);


		// create materials:
		DefineMaterials();

		// create the volumes
		ConstructScintiVolume();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/**
	 *  Constructor to construct replications of a scintillator tile:
	 *  - sets class variables to default values
	 *  - loads the property file
	 *  - sets the variables for the scintillator tile's placement (considering the transformation of the scintillator tile relative to the reference volume and the transformation of the reference volume relative to the mother volume)
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructScintiVolume())
	 */
	G4ScintillatorTile( G4ThreeVector scintillator_dimensions,				/**< dimensions of the G4ScintillatorTile */
			    G4String scintillator_properties_file,				/**< path to the file containing the G4ScintillatorTile%'s properties */
			    G4VPhysicalVolume* mother_volume,					/**< G4ScintillatorTile%'s mother volume */
			    G4String scintillator_name,						/**< name of the G4ScintillatorTile%'s volume (it will be extended to distinguish between different G4ScintillatorTile%s and different volumes of one G4ScintillatorTile) */
			    G4Transform3D transf_relative_to_reference,				/**< G4ScintillatorTile%'s transformation relative to the reference volume */
			    G4VPhysicalVolume* reference_volume,					/**< reference volume relative to which the G4ScintillatorTile will be orientated (if 0, the mother volume is the reference volume) */
			    G4Transform3D reference_transformation,				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
			    G4bool cut_dimple,							/**< a dimple is to be cut out out of the G4ScintillatorTile ("true" or "false") */
			    G4double dimple_radius,						/**< radius of the dimple */
			    G4ThreeVector dimple_centre_relative_to_scintillator_centre,	/**< position of the dimple's centre relative to the G4ScintillatorTile%'s centre */
			    G4bool twin_tile,							/**< the G4ScintillatorTile is to be a twin tile ("true" or "false") */
			    G4String twin_tile_properties_file,					/**< path to the file containing the properties of the twin tile's reflective layer */
			    G4bool constructSensitiveDetector,					/**< a sensitive detector is to be constructed ("true" or "false") */
			    G4bool searchOverlaps,						/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
			    PropertyToolsManager * propertyTools,				/**< pointer to the PropertyToolsManager that is to be used */
			    G4OpticalSurfaceModel surfaceModel,					/**< model for the optical surfaces */
			    G4OpticalSurfaceFinish surfaceFinish,				/**< finish for the optical surfaces */
			    GODDeSS_DataStorage * dataStorage,					/**< pointer to the GODDeSS_DataStorage that is to be used */
			    G4Transform3D replication_transf					/**< transformation of each replication, relative to the previous one */
			  )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, SurfaceModel(surfaceModel)
	, SurfaceFinish(surfaceFinish)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, Dimensions(scintillator_dimensions)
	, ScintillatorName(scintillator_name)
	, CutDimple(cut_dimple)
	, DimpleRadius(dimple_radius)
	, DimpleCentrePosition(dimple_centre_relative_to_scintillator_centre)
	, TwinTile(twin_tile)
	{
		if(reference_volume != MotherVolume_physical && reference_volume->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		if(TwinTile && twin_tile_properties_file == "")
		{
			G4cerr << "##########" << "# WARNING:" << "# A scintillator twin tile is to be build, but no properties file is specified. A normal tile will be build." << "##########" << G4endl;
			TwinTile = false;
		}

		// set default values
		SetDefaults();

		// the obligatory parameters
		ScintiProperties.load(scintillator_properties_file);
		if(TwinTile) TwinTileProperties.load(twin_tile_properties_file);


		G4RotationMatrix rotation = MotherVolume_physical->GetObjectRotationValue().inverse() * reference_volume->GetObjectRotationValue() * reference_transformation.getRotation() * transf_relative_to_reference.getRotation() * replication_transf.getRotation();

		// WARNING-NOTE: as transform() changes the vector it is applied to, the following 5 commands are needed instead of a 1 line command
		G4ThreeVector translation = replication_transf.getTranslation();
		translation.transform(transf_relative_to_reference.getRotation());
		translation += transf_relative_to_reference.getTranslation();
		translation.transform(reference_volume->GetObjectRotationValue());
		translation += reference_volume->GetObjectTranslation();
		translation.transform(reference_transformation.getRotation());
		translation += reference_transformation.getTranslation();
		translation -= MotherVolume_physical->GetObjectTranslation();
		translation.transform(MotherVolume_physical->GetObjectRotationValue().inverse());

		Transformation = G4Transform3D(rotation, translation);


		// create materials:
		DefineMaterials();

		// create the volumes
		ConstructScintiVolume();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

	/**
	 *  Destructor (empty)
	 */
	~G4ScintillatorTile()
	{
	}


//   #######################   //
//   #  Getter functions:  #   //
//   #######################   //
// basic name of the object:
	/** @return the G4ScintillatorTile%'s name */
	G4String GetScintillatorName()
	{ return ScintillatorName; }

// dimensions of the object:
	/** @return the G4ScintillatorTile%'s dimensions */
	G4ThreeVector GetScintillatorDimensions()
	{ return Dimensions; }

// position and orientation of the object:
	/** @return the transformation of the G4ScintillatorTile inside the mother volume */
	G4Transform3D GetScintillatorTransformation()
	{ return Transformation; }

// mother volume of the object:
	/** @return pointer to the G4ScintillatorTile%'s mother volume */
	G4VPhysicalVolume * GetMotherVolume_physicalVolume()
	{
		return MotherVolume_physical;
	}

// volumes of the object:
	/** @return pointer to the physical volume of the G4ScintillatorTile */
	G4VPhysicalVolume * GetScintillator_physicalVolume()
	{ return Scintillator_physical; }

	/**
	 *  @return <b> if the G4ScintillatorTile is a twin tile: </b> pointer to the physical volume of the twin tile's reflective layer
	 *  @return <b> else: </b> 0
	 */
	G4VPhysicalVolume * GetTwinTileReflector_physicalVolume()
	{
		if(TwinTileReflector_physical) return TwinTileReflector_physical;
		else return 0;
	}

	/** @return pointer to the logical volume of the G4ScintillatorTile */
	G4LogicalVolume * GetScintillator_logicalVolume()
	{ return GetScintillator_physicalVolume()->GetLogicalVolume(); }

	/**
	 *  @return <b> if the G4ScintillatorTile is a twin tile: </b> pointer to the logical volume of the twin tile's reflective layer
	 *  @return <b> else: </b> 0
	 */
	G4LogicalVolume * GetTwinTileReflector_logicalVolume()
	{
		if(GetTwinTileReflector_physicalVolume()) return GetTwinTileReflector_physicalVolume()->GetLogicalVolume();
		else return 0;
	}

	/** @return pointer to the solid volume of the G4ScintillatorTile */
	G4VSolid * GetScintillator_solidVolume()
	{ return GetScintillator_logicalVolume()->GetSolid(); }

	/** @return pointer to the basic solid volume of the G4ScintillatorTile, i.e. the solid volume before the dimple has been cut of the G4ScintillatorTile */
	G4Box * GetScintillator_basicSolidVolume()
	{ return (G4Box *) GetScintillator_solidVolume()->GetConstituentSolid(0); }

	/**
	 *  @return <b> if the G4ScintillatorTile is a twin tile: </b> pointer to the solid volume of the twin tile's reflective layer
	 *  @return <b> else: </b> 0
	 */
	G4Box * GetTwinTileReflector_solidVolume()
	{
		if(GetTwinTileReflector_logicalVolume()) return (G4Box *) GetTwinTileReflector_logicalVolume()->GetSolid();
		else return 0;
	}

// materials of the object:
	/** @return pointer to the material of the scintillator */
	G4Material * GetScintillatorMaterial()
	{ return Material_Scinti; }

	/**
	 *  @return <b> if the G4ScintillatorTile is a twin tile: </b> pointer to the material of the twin tile's reflective layer
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetTwinTileReflectorMaterial()
	{ return Material_TwinTile; }

	/**
	 *  @return <b> if the G4ScintillatorTile is a twin tile: </b> pointer to the material of the twin tile's air gap
	 *  @return <b> else: </b> 0
	 */
	G4Material * GetTwinTileAirGapMaterial()
	{ return Material_Air; }

// optical surfaces of the object:
	/** @return pointer to the optical surface of the scintillator */
	G4OpticalSurface * GetScintillatorOpticalSurface()
	{ return OptSurf_ScintiMother; }

	/**
	 *  @return <b> if the G4ScintillatorTile is a twin tile: </b> pointer to the optical surface of the twin tile's reflective layer
	 *  @return <b> else: </b> 0
	 */
	G4OpticalSurface * GetTwinTileReflectorOpticalSurface()
	{ return OptSurf_reflective; }

// sensitive detector:
	/** @return G4bool, if the G4ScintillatorTile has a sensitive detector */
	G4bool HasSensitiveDetector()
	{ return ConstructSensitiveDetector; }

	/** @return pointer to the sensitive detector */
	ScintillatorSensitiveDetector * GetSensitiveDetector()
	{ return ((ScintillatorSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("ScintillatorSD", false)); }

// other properties of the object:
	/** @return G4bool, if the G4ScintillatorTile has a dimple */
	G4bool HasDimple()
	{ return CutDimple; }

	/**
	 *  @return <b> if the G4ScintillatorTile has a dimple: </b> radius of the dimple
	 *  @return <b> else: </b> NAN
	 */
	G4double GetDimpleRadius()
	{
		if(HasDimple()) return DimpleRadius;
		else return NAN;
	}

	/**
	 *  @return <b> if the G4ScintillatorTile has a dimple: </b> position of the dimple's centre relative to the G4ScintillatorTile%'s centre
	 *  @return <b> else: </b> NAN
	 */
	G4ThreeVector GetDimpleCentrePosition()
	{
		if(HasDimple()) return DimpleCentrePosition;
		else return G4ThreeVector(NAN, NAN, NAN);
	}

	/** @return G4bool, if the G4ScintillatorTile is a twin tile */
	G4bool IsTwinTile()
	{ return TwinTile; }



private:
	void DefineMaterials();
	void DefineScintillatorMaterialProperties();
	void DefineTwinTileMaterialProperties();
	void DefineAirMaterialProperties();
	void DefineSurfacesProperties();

	void ConstructScintiVolume();
	void ConstructScintiSurface();

	void SetDefaults();



	G4bool CriticalErrorOccured;
	G4bool SearchOverlaps;
	G4bool ConstructSensitiveDetector;
	PropertyToolsManager * PropertyTools;
	G4OpticalSurfaceModel SurfaceModel;
	G4OpticalSurfaceFinish SurfaceFinish;
	GODDeSS_DataStorage * DataStorage;

	G4VPhysicalVolume * MotherVolume_physical;

// Materials & Elements
	G4Material * Material_Scinti;
	G4Material * Material_TwinTile;
	G4Material * Material_Air;

// Scintillator
	GoddessProperties ScintiProperties;
	G4ThreeVector Dimensions;
	G4String ScintillatorName;
	G4Transform3D Transformation;
	G4VPhysicalVolume * Scintillator_physical;
	G4VPhysicalVolume * TwinTileAirGap_physical;
	G4OpticalSurface * OptSurf_ScintiMother;

// Dimple
	G4bool CutDimple;
	G4double DimpleRadius;
	G4ThreeVector DimpleCentrePosition;

// TwinTile
	G4bool TwinTile;
	GoddessProperties TwinTileProperties;
	G4bool TwinTileReflectorIsMetal;
	G4VPhysicalVolume * TwinTileReflector_physical;
	G4OpticalSurface * OptSurf_reflective;
};

#endif
