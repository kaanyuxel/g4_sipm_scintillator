/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef G4WRAPPING_H
#define G4WRAPPING_H

#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VSolid.hh>
#include <G4VPhysicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4SurfaceProperty.hh>
#include <G4Material.hh>
#include <G4SDManager.hh>

#include <GoddessProperties.hh>
#include <PropertyToolsManager.hh>
#include <WrappingSensitiveDetector.hh>

#include <GODDeSS_DataStorage.hh>
#include <G4ScintillatorTile.hh>



// class variables begin with capital letters, local variables with small letters



///  a class generating the materials, volumes, and optical properties needed for a wrapping and allows to access them.
class G4Wrapping
{
public:

	/**
	 *  Constructor to construct a wrapping volume around a G4ScintillatorTile:
	 *  - sets class variables to default values
	 *  - loads the property file
	 *  - creates materials (DefineWrappingMaterials())
	 *  - constructs the volumes (ConstructWrappingVolumes())
	 */
	G4Wrapping( G4ScintillatorTile* scinti_to_be_wrapped,		/**< G4ScintillatorTile that is to be wrapped with the G4Wrapping */
		    G4String wrapping_property_file,			/**< path to the file containing the G4Wrapping%'s properties */
		    G4VPhysicalVolume* mother_volume,			/**< G4Wrapping%'s mother volume */
		    G4String wrapping_name,				/**< name of the G4Wrapping%'s volume (it will be extended to distinguish between different G4Wrapping%s and different volumes of one G4Wrapping) */
		    std::vector<G4VPhysicalVolume *> cut_volumes,	/**< volumes that is to be cut out of the G4Wrapping */
		    G4bool constructSensitiveDetector,			/**< a sensitive detector is to be constructed ("true" or "false") */
		    G4bool searchOverlaps,				/**< Geant should search for overlaps when placing the physical volumes of G4Fibre%'s ("true" or "false") */
		    PropertyToolsManager * propertyTools,		/**< pointer to the PropertyToolsManager that is to be used */
		    G4OpticalSurfaceModel surfaceModel,			/**< model for the optical surfaces */
		    G4OpticalSurfaceFinish surfaceFinish,		/**< finish for the optical surfaces */
		    GODDeSS_DataStorage * dataStorage			/**< pointer to the GODDeSS_DataStorage that is to be used */
		  )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, ConstructSensitiveDetector(constructSensitiveDetector)
	, PropertyTools(propertyTools)
	, SurfaceModel(surfaceModel)
	, SurfaceFinish(surfaceFinish)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ScintiToBeWrapped_physical(scinti_to_be_wrapped)
	, WrappingName(wrapping_name)
	, CutVolumes(cut_volumes)
	{
		// set default values
		SetDefaults();

		if(!MotherVolume_physical) MotherVolume_physical = ScintiToBeWrapped_physical->GetMotherVolume_physicalVolume();

		WrappingProperties.load(wrapping_property_file);
		Dimensions = ScintiToBeWrapped_physical->GetScintillatorDimensions();
		Transformation = ScintiToBeWrapped_physical->GetScintillatorTransformation();
		AirGapThickness = WrappingProperties.getNumber("thick_air");
		if(isnan(AirGapThickness)) AirGapThickness = 0.;

		// does a second wrapping layer exist?
		if(WrappingProperties.containsNumber("thick_layer_2"))
		{
			Wrapping2Exists = true;
			Wrapping2Thickness = WrappingProperties.getNumber("thick_layer_2");
			if(WrappingProperties.containsString("is_metal_layer_2")) Wrapping2IsMetal = (WrappingProperties.getString("is_metal_layer_2") == "true");
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# WARNING:" <<
				std::endl << "# G4Wrapping():" <<
				std::endl << "# It is not specified if the second wrapping layer is metal-like. It will be set to false." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				Wrapping2IsMetal = false;
			}
		}
		else Wrapping2Exists = false;

		// (first) wrapping layer
		if(WrappingProperties.containsNumber("thick_layer_1")) Wrapping1Thickness = WrappingProperties.getNumber("thick_layer_1");
		else if(WrappingProperties.containsNumber("thick")) Wrapping1Thickness = WrappingProperties.getNumber("thick");
		else if(Wrapping2Exists)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping():" <<
			std::endl << "# The first wrapping layer's thickness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping():" <<
			std::endl << "# The wrapping layer's thickness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(WrappingProperties.containsString("is_metal_layer_1")) Wrapping1IsMetal = (WrappingProperties.getString("is_metal_layer_1") == "true");
		else if(WrappingProperties.containsString("is_metal")) Wrapping1IsMetal = (WrappingProperties.getString("is_metal") == "true");
		else
		{
			G4String counterString = "";
			if(Wrapping2Exists) counterString = " first";

			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Wrapping():" <<
			std::endl << "# It is not specified if the" << counterString << " wrapping layer is metal-like. It will be set to false." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			Wrapping1IsMetal = false;
		}

		// create materials:
		DefineWrappingMaterials();

		// create the volumes
		ConstructWrappingVolumes();

		// abort the whole programme, if a critical error occured
		if(CriticalErrorOccured) exit(1);
	}

	/**
	 *  Destructor (empty)
	 */
	~G4Wrapping()
	{
	}



	static void ConstructLUTWrapping( G4ScintillatorTile* scinti_to_be_wrapped,	/**< G4ScintillatorTile that is to be wrapped with the G4Wrapping */
					  G4OpticalSurfaceModel surfaceModel,		/**< model for the optical surfaces */
					  G4OpticalSurfaceFinish surfaceFinish,		/**< finish for the optical surfaces */
					  G4SurfaceType surfaceType,			/**< type of the optical surfaces */
					  G4double roughness,				/**< roughness of the wrapping surface (the parameter defines the standard deviation of the Gaussian distribution of micro-facets' normals around the average surface normal in degree or radian) */
					  G4VPhysicalVolume* surrounding_volume = 0	/**< volume sharing the surface with the G4ScintillatorTile that is to be wrapped (if 0, the whole G4ScintillatorTile%'s surface will be wrapped, otherwise only the shared surface will be wrapped) */
					);


//   #######################   //
//   #  Getter functions:  #   //
//   #######################   //
// basic name of the object:
	/** @return the G4Wrapping%'s name */
	G4String GetWrappingName()
	{ return WrappingName; }

// dimensions of the object:
	/** @return the G4Wrapping%'s dimensions */
	G4ThreeVector GetWrappingDimensions()
	{ return Dimensions; }

// position and orientation of the object:
	/** @return the transformation of the G4Wrapping inside the mother volume */
	G4Transform3D GetWrappingTransformation()
	{ return Transformation; }

// mother volume of the object:
	/** @return pointer to the G4Wrapping%'s mother volume */
	G4VPhysicalVolume * GetMotherVolume_physicalVolume()
	{
		return MotherVolume_physical;
	}

// volumes of the object:
	/**
	 *  @return <b> if an outer wrapping layer exists: </b> pointer to its physical volume
	 *  @return <b> else: </b> pointer to the physical volume of the inner wrapping layer
	 */
	G4VPhysicalVolume * GetWrapping_physicalVolume()
	{
		if(Wrapping2Exists) return Wrapping2_physical;
		else return Wrapping1_physical;
	}

	/** @return pointer to the physical volume of the inner wrapping layer */
	G4VPhysicalVolume * GetInnerWrappingLayer_physicalVolume()
	{ return Wrapping1_physical; }

	/**
	 *  @return <b> if an outer wrapping layer exists: </b> pointer to its logical volume
	 *  @return <b> else: </b> pointer to the logical volume of the inner wrapping layer
	 */
	G4LogicalVolume * GetWrapping_logicalVolume()
	{ return GetWrapping_physicalVolume()->GetLogicalVolume(); }

	/** @return pointer to the logical volume of the inner wrapping layer */
	G4LogicalVolume * GetInnerWrapping_layerLogicalVolume()
	{ return GetInnerWrappingLayer_physicalVolume()->GetLogicalVolume(); }

	/**
	 *  @return <b> if an outer wrapping layer exists: </b> pointer to its solid volume
	 *  @return <b> else: </b> pointer to the solid volume of the inner wrapping layer
	 */
	G4Box * GetWrapping_solidVolume()
	{ return (G4Box *) GetWrapping_logicalVolume()->GetSolid(); }

	/** @return pointer to the solid volume of the inner wrapping layer */
	G4Box * GetInnerWrappingLayer_solidVolume()
	{ return (G4Box *) GetInnerWrapping_layerLogicalVolume()->GetSolid(); }

// materials of the object:
	/**
	 *  @return <b> if an outer wrapping layer exists: </b> pointer to its material
	 *  @return <b> else: </b> pointer to the material of the inner wrapping layer
	 */
	G4Material * GetWrappingMaterial()
	{
		if(Wrapping2Exists) return Material_Wrapping2;
		else return Material_Wrapping1;
	}

	/** @return pointer to the material of the inner wrapping layer */
	G4Material * GetInnerWrappingLayerMaterial()
	{ return Material_Wrapping1; }

// optical surfaces of the object:
	/** @return pointer to the optical surface of the inner wrapping layer */
	G4OpticalSurface * GetWrappingInnerOpticalSurface()
	{ return OptSurf_MotherWrapping1; }

	/**
	 *  @return <b> if an outer wrapping layer exists: </b> pointer to its optical surface
	 *  @return <b> else: </b> pointer to the optical surface of the inner wrapping layer
	 */
	G4OpticalSurface * GetWrappingOuterOpticalSurface()
	{
		if(Wrapping2Exists) return OptSurf_MotherWrapping2;
		else return OptSurf_MotherWrapping1;
	}

	/**
	 *  @return <b> if an outer wrapping layer exists: </b> pointer to the optical surface between both wrapping layers
	 *  @return <b> else: </b> 0
	 */
	G4OpticalSurface * GetWrappingInternalOpticalSurface()
	{ return OptSurf_Wrapping1Wrapping2; }

// sensitive detector:
	/** @return G4bool, if the G4Wrapping has a sensitive detector */
	G4bool HasSensitiveDetector()
	{ return ConstructSensitiveDetector; }

	/** @return pointer to the sensitive detector */
	WrappingSensitiveDetector * GetSensitiveDetector()
	{ return ((WrappingSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("WrappingSD", false)); }

// other properties of the object:
	/** @return thickness of the air gap */
	G4double GetAirGapThickness()
	{ return AirGapThickness; }

	/**
	 *  @return <b> if an outer wrapping layer exists: </b> its thickness
	 *  @return <b> else: </b> thickness of the inner wrapping layer
	 */
	G4double GetWrappingThickness()
	{
		if(Wrapping2Exists) return Wrapping2Thickness;
		else return Wrapping1Thickness;
	}

	/** @return thickness of the inner wrapping layer */
	G4double GetInnerWrappingLayerThickness()
	{ return Wrapping1Thickness; }

private:
	void DefineWrappingMaterials();
	void DefineWrapping1MaterialProperties();
	void DefineWrapping2MaterialProperties();
	void DefineWrappingSurfacesProperties();

	void ConstructWrappingVolumes();
	void ConstructWrappingSurface();

	void SetDefaults();



	G4bool CriticalErrorOccured;
	G4bool SearchOverlaps;
	G4bool ConstructSensitiveDetector;
	PropertyToolsManager * PropertyTools;
	G4OpticalSurfaceModel SurfaceModel;
	G4OpticalSurfaceFinish SurfaceFinish;
	GODDeSS_DataStorage * DataStorage;

	G4VPhysicalVolume * MotherVolume_physical;
	G4ScintillatorTile * ScintiToBeWrapped_physical;

// Materials & Elements
	G4Material* Material_Wrapping1;
	G4Material* Material_Wrapping2;

// Wrapping
	G4bool Wrapping2Exists;
	G4bool Wrapping1IsMetal;
	G4bool Wrapping2IsMetal;
	G4double AirGapThickness;
	G4double Wrapping1Thickness;
	G4double Wrapping2Thickness;
	GoddessProperties WrappingProperties;
	G4ThreeVector Dimensions;
	G4Transform3D Transformation;
	G4String WrappingName;
	std::vector<G4VPhysicalVolume *> CutVolumes;
	G4VPhysicalVolume * Wrapping1_physical;
	G4VPhysicalVolume * Wrapping2_physical;
	G4OpticalSurface * OptSurf_MotherWrapping1;
	G4OpticalSurface * OptSurf_Wrapping1Wrapping2;
	G4OpticalSurface * OptSurf_MotherWrapping2;
};

#endif
