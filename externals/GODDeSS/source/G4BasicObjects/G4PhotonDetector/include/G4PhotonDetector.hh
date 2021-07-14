/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef G4PHOTONDETECTOR_H
#define G4PHOTONDETECTOR_H

#include <G4Box.hh>
#include <G4Transform3D.hh>
#include <G4OpticalSurface.hh>

#include <PropertyToolsManager.hh>

#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class generating the materials, volumes, and optical properties needed for a photon detector and allows to access them.
class G4PhotonDetector
{
public:

	/**
	 *  Constructor to construct a photon detector:
	 *  - sets class variables to default values
	 *  - sets the variables for the photon detector's placement (considering the transformation of the photon detector relative to the reference volume and the transformation of the reference volume relative to the mother volume)
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructScintiVolume())
	 */
	G4PhotonDetector( G4double edge_length,							/**< edge length of the G4PhotonDetector%'s sensitive area */
			  G4VPhysicalVolume* mother_volume,					/**< G4PhotonDetector%'s mother volume */
			  G4String photon_detector_name,					/**< name of the G4PhotonDetector%'s volume (it will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector) */
			  G4ThreeVector sensitive_surface_normal_relative_to_reference,		/**< surface normal of the G4PhotonDetector%'s sensitive area relative to the reference volume */
			  G4ThreeVector sensitive_surface_position_relative_to_reference,	/**< position of the G4PhotonDetector%'s sensitive area relative to the reference volume */
			  G4VPhysicalVolume* reference_volume,					/**< reference volume relative to which the G4PhotonDetector will be orientated (if 0, the mother volume is the reference volume) */
			  G4Transform3D reference_transformation,				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
			  G4bool searchOverlaps,							/**< Geant should search for overlaps when placing the physical volumes of G4PhotonDetector%'s ("true" or "false") */
			  PropertyToolsManager * propertyTools,					/**< pointer to the PropertyToolsManager that is to be used */
			  GODDeSS_DataStorage * dataStorage					/**< pointer to the GODDeSS_DataStorage that is to be used */
			)
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, EdgeLength(edge_length)
	, PhotonDetectorName(photon_detector_name)
	, SurfaceNormal_rel(sensitive_surface_normal_relative_to_reference / sensitive_surface_normal_relative_to_reference.mag())
	, SurfacePos_rel(sensitive_surface_position_relative_to_reference)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		// the obligatory parameters
		Depth = 0.5 * CLHEP::mm;
		CoatingWidth = 0.005 * CLHEP::mm;

		// Determine rotation parameters for given surface normal.
		G4ThreeVector originalSurfaceNormal(0, 0, 1);
		G4ThreeVector rotationAxis = originalSurfaceNormal.cross(SurfaceNormal_rel);
		G4double rotationAngle = 0.;
		if(fabs(rotationAxis.mag()) > 1e-12)
		{
			rotationAngle = acos(originalSurfaceNormal * SurfaceNormal_rel / (originalSurfaceNormal.mag() * SurfaceNormal_rel.mag()));
		}
		else if(SurfaceNormal_rel == -originalSurfaceNormal)
		{
			rotationAngle = 180. * CLHEP::deg;
			rotationAxis = G4ThreeVector(1, 0, 0);
		}

		// calculate the transformation
		G4RotationMatrix photonDetectorRotation = G4RotationMatrix();
		photonDetectorRotation.rotate(rotationAngle, rotationAxis);

		// considering the transformation of the reference volume
		photonDetectorRotation = MotherVolume_physical->GetObjectRotationValue().inverse() * ReferenceVolume_physical->GetObjectRotationValue() * reference_transformation.getRotation() * photonDetectorRotation;

		// WARNING-NOTE: as transform() changes the vector it is applied to, the following 7 commands are needed instead of a 1 line command
		G4ThreeVector photonDetectorTranslation = SurfacePos_rel - SurfaceNormal_rel * (Depth + CoatingWidth) / 2.;
		photonDetectorTranslation.transform(ReferenceVolume_physical->GetObjectRotationValue());
		photonDetectorTranslation += ReferenceVolume_physical->GetObjectTranslation();
		photonDetectorTranslation.transform(reference_transformation.getRotation());
		photonDetectorTranslation += reference_transformation.getTranslation();
		photonDetectorTranslation -= MotherVolume_physical->GetObjectTranslation();
		photonDetectorTranslation.transform(MotherVolume_physical->GetObjectRotationValue().inverse());

		Transformation = G4Transform3D(photonDetectorRotation, photonDetectorTranslation);


		// create materials:
		DefineMaterials();

		// create the volumes
		ConstructVolumes();
	}

// option to create replications: the first object will be placed at transf_relative_to_reference, the other number_of_replications - 1 objects will be placed at replication_transf relative to the previous one
	/**
	 *  Constructor to construct replications of a photon detector:
	 *  - sets class variables to default values
	 *  - sets the variables for the photon detector's placement (considering the transformation of the photon detector relative to the reference volume and the transformation of the reference volume relative to the mother volume)
	 *  - creates materials (DefineMaterials())
	 *  - constructs the volumes (ConstructScintiVolume())
	 */
	G4PhotonDetector( G4double edge_length,							/**< edge length of the G4PhotonDetector%'s sensitive area */
			  G4VPhysicalVolume* mother_volume,					/**< G4PhotonDetector%'s mother volume */
			  G4String photon_detector_name,					/**< name of the G4PhotonDetector%'s volume (it will be extended to distinguish between different G4PhotonDetector%s and different volumes of one G4PhotonDetector) */
			  G4ThreeVector sensitive_surface_normal_relative_to_reference,		/**< surface normal of the G4PhotonDetector%'s sensitive area relative to the reference volume */
			  G4ThreeVector sensitive_surface_position_relative_to_reference,	/**< position of the G4PhotonDetector%'s sensitive area relative to the reference volume */
			  G4VPhysicalVolume* reference_volume,					/**< reference volume relative to which the G4PhotonDetector will be orientated (if 0, the mother volume is the reference volume) */
			  G4Transform3D reference_transformation,				/**< transformation between mother volume and reference volume's mother volume (necessary if the reference volume is located in a different volume) */
			  G4bool searchOverlaps,							/**< Geant should search for overlaps when placing the physical volumes of G4PhotonDetector%'s ("true" or "false") */
			  PropertyToolsManager * propertyTools,					/**< pointer to the PropertyToolsManager that is to be used */
			  GODDeSS_DataStorage * dataStorage,					/**< pointer to the GODDeSS_DataStorage that is to be used */
			  G4Transform3D replication_transf					/**< transformation of each replication, relative to the previous one */
			)
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: SearchOverlaps(searchOverlaps)
	, PropertyTools(propertyTools)
	, DataStorage(dataStorage)
	, MotherVolume_physical(mother_volume)
	, ReferenceVolume_physical(reference_volume)
	, EdgeLength(edge_length)
	, PhotonDetectorName(photon_detector_name)
	, SurfaceNormal_rel(sensitive_surface_normal_relative_to_reference / sensitive_surface_normal_relative_to_reference.mag())
	, SurfacePos_rel(sensitive_surface_position_relative_to_reference)
	{
		if(ReferenceVolume_physical != MotherVolume_physical && ReferenceVolume_physical->GetMotherLogical() != MotherVolume_physical->GetLogicalVolume() && reference_transformation == G4Transform3D()) G4cerr << G4endl << "The reference volume is neither located in the mother volume nor is it the mother volume. This might lead to geometry problems, if no transformation between mother volume and reference volume's mother volume is defined!" << G4endl << G4endl;

		// set default values
		SetDefaults();

		// the obligatory parameters
		Depth = 0.5 * CLHEP::mm;
		CoatingWidth = 0.005 * CLHEP::mm;

		// Determine rotation parameters for given surface normal.
		G4ThreeVector originalSurfaceNormal(0, 0, 1);
		G4ThreeVector rotationAxis = originalSurfaceNormal.cross(SurfaceNormal_rel);
		G4double rotationAngle = 0.;
		if(fabs(rotationAxis.mag()) > 1e-12)
		{
			rotationAngle = acos(originalSurfaceNormal * SurfaceNormal_rel / (originalSurfaceNormal.mag() * SurfaceNormal_rel.mag()));
		}
		else if(SurfaceNormal_rel == -originalSurfaceNormal)
		{
			rotationAngle = 180. * CLHEP::deg;
			rotationAxis = G4ThreeVector(1, 0, 0);
		}

		// calculate the transformation
		G4RotationMatrix photonDetectorRotation = G4RotationMatrix();
		photonDetectorRotation.rotate(rotationAngle, rotationAxis);

		// considering the transformation of the reference volume
		photonDetectorRotation = MotherVolume_physical->GetObjectRotationValue().inverse() * ReferenceVolume_physical->GetObjectRotationValue() * reference_transformation.getRotation() * photonDetectorRotation * replication_transf.getRotation();

		// WARNING-NOTE: as transform() changes the vector it is applied to, the following 9 commands are needed instead of a 1 line command
		G4ThreeVector photonDetectorTranslation = replication_transf.getTranslation();
		photonDetectorTranslation.transform(photonDetectorRotation);
		photonDetectorTranslation += SurfacePos_rel - SurfaceNormal_rel * (Depth + CoatingWidth) / 2.;
		photonDetectorTranslation.transform(ReferenceVolume_physical->GetObjectRotationValue());
		photonDetectorTranslation += ReferenceVolume_physical->GetObjectTranslation();
		photonDetectorTranslation.transform(reference_transformation.getRotation());
		photonDetectorTranslation += reference_transformation.getTranslation();
		photonDetectorTranslation -= MotherVolume_physical->GetObjectTranslation();
		photonDetectorTranslation.transform(MotherVolume_physical->GetObjectRotationValue().inverse());

		Transformation = G4Transform3D(photonDetectorRotation, photonDetectorTranslation);


		// create materials:
		DefineMaterials();

		// create the volumes
		ConstructVolumes();
	}

	/**
	 *  Destructor (empty)
	 */
	~G4PhotonDetector()
	{
	}


//   #######################   //
//   #  Getter functions:  #   //
//   #######################   //
// basic name of the object:
	/** @return the G4PhotonDetector%'s name */
	G4String GetPhotonDetectorName()
	{ return PhotonDetectorName; }

// dimensions of the object:
	/** @return the dimensions of the G4PhotonDetector%'s coating */
	G4ThreeVector GetCoatingDimensions()
	{ return G4ThreeVector(EdgeLength + 2. * CoatingWidth, EdgeLength + 2. * CoatingWidth, Depth + CoatingWidth); }

	/** @return the dimensions of the G4PhotonDetector%'s sensitive volume */
	G4ThreeVector GetSensitiveVolumeDimensions()
	{ return G4ThreeVector(EdgeLength, EdgeLength, Depth); }

// position and orientation of the object:
	/** @return the transformation of the G4PhotonDetector inside the mother volume */
	G4Transform3D GetPhotonDetectorTransformation()
	{ return Transformation; }

	/** @return the surface normal of the G4PhotonDetector%'s sensitive area relative to the base volume */
	G4ThreeVector GetSensitiveVolumeRelativeSurfaceNormal()
	{ return SurfaceNormal_rel; }

	/** @return the G4PhotonDetector%'s sensitive area position relative to the base volume */
	G4ThreeVector GetSensitiveVolumeRelativeSurfacePosition()
	{ return SurfacePos_rel; }

// mother volume of the object:
	/** @return pointer to the G4PhotonDetector%'s mother volume */
	G4VPhysicalVolume * GetMotherVolume_physicalVolume()
	{
		return MotherVolume_physical;
	}

// volumes of the object:
	/** @return pointer to the physical volume of the G4PhotonDetector%'s coating */
	G4VPhysicalVolume * GetCoating_physicalVolume()
	{ return PhotonDetector_physical; }

	/** @return pointer to the physical volume of the G4PhotonDetector%'s sensitive area */
	G4VPhysicalVolume * GetSensitiveVolume_physicalVolume()
	{ return SensitiveVolume_physical; }

	/** @return pointer to the logical volume of the G4PhotonDetector%'s coating */
	G4LogicalVolume * GetCoating_logicalVolume()
	{ return GetCoating_physicalVolume()->GetLogicalVolume(); }

	/** @return pointer to the logical volume of the G4PhotonDetector%'s sensitive area */
	G4LogicalVolume * GetSensitiveVolume_logicalVolume()
	{ return GetSensitiveVolume_physicalVolume()->GetLogicalVolume(); }

	/** @return pointer to the solid volume of the G4PhotonDetector%'s coating */
	G4Box * GetCoating_solidVolume()
	{ return (G4Box *) GetCoating_logicalVolume()->GetSolid(); }

	/** @return pointer to the solid volume of the G4PhotonDetector%'s sensitive area */
	G4Box * GetSensitiveVolume_solidVolume()
	{ return (G4Box *) GetSensitiveVolume_logicalVolume()->GetSolid(); }
// materials of the object:
	/** @return pointer to the material of the G4PhotonDetector%'s coating */
	G4Material * GetPhotonDetectorMaterial()
	{ return Material_PhotonDetector; }

// optical surfaces of the object:
	/** @return pointer to the optical surface of the G4PhotonDetector%'s coating */
	G4OpticalSurface * GetCoatingOpticalSurface()
	{ return OptSurf_photonDetector; }

	/** @return pointer to the optical surface of the G4PhotonDetector%'s sensitive area */
	G4OpticalSurface * GetSensitiveVolumeOpticalSurface()
	{ return OptSurf_sensitiveVolume; }

private:
	void DefineMaterials();
	void DefineMaterialProperties();
	void DefineSurfacesProperties();

	void ConstructVolumes();
	void ConstructSurface();

	void SetDefaults();



	G4bool SearchOverlaps;
	PropertyToolsManager * PropertyTools;
	GODDeSS_DataStorage * DataStorage;
	G4VPhysicalVolume * MotherVolume_physical;
	G4VPhysicalVolume * ReferenceVolume_physical;

// Materials & Elements
	G4Material* Material_PhotonDetector;

// Photon Detector
	G4double EdgeLength;
	G4double Depth;
	G4double CoatingWidth;
	G4String PhotonDetectorName;
	G4ThreeVector SurfaceNormal_rel;
	G4ThreeVector SurfacePos_rel;
	G4Transform3D Transformation;
	G4VPhysicalVolume * PhotonDetector_physical;
	G4VPhysicalVolume * SensitiveVolume_physical;
	G4OpticalSurface * OptSurf_photonDetector;
	G4OpticalSurface * OptSurf_sensitiveVolume;
};

#endif
