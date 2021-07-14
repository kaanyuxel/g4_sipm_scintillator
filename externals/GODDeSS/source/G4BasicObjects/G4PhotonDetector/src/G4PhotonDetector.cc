/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4PVPlacement.hh>
#include <G4VisAttributes.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4SDManager.hh>

#include <PropertyToolsManager.hh>
#include <PhotonSensitiveDetector.hh>

#include <G4PhotonDetector.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to construct all volumes and surfaces for photon detectors.
 */
void G4PhotonDetector::ConstructVolumes()
{
// solids (dimensions):
	G4String photon_detector_solid_name = PhotonDetectorName + "_coating_solid";
	G4VSolid * photonDetector_solid = new G4Box(photon_detector_solid_name, EdgeLength / 2. + CoatingWidth, EdgeLength / 2. + CoatingWidth, (Depth + CoatingWidth) / 2.);

	G4String sensitive_volume_solid_name = PhotonDetectorName + "_solid";
	G4VSolid * sensitiveVolume_solid = new G4Box(sensitive_volume_solid_name, EdgeLength / 2., EdgeLength / 2., Depth / 2.);

// logical volumes (material):
	G4String photon_detector_logical_name = PhotonDetectorName + "_coating_logical";
	G4LogicalVolume * photonDetector_logical = new G4LogicalVolume(photonDetector_solid, Material_PhotonDetector, photon_detector_logical_name, 0, 0, 0);

	G4VisAttributes photonDetectorVisAtt(G4Colour::Black());
	photonDetectorVisAtt.SetForceWireframe(true);
	photonDetectorVisAtt.SetLineWidth(3.);
	photonDetector_logical->SetVisAttributes(photonDetectorVisAtt);

	G4String sensitive_volume_logical_name = PhotonDetectorName + "_logical";
	G4LogicalVolume * sensitiveVolume_logical = new G4LogicalVolume(sensitiveVolume_solid, Material_PhotonDetector, sensitive_volume_logical_name, 0, 0, 0);

	G4VisAttributes sensitiveVolumeVisAtt(G4Colour::Yellow());
	sensitiveVolumeVisAtt.SetForceSolid(true);
	sensitiveVolume_logical->SetVisAttributes(sensitiveVolumeVisAtt);

	// try to find an already existing senitive detector
	PhotonSensitiveDetector * sensitiveDetector = (PhotonSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("PhotonDetectorSD", false);
	if(!sensitiveDetector)
	{
		// create a new senitive detector
		sensitiveDetector = new PhotonSensitiveDetector("PhotonDetectorSD", DataStorage);
		// register it to Geant4
		G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDetector);
	}
	// assign it to detector parts
	sensitiveVolume_logical->SetSensitiveDetector(sensitiveDetector);
	DataStorage->AddPhotonSensitiveDetectorVolumeName(PhotonDetectorName);

// physical volumes (placement):
//NOTE: G4PVPlacement puts the volume to be placed into EVERY physical volume emanating from the same logical volume (no matter whether the logical or physical volume is specified as mother volume)!
//      => If G4PVPlacement is to be able to distinguish physical volumes, for each physical volume a separate logical volume has to be created.
//      => rule of thumb: For each volume that might become a mother volume, a separate logical volume should to be created.
	G4String photon_detector_physical_name = PhotonDetectorName + "_coating";
	PhotonDetector_physical = new G4PVPlacement(Transformation, photon_detector_physical_name, photonDetector_logical, MotherVolume_physical, false, 0, SearchOverlaps);

	SensitiveVolume_physical = new G4PVPlacement(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., CoatingWidth / 2.)), PhotonDetectorName, sensitiveVolume_logical, PhotonDetector_physical, false, 0, SearchOverlaps);

// surfaces:
//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.
	ConstructSurface();
}



/**
 *  Function to define materials (compounds, alloys) and their properties:
 *  - creates new materials according to the specifications from the property file(s)
 *  - if the same materials already exist, the newly created ones are deleted
 */
void G4PhotonDetector::DefineMaterials()
{
	// NOTE The materials which have already been define (e.g. in the DetectorConstruction) are used here

	// G4Material(const G4String &name, G4double z, G4double a, G4double density, G4State state=kStateUndefined, G4double temp=STP_Temperature, G4double pressure=STP_Pressure)
	Material_PhotonDetector = new G4Material("tempName", 1., 1.01 * CLHEP::g / CLHEP::mole, CLHEP::universe_mean_density, kStateGas, 0.1 * CLHEP::kelvin, 1.e-19 * CLHEP::pascal);		// definition from Geant

	DefineMaterialProperties();

	PropertyTools->checkIfMaterialAlreadyExists("PhotonDetector", Material_PhotonDetector);
}



//NOTE: possible Properties:
//         "RINDEX":			(spectrum (in dependence of the photon energy))		(obligatory property!)
//		defines the refraction index of the material, used for boundary processes, Cerenkov radiation and Rayleigh scattering
//         "ABSLENGTH":			(spectrum (in dependence of the photon energy))
//		defines the absorption length (absorption spectrum) of the material, used for the "normal" absorption of optical photons (default is infinity, i.e. no absorption)
//		(the absorption length for the WLS process of WLS materials is specified by "WLSABSLENGTH", "ABSLENGTH" can be specified additionally to simulate a non-WLS absorption fraction)
//         "RAYLEIGH":			(spectrum (in dependence of the photon energy))
//		defines the absorption length of the material, used for the rayleigh scattering of optical photons (default is infinity, i.e. no scattering)
//
//         "SCINTILLATIONYIELD":	(constant value (energy independent))			(obligatory property for scintillator materials!)
//		defines the mean number of photons, emitted per MeV energy deposition in the scintillator material (the real number is Poisson/Gauss distributed)
//		(can also be specified separately for different particles by putting "ELECTRON...", "PROTON...", "DEUTERON...", "TRITON...", "ALPHA...", "ION..." infront of "SCINTILLATIONYIELD")
//		(default is 0, i.e. no scintillation process)
//         "RESOLUTIONSCALE":		(constant value (energy independent))
//		defines the intrinsic resolution of the scintillator material, used for the statistical distribution of the number of generated photons in the scintillation process
//		(values > 1 result in a wider distribution, values < 1 result in a narrower distribution -> 1 is to be chosen as default)
//		(default is 0)
//         "FASTCOMPONENT":		(spectrum (in dependence of the photon energy))		(at least one "...COMPONENT" is obligatory for scintillator materials!)
//		defines the emission spectrum of the material, used for the fast scintillation process	NOTE: emission spectra are NOT linearly extrapolated between two given points!
//         "SLOWCOMPONENT":		(spectrum (in dependence of the photon energy))		(at least one "...COMPONENT" is obligatory for scintillator materials!)
//		defines the emission spectrum of the material, used for the slow scintillation process	NOTE: emission spectra are NOT linearly extrapolated between two given points!
//         "FASTTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between energy deposition and photon emission), used for the fast scintillation process (default is 0)
//         "SLOWTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between energy deposition and photon emission), used for the slow scintillation process (default is 0)
//         "FASTSCINTILLATIONRISETIME":	(constant value (energy independent))
//		defines the rise time (time between the start of the emission and the emission peak), used for the fast scintillation process (default is 0)
//         "SLOWSCINTILLATIONRISETIME":	(constant value (energy independent))
//		defines the rise time (time between the start of the emission and the emission peak), used for the slow scintillation process (default is 0)
//         "YIELDRATIO":		(constant value (energy independent))			(obligatory property for scintillator materials, if both "...COMPONENT"s are specified!)
//		defines relative strength of the fast scintillation process as a fraction of total scintillation yield (default is 0)
//
//         "WLSABSLENGTH":		(spectrum (in dependence of the photon energy))		(obligatory property for WLS materials!)
//		defines the absorption length (absorption spectrum) of the material, used for the WLS process (default is infinity, i.e. no WLS process)
//         "WLSCOMPONENT":		(spectrum (in dependence of the photon energy))		(obligatory property for WLS materials!)
//		defines the emission spectrum of the material, used for the WLS process	NOTE: emission spectra are NOT linearly extrapolated between two given points!
//         "WLSTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between absorption and emission), used for the WLS process (default is 0)
//         "WLSMEANNUMBERPHOTONS":	(constant value (energy independent))
//		defines the mean number of photons, emitted for each photon that was absorbed by the WLS material
//		(if specified, the real number of emitted photons is Poisson distributed, else the real number of emitted photons is 1)

//**
// *  Function to define the following material properties of the photon detector (according to the specifications from the property file):
// *  - refractive index spectrum
// */
void G4PhotonDetector::DefineMaterialProperties()
{
	// RINDEX has to be set and RINDEX = 1 should ensure that there is no Cerenkov radiation
	G4MaterialPropertyVector * refractiveIndex_PhotonDetector = PropertyTools->GetPropertyDistribution(1.);

	G4MaterialPropertiesTable * MPT_PhotonDetector = new G4MaterialPropertiesTable();
	MPT_PhotonDetector->AddProperty("RINDEX", refractiveIndex_PhotonDetector);
	Material_PhotonDetector->SetMaterialPropertiesTable(MPT_PhotonDetector);
}



/**
 *  Function to construct surfaces between different volumes and define their mechanical and optical properties (according to the specifications from the property file(s)). \n
 *  The following properties are defined:
 *  - surface type ("dielectric_metal" or "dielectric_dielectric")
 *  - reflectivity (DefineSurfacesProperties())
 */
void G4PhotonDetector::ConstructSurface()
{
//NOTE: !!!!!! In order to understand simulation of optical surface properties (and esspecially the relations between the different options and parameters),
//             you should definitely start with reading the diagram to be found at:
//             http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2012/05/23/23.20-70533-ces_in_geant4_revised.png
//             (in http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/445.html)
//             The following comments are only a excerpt and CANNOT compete with the precision of this diagram !!!!!!
//
//      possible reflection models:
//         glisur:  original and obsolete GEANT3 model
//         unified: unified model (provides a range of different reflection mechanisms, cf. "possible Properties" below)
//         LUT:     using the Look-Up-Table for the surface simulation
//
//      possible surface types:
//         dielectric_dielectric: if both materials are dielectric, i.e. the reflection and refraction probabilities are defined by the refractive indices via the Fresnel equations
//                                (the reflection probability CANNOT be defined by the material's reflectivity, cf. "REFLECTIVITY" below)
//         dielectric_metal:      if one material is dielectric and the other is metal-like (i.e. the photons can only be reflected or absorbed but not refracted)
//                                this surface type has to be used if the reflection probability is to be defined by the material's reflectivity
//                                (this does not work for dielectric_dielectric, cf. "REFLECTIVITY" below)
//         dielectric_LUT:        if Look-Up-Tables is to be used
//         firsov:                for Firsov Process (O.B. Firsov, “Reflection of fast ions from a dense medium at glancing angles”, Sov. Phys.-Docklady, vol. 11, No. 8, pp.732-733, 1967)
//         x_ray:                 for x-ray mirror process
//
//      possible surface finishs: (cf. http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2012/05/23/23.20-70533-ces_in_geant4_revised.png !!!)
//         for dielectric_metal:
//            polished:             perfectly polished surface
//                                  -> specular spike reflection
//                                     (requiring the reflectivity of the metal)
//            ground:               rough surface
//                                  -> specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the reflectivity of the metal and the angle alpha between a micro-facet normal and the average surface normal)
//
//         for dielectric_dielectric:
//            polished:             perfectly polished surface
//                                  -> specular spike reflection
//                                     (requiring the refraction indices of both materials)
//            ground:               rough surface
//                                  -> specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of both materials and the angle alpha between a micro-facet normal and the average surface normal)
//            polishedfrontpainted: the volume has a perfectly polished surface and an absorbing paint without air gap
//                                  -> specular spike reflection
//                                     (requiring the reflectivity of the paint)
//            groundfrontpainted:   the volume has a rough surface and an absorbing paint without air gap
//                                  -> (diffuse) Lambertian reflection
//                                     (requiring the reflectivity of the paint)
//            polishedbackpainted:  the volume has a rough surface and a perfectly polished, absorbing paint with air gap
//                                  -> volume-air surface: specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of the volume material and the surface as well as the angle alpha between a micro-facet normal and the average surface normal)
//                                  -> air-paint surface: specular spike reflection
//                                     (requiring the reflectivity of the paint)
//            groundbackpainted:    the volume has a rough surface and a rough, absorbing paint with air gap
//                                  -> volume-air surface: specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of the volume material and the surface as well as the angle alpha between a micro-facet normal and the average surface normal)
//                                  -> air-paint surface: (diffuse) Lambertian reflection
//                                     (requiring the reflectivity of the paint)
//
//         for the Look-Up-Table model:
//            polishedlumirrorair:
//            polishedlumirrorglue:
//            polishedair:
//            polishedteflonair:
//            polishedtioair:
//            polishedtyvekair:
//            polishedvm2000air:
//            polishedvm2000glue:
//            etchedlumirrorair:
//            etchedlumirrorglue:
//            etchedair:
//            etchedteflonair:
//            etchedtioair:
//            etchedtyvekair:
//            etchedvm2000air:
//            etchedvm2000glue:
//            groundlumirrorair:
//            groundlumirrorglue:
//            groundair:
//            groundteflonair:
//            groundtioair:
//            groundtyvekair:
//            groundvm2000air:
//            groundvm2000glue:
//
//      sigma alpha: defines the roughness of the surface (sigma alpha specifies the standard deviation of the distribution of the micro-facet normals in [rad])

	// create the surface an define its properties
	OptSurf_photonDetector = new G4OpticalSurface("coating surface", unified, polished, dielectric_metal, 0.);	// the 0. defines the sigma alpha
	OptSurf_sensitiveVolume = new G4OpticalSurface("sensitive volume surface", unified, polished, dielectric_dielectric, 0.);	// the 0. defines the sigma alpha

	// define the optical surface properties
	DefineSurfacesProperties();

//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.

	// specify between which volumes the previously defined surface is to be applied
	// the first volume is where the photon comes from, the second volume is where the photon heads for.
	new G4LogicalSkinSurface("photon detector surface", PhotonDetector_physical->GetLogicalVolume(), OptSurf_photonDetector);
	new G4LogicalSkinSurface("photon detector surface", SensitiveVolume_physical->GetLogicalVolume(), OptSurf_sensitiveVolume);
}

void G4PhotonDetector::DefineSurfacesProperties()
{
//NOTE: possible Properties (for reflection model "unified"):
//         "RINDEX":				(spectrum (in dependence of the photon energy))			(obligatory property for surfaces with surface finish "...backpainted"!)
//		in case of surface finish "...backpainted", defines the refraction index of the thin material layer between the surface and the paint
//		(only used for boundary processes and only in case of surface finish "...backpainted")
//         "REFLECTIVITY":			(spectrum (in dependence of the photon energy))
//		does NOT IN ANY CASE define the reflectivity of the surface but defines (1 - absorption coefficient), which may make a difference:
//			for dielectric_metal:      it is really the reflectivity of the metal, i.e. the probability that the photon is reflected at all
//						   (every photon not absorbed by the metal will be reflected)
//			for dielectric_dielectric: it is NOT the reflectivity but defines the absorption coefficient (absorption coefficient = 1 - "reflectivity") of a dirty surface
//						   (every photon not absorbed by dirt will be reflected or refracted as normal for "dielectric_dielectric" surfaces and the reflectivity
//						   for this process can be specified via "TRANSMITTANCE" or is calculated via the Fresnel equations, cf. "TRANSMITTANCE")
//		(default is 1., i.e. nothing is absorbed)
//         "REALRINDEX":			(spectrum (in dependence of the photon energy))
//		defines the real part of the complex refraction index of the surface   //FIXME surface <-> material ???
//		(in case that "REFLECTIVITY" is not specified and if "REALRINDEX" and "IMAGINARYRINDEX" are both specified,
//		the refectivity is calculated from the complex refraction index via the Fresnel equations   //FIXME shouldn't it be two refraction indices???
//		-> therefore, the complex refraction index should not be specified for surface type "dielectric_dielectric", cf. "REFLECTIVITY"
//		   (if one wants to define the reflectivity of a "dielectric_dielectric" surface from the complex refraction index,
//		   one has to calculate the transmittance via the Fresnel equations and specify it with "TRANSMITTANCE"))
//         "IMAGINARYRINDEX":			(spectrum (in dependence of the photon energy))
//		defines the imaginary part of the complex refraction index of the surface   //FIXME surface <-> material ???
//		(in case that "REFLECTIVITY" is not specified and if "REALRINDEX" and "IMAGINARYRINDEX" are both specified,
//		the refectivity is calculated from the complex refraction index via the Fresnel equations   //FIXME shouldn't it be two refraction indices???
//		-> therefore, the complex refraction index should not be specified for surface type "dielectric_dielectric", cf. "REFLECTIVITY"
//		   (if one wants to define the reflectivity of a "dielectric_dielectric" surface from the complex refraction index,
//		   one has to calculate the transmittance via the Fresnel equations and specify it with "TRANSMITTANCE"))
//         "EFFICIENCY":			(spectrum (in dependence of the photon energy))
//		defines the detection efficiency of absorbed photons (default is 0)
//         "TRANSMITTANCE":			(spectrum (in dependence of the photon energy))
//		in case of "dielectric_dielectric" surfaces, defines the fraction of photons (reaching the surface despite of dirt (cf. "REFLECTIVITY"))
//		which are refracted by the surface instead of being reflected
//		(if "TRANSMITTANCE" is not specified, the transmittance is calculated from the (real) refraction indices of the two materials forming the surface via the Fresnel equations)
//		(only used for boundary processes and only in case of surface type "dielectric_dielectric")
//         "SPECULARSPIKECONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability for reflection at the average surface
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//         "SPECULARLOBECONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability for reflection at a micro facet surface, i.e. the direction is smeared around the direction of the specular spike reflection
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//         "BACKSCATTERCONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability of back scattering, caused by several reflections within a deep grove
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//	    !!! the probability of internal (diffuse) Lambertian reflection can not be specified directly but is defined via
//		    100% = "SPECULARSPIKECONSTANT" + "SPECULARLOBECONSTANT" + "BACKSCATTERCONSTANT" + Lambertian (-> default is 1) !!!

	// define the optical surface properties
	G4MaterialPropertyVector * reflectivity_photonDetector = PropertyTools->GetPropertyDistribution(1.);   // metallic surface with 100% reflection => no absorption
// 	G4MaterialPropertyVector * reflectivity_sensitiveVolume = PropertyTools->GetPropertyDistribution(0.);  // surface with 100% reflection => no transmission

	G4MaterialPropertiesTable * MPT_OptSurf_photonDetector = new G4MaterialPropertiesTable();
	MPT_OptSurf_photonDetector->AddProperty("REFLECTIVITY", reflectivity_photonDetector);
	MPT_OptSurf_photonDetector->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(0.));
	MPT_OptSurf_photonDetector->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(1.));
	MPT_OptSurf_photonDetector->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
	OptSurf_photonDetector->SetMaterialPropertiesTable(MPT_OptSurf_photonDetector);

	G4MaterialPropertiesTable * MPT_OptSurf_sensitiveVolume = new G4MaterialPropertiesTable();
// 	MPT_OptSurf_sensitiveVolume->AddProperty("TRANSMITTANCE", reflectivity_sensitiveVolume);
	MPT_OptSurf_sensitiveVolume->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(0.));
	MPT_OptSurf_sensitiveVolume->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(1.));
	MPT_OptSurf_sensitiveVolume->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
	OptSurf_sensitiveVolume->SetMaterialPropertiesTable(MPT_OptSurf_sensitiveVolume);
}



void G4PhotonDetector::SetDefaults()
{
// Materials & Elements
	Material_PhotonDetector = 0;

// Photon Detector
	Depth = NAN;
	CoatingWidth = NAN;
	Transformation = G4Transform3D();
	PhotonDetector_physical = 0;
	SensitiveVolume_physical = 0;
	OptSurf_photonDetector = 0;
	OptSurf_sensitiveVolume = 0;
}
