/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4VisAttributes.hh>
#include <G4PVPlacement.hh>
#include <G4GeometryManager.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4Element.hh>
#include <G4MaterialTable.hh>

#include <DetectorConstruction.hh>

using namespace CLHEP;   //for mathematics (e.g. CLHEP::sqrt, CLHEP::pi, CLHEP::c_light, CLHEP::h_Planck,...)



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to create the simulation setup and register it to Geant4:
 *  - will be called by Geant4 in the initialisation process and <b> has to create materials and volumes as well as to return the pointer to the world volume (the volume that contains all other volumes) </b>
 *  - calls DefineVariables(), DefineElements(), DefineMaterials(), and DefineMaterialProperties().
 *  - creates the setup that is defined inside and registers it to Geant4
 */
G4VPhysicalVolume* DetectorConstruction::Construct()
{
	DefineVariables();
	DefineElements();
	DefineMaterials();
	DefineMaterialProperties();

//   ######  solids (dimensions):  ######   //

	//---------- experimental hall (world volume) ----------//
	G4Box * world_solid = new G4Box("world_solid", WorldDimensions.x(), WorldDimensions.y(), WorldDimensions.z());

//   ######  logical volumes (material):  ######   //

	//---------- experimental hall (world volume) ----------//
	G4LogicalVolume * world_logical = new G4LogicalVolume(world_solid, Material_Vacuum, "world_logical", 0, 0, 0);
	world_logical->SetVisAttributes(G4VisAttributes::Invisible);

//   ######  physical volumes (placement):  ######   //
//NOTE: G4PVPlacement puts the volume to be placed into EVERY physical volume emanating from the same logical volume (no matter whether the logical or physical volume is specified as mother volume)!
//      => If G4PVPlacement should be able to distinguish physical volumes, for each physical volume a separate logical volume has to be created.
//      => rule of thumb: For each volume that might become a mother volume, a separate logical volume should to be created.

	//---------- experimental hall (world volume) ----------//
	G4VPhysicalVolume * world_physical = new G4PVPlacement(G4Transform3D(), world_logical, "world", 0, false, 0);





	G4double edgeLength = 0.;
	G4ThreeVector sensitiveSurfaceNormalRelativeToReferenceVolume = G4ThreeVector();
	G4ThreeVector sensitiveSurfacePositionRelativeToReferenceVolume = G4ThreeVector();

// >>> Setup: 1 SiPM at the side <<< //

	ScintiTransform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0., 0.));
	edgeLength = 3. * CLHEP::mm;
	sensitiveSurfaceNormalRelativeToReferenceVolume = G4ThreeVector(0., 0., -1.);
	sensitiveSurfacePositionRelativeToReferenceVolume = G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, ScintiDimensions[2] / 2.);

	//---------- scintillator tile ----------//
	STConstructor->SetScintillatorTransformation(ScintiTransform);
	STConstructor->SetScintillatorName("scintillator 1");
	STConstructor->ConstructASensitiveDetector();
	G4ScintillatorTile * scintillator_1 = STConstructor->ConstructScintillator(ScintiDimensions, ScintiPropertyFile, world_physical);
	G4VPhysicalVolume * scintillator_1_physical = scintillator_1->GetScintillator_physicalVolume();

	//---------- SiPM ----------//
	PDConstructor->SetPhotonDetectorName("SiPM 1");
	PDConstructor->SetPhotonDetectorReferenceVolume(scintillator_1_physical);
	PDConstructor->SetSensitiveSurfaceNormalRelativeToReferenceVolume(sensitiveSurfaceNormalRelativeToReferenceVolume);
	PDConstructor->SetSensitiveSurfacePositionRelativeToReferenceVolume(sensitiveSurfacePositionRelativeToReferenceVolume);
	PDConstructor->ConstructPhotonDetector(edgeLength, world_physical);

	//---------- external reflector / wrapping ----------//
	STConstructor->SetWrappingName("wrapping 1");
	STConstructor->ConstructWrapping(scintillator_1, WrappingPropertyFile);


	// >>> Setup: 1 glued fibre with 1 SiPM and a reflective end <<< //

	ScintiTransform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 50. * CLHEP::mm, 0.));
	FibreEndPoint = G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, ScintiDimensions[2] / 2. + 5. * CLHEP::mm);
	FibreStartPoint = G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, -ScintiDimensions[2] / 2.);
	edgeLength = 1. * CLHEP::mm;

	//---------- scintillator tile ----------//
	STConstructor->SetScintillatorTransformation(ScintiTransform);
	STConstructor->SetScintillatorName("scintillator 2");
	G4ScintillatorTile * scintillator_2 = STConstructor->ConstructScintillator(ScintiDimensions, ScintiPropertyFile, world_physical);

	//---------- fibre ----------//
	FConstructor->SetFibreStartPointReflectivity(FibreEndReflectivity);
	G4Fibre * fibre_2 = FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_2->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint);

	G4VPhysicalVolume * fibre_2_physical = fibre_2->GetOutermostVolumeOutsideMother_physicalVolume();

	//---------- SiPM ----------//
	G4double fibreLength = fibre_2->GetFibreLength();
	sensitiveSurfaceNormalRelativeToReferenceVolume = G4ThreeVector(0., 0., -1.);
	sensitiveSurfacePositionRelativeToReferenceVolume = G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, fibreLength / 2.);

	PDConstructor->SetPhotonDetectorName("SiPM 2");
	PDConstructor->SetPhotonDetectorReferenceVolume(fibre_2_physical);
	PDConstructor->SetSensitiveSurfaceNormalRelativeToReferenceVolume(sensitiveSurfaceNormalRelativeToReferenceVolume);
	PDConstructor->SetSensitiveSurfacePositionRelativeToReferenceVolume(sensitiveSurfacePositionRelativeToReferenceVolume);
	PDConstructor->ConstructPhotonDetector(edgeLength, world_physical);

	//---------- external reflector / wrapping ----------//
	STConstructor->SetWrappingName("wrapping 2");
	STConstructor->ConstructWrapping(scintillator_2, WrappingPropertyFile);


	// >>> Setup: 1 glued sigma fibre with 1 SiPM and a reflective end <<< //

	ScintiTransform = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., -50. * CLHEP::mm, 0.));
	G4double fibreDistanceToScintillatorEdge = 50. * CLHEP::mm;
	G4double fibreBendigRadius = 50. * CLHEP::mm;
	if(fabs(ScintiDimensions[0] - 100. * CLHEP::mm) < 1e-12)
	{
		fibreDistanceToScintillatorEdge = 15. * CLHEP::mm;
		fibreBendigRadius = 15. * CLHEP::mm;
	}
	edgeLength = 1. * CLHEP::mm;

	//---------- scintillator tile ----------//
	STConstructor->SetScintillatorTransformation(ScintiTransform);
	STConstructor->SetScintillatorName("scintillator 3");
	G4ScintillatorTile * scintillator_3 = STConstructor->ConstructScintillator(ScintiDimensions, ScintiPropertyFile, world_physical);

	//---------- fibre ----------//
	G4double straightFibreLength_x = 2. * (ScintiDimensions[0] / 2. - (fibreDistanceToScintillatorEdge + fibreBendigRadius));
	G4double straightFibreLength_z = 2. * (ScintiDimensions[2] / 2. - (fibreDistanceToScintillatorEdge + fibreBendigRadius));

	FibreStartPoint = G4ThreeVector(ScintiDimensions[0] / 2. - fibreDistanceToScintillatorEdge - 5. * CLHEP::mm, 2. * CLHEP::mm, ScintiDimensions[2] / 2. - fibreDistanceToScintillatorEdge);
	FibreEndPoint = FibreStartPoint + G4ThreeVector(-straightFibreLength_x, 0. * CLHEP::mm, 0. * CLHEP::mm) + G4ThreeVector(-(fibreBendigRadius - 5. * CLHEP::mm), 0. * CLHEP::mm, 0. * CLHEP::mm);
	FConstructor->SetFibreStartPointReflectivity(FibreEndReflectivity);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "E");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint);

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(-fibreBendigRadius, 0. * CLHEP::mm, -fibreBendigRadius);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "SE");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint, 90. * deg, G4ThreeVector(0, -1, 0));

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, -straightFibreLength_z);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "SE");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint);

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(fibreBendigRadius, 0. * CLHEP::mm, -fibreBendigRadius);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "SE");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint, 90. * deg, G4ThreeVector(0, -1, 0));

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(straightFibreLength_x, 0. * CLHEP::mm, 0. * CLHEP::mm);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "SE");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint);

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(fibreBendigRadius, 0. * CLHEP::mm, fibreBendigRadius);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "SE");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint, 90. * deg, G4ThreeVector(0, -1, 0));

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, straightFibreLength_z) + G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, fibreDistanceToScintillatorEdge + fibreBendigRadius);
	FConstructor->SetFibreGlued(OpticalCementPropertyFile, "quadratic", 0, "SE");
	FConstructor->ConstructFibre(WLSFibrePropertyFile, scintillator_3->GetScintillator_physicalVolume(), FibreStartPoint, FibreEndPoint);

	FibreStartPoint = FibreEndPoint;
	FibreEndPoint = FibreStartPoint + G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, 2000. * CLHEP::mm);
	G4Fibre * fibre_3 = FConstructor->ConstructFibre(LightGuidingFibrePropertyFile, world_physical, FibreStartPoint, FibreEndPoint, scintillator_3->GetScintillator_physicalVolume());

	G4VPhysicalVolume * fibre_3_physical = fibre_3->GetOutermostVolumeInsideMother_physicalVolume();

	///---------- SiPM ----------///
	fibreLength = fibre_3->GetFibreLength();
	sensitiveSurfaceNormalRelativeToReferenceVolume = G4ThreeVector(0., 0., -1.);
	sensitiveSurfacePositionRelativeToReferenceVolume = G4ThreeVector(0. * CLHEP::mm, 0. * CLHEP::mm, fibreLength / 2.);

	PDConstructor->SetPhotonDetectorName("SiPM");
	PDConstructor->SetPhotonDetectorReferenceVolume(fibre_3_physical);
	PDConstructor->SetSensitiveSurfaceNormalRelativeToReferenceVolume(sensitiveSurfaceNormalRelativeToReferenceVolume);
	PDConstructor->SetSensitiveSurfacePositionRelativeToReferenceVolume(sensitiveSurfacePositionRelativeToReferenceVolume);

	///---------- external reflector / wrapping ----------///
	STConstructor->SetWrappingName("wrapping 3");
	STConstructor->ConstructWrapping(scintillator_3, WrappingPropertyFile);
	PDConstructor->ConstructPhotonDetector(edgeLength, world_physical);




	return world_physical;   //return experimentalHall-volume with all volumes placed inside
}



/**
 *  Function for the initialisation of variables.
 */
void DetectorConstruction::DefineVariables()
{
//---------- experimental hall (world volume) ----------//
	WorldDimensions = G4ThreeVector(10. * m, 10. * m, 10. * m);

//---------- scintillator tile ----------//

	ScintiDimensions = Messenger->GetScintillatorDimensions();
	if(!(!isnan(ScintiDimensions.x()) && !isnan(ScintiDimensions.y()) && !isnan(ScintiDimensions.z())))
	{
		ScintiDimensions = G4ThreeVector(100. * CLHEP::mm, 10. * CLHEP::mm, 100. * CLHEP::mm);
	}

//---------- fibre ----------//

	FibreEndReflectivity = 0.9;
}



/**
 *  Function to create chemical elements and register them to Geant4:
 *  - creates hydrogen, carbon, nitrogen, oxygen, fluorine, aluminum, titanium, and lead
 *  - <b> if other chemical elements are needed for the simulation, they have to be defined here </b>
 */
void DetectorConstruction::DefineElements()
{
// http://pdg.lbl.gov/2009/AtomicNuclearProperties/index.html
	new G4Element("Hydrogen", "H", 1., 1.00794 * g/mole);
	new G4Element("Carbon", "C", 6., 12.0107 * g/mole);
	new G4Element("Nitrogen", "N", 7., 14.0067 * g/mole);
	new G4Element("Oxygen", "O", 8., 15.9994 * g/mole);
	new G4Element("Fluorine", "F", 9., 18.9984032 * g/mole);
	new G4Element("Aluminum", "Al", 13., 26.9815386 * g/mole);
	new G4Element("Titanium", "Ti", 22., 47.867 * g/mole);
	new G4Element("Lead", "Pb", 82., 207.2 * g/mole);
}



/**
 *  Function to create materials and register them to Geant4.
 */
void DetectorConstruction::DefineMaterials()
{
//---------- vacuum ----------//
	// G4Material(const G4String &name, G4double z, G4double a, G4double density, G4State state=kStateUndefined, G4double temp=STP_Temperature, G4double pressure=STP_Pressure)
	Material_Vacuum = new G4Material("Vacuum", 1., 1.01 * g/mole, universe_mean_density, kStateGas, 0.1 * kelvin, 1.e-19 * pascal);		// definition from Geant

//---------- air ----------//
	// G4Material(const G4String &name, G4double density, G4int nComponents, G4State state=kStateUndefined, G4double temp=STP_Temperature, G4double pressure=STP_Pressure)
	Material_Air = new G4Material("Air", 1.293 * kg/m3, 2);
	Material_Air->AddElement(G4Element::GetElement("Nitrogen"), 70 * perCent);
	Material_Air->AddElement(G4Element::GetElement("Oxygen"), 30 * perCent);
}



//   ######  material properties:  ######   //
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
//		(values > 1 result in a wider distribution, values < 1 result in a narrower distribution -> 1 should be chosen as default)
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

/**
 *  Function to set material properties and register them to Geant4.
 */
void DetectorConstruction::DefineMaterialProperties()
{
//---------- experimental hall (world volume) -> vacuum or air ----------//
	// by now, GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Vacuum = PropertyTools->GetPropertyDistribution(1.);

	G4MaterialPropertiesTable * mpt_Vacuum = new G4MaterialPropertiesTable();
	mpt_Vacuum->AddProperty("RINDEX", refractiveIndex_Vacuum);
	Material_Vacuum->SetMaterialPropertiesTable(mpt_Vacuum);


	// by now, GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Air = PropertyTools->GetPropertyDistribution(1.);
	// the refractive index of air can approximated by 1 as it only varies from 1.000308 (230nm) to 1.00027417 (1000nm)
	// http://refractiveindex.info/?group=GASES&material=Air

	// absorption and rayleigh scattering can be neglected, as only relatively thin layers are used in this simulation

	G4MaterialPropertiesTable * mpt_Air = new G4MaterialPropertiesTable();
	mpt_Air->AddProperty("RINDEX", refractiveIndex_Air);
	Material_Air->SetMaterialPropertiesTable(mpt_Air);
}



void DetectorConstruction::CleanUp()
{
	G4GeometryManager::GetInstance()->OpenGeometry();

	G4LogicalSkinSurface::CleanSurfaceTable();
	G4LogicalBorderSurface::CleanSurfaceTable();

	G4PhysicalVolumeStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4SolidStore::GetInstance()->Clean();
}



void DetectorConstruction::DeleteMaterialPropertiesTables()
{
	G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();   // getting the table of materials from Geant
	for(size_t i = 0; i < matTable->size(); i++) delete (*(matTable))[i]->GetMaterialPropertiesTable();
}



void DetectorConstruction::DeleteMaterials()
{
	G4MaterialTable* matTable = (G4MaterialTable*) G4Material::GetMaterialTable();   // getting the table of materials from Geant
	for(size_t i = 0; i < matTable->size(); i++) delete (*(matTable))[i];
	matTable->clear();
}
