/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef DETECTORCONSTRUCTION_HH_
#define DETECTORCONSTRUCTION_HH_

#include <G4VUserDetectorConstruction.hh>
#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Cons.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4Material.hh>
#include <G4OpticalSurface.hh>

#include <SimulationMessenger.hh>
#include <GoddessProperties.hh>
#include <ScintillatorTileConstructor.hh>
#include <FibreConstructor.hh>
#include <PhotonDetectorConstructor.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Builds the simulated setup, using functions inherited from G4VUserDetectorConstruction.
class DetectorConstruction: public G4VUserDetectorConstruction
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	DetectorConstruction( SimulationMessenger * simulationMessenger   /**< class storing and providing variables which are needed in different parts of the simulation. */
			    )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: Messenger(simulationMessenger)
	, PropertyTools(Messenger->GetGoddessMessenger()->GetPropertyToolsManager())
	, SearchOverlaps(Messenger->GetSearchOverlaps())
	, STConstructor(Messenger->GetGoddessMessenger()->GetScintillatorTileConstructor())
	, FConstructor(Messenger->GetGoddessMessenger()->GetFibreConstructor())
	, PDConstructor(Messenger->GetGoddessMessenger()->GetPhotonDetectorConstructor())
	, Material_Vacuum(0)
	, Material_Air(0)
	, WorldDimensions(G4ThreeVector(NAN, NAN, NAN))
	, ScintiPropertyFile(Messenger->GetScintillatorPropertyFile())
	, WrappingPropertyFile(Messenger->GetWrappingPropertyFile())
	, ScintiDimensions(G4ThreeVector(NAN, NAN, NAN))
	, ScintiTransform(G4Transform3D())
	, LightGuidingFibrePropertyFile(Messenger->GetLightGuidingFibrePropertyFile())
	, WLSFibrePropertyFile(Messenger->GetWLSFibrePropertyFile())
	, ScintiFibrePropertyFile(Messenger->GetScintiFibrePropertyFile())
	, OpticalCementPropertyFile(Messenger->GetOpticalCementPropertyFile())
	, FibreStartPoint(G4ThreeVector(NAN, NAN, NAN))
	, FibreEndPoint(G4ThreeVector(NAN, NAN, NAN))
	{
	}

	/**
	 *  Destructor:
	 *  - deletes the objects, that have been created and will NOT automatically be created by Geant4
	 */
	~DetectorConstruction()
	{
		CleanUp();
		DeleteMaterialPropertiesTables();
		DeleteMaterials();
	}



	// the following function will be called by Geant4 in the initialisation process:
	G4VPhysicalVolume* Construct();

private:
	void DefineVariables();
	void DefineElements();
	void DefineMaterials();
	void DefineMaterialProperties();

	void CleanUp();
	void DeleteMaterials();
	void DeleteMaterialPropertiesTables();



	SimulationMessenger * Messenger;
	PropertyToolsManager * PropertyTools;
	G4bool SearchOverlaps;
	ScintillatorTileConstructor * STConstructor;
	FibreConstructor * FConstructor;
	PhotonDetectorConstructor * PDConstructor;

// Materials
	G4Material * Material_Vacuum;
	G4Material * Material_Air;

// the world
	G4ThreeVector WorldDimensions;

// Scintillator
	G4String ScintiPropertyFile;
	G4String WrappingPropertyFile;
	G4ThreeVector ScintiDimensions;
	G4Transform3D ScintiTransform;

// fibre
	G4String LightGuidingFibrePropertyFile;
	G4String WLSFibrePropertyFile;
	G4String ScintiFibrePropertyFile;
	G4String OpticalCementPropertyFile;
	G4ThreeVector FibreStartPoint;
	G4ThreeVector FibreEndPoint;
	G4double FibreEndReflectivity;
};

#endif
