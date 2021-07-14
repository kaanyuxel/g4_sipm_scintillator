#include "DetectorConstruction.hh"
#include "G4Sipm.hh"
#include "MaterialFactory.hh"
#include "model/G4SipmModelFactory.hh"
#include "housing/G4SipmHousing.hh"
#include "housing/impl/HamamatsuCeramicHousing.hh"
#include "housing/impl/HamamatsuSmdHousing.hh"

#include <PropertyToolsManager.hh>
#include <Properties.hh>

#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4NistManager.hh>
#include <G4SystemOfUnits.hh>
#include <G4VisAttributes.hh>
#include <G4Box.hh>
#include <G4Orb.hh>
#include <G4SDManager.hh>
#include <G4GlobalMagFieldMessenger.hh>

#include <sstream>

using namespace std;

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    G4NistManager* nist = G4NistManager::Instance();
    G4double worldSizeX = 5 * m;
    G4double worldSizeY = 5 * m;
    G4double worldSizeZ = 5 * m;

    DefineElements();

    // 1) Solid
    G4VSolid* worldBox = new G4Box("world", worldSizeX / 2, worldSizeY / 2, worldSizeZ / 2);

    // 2) Logical volume
    G4LogicalVolume* worldLog = new G4LogicalVolume(worldBox, nist->FindOrBuildMaterial("G4_AIR"), "world");
    G4VisAttributes* visAttr = new G4VisAttributes();
    visAttr->SetVisibility(false);
    worldLog->SetVisAttributes(visAttr);

    // 3) Physical volume
    G4VPhysicalVolume* worldPhys = new G4PVPlacement(nullptr, {}, worldLog, "world", nullptr, false, 0);

    ScintillatorTileConstructor * STConstructor1 = goddessMessenger->GetScintillatorTileConstructor();
    G4ThreeVector ScintiDimensions1 = G4ThreeVector(200. * mm, 10. * mm, 200. * mm);
    STConstructor1->SetScintillatorTransformation(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., -0.8*m, 0.)));
    STConstructor1->SetScintillatorName("scintillator1");
    STConstructor1->ConstructASensitiveDetector();
    
    G4ScintillatorTile * scintillator1 = STConstructor1->ConstructScintillator(ScintiDimensions1, "externals/GODDeSS/MaterialProperties/Scintillator/Saint-Gobain_BC-408.properties", worldPhys);    
    STConstructor1->SetWrappingName("wrapping");
    STConstructor1->ConstructWrapping(scintillator1, "externals/GODDeSS/MaterialProperties/Scintillator/Wrapping_Teflon.properties");    
   
    ScintillatorTileConstructor * STConstructor2 = goddessMessenger->GetScintillatorTileConstructor();   
    G4ThreeVector ScintiDimensions2 = G4ThreeVector(200. * mm, 10. * mm, 200. * mm);
    STConstructor2->SetScintillatorTransformation(G4Transform3D(G4RotationMatrix(), G4ThreeVector(0., 0.8* m, 0.)));
    STConstructor2->SetScintillatorName("scintillator2");
    STConstructor2->ConstructASensitiveDetector();
    
    G4ScintillatorTile * scintillator2 = STConstructor2->ConstructScintillator(ScintiDimensions2, "externals/GODDeSS/MaterialProperties/Scintillator/Saint-Gobain_BC-408.properties", worldPhys);    
    STConstructor2->SetWrappingName("wrapping");
    STConstructor2->ConstructWrapping(scintillator2, "externals/GODDeSS/MaterialProperties/Scintillator/Wrapping_Teflon.properties");    
 
    G4SipmModel* model = G4SipmModelFactory::getInstance()->createHamamatsuS12651050();
    G4Sipm* sipm1 = new G4Sipm(model);
    G4SipmHousing* house1 = new G4SipmHousing(sipm1);
    house1->setPosition(G4ThreeVector(0., 0., -200. * mm/2));
    house1->buildAndPlace(scintillator1->GetScintillator_physicalVolume());

    G4Sipm* sipm2 = new G4Sipm(model);
    G4SipmHousing* house2 = new G4SipmHousing(sipm2);
    house2->setPosition(G4ThreeVector(0., 0., -200. * mm/2));
    house2->buildAndPlace(scintillator2->GetScintillator_physicalVolume());    
   
    return worldPhys;
}

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

void DetectorConstruction::ConstructSDandField()
{
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    sdManager->SetVerboseLevel(2); 
}