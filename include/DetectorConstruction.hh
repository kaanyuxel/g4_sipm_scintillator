#ifndef DETECTOR_CONSTRUCTION_HH
#define DETECTOR_CONSTRUCTION_HH

#include <G4VUserDetectorConstruction.hh>

#include <ScintillatorTileConstructor.hh>
#include <FibreConstructor.hh>
#include <PhotonDetectorConstructor.hh>
#include <GODDeSS_Messenger.hh>

class G4LogicalVolume;
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

    DetectorConstruction(GODDeSS_Messenger* gMessenger){ goddessMessenger = gMessenger;};
    G4VPhysicalVolume* Construct() override;

    void DefineElements();
    void ConstructSDandField() override;

private:
    GODDeSS_Messenger *goddessMessenger;
};
#endif
