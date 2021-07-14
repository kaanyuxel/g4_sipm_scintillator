#include <vector>

#ifdef G4MULTITHREADED
    #include <G4MTRunManager.hh>
    using RunManager = G4MTRunManager;
#else
    #include <G4RunManager.hh>
    using RunManager = G4RunManager;
#endif

#ifdef G4VIS_USE
    #include <G4VisExecutive.hh>
#endif

#ifdef G4UI_USE
    #include <G4UIExecutive.hh>
#endif

#include <G4String.hh>
#include <G4UImanager.hh>

#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include <GODDeSS_Messenger.hh>
#include <ScintillatorTileConstructor.hh>
#include <FibreConstructor.hh>
#include <PhotonDetectorConstructor.hh>    
#include "G4ScoringManager.hh"

using namespace std;

int main(int argc, char** argv)
{
    std::cout << "Application starting..." << std::endl;

    vector<G4String> macros;
    bool interactive = false;

    // Parse command line arguments
    if  (argc == 1)
    {
        interactive = true;
    }
    else
    {
        for (int i = 1; i < argc; i++)
        {
            G4String arg = argv[i];
            if (arg == "-i" || arg == "--interactive")
            {
                interactive = true;
                continue;
            }
            else
            {
                macros.push_back(arg);
            }
        }
    }

    // Create the run manager (MT or non-MT) and make it a bit verbose.
    auto runManager = new RunManager();
    runManager->SetVerboseLevel(1);

    #ifdef G4VIS_USE
        G4VisManager* visManager = new G4VisExecutive();
        visManager->Initialize();
    #endif

    G4VModularPhysicsList * physicsList = new PhysicsList();
    runManager->SetUserInitialization(physicsList);

    // set minimal and maximal possible photon energy and create a vector containing them:
    G4double energiesMin = 1.0 * CLHEP::eV;   // at least minimal and maximal possible photon energy, necessary to define the refraction index etc. for the full spectrum of possible photon energies
    G4double energiesMax = 5.0 * CLHEP::eV;   // from ~180nm to ~1000nm -> Hamamatsu S13370 (UV-sensitive SiPM)

    // GODDeSS: energy range for the property distributions
    vector<G4double> energyRangeVector;
    G4double energyRangeVectorSize = 10;
    for(int iter = 0; iter < energyRangeVectorSize; iter++) energyRangeVector.push_back( energiesMin + iter * (energiesMax - energiesMin) / (energyRangeVectorSize - 1) );

/// initialise the simulation:
    // GODDeSS Messenger:
    GODDeSS_Messenger * goddessMessenger = new GODDeSS_Messenger(energyRangeVector);
    G4bool SearchOverlaps = false; // set to true to check your geometry
    ScintillatorTileConstructor *scintillatorTileConstructor = new ScintillatorTileConstructor(physicsList, goddessMessenger->GetPropertyToolsManager(), goddessMessenger->GetDataStorage(), SearchOverlaps);
    goddessMessenger->SetScintillatorTileConstructor(scintillatorTileConstructor);

    runManager->SetUserInitialization(new DetectorConstruction(goddessMessenger));
    runManager->SetUserInitialization(new ActionInitialization());

    #ifdef G4UI_USE
        G4UIExecutive* ui = nullptr;
        if (interactive)
        {
            ui = new G4UIExecutive(argc, argv);
        }
    #endif

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    G4ScoringManager::GetScoringManager();

    for (auto macro : macros)
    {
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macro);
    }

    #ifdef G4UI_USE
        if (interactive)
        {
            if (ui->IsGUI()) 
            { 
                UImanager->ApplyCommand("/control/execute macros/vis.mac"); 
            } 
            else 
            { 
                UImanager->ApplyCommand("/run/initialize"); 
            } 
            ui->SessionStart();
            delete ui;
        }
    #endif

    delete runManager;
    delete goddessMessenger;
    delete scintillatorTileConstructor;    

    std::cout << "Application successfully ended.\nBye :-)" << std::endl;

    return 0;
}
