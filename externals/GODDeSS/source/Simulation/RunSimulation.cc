/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <iostream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4UIExecutive.hh>
#include <G4VisExecutive.hh>

#include <PhysicsList.hh>
#include <DetectorConstruction.hh>
#include <GeneralParticleSource.hh>
#include <PhotonListSource.hh>
#include <GeneralParticleSourceMessenger.hh>
#include <RunAction.hh>
#include <SimulationMessenger.hh>
#include <EventAction.hh>
#include <TrackingAction.hh>
#include <SteppingAction.hh>

// GODDeSS
#include <GODDeSS_Messenger.hh>
#include <ScintillatorTileConstructor.hh>
#include <FibreConstructor.hh>
#include <PhotonDetectorConstructor.hh>

G4String assembleOutputFileName(G4String outDir, G4String outFileName, G4String addPhrase);

void ParseCommandLine(int argc, char** argv);
G4bool getBooleanOption(int arg, char** argv, G4String programParameterString, G4bool & variable, G4String coutString = "");
G4bool getStringlikeOption(int & arg, char** argv, G4String programParameterString, G4String & variable, G4String coutStringFragment);
G4bool getIntegerOption(int & arg, char** argv, G4String programParameterString, G4int & variable, G4String coutStringFragment);
void ParseInitFile(G4String initFile);

void Usage();
void Help();

G4String InitFile;
G4String MacroFile;
G4int NumEvents;
G4bool UsePhotonListSource;
G4String PSInputFile;
G4ThreeVector ScintillatorDimensions;
G4bool SearchOverlaps;
G4String ScintiFile;
G4String TwinTileReflectorFile;
G4String WrappingFile;
G4String LightGuidingFibreFile;
G4String WLSFibreFile;
G4String ScintiFibreFile;
G4String CementFile;
G4String OutDir;
G4String SimulationOutputFile;
G4String Phrase;
G4int ControlVerbosity;
G4int RunVerbosity;
G4int EventVerbosity;
G4int TrackingVerbosity;
G4bool DontTrackOpticalPhotons;








int main(int argc, char** argv)
{
/// set variables
	InitFile = "";
	MacroFile = "";
	NumEvents = -1;
	UsePhotonListSource = false;
	PSInputFile = "";
	ScintillatorDimensions = G4ThreeVector(NAN, NAN, NAN);
	SearchOverlaps = false;
	ScintiFile = "";
	TwinTileReflectorFile = "";
	WrappingFile = "";
	LightGuidingFibreFile = "";
	WLSFibreFile = "";
	ScintiFibreFile = "";
	CementFile = "";
	OutDir = "";
	SimulationOutputFile = "outputSimulation.data";
	G4String dataSubOutDir = "Data/";
	G4String controlDataSubOutDir = "ControlData/";
	Phrase = "";
	ControlVerbosity = 0;
	RunVerbosity = 0;
	EventVerbosity = 0;
	TrackingVerbosity = 0;
	DontTrackOpticalPhotons = false;

	ParseCommandLine(argc, argv);

	if(InitFile != "") ParseInitFile(InitFile);

	// get environment variable "SIMDIR":
	G4String simDir = "";
	if(getenv("SIMDIR")) simDir = getenv("SIMDIR");
	else
	{
		G4cerr << G4endl << "The environment variable \"SIMDIR\" (pointing to the directory of your simulation code) has not been specified!" << G4endl << G4endl;
		return 1;
	}

	if(OutDir == "" && getenv("BUILDDIR")) OutDir = getenv("BUILDDIR");

	std::string::reverse_iterator dirIter = OutDir.rbegin();
	if ((G4String) *dirIter != "/") OutDir += "/";

	// get current directory:
// 	G4String pwd = getenv("PWD");

	// change current working directory:
	G4String buildDir = "";
	if(getenv("BUILDDIR")) buildDir = getenv("BUILDDIR");
	else
	{
		G4cerr << G4endl << "The environment variable \"BUILDDIR\" (pointing to the directory of where the executable is build) has not been specified!" << G4endl << G4endl;
		return 1;
	}
	chdir(buildDir.c_str());

	// set minimal and maximal possible photon energy and create a vector containing them:
	G4double energiesMin = 1.0 * CLHEP::eV;   // at least minimal and maximal possible photon energy, necessary to define the refraction index etc. for the full spectrum of possible photon energies
	G4double energiesMax = 7.0 * CLHEP::eV;   // from ~180nm to ~1000nm -> Hamamatsu S13370 (UV-sensitive SiPM)

	// GODDeSS: energy range for the property distributions
	vector<G4double> energyRangeVector;
	G4double energyRangeVectorSize = 10;
	for(int iter = 0; iter < energyRangeVectorSize; iter++) energyRangeVector.push_back( energiesMin + iter * (energiesMax - energiesMin) / (energyRangeVectorSize - 1) );

/// initialise the simulation:
	// GODDeSS Messenger:
	GODDeSS_Messenger * goddessMessenger = new GODDeSS_Messenger(energyRangeVector);

	// create simulation messenger:
	SimulationMessenger * simulationMessenger = new SimulationMessenger(goddessMessenger);
	simulationMessenger->SetUsePhotonListSource(UsePhotonListSource);
	simulationMessenger->SetScintillatorDimensions(ScintillatorDimensions);
	simulationMessenger->SetSearchOverlaps(SearchOverlaps);
	simulationMessenger->SetScintillatorPropertyFile(ScintiFile);
	simulationMessenger->SetTwinTileReflectorPropertyFile(TwinTileReflectorFile);
	simulationMessenger->SetWrappingPropertyFile(WrappingFile);
	simulationMessenger->SetLightGuidingFibrePropertyFile(LightGuidingFibreFile);
	simulationMessenger->SetWLSFibrePropertyFile(WLSFibreFile);
	simulationMessenger->SetScintiFibrePropertyFile(ScintiFibreFile);
	simulationMessenger->SetOpticalCementPropertyFile(CementFile);
	simulationMessenger->SetTrackPhotons(!DontTrackOpticalPhotons);

	// compose the output file names
	G4String simOutFile = assembleOutputFileName(OutDir + dataSubOutDir, SimulationOutputFile, Phrase);
	simulationMessenger->SetDataFileName(simOutFile);

	G4String simControlOutFile = assembleOutputFileName(OutDir + controlDataSubOutDir, SimulationOutputFile, Phrase);
	simulationMessenger->SetControlFileName(simControlOutFile);

	// Construct the (default) run manager:
	G4RunManager * runManager = new G4RunManager;

	// set mandatory initialization classes and create constructors:
	G4VModularPhysicsList * physicsList = new PhysicsList();

	// GODDeSS Constructor:
	ScintillatorTileConstructor * scintillatorTileConstructor = new ScintillatorTileConstructor(physicsList, goddessMessenger->GetPropertyToolsManager(), goddessMessenger->GetDataStorage(), SearchOverlaps);
	goddessMessenger->SetScintillatorTileConstructor(scintillatorTileConstructor);

	FibreConstructor * fibreConstructor = new FibreConstructor(physicsList, goddessMessenger->GetPropertyToolsManager(), goddessMessenger->GetDataStorage(), SearchOverlaps);
	goddessMessenger->SetFibreConstructor(fibreConstructor);

	G4String hitFile = "hittingPhotons.hitData";
	hitFile = assembleOutputFileName(OutDir + dataSubOutDir, hitFile, Phrase);
	goddessMessenger->GetDataStorage()->SetPhotonDetectorHitFile(hitFile);
	PhotonDetectorConstructor * photonDetectorConstructor = new PhotonDetectorConstructor(physicsList, goddessMessenger->GetPropertyToolsManager(), goddessMessenger->GetDataStorage(), SearchOverlaps);
	goddessMessenger->SetPhotonDetectorConstructor(photonDetectorConstructor);

	// initialise the physics list (this has to be done AFTER the physics processes have been defined!!!)
	runManager->SetUserInitialization(physicsList);

	DetectorConstruction * detector = new DetectorConstruction(simulationMessenger);   // this has to be done AFTER the simulationMessenger has been filled!!!
	runManager->SetUserInitialization(detector);

	// set mandatory user action classes:
	GeneralParticleSourceMessenger * gpsm = 0;
	if(UsePhotonListSource)
	{
		PhotonListSource* pls = new PhotonListSource();
		runManager->SetUserAction(pls);
	}
	else
	{
		GeneralParticleSource* gps = new GeneralParticleSource();
		gpsm = gps->getMessenger();
		runManager->SetUserAction(gps);
	}

	// A Run is started by /run/beamOn and may consist of one or more events.
// 	G4UserRunAction* run_action = new RunAction();
	RunAction* run_action = new RunAction(simulationMessenger);
	runManager->SetUserAction(run_action);

	// An Event covers everything happening when the particle gun is fired once (no matter how many particles it fires).
	EventAction* event_action = new EventAction(simulationMessenger);
	runManager->SetUserAction(event_action);

// 	G4UserStackingAction* stacking_action = new ScintiStackingAction(&data);
// 	runManager->SetUserAction(stacking_action);

	// A Track is the full trajectory of a particle, it can be divided into steps (which is the part of the trajectory between two interactions).
	// In Geant4, each particle is described by one G4Track, but the processing of a particle's G4Track is suspended every time the original particle creates new particles (G4Tracks). Thus, for the user (i.e. in the G4UserTrackingAction and the G4UserSteppingAction) some particles appear to own several G4Tracks (all with the same TrackID and ParentID), which can lead to problems, e.g. when trying to determine the energy deposit of a single particle (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/eventtrackmanage/1072.html)!!!
	TrackingAction* tracking_action = new TrackingAction(simulationMessenger);
	runManager->SetUserAction(tracking_action);

	// A Step is a part of the full trajectory of a particle, confined by two interactions. The sum of all steps of one particle gives its full trajectory.
	SteppingAction* stepping_action = new SteppingAction(simulationMessenger);
	runManager->SetUserAction(stepping_action);

/// initialise G4 kernel (i.e. construct detector, define physics,...)
	runManager->Initialize();

	// get the pointer to the UI manager to apply commands
	G4UImanager * UI = G4UImanager::GetUIpointer();

	// Verbosity
	UI->ApplyCommand("/control/verbose " + boost::lexical_cast<G4String>(ControlVerbosity));
	UI->ApplyCommand("/run/verbose " + boost::lexical_cast<G4String>(RunVerbosity));
	UI->ApplyCommand("/event/verbose " + boost::lexical_cast<G4String>(EventVerbosity));
	UI->ApplyCommand("/tracking/verbose " + boost::lexical_cast<G4String>(TrackingVerbosity));

	// load the particle source input file
	if(UsePhotonListSource)
	{
		UI->ApplyCommand("/pls/photonListPath " + PSInputFile);
	}
	else
	{
		UI->ApplyCommand("/control/execute " + PSInputFile);

		G4String gpsOutFile = gpsm->getOutputFileName();
		gpsOutFile = assembleOutputFileName(OutDir + controlDataSubOutDir, gpsOutFile, Phrase);
		gpsm->setOutputFileName(gpsOutFile);
		simulationMessenger->SetGPSFileName(gpsOutFile);
	}

	// get the settings from a macro
	if(MacroFile != "")
	{
		G4cout << "Using macro file : " << MacroFile << G4endl;
		UI->ApplyCommand("/control/execute " + MacroFile);
	}

/// run Geant4
	// batch mode: no visualisation
	if(NumEvents > 0)
	{
		G4cout << "executing command: run/beamOn " << NumEvents << G4endl;
		UI->ApplyCommand("/run/beamOn " + boost::lexical_cast<G4String>(NumEvents));
	}
	// interactive mode: define visualization and UI terminal
	else if(NumEvents < 0)
	{
		/// Visualisation:
		// initialise visualisation manager
		G4VisManager * visManager = new G4VisExecutive;
		visManager->Initialize();

		// create session
		G4UIExecutive * session = new G4UIExecutive(argc, argv);

		// get visualisation commands
		UI->ApplyCommand("/control/execute " + simDir + "/macros/Visualisation.mac");

		// start session
		session->SessionStart();

		/// Job termination:
		// Everything initialised by the G4RunManager initialisation is owned and deleted by the G4RunManager and should not be deleted manually in the main() program!
		delete session;
		delete visManager;
	}

/// Job termination:
	// Everything initialised by the G4RunManager initialisation is owned and deleted by the G4RunManager and should not be deleted manually in the main() program!
	delete runManager;

	delete simulationMessenger;

	// GODDeSS
	delete goddessMessenger;
	delete scintillatorTileConstructor;
	delete fibreConstructor;
	delete photonDetectorConstructor;

// 	system( ("cd " + pwd).c_str() );

	return 0;
}



G4String assembleOutputFileName(G4String outDir, G4String outFileName, G4String addPhrase)
{
	system( ("mkdir -p " + outDir).c_str() );

	boost::match_results<std::string::const_iterator> result;
	if (regex_search(outFileName, result, boost::regex("([^.]+)[.]([^.]+)"), boost::match_all)) outFileName = result[1] + addPhrase + "." + result[2];
	else  outFileName += addPhrase;

	outFileName = outDir + outFileName;

	return outFileName;
}



void ParseCommandLine(int argc, char** argv)
{
	if(argc!=1)
	{
		for(int arg = 1; arg < argc; arg++)   //check if the help-option has been chosen
		{
			if( (strcmp(argv[arg], "--help") == 0) || (strcmp(argv[arg], "-h") == 0) )
			{
				Help();
				exit(-1);
			}
		}

		for(int arg = 1; arg < argc; arg++)
		{
			if(strcmp(argv[arg], "--init") == 0)
			{
				if(argc == 3)
				{
					if(argv[arg+1] && argv[arg+1][0] != '-')
					{
						InitFile = argv[arg+1];
						G4cout << "Used init file: " << InitFile << G4endl;
						break;
					}
					else
					{
						G4cout << "No init file specified!" << G4endl;
						Usage();
						exit(-1);
					}
				}
				else
				{
					G4cout << "Init file AND other options specified!" << G4endl;
					Usage();
					exit(-1);
				}
			}
			if(getStringlikeOption(arg, argv, "--macro", MacroFile, "macro filename")) continue;
			if(getIntegerOption(arg, argv, "--batch", NumEvents, "number of events"))
			{
				G4cout << "Running in batch mode!" << G4endl;
				continue;
			}
			if(getBooleanOption(arg, argv, "--usePhotonList", UsePhotonListSource)) continue;
			if(getStringlikeOption(arg, argv, "--particleSourceInput", PSInputFile, "particle source input file")) continue;
			if(strcmp(argv[arg], "--tileDims") == 0)
			{
				if(argv[arg+1] && argv[arg+1][0] != '-')
				{
					G4double vx = NAN;
					G4double vy = NAN;
					G4double vz = NAN;

					sscanf(argv[arg+1], "(%lf,%lf,%lf)", &vx, &vy, &vz);

					ScintillatorDimensions = G4ThreeVector(vx, vy, vz);

					arg++;

					continue;
				}
				else
				{
					G4cout << "No scintillator dimensions specified!" << G4endl;
					Usage();
					exit(-1);
				}
			}
			if(getBooleanOption(arg, argv, "--overlap", SearchOverlaps, "Searching for overlaps in the setup.")) continue;
			if(getStringlikeOption(arg, argv, "--scinti", ScintiFile, "scinti property file")) continue;
			if(getStringlikeOption(arg, argv, "--twin", TwinTileReflectorFile, "twin tile reflector property file")) continue;
			if(getStringlikeOption(arg, argv, "--wrapping", WrappingFile, "wrapping property file")) continue;
			if(getStringlikeOption(arg, argv, "--lgFibre", LightGuidingFibreFile, "fibre property file")) continue;
			if(getStringlikeOption(arg, argv, "--wlsFibre", WLSFibreFile, "fibre property file")) continue;
			if(getStringlikeOption(arg, argv, "--scintiFibre", ScintiFibreFile, "fibre property file")) continue;
			if(getStringlikeOption(arg, argv, "--cement", CementFile, "cement property file")) continue;
			if(getStringlikeOption(arg, argv, "--outDir", OutDir, "output directory")) continue;
			if(getStringlikeOption(arg, argv, "--outFile", SimulationOutputFile, "simulation output file")) continue;
			if(getStringlikeOption(arg, argv, "--add", Phrase, "phrase to be added to filename")) continue;
			if(getIntegerOption(arg, argv, "--controlVerbose", ControlVerbosity, "control verbosity")) continue;
			if(getIntegerOption(arg, argv, "--runVerbose", RunVerbosity, "run verbosity")) continue;
			if(getIntegerOption(arg, argv, "--eventVerbose", EventVerbosity, "event verbosity")) continue;
			if(getIntegerOption(arg, argv, "--trackingVerbose", TrackingVerbosity, "tracking verbosity")) continue;
			if(getBooleanOption(arg, argv, "--noOpticalPhotonTracking", DontTrackOpticalPhotons, "Don't track optical photons.")) {G4cout << "--------DontTrackOpticalPhotons " << DontTrackOpticalPhotons << G4endl; continue;}

			G4cout << "Unknown option: " << argv[arg] << G4endl;
			Usage();
			exit(-1);
		}
	}
}



G4bool getBooleanOption(int arg, char** argv, G4String programParameterString, G4bool & variable, G4String coutString)
{
	if(strcmp(argv[arg], programParameterString.c_str()) == 0)
	{
		variable = true;

		if(coutString != "") G4cout << G4endl << G4endl << coutString << G4endl << G4endl << G4endl;

		return true;
	}

	return false;
}



G4bool getStringlikeOption(int & arg, char** argv, G4String programParameterString, G4String & variable, G4String coutStringFragment)
{
	if(strcmp(argv[arg], programParameterString.c_str()) == 0)
	{
		if(argv[arg+1] && (argv[arg+1])[0] != '-')
		{
			variable = argv[arg+1];
			arg++;

			return true;
		}

		G4cout << "No " << coutStringFragment << " specified!" << G4endl;
		Usage();
		exit(-1);
	}

	return false;
}



G4bool getIntegerOption(int & arg, char** argv, G4String programParameterString, G4int & variable, G4String coutStringFragment)
{
	if(strcmp(argv[arg], programParameterString.c_str()) == 0)
	{
		if(argv[arg+1] && (argv[arg+1])[0] != '-')
		{
			variable = atoi(argv[arg+1]);
			arg++;

			return true;
		}

		G4cout << "No " << coutStringFragment << " specified!" << G4endl;
		Usage();
		exit(-1);
	}

	return false;
}



void ParseInitFile(G4String initFile)
{
	// Open the file.
	FILE * infile;
	infile = fopen(initFile.c_str(), "r");

	if(!infile)
	{
		std::cerr << "Cannot open file \"" << initFile << "\"" << std::endl;
		exit(-1);
	}
	else
	{
		std::vector<G4String> optionVector;
		optionVector.push_back("ProgramNameDummy");   // the first entry of argv (argv[0]) is the name of the program

		while(!feof(infile))
		{
			char line[1000];

			// get the lines
			fgets(line, 1000, infile);

			if(strncmp(line, "#", strlen("#")) != 0)   // if the line is not commented out
			{
				char rchar1[200];
				char rchar2[200];
				int numVariables = sscanf(line, "%s %s", rchar1, rchar2);
				if (numVariables > 0) optionVector.push_back(rchar1);
				if (numVariables > 1) optionVector.push_back(rchar2);
			}
		}

		int ARGC = optionVector.size();
		char **ARGV = (char**)malloc(ARGC * sizeof *ARGV);

		for(int iter = 0; iter < ARGC; iter++) {
			ARGV[iter] = (char*)malloc(strlen(optionVector[iter]) + 1);
			strcpy(ARGV[iter], optionVector[iter]);
		}

		ParseCommandLine(ARGC, ARGV);
	}

	fclose(infile);
}



void Usage()
{
	G4cout << "Type \"./RunSimulation --help\" for run options!" << G4endl;
}



void Help()
{
	G4cout << G4endl;
	G4cout << "== RunSimulation ========================================================================================================================" << G4endl << G4endl;
	G4cout << "   Options:"                                                                                                                               << G4endl;
	G4cout << "            --help                          : Show this output."                                                                           << G4endl;
	G4cout << "            --init <init file>              : Get run options from <init file>."                                                           << G4endl;
	G4cout << "                                              In this case, no other options (but \"--help\") when calling the program."                   << G4endl;
	G4cout << "            --macro <macroname.mac>         : Run with the settings specified in this macro."                                              << G4endl;
	G4cout << "            --batch <num of events>         : Force program to run in batch mode (without visualisation, i.e. faster)."                    << G4endl;
	G4cout << "                                              The given number of events will be simulated."                                               << G4endl;
	G4cout << "            --usePhotonList                 : Use a list with photons to define the primary particles."                                    << G4endl;
	G4cout << "                                              This list has to be specified with the \"--particleSourceInput\" option."                    << G4endl;
	G4cout << "            --particleSourceInput <filename>: Use <filename> as input for the particle source."                                            << G4endl;
	G4cout << "                                              This can either be a macro for the General Particle Source (if the \"--usePhotonList\""      << G4endl;
	G4cout << "                                              option IS NOT set) or a list with photons (if the \"--usePhotonList\" option IS set)."       << G4endl;
	G4cout << "            --tileDims <(X,Y,Z)>            : Use different dimensions (in mm) for the scintillator tile than the hardcoded ones."         << G4endl;
	G4cout << "            --overlap                       : Search for overlaps in the setup. Default is: \"false\""                                     << G4endl;
	G4cout << "            --scinti <filename.properties>  : Get material properties of the scintillator from <filename.properties>"                      << G4endl;
	G4cout << "            --twin <filename.properties>    : Get material properties of the scintillator twin tile reflector from <filename.properties>"  << G4endl;
	G4cout << "            --wrapping <filename.properties>: Get material properties of the wrapping from <filename.properties>"                          << G4endl;
	G4cout << "            --fibre <filename.properties>   : Get material properties of the fibre from <filename.properties>"                             << G4endl;
	G4cout << "            --cement <filename.properties>  : Get material properties of the optical cement from <filename.properties>"                    << G4endl;
	G4cout << "            --outDir <directory>            : Define output directory. This directory must exisit! Default is: \"./\""                     << G4endl;
	G4cout << "            --outFile <file name>           : Define file name for output files of the simulation. Default is: \"outputSimulation.data\""  << G4endl;
	G4cout << "            --add <phrase>                  : Add <phrase> to output file names."                                                          << G4endl;
	G4cout << "            --controlVerbose <number>       : Specify the GEANT4 verbosity level set using \"/control/verbose\". Default is: \"0\""        << G4endl;
	G4cout << "            --runVerbose <number>           : Specify the GEANT4 verbosity level set using \"/run/verbose\". Default is: \"0\""            << G4endl;
	G4cout << "            --eventVerbose <number>         : Specify the GEANT4 verbosity level set using \"/event/verbose\". Default is: \"0\""          << G4endl;
	G4cout << "            --trackingVerbose <number>      : Specify the GEANT4 verbosity level set using \"/tracking/verbose\". Default is: \"0\""       << G4endl;
	G4cout << "            --noOpticalPhotonTracking       : Do not track optical photons when simulating events."                                        << G4endl << G4endl;
	G4cout << "=========================================================================================================================================" << G4endl << G4endl;
}
