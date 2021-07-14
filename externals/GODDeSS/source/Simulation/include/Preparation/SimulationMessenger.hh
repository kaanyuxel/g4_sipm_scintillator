/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef SimulationMessenger_h
#define SimulationMessenger_h 1

#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>
#include <UserRunInformation.hh>

#include <GODDeSS_Messenger.hh>



// class variables begin with capital letters, local variables with small letters



/// <b> (Part of the example simulation and not belonging to GODDeSS.) </b> Class to store and provide variables which are needed in different parts of the simulation.
class SimulationMessenger: public G4UImessenger
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 *  - creates commands for operating in the interactive mode
	 */
	SimulationMessenger( GODDeSS_Messenger * goddess_messenger   /**< pointer to the GODDeSS_DataStorage that is to be used */
			   )
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: GoddessMessenger(goddess_messenger)
	, UsePhotonListSource(false)
	, SetupIdentificationNumber(-1)
	, ScintillatorDimensions(G4ThreeVector(NAN, NAN, NAN))
	, SearchOverlaps(false)
	, RunInformation(0)
	, ScintillatorPropertyFile("")
	, TwinTileReflectorPropertyFile("")
	, WrappingPropertyFile("")
	, LightGuidingFibrePropertyFile("")
	, WLSFibrePropertyFile("")
	, ScintiFibrePropertyFile("")
	, OpticalCementPropertyFile("")
	, DiffuserPropertyFile("")
	, TrackPhotons(true)
	, GPSFile("outputParticleSource.data")
	, DataFile("outputSimulation.data")
	, ControlFile("outputSimulation.data")
	, HitFile("hittingPhotons.data")
	{
		SimulationDir = new G4UIdirectory("/simulation/");
		SimulationDir->SetGuidance("UI commands of this setup");

		TrackPhotonsCmd = new G4UIcmdWithABool("/simulation/trackOpticalPhotons",this);
		TrackPhotonsCmd->SetGuidance("Switch optical photon tracking on/off. Default is \"true\" ");
		TrackPhotonsCmd->SetDefaultValue(true);

		OutputFileCmd = new G4UIcmdWithAString("/simulation/outFileName", this);
		OutputFileCmd->SetGuidance("Set the output file.");
		OutputFileCmd->SetDefaultValue("outputSimulation.txt");
	}

	/**
	 *  Destructor:
	 *  - deletes the command objects
	 */
	~SimulationMessenger()
	{
		delete OutputFileCmd;
		delete TrackPhotonsCmd;
		delete SimulationDir;
	}



	void SetNewValue(G4UIcommand*, G4String);

private:
	G4UIdirectory* SimulationDir;

	G4UIcmdWithABool* TrackPhotonsCmd;
	G4UIcmdWithAString* OutputFileCmd;

	GODDeSS_Messenger * GoddessMessenger;
	G4bool UsePhotonListSource;
	G4int SetupIdentificationNumber;
	G4ThreeVector ScintillatorDimensions;
	G4bool SearchOverlaps;
	UserRunInformation * RunInformation;
	G4String ScintillatorPropertyFile;
	G4String TwinTileReflectorPropertyFile;
	G4String WrappingPropertyFile;
	G4String LightGuidingFibrePropertyFile;
	G4String WLSFibrePropertyFile;
	G4String ScintiFibrePropertyFile;
	G4String OpticalCementPropertyFile;
	G4String DiffuserPropertyFile;
	G4bool TrackPhotons;
	G4String GPSFile;
	G4String DataFile;
	G4String ControlFile;
	G4String HitFile;

public:
	/**
	 *  @return pointer to the GODDeSS_DataStorage that is to be used
	 */
	GODDeSS_Messenger * GetGoddessMessenger() const
	{ return GoddessMessenger; }

	/**
	 *  Set whether a list of photons or the "normal" particle source is to be used for the simulation.
	 */
	void SetUsePhotonListSource(G4bool usePhotonListSource)
	{ UsePhotonListSource = usePhotonListSource; }
	/**
	 *  @return G4bool, whether a list of photons or the "normal" particle source is to be used for the simulation
	 */
	G4bool GetUsePhotonListSource() const
	{ return UsePhotonListSource; }

	/**
	 *  Set the dimensions of the scintillator tiles used in the simulation. If nothing is set, default values (defined in the DetectorConstruction) will be used.
	 */
	void SetScintillatorDimensions(G4ThreeVector scintillatorDimensions)
	{ ScintillatorDimensions = scintillatorDimensions; }
	/**
	 *  @return dimensions of the scintillator tiles to be used in the simulation (if nothing is set, default values (defined in the DetectorConstruction) will be used)
	 */
	G4ThreeVector GetScintillatorDimensions() const
	{ return ScintillatorDimensions; }

	/**
	 *  Set whether Geant4 should search for overlaps when placing the physical volumes
	 */
	void SetSearchOverlaps(G4bool searchOverlaps)
	{ SearchOverlaps = searchOverlaps; }
	/**
	 *  @return G4bool, whether Geant4 should search for overlaps when placing the physical volumes
	 */
	G4bool GetSearchOverlaps() const
	{ return SearchOverlaps; }

	/**
	 *  Set the pointer to the UserRunInformation.
	 */
	void SetRunInformation(UserRunInformation * runInformation)
	{ RunInformation = runInformation; }
	/**
	 *  @return pointer to the UserRunInformation
	 */
	UserRunInformation * GetRunInformation() const
	{ return RunInformation; }

	/**
	 *  Set the path to the file containing the scintillator tile's properties.
	 */
	void SetScintillatorPropertyFile(G4String scintillatorPropertyFile)
	{ ScintillatorPropertyFile = scintillatorPropertyFile; }
	/**
	 *  @return path to the file containing the scintillator tile's properties
	 */
	G4String GetScintillatorPropertyFile() const
	{ return ScintillatorPropertyFile; }

	/**
	 *  Set the path to the file containing the properties of the scintillator twin tile reflector.
	 */
	void SetTwinTileReflectorPropertyFile(G4String twinTileReflectorPropertyFile)
	{ TwinTileReflectorPropertyFile = twinTileReflectorPropertyFile; }
	/**
	 *  @return path to the file containing the properties of the scintillator twin tile reflector
	 */
	G4String GetTwinTileReflectorPropertyFile() const
	{ return TwinTileReflectorPropertyFile; }

	/**
	 *  Set the path to the file containing the wrapping's properties.
	 */
	void SetWrappingPropertyFile(G4String wrappingPropertyFile)
	{ WrappingPropertyFile = wrappingPropertyFile; }
	/**
	 *  @return path to the file containing the wrapping's properties
	 */
	G4String GetWrappingPropertyFile() const
	{ return WrappingPropertyFile; }

	/**
	 *  Set the path to the file containing the light-guiding fibre's properties.
	 */
	void SetLightGuidingFibrePropertyFile(G4String fibrePropertyFile)
	{ LightGuidingFibrePropertyFile = fibrePropertyFile; }
	/**
	 *  @return path to the file containing the light-guiding fibre's properties.
	 */
	G4String GetLightGuidingFibrePropertyFile() const
	{ return LightGuidingFibrePropertyFile; }

	/**
	 *  Set the path to the file containing the wavelength-shifting fibre's properties.
	 */
	void SetWLSFibrePropertyFile(G4String fibrePropertyFile)
	{ WLSFibrePropertyFile = fibrePropertyFile; }
	/**
	 *  @return path to the file containing the wavelength-shifting fibre's properties
	 */
	G4String GetWLSFibrePropertyFile() const
	{ return WLSFibrePropertyFile; }

	/**
	 *  Set the path to the file containing the scintillating fibre's properties.
	 */
	void SetScintiFibrePropertyFile(G4String fibrePropertyFile)
	{ ScintiFibrePropertyFile = fibrePropertyFile; }
	/**
	 *  @return path to the file containing the scintillating fibre's properties
	 */
	G4String GetScintiFibrePropertyFile() const
	{ return ScintiFibrePropertyFile; }

	/**
	 *  Set the path to the file containing the optical cement's properties.
	 */
	void SetOpticalCementPropertyFile(G4String opticalCementPropertyFile)
	{ OpticalCementPropertyFile = opticalCementPropertyFile; }
	/**
	 *  @return path to the file containing the optical cement's properties
	 */
	G4String GetOpticalCementPropertyFile() const
	{ return OpticalCementPropertyFile; }

	/**
	 *  Set whether optical photons is to be tracked in the simulaton. (If "false", optical photons are killed in their first step.)
	 */
	void SetTrackPhotons(G4bool trackPhotons)
	{ TrackPhotons = trackPhotons; }
	/**
	 *  @return G4bool, whether optical photons is to be tracked in the simulaton (if "false", optical photons are killed in their first step)
	 */
	G4bool GetTrackPhotons() const
	{ return TrackPhotons; }

	/**
	 *  Set path to the file containing the particle source's settings.
	 */
	void SetGPSFileName(G4String gpsFileName)
	{ GPSFile = gpsFileName; }
	/**
	 *  @return path to the file containing the particle source's settings
	 */
	G4String GetGPSFileName() const
	{ return GPSFile; }

	/**
	 *  Set path to the file which the (shortened) simulation results is to be written to.
	 */
	void SetDataFileName(G4String dataFileName)
	{ DataFile = dataFileName; }
	/**
	 *  @return path to the file which the (shortened) simulation results is to be written to
	 */
	G4String GetDataFileName() const
	{ return DataFile; }

	/**
	 *  Set path to the file which the simulation results is to be written to.
	 */
	void SetControlFileName(G4String controlFileName)
	{ ControlFile = controlFileName; }
	/**
	 *  @return path to the file which the simulation results is to be written to
	 */
	G4String GetControlFileName() const
	{ return ControlFile; }
};

#endif
