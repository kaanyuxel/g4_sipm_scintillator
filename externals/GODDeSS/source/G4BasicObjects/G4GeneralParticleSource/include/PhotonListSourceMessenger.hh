/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef PhotonListSourceMessenger_h
#define PhotonListSourceMessenger_h 1

#include <G4UImessenger.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithAString.hh>



// class variables begin with capital letters, local variables with small letters



///  a class for operating the PhotonListSource and providing the necessary information for this (taken from compressed (.gz) text files)
class PhotonListSourceMessenger: public G4UImessenger
{
public:

	/**
	 *  Constructor:
	 *  - creates the commands for operating the PhotonListSource in the interactive mode
	 */
	PhotonListSourceMessenger()
	{
		PhotonListSourceDir = new G4UIdirectory("/pls/");
		PhotonListSourceDir->SetGuidance("UI commands of this setup");

		InfilePathCmd = new G4UIcmdWithAString("/pls/photonListPath", this);
		InfilePathCmd->SetGuidance("Set the path to the compressed (.gz) text file containing the list of photons that is to be processed.");
	}

	/**
	 *  Destructor:
	 *  - deletes the command objects
	 */
	~PhotonListSourceMessenger()
	{
		delete InfilePathCmd;
		delete PhotonListSourceDir;
	}



	void SetNewValue(G4UIcommand*, G4String);



	/**
	 *  @return G4bool, if a list with events is available
	 */
	G4bool IsEventAvailable()
	{ return (G4bool) ParticleInformations.size(); }

	/**
	 *  @return G4bool, if a list with start conditions of optical photons is available for the next event
	 */
	G4bool IsPhotonDataAvailable()
	{
		if(IsEventAvailable()) return (G4bool) (ParticleInformations[0]).size();
		else return false;
	}

	/**
	 *  @return start time of the next optical photon to be generated
	 */
	G4double GetNextPhotonStartTime_ns()
	{ return (((ParticleInformations[0])[0])[0]).x(); }

	/**
	 *  @return start position of the next optical photon to be generated
	 */
	G4ThreeVector GetNextPhotonStartPosition_mm()
	{ return ((ParticleInformations[0])[0])[1]; }

	/**
	 *  @return start momentum of the next optical photon to be generated
	 */
	G4ThreeVector GetNextPhotonStartMomentum_MeV()
	{ return ((ParticleInformations[0])[0])[2]; }

	/**
	 *  @return start polarisation of the next optical photon to be generated
	 */
	G4ThreeVector GetNextPhotonStartPolarisation()
	{ return ((ParticleInformations[0])[0])[3]; }

	/**
	 *  @return remove the data of the next optical photon from the list
	 */
	void RemoveNextPhotonData()
	{ (ParticleInformations[0]).erase((ParticleInformations[0]).begin()); }

	/**
	 *  @return remove the event from the list
	 */
	void RemoveEvent()
	{ ParticleInformations.erase(ParticleInformations.begin()); }


	std::vector< std::vector< std::vector<G4ThreeVector> > > ParticleInformations;
private:
	void ReadInPhotonList( G4String infilePath	/**< path to the compressed (.gz) text file containing the optical photon's start conditions */
			     );

	G4bool getIntValue(const char * & data, G4String keyword, G4int & output);
	G4bool getDoubleValue(const char * & data, G4String keyword, G4double & output);
	G4bool getVectorValue(const char * & data, G4String keyword, G4ThreeVector & output);
	void relocateCharArrayPointer(const char * & charArray, const char * charArrayEnd, G4String relocationDestination);
	void createNewEvent();
	void savePhotonData();

	void setDefaultValues();


	G4UIdirectory* PhotonListSourceDir;
	G4UIcmdWithAString* InfilePathCmd;

	G4String InfilePath;

	G4int NumEvent;
	G4double Start_time_ns;
	G4bool Start_time_ns_set;
	G4ThreeVector Start_position_mm;
	G4bool Start_position_mm_set;
	G4ThreeVector Start_momentum_MeV;
	G4bool Start_momentum_MeV_set;
	G4ThreeVector Start_polarisation;
	G4bool Start_polarisation_set;
};

#endif
