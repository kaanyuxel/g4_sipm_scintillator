/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef PHOTONLISTSOURCE_HH_
#define PHOTONLISTSOURCE_HH_

#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4ParticleGun.hh>

#include <PhotonListSourceMessenger.hh>



// class variables begin with capital letters, local variables with small letters



///  a class making it possible to generate single optical photons with  (specified in a compressed (.gz) text file and read in via PhotonListSourceMessenger)
class PhotonListSource: public G4VUserPrimaryGeneratorAction
{
public:

	/**
	 *  Constructor:
	 *  - creates an object of the G4ParticleGun and the PhotonListSourceMessenger classes
	 */
	PhotonListSource()
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: ParticleGun(new G4ParticleGun(1))
	, Messenger(new PhotonListSourceMessenger())
	{
	}

	/**
	 *  Destructor (empty)
	 */
	~PhotonListSource()
	{
	}



	virtual void GeneratePrimaries(G4Event* anEvent);

private:
	G4ParticleGun* ParticleGun;

	PhotonListSourceMessenger * Messenger;
};

#endif /* PHOTONLISTSOURCE_HH_ */
