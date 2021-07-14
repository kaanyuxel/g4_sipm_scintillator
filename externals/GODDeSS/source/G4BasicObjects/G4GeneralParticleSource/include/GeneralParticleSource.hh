/*
 * GeneralParticleSource.hh
 *
 * @date Jan 18, 2011
 * @author niggemann
 * @copyright GNU General Public License
 */

#ifndef GENERALPARTICLESOURCE_HH_
#define GENERALPARTICLESOURCE_HH_

#include <globals.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4Event.hh>
#include <G4ParticleGun.hh>

#include <fstream>

#include "GeneralParticleSourceMessenger.hh"

/**
 * The general particle source follows the example of its Geant4 pendant. The type, energy-, spacial- and
 * angular distribution of the particles shoot by the source can be manipulated via the GeneralParticleSourceMessenger
 * class.
 */
class GeneralParticleSource: public G4VUserPrimaryGeneratorAction {
private:
	G4ParticleGun* particleGun;
	GeneralParticleSourceMessenger* messenger;
	G4int runIdDeferred;
	std::vector<G4double> dicedTimes;
	std::ofstream OutFile;

	void dicePos(G4double& x, G4double& y);
	void dicePosInCircle(G4double& x, G4double& y);
	void dicePosInRect(G4double& x, G4double& y);
	void dicePosInHexagon(G4double& x, G4double& y);
	void dicePosForDistribution(G4double& x, G4double& y, G4double xMax, G4double yMax);

protected:
	/**
	 * Dices the photon's polarization.
	 */
	void setOptPhotonPolar();

	/**
	 * Sets the photon's polarization.
	 */
	void setOptPhotonPolar(G4double);

	/**
	 * Refreshes the cache of the particle times.
	 */
	void diceTimingsIfNecessary();

public:
	GeneralParticleSource();

	virtual ~GeneralParticleSource();

	virtual void GeneratePrimaries(G4Event* anEvent);

	/**
	 * GeneralParticleSourceMessenger - the messenger instance.
	 */
	GeneralParticleSourceMessenger* getMessenger() const {
		return messenger;
	}

};

#endif /* GENERALPARTICLESOURCE_HH_ */
