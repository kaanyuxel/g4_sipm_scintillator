/*
 * GeneralParticleSource.cc
 *
 * @date Jan 18, 2011
 * @author niggemann
 * @copyright GNU General Public License
 */

#include "GeneralParticleSource.hh"

#include <G4ParticleTable.hh>
#include <G4UImanager.hh>
#include <G4RunManager.hh>
#include <G4Run.hh>
#include <G4OpticalPhoton.hh>

#include <algorithm>

#include <boost/any.hpp>

#include <Randomize.hh>

#include "GeneralParticleSourceMessenger.hh"

using namespace CLHEP;

GeneralParticleSource::GeneralParticleSource() {
	// Init gun.
	particleGun = new G4ParticleGun(1);
	// Init messenger.
	messenger = GeneralParticleSourceMessenger::getInstance();
	//
	runIdDeferred = -1;
	//
	CLHEP::HepRandom::createInstance(); // to force instantiation of static generator before any usage of HepRandom classes!
	CLHEP::HepRandom::setTheSeed(time(NULL));
}

GeneralParticleSource::~GeneralParticleSource() {
	delete messenger;
}

void GeneralParticleSource::setOptPhotonPolar() {
	setOptPhotonPolar(CLHEP::RandFlat::shoot(360. * CLHEP::deg));
}

void GeneralParticleSource::setOptPhotonPolar(G4double angle) {
	// Check particle type.
	if (particleGun->GetParticleDefinition() != G4OpticalPhoton::Definition()) {
		return;
	}
	// Calculate polarization.
	G4ThreeVector normal(1., 0., 0.);
	G4ThreeVector kphoton = particleGun->GetParticleMomentumDirection();
	G4ThreeVector product = normal.cross(kphoton);
	G4double modul2 = product * product;
	G4ThreeVector e_perpend(0., 0., 1.);
	if (modul2 > 0.) {
		e_perpend = (1. / std::sqrt(modul2)) * product;
	}
	G4ThreeVector e_paralle = e_perpend.cross(kphoton);
	G4ThreeVector polar = std::cos(angle) * e_paralle + std::sin(angle) * e_perpend;
	particleGun->SetParticlePolarization(polar);
	//	std::cout << "GeneralParticleSource::setOptPhotonPolar(" << angle / CLHEP::deg << " CLHEP::deg): polarization vector is " << polar
	//			<< "." << std::endl;
}

void GeneralParticleSource::dicePos(G4double& x, G4double& y) {
	if (messenger->getShape() == GeneralParticleSourceMessenger::SHAPE_CIRCLE_NAME) {
		return dicePosInCircle(x, y);
	}
	if (messenger->getShape() == GeneralParticleSourceMessenger::SHAPE_RECT_NAME) {
		return dicePosInRect(x, y);
	}
	if (messenger->getShape() == GeneralParticleSourceMessenger::SHAPE_HEXAGON_NAME) {
		return dicePosInHexagon(x, y);
	}
	G4cerr << "GeneralParticleSource::dicePos(x, y): Shape " << messenger->getShape() << " not supported." << G4endl;
	G4Exception("GeneralParticleSource::dicePos(x, y)", "InvalidSetup", FatalException, "Shape not recognized.");
}

void GeneralParticleSource::dicePosInCircle(G4double& x, G4double& y) {
	x = -1;
	y = -1;
	G4double r = -1;
	while (r > messenger->getRMax() || r < messenger->getRMin()) {
		dicePosForDistribution(x, y, messenger->getRMax(), messenger->getRMax());
		r = sqrt(x * x + y * y);
	}
}

void GeneralParticleSource::dicePosInRect(G4double& x, G4double& y) {
	dicePosForDistribution(x, y, messenger->getA() / 2., messenger->getB() / 2.);
}

void GeneralParticleSource::dicePosInHexagon(G4double& x, G4double& y) {
	G4double a = messenger->getA() / 2.;
	G4double ri = sqrt(3.) / 2. * a;
	// Y intercept for a straight line with slope -60 deg.
	G4double m = tan(CLHEP::pi / 3.) * a;
	x = a;
	y = ri;
	// Check if the point is within the hexagon.
	while (-tan(CLHEP::pi / 3.) * fabs(x) + m < fabs(y)) {
		dicePosForDistribution(x, y, a, ri);
	}
}

void GeneralParticleSource::dicePosForDistribution(G4double& x, G4double& y, G4double xMax, G4double yMax) {
	if (messenger->getPosDist() == GeneralParticleSourceMessenger::POS_DIST_UNIFORM) {
		// Dice uniform.
		x = CLHEP::RandFlat::shoot(-xMax, xMax);
		y = CLHEP::RandFlat::shoot(-yMax, yMax);
		return;
	} else if (messenger->getPosDist() == GeneralParticleSourceMessenger::POS_DIST_GRID) {
		// Dice for grid distribution.
		const G4double a = messenger->getPosGridDistA();
		// Determine number of grid points.
		int nGridX = (int) floor(2. * xMax / a);
		int nGridY = (int) floor(2. * yMax / a);
		// Dice grid position.
		int gridX = (int) round(CLHEP::RandFlat::shoot(nGridX));
		int gridY = (int) round(CLHEP::RandFlat::shoot(nGridY));
		// Calculate coordinates.
		x = gridX * a - xMax;
		y = gridY * a - yMax;
		return;
	}
	G4cerr << "GeneralParticleSource::dicePosForDistribution(x, y, xMax, yMax): Distribution "
			<< messenger->getPosDist() << " not supported." << G4endl;
	G4Exception("GeneralParticleSource::dicePosForDistribution(x, y, xMax, yMax)", "InvalidSetup", FatalException,
			"Distribution not recognized.");
}

void GeneralParticleSource::diceTimingsIfNecessary() {
	// Check whether run id has changed.
	const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
	if (runIdDeferred == run->GetRunID() && dicedTimes.size() >= (unsigned) messenger->getNParticles()) {
		return;
	}
	runIdDeferred = run->GetRunID();
	dicedTimes.clear();
	// Get number of events to be fired.
	int nEvents = messenger->getNParticles();
	// Dice nEvents times.
	const double tMax = messenger->getTMax();
	const double tMin = messenger->getTMin();
	for (int i = 0; i < nEvents; i++) {
		dicedTimes.push_back(CLHEP::RandFlat::shoot(tMin, tMax));
	}
	// Sort by value to have them ordered.
	std::sort(dicedTimes.begin(), dicedTimes.end());
	// Reverse order.
	std::reverse(dicedTimes.begin(), dicedTimes.end());
}

void GeneralParticleSource::GeneratePrimaries(G4Event* event) {
	G4String outFileName = messenger->getOutputFileName();
	OutFile.open(outFileName.c_str(), std::ios_base::out | std::ios_base::app);

	// Persist current settings.
	const int nParticles = messenger->getNParticles();
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

	diceTimingsIfNecessary();

	if (event->GetEventID() == 0)
	{
		OutFile << "GPS data\n" << "--------------------\n" << "\n";
		OutFile << "RunID:\t\t" << runIdDeferred << "\n";
	}
	OutFile << "\n" << "EventID:\t" << event->GetEventID() << "\n" << "\n";
	//
	// Loop.
	for (int i = 0; i < nParticles; i++) {
		// Set particle type.
		G4String particleName = messenger->getParticleName();
		// If particle name does not exist ...
		if (!particleTable->contains(particleName)) {
			// ... the particle name does not contain a charge specification ...
			if (!(particleName.contains("+") || particleName.contains("-"))) {
				// ... but particle name + charge specification exists ...
				if (particleTable->contains(particleName + "+") || particleTable->contains(particleName + "-")) {
					// ... alternate the charges.
					bool hasEvenEventID = event->GetEventID() % 2;
					if (hasEvenEventID) {
						particleName += "+";
					} else {
						particleName += "-";
					}
				}
			}
		}
		OutFile << "particleName:\t\t\t\t" << particleName << "\n";
		OutFile << "particleID:\t\t\t\t" << particleTable->FindParticle(particleName)->GetPDGEncoding() << "\n";
		particleGun->SetParticleDefinition(particleTable->FindParticle(particleName));
		// Dice energy.
		G4double eMin = messenger->getEMin();
		G4double eMax = messenger->getEMax();
		G4double eMin_default = 1 * eV;
		G4double eMax_default = 10 * eV;
		G4bool diceBetaGamma = messenger->getDiceBetaGamma();
		G4double betaGammaMin = messenger->getBetaGammaMin();
		G4double betaGammaMax = messenger->getBetaGammaMax();
		G4double betaGammaMin_default = 3;
		G4double betaGammaMax_default = 4;
		G4double mass = particleGun->GetParticleDefinition()->GetPDGMass();
		if (fabs(mass) > 1e-12) {
			if(diceBetaGamma) {
				if (isnan(betaGammaMin)) {
					if (!isnan(eMin)) betaGammaMin = sqrt( (eMin / mass + 1.) * (eMin / mass + 1.) - 1. );
					else betaGammaMin = betaGammaMin_default;
				}
				if (isnan(betaGammaMax)) {
					if (!isnan(eMax)) betaGammaMax = sqrt( (eMax / mass + 1.) * (eMax / mass + 1.) - 1. );
					else betaGammaMax = betaGammaMax_default;
				}
			}
			else {
				if (isnan(eMin)) {
					if (!isnan(betaGammaMin)) eMin = (sqrt(betaGammaMin * betaGammaMin + 1.) - 1.) * mass;
					else eMin = eMin_default;
				}
				if (isnan(eMax)) {
					if (!isnan(betaGammaMax)) eMax = (sqrt(betaGammaMax * betaGammaMax + 1.) - 1.) * mass;
					else eMax = eMax_default;
				}
			}
		}
		else {
			diceBetaGamma = false;
			if (isnan(eMin)) eMin = eMin_default;
			if (isnan(eMax)) eMax = eMax_default;
		}
		G4double energy = 0.;
		// If it is a massiv particle and beta * gamma is to be diced ...
		if(fabs(mass) > 1e-12 && diceBetaGamma) {
			//... dice beta * gamma ...
			G4double betaGamma = CLHEP::RandFlat::shoot(betaGammaMin, betaGammaMax);
			// ... determine the energy by beta * gamma.
			energy = ( sqrt( pow(betaGamma, 2) + 1.) - 1.) * mass;
			OutFile << "initial_energy_/_MeV:\t\t\t" << energy / CLHEP::MeV << "\n";
			OutFile << "initial_beta*gamma:\t\t\t" << betaGamma << "\n";
		}
		// Else ...
		else {
			//... dice the energy.
			energy = CLHEP::RandFlat::shoot(eMin, eMax);
			OutFile << "initial_energy_/_MeV:\t\t\t" << energy / CLHEP::MeV << "\n";
			if (fabs(mass) > 1e-12) {
				G4double betaGamma = sqrt( pow(energy / mass + 1., 2) - 1. );
				OutFile << "initial_beta*gamma:\t\t\t" << betaGamma << "\n";
			}
		}
		particleGun->SetParticleEnergy(energy);
		// The particle position is diced in the x-y-plane by default, cf. below.
		G4ThreeVector dicingPlaneNormal(0, 0, 1);
		// Determine rotation parameters for given surface normal.
		G4ThreeVector sourcePlaneNormal = messenger->getSourceSurfaceNormal();
		G4ThreeVector rotationAxis = dicingPlaneNormal.cross(sourcePlaneNormal);
		G4double rotationAngle = 0.;
		if (fabs(rotationAxis.mag()) > 1e-12) {
			rotationAngle = acos(
					dicingPlaneNormal * sourcePlaneNormal / (dicingPlaneNormal.mag() * sourcePlaneNormal.mag()));
		} else if (sourcePlaneNormal == -dicingPlaneNormal) {
			rotationAngle = 180. * CLHEP::deg;
		}
		// Dice position.
		G4double x;
		G4double y;
		dicePos(x, y);
		G4double z = 0.;
		if(messenger->getShift() != 0.){
			z = CLHEP::RandFlat::shoot(-messenger->getShift(), messenger->getShift());
		}

		OutFile << "x_position_in_source_plain_/_mm:\t" << x / CLHEP::mm << "\n";
		OutFile << "y_position_in_source_plain_/_mm:\t" << y / CLHEP::mm << "\n";
		// Rotate particle position.
		G4ThreeVector particlePos(x, y, z);
		if (rotationAngle == 180. * CLHEP::deg) {
			particlePos *= -1;
		} else if (fabs(rotationAngle) > 1e-12) {
			particlePos.rotate(rotationAxis, rotationAngle);
		}
		particlePos += messenger->getPos();
		OutFile << "global_initial_particle_position_/_mm:\t" << particlePos / CLHEP::mm << "\n";
		particleGun->SetParticlePosition(particlePos);
		// Dice momentum regarding the cosine law.
		G4double phi = CLHEP::RandFlat::shoot(messenger->getPhiMin(), messenger->getPhiMax());
		OutFile << "momentum_direction_phi_/_deg:\t\t" << phi / CLHEP::deg << "\n";

		G4ThreeVector momentumDirection;
		if(messenger->getCenterOrientated() && messenger->getShape() == GeneralParticleSourceMessenger::SHAPE_CIRCLE_NAME){
			G4double theta = CLHEP::RandFlat::shoot(messenger->getThetaMin(), messenger->getThetaMax());

			G4ThreeVector vectorTowardsCenter = G4ThreeVector(0., 0., z);
			if (rotationAngle == 180. * CLHEP::deg) {
				vectorTowardsCenter *= -1;
			} else if (fabs(rotationAngle) > 1e-12) {
				vectorTowardsCenter.rotate(rotationAxis, rotationAngle);
			}
			vectorTowardsCenter += messenger->getPos() - particlePos;
			vectorTowardsCenter /= vectorTowardsCenter.mag();
			vectorTowardsCenter.rotate(sourcePlaneNormal, theta);
			vectorTowardsCenter.rotate(vectorTowardsCenter.cross(sourcePlaneNormal), phi);
			momentumDirection = vectorTowardsCenter;
			if (rotationAngle == 180. * CLHEP::deg) {
				vectorTowardsCenter *= -1;
			} else if (fabs(rotationAngle) > 1e-12) {
				vectorTowardsCenter.rotate(-rotationAxis, rotationAngle);
			}
			OutFile << "momentum_direction:\t\t\t" << vectorTowardsCenter << "\n";
			OutFile << "global_momentum_direction:\t\t" << momentumDirection << "\n";
		}else{
			G4double cosTheta = CLHEP::RandFlat::shoot(cos(messenger->getThetaMin()), cos(messenger->getThetaMax()));
			OutFile << "momentum_direction_cosTheta:\t\t" << cosTheta << "\n";
			G4double sinTheta = sqrt(1. - cosTheta * cosTheta);

			momentumDirection.set(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
			OutFile << "momentum_direction:\t\t\t" << momentumDirection << "\n";
			// Rotate momentum direction.
			if (rotationAngle == 180. * CLHEP::deg) {
				momentumDirection *= -1;
			} else if (fabs(rotationAngle) > 1e-12) {
				momentumDirection.rotate(rotationAxis, rotationAngle);
			}
			OutFile << "global_momentum_direction:\t\t" << momentumDirection << "\n";
		}
		particleGun->SetParticleMomentumDirection(momentumDirection);
		// Set polarization.
		if (messenger->getPolar() >= 0) {
			setOptPhotonPolar(messenger->getPolar());
		} else {
			// Dice polarization.
			setOptPhotonPolar();
			// particleGun->SetParticlePolarization(G4ThreeVector(0, 0, 0));
		}
		// The last element is the smallest one.
		particleGun->SetParticleTime(dicedTimes.back());
		dicedTimes.pop_back();
		// Fire!
		particleGun->SetNumberOfParticles(1);
		particleGun->GeneratePrimaryVertex(event);

		OutFile << "\n";
	}

	OutFile.close();
}
