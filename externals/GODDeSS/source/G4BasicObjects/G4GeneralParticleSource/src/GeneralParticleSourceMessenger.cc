/*
 * GeneralParticleSourceMessenger.cc
 *
 * @date Jan 18, 2011
 * @author niggemann
 * @copyright GNU General Public License
 */

#include "GeneralParticleSourceMessenger.hh"

#include <G4SystemOfUnits.hh>
#include <G4UImanager.hh>
#include <G4OpticalPhoton.hh>

#include "G4UiMessengerUtil.hh"
#include "GeneralParticleSource.hh"

GeneralParticleSourceMessenger* GeneralParticleSourceMessenger::INSTANCE = NULL;
const std::string GeneralParticleSourceMessenger::MACRO_FILE_NAME = "gps.mac";
const G4String GeneralParticleSourceMessenger::SHAPE_CIRCLE_NAME = "circle";
const G4String GeneralParticleSourceMessenger::SHAPE_RECT_NAME = "rect";
const G4String GeneralParticleSourceMessenger::SHAPE_HEXAGON_NAME = "hexagon";
const G4String GeneralParticleSourceMessenger::POS_DIST_UNIFORM = "uniform";
const G4String GeneralParticleSourceMessenger::POS_DIST_GRID = "grid";


GeneralParticleSourceMessenger::GeneralParticleSourceMessenger() :
		G4UImessenger() {
	setDefaultValues();
	// Create directories.
	new G4UIdirectory("/gps/");
	new G4UIdirectory("/gps/energy/");
	new G4UIdirectory("/gps/angle/");
	new G4UIdirectory("/gps/plane/");
	// Create basic commands.
	verboseCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "verbose", verbose);
	outFileNameCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "outFile", outFileName);
	particleNameCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "particle", particleName);
	nParticlesCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "nParticles", nParticles);
	polarCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "polar", polar, "deg");
	tMinCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "tMin", tMin, "s");
	tMaxCmd = G4UiMessengerUtil::createCmd(this, "/gps/", "tMax", tMax, "s");
	// Create energy commands.
	eMinCmd = G4UiMessengerUtil::createCmd(this, "/gps/energy/", "eMin", eMin, "eV");
	eMaxCmd = G4UiMessengerUtil::createCmd(this, "/gps/energy/", "eMax", eMax, "eV");
	diceBetaGammaCmd = G4UiMessengerUtil::createCmd(this, "/gps/energy/", "diceBetaGamma", diceBetaGamma);
	betaGammaMinCmd = G4UiMessengerUtil::createCmd(this, "/gps/energy/", "betaGammaMin", betaGammaMin);
	betaGammaMaxCmd = G4UiMessengerUtil::createCmd(this, "/gps/energy/", "betaGammaMax", betaGammaMax);
	// Create angle commands.
	phiMinCmd = G4UiMessengerUtil::createCmd(this, "/gps/angle/", "phiMin", phiMin, "deg");
	phiMaxCmd = G4UiMessengerUtil::createCmd(this, "/gps/angle/", "phiMax", phiMax, "deg");
	thetaMinCmd = G4UiMessengerUtil::createCmd(this, "/gps/angle/", "thetaMin", thetaMin, "deg");
	thetaMaxCmd = G4UiMessengerUtil::createCmd(this, "/gps/angle/", "thetaMax", thetaMax, "deg");
	// Create plane commands.
	surfaceNormalCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "surfaceNormal", surfaceNormal, "mm");
	shiftCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "shift", shift, "mm");
	shapeCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "shape", shape);
	centerOrientatedCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "centerOrientated", centerOrientated);
	rMinCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "rMin", rMin, "mm");
	rMaxCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "rMax", rMax, "mm");
	aCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "a", a, "mm");
	bCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "b", b, "mm");
	posCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "pos", pos, "mm");
	posDistCmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "posDist", posDist);
	// Set guidance.
	polarCmd->SetGuidance("Set linear polarization angle w.r.t. (k,n) plane. "
			"Negative values will trigger a random polarization.");
	// Set parameter boundaries.
	verboseCmd->SetParameterName("verbose", false);
	verboseCmd->SetRange("verbose >= 0 && verbose <= 2");
	nParticlesCmd->SetParameterName("nParticles", false);
	nParticlesCmd->SetRange("nParticles > 0");
	betaGammaMinCmd->SetParameterName("betaGammaMin", false);
	betaGammaMinCmd->SetRange("betaGammaMin == -1 || betaGammaMin >= 0");
	betaGammaMaxCmd->SetParameterName("betaGammaMax", false);
	betaGammaMaxCmd->SetRange("betaGammaMax == -1 || betaGammaMax >= 0");
	// Set parameter candidates.
	shapeCmd->SetCandidates(
			(SHAPE_CIRCLE_NAME + G4String(" ") + SHAPE_RECT_NAME + G4String(" ") + SHAPE_HEXAGON_NAME).data());
	posGridDistACmd = G4UiMessengerUtil::createCmd(this, "/gps/plane/", "posGridDistA", posGridDistA, "mm");
	posDistCmd->SetCandidates((POS_DIST_UNIFORM + G4String(" ") + POS_DIST_GRID).data());
	// Initialize from macro.
	G4UiMessengerUtil::executeMacro(MACRO_FILE_NAME, true);
}

GeneralParticleSourceMessenger::~GeneralParticleSourceMessenger() {
	delete verboseCmd;
	delete outFileNameCmd;
	delete particleNameCmd;
	delete nParticlesCmd;
	delete polarCmd;
	delete eMinCmd;
	delete eMaxCmd;
	delete betaGammaMinCmd;
	delete betaGammaMaxCmd;
	delete diceBetaGammaCmd;
	delete phiMinCmd;
	delete phiMaxCmd;
	delete thetaMinCmd;
	delete thetaMaxCmd;
	delete surfaceNormalCmd;
	delete shiftCmd;
	delete shapeCmd;
	delete centerOrientatedCmd;
	delete rMinCmd;
	delete rMaxCmd;
	delete aCmd;
	delete bCmd;
	delete posCmd;
	delete posDistCmd;
	delete posGridDistACmd;
	delete tMinCmd;
	delete tMaxCmd;
}

void GeneralParticleSourceMessenger::setDefaultValues() {
	verbose = 0;
	particleName = G4OpticalPhoton::Definition()->GetParticleName();
	outFileName = "outputParticleSource.data";
	nParticles = 1;
	polar = -360 * CLHEP::deg;
	eMin = NAN;
	eMax = NAN;
	betaGammaMin = NAN;
	betaGammaMax = NAN;
	diceBetaGamma = false;
	phiMin = 0;
	phiMax = 360 * CLHEP::deg;
	thetaMin = 0;
	thetaMax = 90 * CLHEP::deg;
	surfaceNormal = G4ThreeVector(0., 0., -1.);
	shift = 0. * mm;
	shape = SHAPE_CIRCLE_NAME;
	centerOrientated = false;
	rMin = 0;
	rMax = 0.5 * CLHEP::mm;
	a = 1 * CLHEP::mm;
	b = 1 * CLHEP::mm;
	pos = G4ThreeVector(0, 0, 5 * cm);
	posDist = POS_DIST_UNIFORM;
	posGridDistA = 1 * CLHEP::mm;
	tMin = 0;
	tMax = 0;
}

void GeneralParticleSourceMessenger::SetNewValue(G4UIcommand* cmd, G4String value) {
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, verboseCmd, value, &verbose);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, outFileNameCmd, value, &outFileName);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, particleNameCmd, value, &particleName);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, nParticlesCmd, value, &nParticles);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, polarCmd, value, &polar);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, tMinCmd, value, &tMin);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, tMaxCmd, value, &tMax);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, eMinCmd, value, &eMin);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, eMaxCmd, value, &eMax);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, betaGammaMinCmd, value, &betaGammaMin);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, betaGammaMaxCmd, value, &betaGammaMax);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, diceBetaGammaCmd, value, &diceBetaGamma);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, phiMinCmd, value, &phiMin);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, phiMaxCmd, value, &phiMax);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, thetaMinCmd, value, &thetaMin);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, thetaMaxCmd, value, &thetaMax);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, surfaceNormalCmd, value, &surfaceNormal);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, shiftCmd, value, &shift);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, shapeCmd, value, &shape);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, centerOrientatedCmd, value, &centerOrientated);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, rMinCmd, value, &rMin);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, rMaxCmd, value, &rMax);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, aCmd, value, &a);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, bCmd, value, &b);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, posCmd, value, &pos);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, posDistCmd, value, &posDist);
	G4UiMessengerUtil::setNewValueIfCmdMatches(cmd, posGridDistACmd, value, &posGridDistA);
}
