/*
 * GeneralParticleSourceMessenger.hh
 *
 * @date Jan 18, 2011
 * @author niggemann
 * @copyright GNU General Public License
 */

#ifndef GENERALPARTICLESOURCEMESSENGER_HH_
#define GENERALPARTICLESOURCEMESSENGER_HH_

#include <globals.hh>
#include <G4UImessenger.hh>
#include <G4UIcommand.hh>
#include <G4UIdirectory.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include "G4UIcmdWithADouble.hh"
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWith3Vector.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithABool.hh>

/**
 * A messenger for configuring the GeneralParticleSource.
 *
 * Commands:
 * /gps/verbose
 * /gps/particleName
 * /gps/nParticles
 * /gps/polar
 * /gps/tMin
 * /gps/tMax
 * /gps/energy/eMin
 * /gps/energy/eMax
 * /gps/energy/diceBetaGamma
 * /gps/energy/betaGammaMin
 * /gps/energy/betaGammaMax
 * /gps/angle/phiMin
 * /gps/angle/phiMax
 * /gps/angle/thetaMin
 * /gps/angle/thetaMax
 * /gps/plane/shape
 * /gps/plane/centerOrientated
 * /gps/plane/surfaceNormal
 * /gps/plane/shift
 * /gps/plane/rMin
 * /gps/plane/rMax
 * /gps/plane/a
 * /gps/plane/b
 * /gps/plane/pos
 * /gps/plane/posDist
 * /gps/plane/posGridDistA
 */
class GeneralParticleSourceMessenger: public G4UImessenger {
private:
	// Plain values
	int verbose;
	G4String particleName;
	G4String outFileName;
	int nParticles;
	double polar;
	double tMin;
	double tMax;
	double eMin;
	double eMax;
	double betaGammaMin;
	double betaGammaMax;
	bool diceBetaGamma;
	double phiMin;
	double phiMax;
	double thetaMin;
	double thetaMax;
	G4ThreeVector surfaceNormal;
	double shift;
	G4String shape;
	bool centerOrientated;
	double rMin;
	double rMax;
	double a;
	double b;
	G4ThreeVector pos;
	G4String posDist;
	double posGridDistA;
	// Commands
	G4UIcmdWithAnInteger* verboseCmd;
	G4UIcmdWithAString* outFileNameCmd;
	G4UIcmdWithAString* particleNameCmd;
	G4UIcmdWithAnInteger* nParticlesCmd;
	G4UIcmdWithADoubleAndUnit* polarCmd;
	G4UIcmdWithADoubleAndUnit* eMinCmd;
	G4UIcmdWithADoubleAndUnit* eMaxCmd;
	G4UIcmdWithADouble* betaGammaMinCmd;
	G4UIcmdWithADouble* betaGammaMaxCmd;
	G4UIcmdWithABool* diceBetaGammaCmd;
	G4UIcmdWithADoubleAndUnit* phiMinCmd;
	G4UIcmdWithADoubleAndUnit* phiMaxCmd;
	G4UIcmdWithADoubleAndUnit* thetaMinCmd;
	G4UIcmdWithADoubleAndUnit* thetaMaxCmd;
	G4UIcmdWith3VectorAndUnit* surfaceNormalCmd;
	G4UIcmdWithADoubleAndUnit* shiftCmd;
	G4UIcmdWithAString* shapeCmd;
	G4UIcmdWithABool* centerOrientatedCmd;
	G4UIcmdWithADoubleAndUnit* rMinCmd;
	G4UIcmdWithADoubleAndUnit* rMaxCmd;
	G4UIcmdWithADoubleAndUnit* aCmd;
	G4UIcmdWithADoubleAndUnit* bCmd;
	G4UIcmdWith3VectorAndUnit* posCmd;
	G4UIcmdWithAString* posDistCmd;
	G4UIcmdWithADoubleAndUnit* posGridDistACmd;
	G4UIcmdWithADoubleAndUnit* tMinCmd;
	G4UIcmdWithADoubleAndUnit* tMaxCmd;

	/**
	 * Singleton.
	 */
	static GeneralParticleSourceMessenger* INSTANCE;

	GeneralParticleSourceMessenger();

public:
	/**
	 * The file name of the macro file.
	 */
	static const std::string MACRO_FILE_NAME;
	/**
	 * String constant for the circular source shape.
	 */
	static const G4String SHAPE_CIRCLE_NAME;
	/**
	 * String constant for the rectangular source shape.
	 */
	static const G4String SHAPE_RECT_NAME;
	/**
	 * String constant for the hexagonal source shape.
	 */
	static const G4String SHAPE_HEXAGON_NAME;
	/**
	 * String constant for the uniform particle distribution
	 */
	static const G4String POS_DIST_UNIFORM;
	/**
	 * String constant for the grid point particle distribution.
	 */
	static const G4String POS_DIST_GRID;

	/**
	 * @return GeneralParticleSourceMessenger - the singleton.
	 */
	static GeneralParticleSourceMessenger* getInstance() {
		if (INSTANCE == NULL) {
			INSTANCE = new GeneralParticleSourceMessenger;
		}
		return INSTANCE;
	}

	virtual ~GeneralParticleSourceMessenger();

	void SetNewValue(G4UIcommand*, G4String);

	void setDefaultValues();

	/**
	 * @return G4String - the position distribution.
	 */
	G4String getPosDist() const {
		return posDist;
	}

	/**
	 * @param _posDist - the posDist to set.
	 */
	void setPosDist(G4String _posDist) {
		posDist = _posDist;
	}

	/**
	 * @return double - the grid pitch of the grid particle distribution.
	 */
	double getPosGridDistA() const {
		return posGridDistA;
	}

	/**
	 * @param _posGridDistA - the posGridDistA to set.
	 */
	void setPosGridDistA(double _posGridDistA) {
		posGridDistA = _posGridDistA;
	}

	/**
	 * @return int - the verbosity.
	 */
	int getVerbose() const {
		return verbose;
	}

	/**
	 * @param _verbose - the verbosity to set.
	 */
	void setVerbose(int _verbose) {
		verbose = _verbose;
	}

	/**
	 * @return G4String - the particle name.
	 */
	G4String getParticleName() const {
		return particleName;
	}

	/**
	 * @param _particleName - the particleName to set.
	 */
	void setParticleName(G4String _particleName) {
		particleName = _particleName;
	}

	G4String getOutputFileName() const {
		return outFileName;
	}

	void setOutputFileName(G4String _outFileName) {
		outFileName = _outFileName;
	}

	/**
	 * @return int - the number of particles per beamOn.
	 */
	int getNParticles() const {
		return nParticles;
	}

	/**
	 * @param _nParticles - the nParticles to set.
	 */
	void setNParticles(int _nParticles) {
		nParticles = _nParticles;
	}

	/**
	 * @return double - the polar angle of photons.
	 */
	double getPolar() const {
		return polar;
	}

	/**
	 * @param _polar - the polar to set.
	 */
	void setPolar(double _polar) {
		polar = _polar;
	}

	/**
	 * @return double - the eMax.
	 */
	double getEMax() const {
		return eMax;
	}

	/**
	 * @param _eMax - the eMax to set.
	 */
	void setEMax(double _eMax) {
		eMax = _eMax;
	}

	/**
	 * @return double - the eMin.
	 */
	double getEMin() const {
		return eMin;
	}

	/**
	 * @param _eMin - the eMin to set.
	 */
	void setEMin(double _eMin) {
		eMin = _eMin;
	}

	/**
	 * @return double - the beta * gamma minimum.
	 */
	double getBetaGammaMin() const {
		return betaGammaMin;
	}

	/**
	 * @param _betaGammaMin - the betaGammaMin to set.
	 */
	void setBetaGammaMin(double _betaGammaMin) {
		betaGammaMin = _betaGammaMin;
	}

	/**
	 * @return double - the beta * gamma maximum.
	 */
	double getBetaGammaMax() const {
		return betaGammaMax;
	}

	/**
	 * @param _betaGammaMax - the betaGammaMax to set.
	 */
	void setBetaGammaMax(double _betaGammaMax) {
		betaGammaMax = _betaGammaMax;
	}

	/**
	 * @return bool - true if beta * gamma is to be diced instead of the energy.
	 */
	bool getDiceBetaGamma() const {
		return diceBetaGamma;
	}

	/**
	 * @param _diceBetaGamma - the diceBetaGamma to set.
	 */
	void setDiceBetaGamma(bool _diceBetaGamma) {
		diceBetaGamma = _diceBetaGamma;
	}

	/**
	 * @return double - the phiMax.
	 */
	double getPhiMax() const {
		return phiMax;
	}

	/**
	 * @param _phiMax - the phiMax to set.
	 */
	void setPhiMax(double _phiMax) {
		phiMax = _phiMax;
	}

	/**
	 * @return double - the phiMin.
	 */
	double getPhiMin() const {
		return phiMin;
	}

	/**
	 * @param _phiMin - the phiMin to set.
	 */
	void setPhiMin(double _phiMin) {
		phiMin = _phiMin;
	}

	/**
	 * @return double - the thetaMax.
	 */
	double getThetaMax() const {
		return thetaMax;
	}

	/**
	 * @param _thetaMax - the thetaMax to set.
	 */
	void setThetaMax(double _thetaMax) {
		thetaMax = _thetaMax;
	}

	/**
	 * @return double - the thetaMin.
	 */
	double getThetaMin() const {
		return thetaMin;
	}

	/**
	 * @param _thetaMin - the thetaMin to set.
	 */
	void setThetaMin(double _thetaMin) {
		thetaMin = _thetaMin;
	}

	/**
	 * @return G4ThreeVector - the position of the particle plane.
	 */
	G4ThreeVector getPos() const {
		return pos;
	}

	/**
	 * @param _pos - the pos to set.
	 */
	void setPos(G4ThreeVector _pos) {
		pos = _pos;
	}

	/**
	 * @return G4ThreeVector - the surface normal.
	 */
	G4ThreeVector getSourceSurfaceNormal() const {
		return surfaceNormal;
	}

	/**
	 * @param _sourceSurfaceNormal - the sourceSurfaceNormal to set.
	 */
	void setSourceSurfaceNormal(G4ThreeVector _sourceSurfaceNormal) {
		surfaceNormal = _sourceSurfaceNormal;
	}

	/**
	 * @return double - the shift.
	 */
	double getShift() const {
		return shift;
	}

	/**
	 * @param _shift - the shift to set.
	 */
	void setShift(double _shift) {
		shift = _shift;
	}

	/**
	 * @return G4String - the shape type name.
	 */
	G4String getShape() const {
		return shape;
	}

	/**
	 * @param _shape - the shape to set.
	 */
	void setShape(G4String _shape) {
		shape = _shape;
	}

	/**
	 * @return bool - true if particles are to be orientated towards the center of the shape
	 */
	bool getCenterOrientated() const {
		return centerOrientated;
	}

	/**
	 * @param _centerOrientated - the centerOrientated to set.
	 */
	void setCenterOrientated(bool _centerOrientated) {
		centerOrientated = _centerOrientated;
	}

	/**
	 * @return double - the rMax.
	 */
	double getRMax() const {
		return rMax;
	}

	/**
	 * @param _rMax - the rMax to set.
	 */
	void setRMax(double _rMax) {
		rMax = _rMax;
	}

	/**
	 * @return double - the rMin.
	 */
	double getRMin() const {
		return rMin;
	}

	/**
	 * @param _rMin - the rMin to set.
	 */
	void setRMin(double _rMin) {
		rMin = _rMin;
	}

	/**
	 * @return double - the length of the rectangular source plane.
	 */
	double getA() const {
		return a;
	}

	/**
	 * @param _a - the a to set.
	 */
	void setA(double _a) {
		a = _a;
	}

	/**
	 * @return double - the width of the rectangular source plane.
	 */
	double getB() const {
		return b;
	}

	/**
	 * @param _b - the b to set.
	 */
	void setB(double _b) {
		b = _b;
	}

	/**
	 * @return double - the minimum particle time.
	 */
	double getTMin() {
		return tMin;
	}

	/**
	 * @param _tMin - the tMin to set.
	 */
	void setTMin(double _tMin) {
		tMin = _tMin;
	}

	/**
	 * @return double - the maximum particle time.
	 */
	double getTMax() {
		return tMax;
	}

	/**
	 * @param _tMax - the tMax to set.
	 */
	void setTMax(double _tMax) {
		tMax = _tMax;
	}
};

#endif /* GENERALPARTICLESOURCEMESSENGER_HH_ */
