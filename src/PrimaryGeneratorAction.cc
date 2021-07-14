#include "PrimaryGeneratorAction.hh"

#include <G4ParticleTable.hh>
#include <G4Event.hh>
#include <G4SystemOfUnits.hh>
#include <G4ParticleGun.hh>
#include <Randomize.hh>

#include <G4GeneralParticleSource.hh>

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
    fGun = new G4ParticleGun();
    G4ParticleDefinition* electron = G4ParticleTable::GetParticleTable()->FindParticle("mu-");

    fGun->SetParticleDefinition(electron);
    fGun->SetParticleEnergy(5. * GeV); 
    fGun->SetParticlePosition(G4ThreeVector(0., 3.*m, 0.));  // along y
    fGun->SetParticleMomentumDirection(G4ThreeVector(0., -1., 0.));  // along y
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
    delete fGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    fGun->GeneratePrimaryVertex(anEvent);
}
