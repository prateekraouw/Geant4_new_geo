#include "PrimaryGeneratorAction.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(nullptr)
{
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);

  // Define the beam as protons
  G4ParticleDefinition* particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fParticleGun->SetParticleDefinition(particleDefinition);

  // Set initial energy of proton beam (8 GeV)
  fParticleGun->SetParticleEnergy(8.0*GeV);
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Set beam direction along the z-axis (towards the tungsten block)
  G4double angle = 0.0*rad;
  G4ThreeVector dir(0., std::sin(angle), std::cos(angle));
  fParticleGun->SetParticleMomentumDirection(dir);
  
  
  // Gaussian beam defination
  G4double fSigma = 5*mm;
  G4double x0 = CLHEP::RandGauss::shoot(0.0, fSigma);
  G4double y0 = CLHEP::RandGauss::shoot(0.0, fSigma);
  G4double z0 = -1.8*m;
  // Position the beam 1 cm before the tungsten block
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0,z0));

  // Generate the primary vertex
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
