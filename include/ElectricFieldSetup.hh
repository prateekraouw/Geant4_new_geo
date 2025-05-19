// ElectricFieldSetup.hh
#ifndef ElectricFieldSetup_h
#define ElectricFieldSetup_h 1

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"  // Changed from G4EqMagElectricField
#include "G4MagIntegratorStepper.hh"
#include "G4VIntegrationDriver.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "globals.hh"

class ElectricFieldSetup
{
public:
  ElectricFieldSetup();
  ~ElectricFieldSetup();
  
  void SetMagneticField(G4ThreeVector);
  G4FieldManager* GetFieldManager() { return fFieldManager; }
  
private:
  G4FieldManager*            fFieldManager;
  G4ChordFinder*             fChordFinder;
  G4Mag_UsualEqRhs*          fEquation;  // Changed from G4EqMagElectricField
  G4MagneticField*           fMagneticField;
  G4MagIntegratorStepper*    fStepper;
  G4VIntegrationDriver*      fIntgrDriver;
  
  G4double                   fMinStep;
};

#endif