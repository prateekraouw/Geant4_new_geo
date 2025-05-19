// Include necessary headers for magnetic field
#include "ElectricFieldSetup.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagneticField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SystemOfUnits.hh"  // This will include the units

// Define a custom limited-region magnetic field class
class LimitedRegionField : public G4MagneticField {
public:
    LimitedRegionField(G4double zMax, G4double zMin) : fZMax(zMax), fZMin(zMin) {}
    
    virtual void GetFieldValue(const G4double point[4], G4double* field) const override
    {
        // Initialize field to zero
        field[0] = 0.0;
        field[1] = 0.0;
        field[2] = 0.0;
        field[3] = 0.0;
        field[4] = 0.0;
        field[5] = 0.0;
        
        // Check if the point is within our defined region (zMin to zMax)
        if (point[2] <= fZMax && point[2] >= fZMin) {
            field[2] = 7.0 * CLHEP::tesla; // 7 Tesla in +z direction
        }
    }
    
private:
    G4double fZMax; // Upper Z limit of the field
    G4double fZMin; // Lower Z limit of the field
};

// Constructor for ElectricFieldSetup
ElectricFieldSetup::ElectricFieldSetup() {
    // Create a limited region field that extends from z=-1 to z=10
    G4double zMax = 10.0 * CLHEP::m;  // Upper limit (5 meters)
    G4double zMin = -1.0 * CLHEP::m; // Lower limit (-5 meters)
    fMagneticField = new LimitedRegionField(zMax, zMin);
    
    // Get the global field manager
    fFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
    // Create equation of motion for this field
    G4Mag_UsualEqRhs* equation = new G4Mag_UsualEqRhs(fMagneticField);
    
    // Create stepper with higher precision
    fStepper = new G4ClassicalRK4(equation);
    
    // Create the chord finder with a smaller min step for better precision at boundary
    G4double minStep = 0.005 * CLHEP::mm; // Smaller minimum step for better precision
    G4MagInt_Driver* driver = new G4MagInt_Driver(minStep,
                                                 fStepper,
                                                 fStepper->GetNumberOfVariables());
    fChordFinder = new G4ChordFinder(driver);
    
    // Increase accuracy requirements for better tracking at field boundary
    fFieldManager->SetDeltaOneStep(0.001 * CLHEP::mm);  // Precision in one step
    fFieldManager->SetDeltaIntersection(0.001 * CLHEP::mm);  // Precision in intersection
    
    // Set the field manager's parameters
    fFieldManager->SetChordFinder(fChordFinder);
    fFieldManager->SetDetectorField(fMagneticField);
    
    G4cout << "Magnetic field of 7 Tesla in +z direction created, limited to region from "
           << zMin/CLHEP::m << " to " << zMax/CLHEP::m << " meters along z-axis" << G4endl;
    G4cout << "Particles will maintain their direction after leaving the field region" << G4endl;
}

// Destructor for ElectricFieldSetup
ElectricFieldSetup::~ElectricFieldSetup() {
    delete fChordFinder;
    delete fStepper;
    delete fMagneticField;
}

// Method to update the magnetic field
void ElectricFieldSetup::SetMagneticField(G4ThreeVector fieldVector) {
    // This method is kept for potential future field value changes
    // It's not needed for the current implementation
    
    G4cout << "Magnetic field set to ("
          << fieldVector.x()/CLHEP::tesla << " "
          << fieldVector.y()/CLHEP::tesla << " "
          << fieldVector.z()/CLHEP::tesla << ") Tesla" << G4endl;
}