#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class ElectricFieldSetup;  // Rename as needed but keep using this for magnetic field
class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction() override;
    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;  // Add this declaration
    
    // Methods to get the scoring volumes
    G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
    G4LogicalVolume* GetDetector1Volume() const { return fDetector1Volume; }
    G4LogicalVolume* GetDetector2Volume() const { return fDetector2Volume; }
    G4LogicalVolume* GetDetector3Volume() const { return fDetector3Volume; }
    G4LogicalVolume* GetDetector4Volume() const { return fDetector4Volume; }
    G4LogicalVolume* GetRFCavityVolume() const { return fRFCavityVolume; } 
    
    // Position getter methods
    G4ThreeVector GetDetector1Position() const { return fDetector1Position; }
    G4ThreeVector GetDetector2Position() const { return fDetector2Position; }
    G4ThreeVector GetDetector3Position() const { return fDetector3Position; }
    G4ThreeVector GetDetector4Position() const { return fDetector4Position; }
    
  private:
    G4LogicalVolume* fScoringVolume;
    G4LogicalVolume* fDetector1Volume;
    G4LogicalVolume* fDetector2Volume;
    G4LogicalVolume* fDetector3Volume;
    G4LogicalVolume* fDetector4Volume;
    G4LogicalVolume* fRFCavityVolume;  // Add this
    
    ElectricFieldSetup* fElectricFieldSetup;  // We'll keep the same name for simplicity
     
    
    
    // Detector positions
    G4ThreeVector fDetector1Position;
    G4ThreeVector fDetector2Position;
    G4ThreeVector fDetector3Position;
    G4ThreeVector fDetector4Position;
    
};
#endif