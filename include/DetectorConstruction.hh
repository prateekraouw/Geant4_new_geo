#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4MagneticField.hh"  
#include "G4FieldManager.hh" 

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

    void CreateChicaneMagnets();
    void SetupChicaneFields();

    
    
  private:
    G4LogicalVolume* fWorldLogical;
    G4LogicalVolume* fScoringVolume;
    G4LogicalVolume* fDetector1Volume;
    G4LogicalVolume* fDetector2Volume;
    G4LogicalVolume* fDetector3Volume;
    G4LogicalVolume* fDetector4Volume;
    G4LogicalVolume* fRFCavityVolume;  // Add this
    
    //ElectricFieldSetup* fElectricFieldSetup;  // We'll keep the same name for simplicity

    // Vector to store magnetic field volumes
    std::vector<G4LogicalVolume*> fMagFieldVolumes;
     
    // Field-related member variables
    std::vector<G4MagneticField*> fMagFields;
    std::vector<G4FieldManager*> fLocalFieldManagers;
    
    // Detector positions
    G4ThreeVector fDetector1Position;
    G4ThreeVector fDetector2Position;
    G4ThreeVector fDetector3Position;
    G4ThreeVector fDetector4Position;

    // Chicane magnet volumes
    G4LogicalVolume* fMagnet1Volume;
    G4LogicalVolume* fMagnet2Volume;
    G4LogicalVolume* fMagnet3Volume;
    G4LogicalVolume* fMagnet4Volume;
    
    // Chicane magnetic fields
    G4MagneticField* fMagField1;
    G4MagneticField* fMagField2;
    G4MagneticField* fMagField3;
    G4MagneticField* fMagField4;
    
    // Chicane field managers
    G4FieldManager* fFieldManager1;
    G4FieldManager* fFieldManager2;
    G4FieldManager* fFieldManager3;
    G4FieldManager* fFieldManager4;
    
    // Chicane parameters
    G4double fMagnetLength;
    G4double fMagnetWidth;
    G4double fMagnetHeight;
    G4double fMagnetSeparation;
    G4double fFieldStrength;
    
};

// Define uniform magnetic field class
class UniformMagField : public G4MagneticField
{
public:
  UniformMagField(const G4ThreeVector& fieldVector) : fFieldValue(fieldVector) {}
  virtual ~UniformMagField() {}
  
  virtual void GetFieldValue(const G4double[4], G4double* field) const override {
    field[0] = fFieldValue.x();
    field[1] = fFieldValue.y();
    field[2] = fFieldValue.z();
    field[3] = 0.0;
    field[4] = 0.0;
    field[5] = 0.0;
  }
  
private:
  G4ThreeVector fFieldValue;
};
#endif