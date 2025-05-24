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

    
    
  private:
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

    G4LogicalVolume* fChicaneMagnet1Volume;
    G4LogicalVolume* fChicaneMagnet2Volume;
    G4LogicalVolume* fChicaneMagnet3Volume;
    G4LogicalVolume* fChicaneMagnet4Volume;
  
  // Individual magnetic fields
    G4MagneticField* fChicaneMagField1;
    G4MagneticField* fChicaneMagField2;
    G4MagneticField* fChicaneMagField3;
    G4MagneticField* fChicaneMagField4;
  
  // Individual field managers
    G4FieldManager* fChicaneFieldManager1;
    G4FieldManager* fChicaneFieldManager2;
    G4FieldManager* fChicaneFieldManager3;
    G4FieldManager* fChicaneFieldManager4;
  
  // Methods to create individual magnets
    void CreateChicaneMagnet1(G4double start_z);
    void CreateChicaneMagnet2(G4double start_z);
    void CreateChicaneMagnet3(G4double start_z);
    void CreateChicaneMagnet4(G4double start_z);
  
  // Methods to setup individual fields
    void SetupChicaneMagnet1Field();
    void SetupChicaneMagnet2Field();
    void SetupChicaneMagnet3Field();
    void SetupChicaneMagnet4Field();
    
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