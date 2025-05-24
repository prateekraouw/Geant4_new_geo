#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "ElectricFieldSetup.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4MagneticField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4Polycone.hh"

DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(nullptr),
  fDetector1Volume(nullptr),
  fDetector2Volume(nullptr),
  fDetector3Volume(nullptr),
  fDetector4Volume(nullptr)
{
}

DetectorConstruction::~DetectorConstruction()
{
  delete fChicaneMagField1;
  delete fChicaneMagField2;
  delete fChicaneMagField3;
  delete fChicaneMagField4;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  G4NistManager* nist = G4NistManager::Instance();
  
  // World material: Air
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  // Tungsten material
  G4Material* tungsten_mat = nist->FindOrBuildMaterial("G4_W");
  
  // Scintillator material for detectors (plastic scintillator)
  G4Material* scintillator_mat = nist->FindOrBuildMaterial("G4_Ar");

   // Helium material
  G4Material* helium_mat = nist->FindOrBuildMaterial("G4_He");


  // Vaccum for RF cavity
  G4Material* rfcavity_mat = nist->FindOrBuildMaterial("G4_Galactic"); // Use vacuum for RF cavity

  // World volume parameters
  G4double world_size = 10000*cm;

  // Tungsten block parameters - 10×10×30 cm
  G4double tungsten_x = 5*cm;
  G4double tungsten_y = 5*cm; 
  G4double tungsten_z = 75*cm;
  
  // Detector parameters - circular discs with 30 cm diameter and 1 cm thickness
  G4double detector_radius = 75*cm;  // 30 cm diameter
  G4double detector_thickness = 0.1*cm;
  
  // Detector positions (distance from the end of the tungsten block)
  G4double detector1_position = tungsten_z/2 + 60*cm;  // 10 cm from tungsten
  G4double detector2_position = tungsten_z/2 + 700*cm; // m from the block
  G4double detector3_position = 420*cm;
  G4double detector4_position = 630*cm;

  G4double helium_start = 100*cm;
  G4double helium_end = 110*cm;
  G4double helium_thickness = helium_end - helium_start; // 10 cm
  G4double helium_radius = 80*cm; // Make it wide enough to interact with particles

  // World volume - cylindrical
G4double world_radius = 0.5*world_size;  // Radius matching the box half-width
G4double world_length = 20*world_size;   // Length matching the box z-dimension

G4double rfcavity_radius = 100*cm;  // Same as detector radius for simplicity
G4double rfcavity_length = 20*cm;   // 20 cm long as requested
G4double rfcavity_position = 150*cm;  





G4Tubs* solidWorld = 
  new G4Tubs("World",
             0,                // inner radius
             world_radius,     // outer radius
             0.5*world_length, // half-length in z
             0*deg,            // starting angle
             360*deg);         // segment angle
  
G4LogicalVolume* logicWorld = 
  new G4LogicalVolume(solidWorld, world_mat, "World");
  
G4VPhysicalVolume* physWorld = 
  new G4PVPlacement(nullptr,               // no rotation
                    G4ThreeVector(),       // at (0,0,0)
                    logicWorld,            // its logical volume
                    "World",               // its name
                    nullptr,               // its mother volume
                    false,                 // no boolean operation
                    0,                     // copy number
                    true);                 // checking overlaps

  // Tungsten block
  G4Box* solidTungsten = 
    new G4Box("Tungsten", 0.5*tungsten_x, 0.5*tungsten_y, 0.5*tungsten_z);
  
  G4LogicalVolume* logicTungsten = 
    new G4LogicalVolume(solidTungsten, tungsten_mat, "Tungsten");

  G4RotationMatrix* rotation = new G4RotationMatrix();
  rotation->rotateX(-20*deg);
  
  new G4PVPlacement(rotation,                // no rotation
                    G4ThreeVector(0, 0, 0), // at (0,0,0)
                    logicTungsten,          // its logical volume
                    "Tungsten",             // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    true);                  // checking overlaps


G4Tubs* solidHelium = 
    new G4Tubs("HeliumCloud",
              0*cm,                   // inner radius
              helium_radius,          // outer radius
              0.5*helium_thickness,   // half-length in z
              0*deg,                  // start angle
              360*deg);               // spanning angle
              
  G4LogicalVolume* logicHelium = 
    new G4LogicalVolume(solidHelium, helium_mat, "HeliumCloud");
    
  // Place helium cloud at the center position between start and end
  G4double helium_center = helium_start + 0.5*helium_thickness;
  
  new G4PVPlacement(nullptr,               // no rotation
                    G4ThreeVector(0, 0, helium_center), // at center of cloud
                    logicHelium,           // its logical volume
                    "HeliumCloud",         // its name
                    logicWorld,            // its mother volume
                    false,                 // no boolean operation
                    0,                     // copy number
                    true);                 // checking overlaps


  // Create circular detectors (discs)
  G4Tubs* solidDetector1 = 
    new G4Tubs("Detector1", 
              0*cm,                   // inner radius
              detector_radius,        // outer radius
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);              // spanning angle             // spanning angle
  
  // Detector 1 (2 cm from tungsten)
  G4LogicalVolume* logicDetector1 = 
    new G4LogicalVolume(solidDetector1, scintillator_mat, "Detector1");
  
  new G4PVPlacement(nullptr,                // no rotation
                    G4ThreeVector(0, 0, detector1_position), // position
                    logicDetector1,         // its logical volume
                    "Detector1",            // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    false);                  // checking overlaps
  
  G4Tubs* solidDetector2 = 
    new G4Tubs("Detector2", 
              0*cm,                   // inner radius
              detector_radius,        // outer radius
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);              // spanning angle             // spanning angle
  G4ThreeVector detector1Pos = G4ThreeVector(0, 0, detector1_position);
  new G4PVPlacement(nullptr, detector1Pos, logicDetector1, "Detector1", logicWorld, false, 0, false);
  fDetector1Position = detector1Pos;


  // Detector 2 (10 m from tungsten)
  G4LogicalVolume* logicDetector2 = 
    new G4LogicalVolume(solidDetector2, scintillator_mat, "Detector2");
  
  new G4PVPlacement(nullptr,                // no rotation
                    G4ThreeVector(0, 0, detector2_position), // position
                    logicDetector2,         // its logical volume
                    "Detector2",            // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    false);                  // checking overlaps
  
  

  G4ThreeVector detector2Pos = G4ThreeVector(0, 0, detector2_position);
  new G4PVPlacement(nullptr, detector2Pos, logicDetector2, "Detector2", logicWorld, false, 0, false);
  fDetector2Position = detector2Pos;


    // Detector 3 (at 375 cm mark)
  G4Tubs* solidDetector3 = 
    new G4Tubs("Detector3", 
              0*cm,                   // inner radius
              detector_radius,        // outer radius
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);               // spanning angle
              
  G4LogicalVolume* logicDetector3 = 
    new G4LogicalVolume(solidDetector3, scintillator_mat, "Detector3");
  
  G4ThreeVector detector3Pos = G4ThreeVector(0, 0, detector3_position);
  new G4PVPlacement(nullptr, detector3Pos, logicDetector3, "Detector3", logicWorld, false, 0, true);
  fDetector3Position = detector3Pos;


  // Detector 3 (at 500 cm mark)
  G4Tubs* solidDetector4 = 
    new G4Tubs("Detector4", 
              0*cm,                   // inner radius
              30*cm,        // outer radius
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);               // spanning angle
              
  G4LogicalVolume* logicDetector4 = 
    new G4LogicalVolume(solidDetector4, scintillator_mat, "Detector4");
  
  G4ThreeVector detector4Pos = G4ThreeVector(0, 0, detector4_position);
  new G4PVPlacement(nullptr, detector4Pos, logicDetector4, "Detector4", logicWorld, false, 0, true);
  fDetector4Position = detector4Pos;



  // Visualize the ParticleGun as a block
  G4Box* gunBox = new G4Box("GunBox", 2*cm, 2*cm, 3*cm); // 2x2x2 cm block

  G4LogicalVolume* gunLog = new G4LogicalVolume(
    gunBox,             // its solid
    tungsten_mat,       // material
    "GunLogical"        // name
  );

  G4RotationMatrix* gunRot = new G4RotationMatrix();
  gunRot->rotateX(-20.*deg);  // Align with beam direction

  // Position 5 cm before tungsten (as per ParticleGun)
  G4ThreeVector gunPos = G4ThreeVector(0, 12.*cm, -40.*cm);
  new G4PVPlacement(
    gunRot,             // rotation
    gunPos,             // position
    gunLog,             // logical volume
    "GunBlock",         // name
    logicWorld,         // mother volume
    false,              // no boolean operation
    0,                  // copy number 
    true                // overlaps
  );

  G4Tubs* solidRFCavity = 
  new G4Tubs("RFCavity",
            0*cm,                   // inner radius
            rfcavity_radius,        // outer radius
            0.5*rfcavity_length,    // half-length in z
            0*deg,                  // start angle
            360*deg);               // spanning angle

G4LogicalVolume* logicRFCavity = 
  new G4LogicalVolume(solidRFCavity, rfcavity_mat, "RFCavity");

G4ThreeVector rfcavityPos = G4ThreeVector(0, 0, rfcavity_position);
new G4PVPlacement(nullptr, rfcavityPos, logicRFCavity, "RFCavity", logicWorld, false, 0, true);



G4double field_start = -40*cm;
G4double helium_z = 340*cm;  // Use the start of the helium cloud

// Calculate total length of the solenoid
G4double total_length = helium_z - field_start - 10*cm;  // Leave 10cm buffer on each side

// Number of segments to create the tapered shape
const G4int numZPlanes = 16;  // 7 segments + 1 (need n+1 points for n segments)

// Define the z-positions of each plane
G4double zPlane[numZPlanes];
for (G4int i = 0; i < numZPlanes; i++) {
    zPlane[i] = field_start + 10*cm + i * (total_length / 16);
}

// Define the inner and outer radii at each z-position
G4double rInner[numZPlanes];
G4double rOuter[numZPlanes];
G4double max_radius = 80*cm;
G4double min_radius = 20*cm;
G4double radius_step = (max_radius - min_radius) / 16;  // For 7 segments

for (G4int i = 0; i < numZPlanes; i++) {
    rInner[i] = 0*cm;  // No inner hole
    rOuter[i] = max_radius - i * radius_step;  // Decreasing outer radius
}

// Create the polycone solid
G4Polycone* solidMagField = 
    new G4Polycone("MagneticFieldSolenoid",
                   0*deg,                 // starting phi angle
                   360*deg,               // delta phi angle
                   numZPlanes,            // number of z planes
                   zPlane,                // z plane positions
                   rInner,                // inner radii
                   rOuter);               // outer radii

// Use air as the material for the field region
G4Material* field_mat = nist->FindOrBuildMaterial("G4_AIR");

// Create logical volume
G4LogicalVolume* logicMagField = 
    new G4LogicalVolume(solidMagField, field_mat, "MagneticFieldSolenoid");

// Place the volume
new G4PVPlacement(nullptr,               // no rotation
                  G4ThreeVector(0, 0, 0), // position at origin (polycone already positioned via z-planes)
                  logicMagField,          // its logical volume
                  "MagneticFieldSolenoid", // its name
                  logicWorld,             // its mother volume
                  false,                  // no boolean operation
                  0,                      // copy number
                  true);                  // checking overlaps


G4double field_start2 = 425*cm;
G4double z2 = 625*cm;  // Use the start of the helium cloud

// Calculate total length of the solenoid
G4double total_length2 = z2 - field_start2 - 10*cm;  // Leave 10cm buffer on each side

// Number of segments to create the tapered shape
const G4int numZPlanes2 = 8;  // 7 segments + 1 (need n+1 points for n segments)

// Define the z-positions of each plane
G4double zPlane2[numZPlanes2];
for (G4int i = 0; i < numZPlanes2; i++) {
    zPlane2[i] = field_start2 + 10*cm + i * (total_length2 / 8);
}

// Define the inner and outer radii at each z-position
G4double rInner2[numZPlanes2];
G4double rOuter2[numZPlanes2];
G4double max_radius2 = 40*cm;
G4double min_radius2 = 5*cm;
G4double radius_step2 = (max_radius2 - min_radius2) / 8;  // For 7 segments

for (G4int i = 0; i < numZPlanes2; i++) {
    rInner2[i] = 0*cm;  // No inner hole
    rOuter2[i] = max_radius2 - i * radius_step2;  // Decreasing outer radius
}

// Create the polycone solid
G4Polycone* solidMagField2 = 
    new G4Polycone("MagneticFieldSolenoid2",
                   0*deg,                 // starting phi angle
                   360*deg,               // delta phi angle
                   numZPlanes2,            // number of z planes
                   zPlane2,                // z plane positions
                   rInner2,                // inner radii
                   rOuter2);               // outer radii

// Create logical volume
G4LogicalVolume* logicMagField2 = 
    new G4LogicalVolume(solidMagField2, field_mat, "MagneticFieldSolenoid2");

// Place the volume
new G4PVPlacement(nullptr,               // no rotation
                  G4ThreeVector(0, 0, 0), // position at origin (polycone already positioned via z-planes)
                  logicMagField2,          // its logical volume
                  "MagneticFieldSolenoid2", // its name
                  logicWorld,             // its mother volume
                  false,                  // no boolean operation
                  0,                      // copy number
                  true);                  // checking overlaps


  // Visual attributes
  G4VisAttributes* tungsten_vis_att = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5)); // Grey
  logicTungsten->SetVisAttributes(tungsten_vis_att);
  
  G4VisAttributes* detector1_vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue
  detector1_vis_att->SetVisibility(true);
  logicDetector1->SetVisAttributes(detector1_vis_att);

  G4VisAttributes* detector2_vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // Blue
  detector1_vis_att->SetVisibility(true);
  logicDetector1->SetVisAttributes(detector2_vis_att);


  G4VisAttributes* detector3_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red for detector 3
  detector3_vis_att->SetVisibility(true);
  logicDetector3->SetVisAttributes(detector3_vis_att);

   G4VisAttributes* detector4_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red for detector 3
  detector3_vis_att->SetVisibility(true);
  logicDetector3->SetVisAttributes(detector4_vis_att);
  
  // Helium cloud visualization - light blue and semi-transparent
  G4VisAttributes* helium_vis_att = new G4VisAttributes(G4Colour(0.6, 0.8, 1.0, 0.3)); // Light blue, semi-transparent
  helium_vis_att->SetVisibility(true);
  helium_vis_att->SetForceSolid(true);
  logicHelium->SetVisAttributes(helium_vis_att);
  
  G4VisAttributes* gunVis = new G4VisAttributes(G4Colour::Red());
  gunVis->SetVisibility(true);
  gunVis->SetForceSolid(true);
  gunLog->SetVisAttributes(gunVis);

  // RF Cavity visualization attributes - make it distinct
  G4VisAttributes* rfcavity_vis_att = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0, 0.7)); // Orange, semi-transparent
  rfcavity_vis_att->SetVisibility(true);
  rfcavity_vis_att->SetForceSolid(true);
  logicRFCavity->SetVisAttributes(rfcavity_vis_att);

  // Create visualization attributes - gradient blue
  G4VisAttributes* field_vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 0.8, 0.4));  // Blue, semi-transparent
  field_vis_att->SetVisibility(true);
  field_vis_att->SetForceSolid(true);
  logicMagField->SetVisAttributes(field_vis_att);

  // Create visualization attributes - gradient blue
  G4VisAttributes* field_vis_att2 = new G4VisAttributes(G4Colour(0.0, 0.0, 0.8, 0.4));  // Blue, semi-transparent
  field_vis_att2->SetVisibility(true);
  field_vis_att2->SetForceSolid(true);
  logicMagField2->SetVisAttributes(field_vis_att2);


  // Make the world volume transparent
  G4VisAttributes* world_vis_att = new G4VisAttributes(G4Colour(1.5, 1.5, 1.5, 1.0)); // Transparent
  world_vis_att->SetVisibility(true);
  world_vis_att->SetForceWireframe(true);
  logicWorld->SetVisAttributes(world_vis_att);

  // In DetectorConstruction::Construct()
  // Set scoring volumes
  fScoringVolume = logicTungsten;
  fDetector1Volume = logicDetector1;
  fDetector2Volume = logicDetector2;
  fDetector3Volume = logicDetector3;
  fDetector4Volume = logicDetector4;
  fRFCavityVolume = logicRFCavity;
  fMagFieldVolumes.push_back(logicMagField);
  fMagFieldVolumes.push_back(logicMagField2);

  return physWorld;
}



void DetectorConstruction::ConstructSDandField()
{
    if (fMagFieldVolumes.empty()) {
    G4cerr << "Error: No magnetic field volumes created!" << G4endl;
    return;
  }

  if (fMagFieldVolumes.size() < 2) {
        G4cerr << "Error: Expected at least two magnetic field volumes!" << G4endl;
        return;
    }
  
  // We have just one volume (the polycone solenoid) in the vector
  G4LogicalVolume* solenoidVolume = fMagFieldVolumes[0];
  
  // Create a 7 Tesla magnetic field in the Z direction
  G4ThreeVector fieldValue = G4ThreeVector(0.0, 0.0, 7.0*tesla);
  G4MagneticField* magField = new UniformMagField(fieldValue);
  fMagFields.push_back(magField);
  
  // Create a local field manager
  G4FieldManager* localFieldManager = new G4FieldManager(magField);
  fLocalFieldManagers.push_back(localFieldManager);
  
  // Create equation of motion for this field
  G4Mag_UsualEqRhs* equation = new G4Mag_UsualEqRhs(magField);
  
  // Create stepper with higher precision for accurate field tracking
  G4ClassicalRK4* stepper = new G4ClassicalRK4(equation);
  
  // Create the chord finder with appropriate step size
  G4double minStep = 0.005*CLHEP::mm;
  G4MagInt_Driver* driver = new G4MagInt_Driver(minStep, stepper, stepper->GetNumberOfVariables());
  G4ChordFinder* chordFinder = new G4ChordFinder(driver);
  
  // Set field manager parameters
  localFieldManager->SetChordFinder(chordFinder);
  localFieldManager->SetDetectorField(magField);
  localFieldManager->SetDeltaOneStep(0.001*CLHEP::mm);
  localFieldManager->SetDeltaIntersection(0.001*CLHEP::mm);
  
  // Apply field manager to the solenoid logical volume
  solenoidVolume->SetFieldManager(localFieldManager, true);
  
  G4cout << "\n-----------------------------------------------------------" << G4endl;
  G4cout << " Magnetic Field of 7 Tesla applied to continuous solenoid" << G4endl;
  G4cout << " Field value: (0, 0, 7) Tesla" << G4endl;
  G4cout << " Solenoid tapers from radius 100 cm to 50 cm" << G4endl;
  G4cout << "-----------------------------------------------------------\n" << G4endl;

  // Get the second solenoid logical volume
  G4LogicalVolume* solenoidVolume2 = fMagFieldVolumes[1];

  // Create a local field manager
  G4FieldManager* localFieldManager2 = new G4FieldManager(magField);
  fLocalFieldManagers.push_back(localFieldManager);

  // Equation of motion
  G4Mag_UsualEqRhs* equation2 = new G4Mag_UsualEqRhs(magField);

  // Stepper
  G4ClassicalRK4* stepper2 = new G4ClassicalRK4(equation2);

  // Chord finder
  G4double minStep2 = 0.005*CLHEP::mm;
  G4MagInt_Driver* driver2 = new G4MagInt_Driver(minStep2, stepper2, stepper2->GetNumberOfVariables());
  G4ChordFinder* chordFinder2 = new G4ChordFinder(driver2);

  // Set parameters on field manager
  localFieldManager2->SetChordFinder(chordFinder2);
  localFieldManager2->SetDetectorField(magField);
  localFieldManager2->SetDeltaOneStep(0.001*CLHEP::mm);
  localFieldManager2->SetDeltaIntersection(0.001*CLHEP::mm);

  // Apply this field manager to the second solenoid volume
  solenoidVolume2->SetFieldManager(localFieldManager2, true);

  G4cout << "\n-----------------------------------------------------------" << G4endl;
  G4cout << " Second Magnetic Field of 7 Tesla applied to second solenoid" << G4endl;
  G4cout << " Field value: (0, 0, 7) Tesla" << G4endl;
  G4cout << " Second Solenoid tapers from radius 30 cm to 5 cm" << G4endl;
  G4cout << "-----------------------------------------------------------\n" << G4endl;


}
