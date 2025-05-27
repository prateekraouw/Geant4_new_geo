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
  fWorldLogical(nullptr),
  fDetector1Volume(nullptr),
  fDetector2Volume(nullptr),
  fDetector3Volume(nullptr),
  fDetector4Volume(nullptr),
  fMagnet1Volume(nullptr),
  fMagnet2Volume(nullptr),
  fMagnet3Volume(nullptr),
  fMagnet4Volume(nullptr),
  fMagField1(nullptr),
  fMagField2(nullptr),
  fMagField3(nullptr),
  fMagField4(nullptr),
  fFieldManager1(nullptr),
  fFieldManager2(nullptr),
  fFieldManager3(nullptr),
  fFieldManager4(nullptr)
{
  // Initialize chicane parameters
  fMagnetLength = 20*cm;
  fMagnetWidth = 50*cm;
  fMagnetHeight = 50*cm;
  fMagnetSeparation = 5*cm;
  fFieldStrength = 0.5*tesla;
}

DetectorConstruction::~DetectorConstruction()
{
  // Clean up chicane magnetic fields
  delete fMagField1;
  delete fMagField2;
  delete fMagField3;
  delete fMagField4;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  G4NistManager* nist = G4NistManager::Instance();
  
  // World material: Air
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  // Tungsten material
  G4Material* graphite_mat = nist->FindOrBuildMaterial("G4_GRAPHITE");
  
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
  G4double detector1_position = 300*cm;  // 10 cm from tungsten
  G4double detector2_position = 543*cm; // m from the block
  G4double detector3_position = 923*cm;
  G4double detector4_position = 1035*cm;

  G4double helium_start = 545*cm;
  G4double helium_end = 555*cm;
  G4double helium_thickness = helium_end - helium_start; // 10 cm
  G4double helium_radius = 75*cm; // Make it wide enough to interact with particles

  // World volume - cylindrical
G4double world_radius = 0.5*world_size;  // Radius matching the box half-width
G4double world_length = 20*world_size;   // Length matching the box z-dimension

G4double rfcavity_radius = 30*cm;  // Same as detector radius for simplicity
G4double rfcavity_length = 20*cm;   // 20 cm long as requested
G4double rfcavity_position = 1250*cm;  

G4Tubs* solidWorld = 
  new G4Tubs("World",
             0,                // inner radius
             world_radius,     // outer radius
             0.5*world_length, // half-length in z
             0*deg,            // starting angle
             360*deg);         // segment angle
  
G4LogicalVolume* logicWorld = 
  new G4LogicalVolume(solidWorld, world_mat, "World");
fWorldLogical = logicWorld;  // Store it as member variable
  
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
  G4Box* solidGraphite = 
    new G4Box("Graphite", 0.5*tungsten_x, 0.5*tungsten_y, 0.5*tungsten_z);
  
  G4LogicalVolume* logicTungsten = 
    new G4LogicalVolume(solidGraphite, graphite_mat, "Graphite");

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
              360*deg);               // spanning angle
  
  // Detector 1 - create logical volume
  G4LogicalVolume* logicDetector1 = 
    new G4LogicalVolume(solidDetector1, scintillator_mat, "Detector1");
    
  // Place Detector 1 and store position
  G4ThreeVector detector1Pos = G4ThreeVector(0, 0, detector1_position);
  new G4PVPlacement(nullptr,                // no rotation
                    detector1Pos,           // position
                    logicDetector1,         // its logical volume
                    "Detector1",            // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    false);                 // checking overlaps
  fDetector1Position = detector1Pos;

  // Detector 2 - create solid and logical volume
  G4Tubs* solidDetector2 = 
    new G4Tubs("Detector2", 
              0*cm,                   // inner radius
              75*cm,        // outer radius
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);               // spanning angle

  G4LogicalVolume* logicDetector2 = 
    new G4LogicalVolume(solidDetector2, scintillator_mat, "Detector2");
    
  // Place Detector 2 and store position
  G4ThreeVector detector2Pos = G4ThreeVector(0, 0, detector2_position);
  new G4PVPlacement(nullptr,                // no rotation
                    detector2Pos,           // position
                    logicDetector2,         // its logical volume
                    "Detector2",            // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    false);                 // checking overlaps
  fDetector2Position = detector2Pos;

  // Detector 3 - create solid and logical volume
  G4Tubs* solidDetector3 = 
    new G4Tubs("Detector3", 
              0*cm,                   // inner radius
              10*cm,       // outer radius
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);               // spanning angle
              
  G4LogicalVolume* logicDetector3 = 
    new G4LogicalVolume(solidDetector3, scintillator_mat, "Detector3");
    
  // Place Detector 3 and store position
  G4ThreeVector detector3Pos = G4ThreeVector(0, 0, detector3_position);
  new G4PVPlacement(nullptr,                // no rotation
                    detector3Pos,           // position
                    logicDetector3,         // its logical volume
                    "Detector3",            // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    false);                  // checking overlaps
  fDetector3Position = detector3Pos;

  // Detector 4 - create solid and logical volume (note: smaller radius of 30cm)
  G4Tubs* solidDetector4 = 
    new G4Tubs("Detector4", 
              0*cm,                   // inner radius
              10*cm,                  // outer radius (smaller than others)
              0.5*detector_thickness, // half-length in z
              0*deg,                  // start angle
              360*deg);               // spanning angle
              
  G4LogicalVolume* logicDetector4 = 
    new G4LogicalVolume(solidDetector4, scintillator_mat, "Detector4");
    
  // Place Detector 4 and store position
  G4ThreeVector detector4Pos = G4ThreeVector(0, 0, detector4_position);
  new G4PVPlacement(nullptr,                // no rotation
                    detector4Pos,           // position
                    logicDetector4,         // its logical volume
                    "Detector4",            // its name
                    logicWorld,             // its mother volume
                    false,                  // no boolean operation
                    0,                      // copy number
                    false);                  // checking overlaps
  fDetector4Position = detector4Pos;

  // Visualize the ParticleGun as a block
  G4Box* gunBox = new G4Box("GunBox", 2*cm, 2*cm, 3*cm); // 2x2x2 cm block

  G4LogicalVolume* gunLog = new G4LogicalVolume(
    gunBox,             // its solid
    graphite_mat,       // material
    "GunLogical"        // name
  );

  G4RotationMatrix* gunRot = new G4RotationMatrix();
  gunRot->rotateX(-20.*deg);  // Align with beam direction

  // Position 5 cm before tungsten (as per ParticleGun)
  G4ThreeVector gunPos = G4ThreeVector(0, 5.*cm, -40.*cm);
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
new G4PVPlacement(nullptr, rfcavityPos, logicRFCavity, "RFCavity", logicWorld, false, 0, false);

G4double field_start = -30*cm;
G4double helium_z = 570*cm;  // Use the start of the helium cloud

// Calculate total length of the solenoid
G4double total_length = helium_z - field_start;  // Leave 10cm buffer on each side

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
G4double min_radius = 80*cm;
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


G4double field_start2 = 535*cm;
G4double z2 = 950*cm;  // Use the start of the helium cloud

// Calculate total length of the solenoid
G4double total_length2 = z2 - field_start2;  // Leave 10cm buffer on each side

// Number of segments to create the tapered shape
const G4int numZPlanes2 = 11;  // 7 segments + 1 (need n+1 points for n segments)

// Define the z-positions of each plane
G4double zPlane2[numZPlanes2];
for (G4int i = 0; i < numZPlanes2; i++) {
    zPlane2[i] = field_start2 + 10*cm + i * (total_length2 / 11);
}

// Define the inner and outer radii at each z-position
G4double rInner2[numZPlanes2];
G4double rOuter2[numZPlanes2];
G4double max_radius2 = 80*cm;
G4double min_radius2 = 10*cm;
G4double radius_step2 = (max_radius2 - min_radius2) / 11;  // For 7 segments

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
  
  // Detector 1 visualization - Blue
  G4VisAttributes* detector1_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Blue
  detector1_vis_att->SetVisibility(true);
  logicDetector1->SetVisAttributes(detector1_vis_att);

  // Detector 2 visualization - Blue
  G4VisAttributes* detector2_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Blue
  detector2_vis_att->SetVisibility(true);
  logicDetector2->SetVisAttributes(detector2_vis_att);

  // Detector 3 visualization - Red
  G4VisAttributes* detector3_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0)); // Red
  detector3_vis_att->SetVisibility(true);
  logicDetector3->SetVisAttributes(detector3_vis_att);

  // Detector 4 visualization - Red
  G4VisAttributes* detector4_vis_att = new G4VisAttributes(G4Colour::Blue()); // Blue
  detector4_vis_att->SetVisibility(true);
  logicDetector4->SetVisAttributes(detector4_vis_att);
  
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

  // Create chicane magnets
  CreateChicaneMagnets();

  return physWorld;
}

void DetectorConstruction::CreateChicaneMagnets()
{
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* magnet_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  // Calculate chicane position - place it after detector 2
  G4double chicane_start_z = 930*cm;
  
  // Create all four chicane magnets
  for (int i = 0; i < 4; i++) {
    G4double magnet_z = chicane_start_z + i * (fMagnetLength + fMagnetSeparation) + fMagnetLength/2.0;
    
    // Create magnet solid
    G4Box* solidMagnet = new G4Box("ChicaneMagnet" + std::to_string(i+1),
                                   0.5*fMagnetWidth,
                                   0.5*fMagnetHeight,
                                   0.5*fMagnetLength);
    
    // Create logical volume
    G4LogicalVolume* magnetVolume = new G4LogicalVolume(solidMagnet, magnet_mat, "ChicaneMagnet" + std::to_string(i+1));
    
    // Place magnet
    new G4PVPlacement(nullptr,
                      G4ThreeVector(0, 0, magnet_z),
                      magnetVolume,
                      "ChicaneMagnet" + std::to_string(i+1),
                      fWorldLogical,  // Use the member variable instead of logicWorld
                      false,
                      i,
                      false);
    
    // Store magnet volumes
    switch(i) {
      case 0: fMagnet1Volume = magnetVolume; break;
      case 1: fMagnet2Volume = magnetVolume; break;
      case 2: fMagnet3Volume = magnetVolume; break;
      case 3: fMagnet4Volume = magnetVolume; break;
    }
    
    // Set visualization attributes
    G4VisAttributes* magnet_vis_att;
    if (i == 0 || i == 3) {
      // Magnets 1 and 4: +X field (Red)
      magnet_vis_att = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.6));
    } else {
      // Magnets 2 and 3: -X field (Blue)
      magnet_vis_att = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.6));
    }
    magnet_vis_att->SetVisibility(true);
    magnet_vis_att->SetForceSolid(true);
    magnetVolume->SetVisAttributes(magnet_vis_att);
  }
  
  G4cout << "Chicane magnets created successfully" << G4endl;
}

void DetectorConstruction::SetupChicaneFields()
{
  if (!fMagnet1Volume || !fMagnet2Volume || !fMagnet3Volume || !fMagnet4Volume) {
    G4cerr << "Error: Chicane magnet volumes not created!" << G4endl;
    return;
  }
  
  // Setup each magnet field individually
  for (int i = 0; i < 4; i++) {
    G4ThreeVector fieldValue;
    G4LogicalVolume* magnetVolume;
    G4MagneticField** magField;
    G4FieldManager** fieldManager;
    
    // Determine field direction and get references
    switch(i) {
      case 0: // Magnet 1: +X field
        fieldValue = G4ThreeVector(+fFieldStrength, 0.0, 0.0);
        magnetVolume = fMagnet1Volume;
        magField = &fMagField1;
        fieldManager = &fFieldManager1;
        break;
      case 1: // Magnet 2: -X field
        fieldValue = G4ThreeVector(-fFieldStrength, 0.0, 0.0);
        magnetVolume = fMagnet2Volume;
        magField = &fMagField2;
        fieldManager = &fFieldManager2;
        break;
      case 2: // Magnet 3: -X field
        fieldValue = G4ThreeVector(-fFieldStrength, 0.0, 0.0);
        magnetVolume = fMagnet3Volume;
        magField = &fMagField3;
        fieldManager = &fFieldManager3;
        break;
      case 3: // Magnet 4: +X field
        fieldValue = G4ThreeVector(+fFieldStrength, 0.0, 0.0);
        magnetVolume = fMagnet4Volume;
        magField = &fMagField4;
        fieldManager = &fFieldManager4;
        break;
    }
    
    // Create magnetic field
    *magField = new UniformMagField(fieldValue);
    
    // Create field manager
    *fieldManager = new G4FieldManager(*magField);
    
    // Create equation of motion
    G4Mag_UsualEqRhs* equation = new G4Mag_UsualEqRhs(*magField);
    
    // Create stepper
    G4ClassicalRK4* stepper = new G4ClassicalRK4(equation);
    
    // Create chord finder
    G4double minStep = 0.005*mm;
    G4MagInt_Driver* driver = new G4MagInt_Driver(minStep, stepper, stepper->GetNumberOfVariables());
    G4ChordFinder* chordFinder = new G4ChordFinder(driver);
    
    // Set field manager parameters
    (*fieldManager)->SetChordFinder(chordFinder);
    (*fieldManager)->SetDetectorField(*magField);
    (*fieldManager)->SetDeltaOneStep(0.001*mm);
    (*fieldManager)->SetDeltaIntersection(0.001*mm);
    
    // Apply field manager to magnet
    magnetVolume->SetFieldManager(*fieldManager, true);
  }
  
  G4cout << "\n-----------------------------------------------------------" << G4endl;
  G4cout << " Chicane Magnetic Fields Setup Complete!" << G4endl;
  G4cout << " Field strength: " << fFieldStrength/tesla << " Tesla" << G4endl;
  G4cout << " Field pattern: +X, -X, -X, +X (standard chicane)" << G4endl;
  G4cout << "-----------------------------------------------------------\n" << G4endl;
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
  
    // First solenoid setup
    G4LogicalVolume* solenoidVolume1 = fMagFieldVolumes[0];
  
    // Create a 7 Tesla magnetic field in the Z direction for first solenoid
    G4ThreeVector fieldValue1 = G4ThreeVector(0.0, 0.0, 7.0*tesla);
    G4MagneticField* magField1 = new UniformMagField(fieldValue1);
    fMagFields.push_back(magField1);
  
    // Create a local field manager for first solenoid
    G4FieldManager* localFieldManager1 = new G4FieldManager(magField1);
    fLocalFieldManagers.push_back(localFieldManager1);
  
    // Create equation of motion for first field
    G4Mag_UsualEqRhs* equation1 = new G4Mag_UsualEqRhs(magField1);
  
    // Create stepper with higher precision for accurate field tracking
    G4ClassicalRK4* stepper1 = new G4ClassicalRK4(equation1);
  
    // Create the chord finder with appropriate step size
    G4double minStep1 = 0.005*CLHEP::mm;
    G4MagInt_Driver* driver1 = new G4MagInt_Driver(minStep1, stepper1, stepper1->GetNumberOfVariables());
    G4ChordFinder* chordFinder1 = new G4ChordFinder(driver1);
  
    // Set field manager parameters for first solenoid
    localFieldManager1->SetChordFinder(chordFinder1);
    localFieldManager1->SetDetectorField(magField1);
    localFieldManager1->SetDeltaOneStep(0.001*CLHEP::mm);
    localFieldManager1->SetDeltaIntersection(0.001*CLHEP::mm);
  
    // Apply field manager to the first solenoid logical volume
    solenoidVolume1->SetFieldManager(localFieldManager1, true);
  
    G4cout << "\n-----------------------------------------------------------" << G4endl;
    G4cout << " First Magnetic Field of 7 Tesla applied to first solenoid" << G4endl;
    G4cout << " Field value: (0, 0, 7) Tesla" << G4endl;
    G4cout << " First Solenoid tapers from radius 80 cm to 20 cm" << G4endl;
    G4cout << "-----------------------------------------------------------\n" << G4endl;

    // Second solenoid setup
    G4LogicalVolume* solenoidVolume2 = fMagFieldVolumes[1];

    // Create a separate 7 Tesla magnetic field for the second solenoid
    G4ThreeVector fieldValue2 = G4ThreeVector(0.0, 0.0, 7.0*tesla);
    G4MagneticField* magField2 = new UniformMagField(fieldValue2);
    fMagFields.push_back(magField2);

    // Create a separate local field manager for second solenoid
    G4FieldManager* localFieldManager2 = new G4FieldManager(magField2);
    fLocalFieldManagers.push_back(localFieldManager2);

    // Create separate equation of motion for second field
    G4Mag_UsualEqRhs* equation2 = new G4Mag_UsualEqRhs(magField2);

    // Create separate stepper for second field
    G4ClassicalRK4* stepper2 = new G4ClassicalRK4(equation2);

    // Create separate chord finder for second field
    G4double minStep2 = 0.005*CLHEP::mm;
    G4MagInt_Driver* driver2 = new G4MagInt_Driver(minStep2, stepper2, stepper2->GetNumberOfVariables());
    G4ChordFinder* chordFinder2 = new G4ChordFinder(driver2);

    // Set parameters on second field manager
    localFieldManager2->SetChordFinder(chordFinder2);
    localFieldManager2->SetDetectorField(magField2);
    localFieldManager2->SetDeltaOneStep(0.001*CLHEP::mm);
    localFieldManager2->SetDeltaIntersection(0.001*CLHEP::mm);

    // Apply this field manager to the second solenoid volume
    solenoidVolume2->SetFieldManager(localFieldManager2, true);

    G4cout << "\n-----------------------------------------------------------" << G4endl;
    G4cout << " Second Magnetic Field of 7 Tesla applied to second solenoid" << G4endl;
    G4cout << " Field value: (0, 0, 7) Tesla" << G4endl;
    G4cout << " Second Solenoid tapers from radius 40 cm to 5 cm" << G4endl;
    G4cout << "-----------------------------------------------------------\n" << G4endl;

    // Setup chicane magnetic fields
    SetupChicaneFields();
}