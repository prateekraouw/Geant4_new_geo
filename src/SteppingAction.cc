#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
// NEW: Additional includes for field logging
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4TransportationManager.hh"
#include <iomanip>

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(nullptr),
  fDetector1Volume(nullptr),
  fDetector2Volume(nullptr),
  fDetector3Volume(nullptr),
  fDetector4Volume(nullptr),
  fRFCavityVolume(nullptr),
  // NEW: Initialize field logging members
  fLogField(false),  // Enable by default
  fLogFileName("magnetic_field_continuous.csv"),
  fFieldStepCounter(0),
  fLogInterval(1),  // Log every step for continuous tracking
  fFileInitialized(false)
{}

SteppingAction::~SteppingAction()
{
  // NEW: Close field logging file
  if (fLogFile.is_open()) {
    fLogFile.close();
    G4cout << "Magnetic field logging file closed: " << fLogFileName << G4endl;
  }

  // Your existing destructor code unchanged
  // Print out the particle count at the end, filtering for muons and pions only
  G4cout << "=== Muons and Pions Generated ===" << G4endl;
  for (auto const& pair : fParticleCounter) {
    // Only include muons and pions in the summary
    if (pair.first == "mu+" || pair.first == "mu-" || 
        pair.first == "pi+" || pair.first == "pi-") {
      G4cout << pair.first << ": " << pair.second << G4endl;
    }
  }
  G4cout << "===================================" << G4endl;
  
  // Add a section for particles detected at detector 1
  G4cout << "\n=== Muons and Pions Detected at Detector 1 ===" << G4endl;
  for (auto const& pair : fDetector1Particles) {
    G4cout << pair.first << ": " << pair.second << G4endl;
  }
  G4cout << "==========================================" << G4endl;
  
  // Add a section for particles detected at detector 2
  G4cout << "\n=== Muons and Pions Detected at Detector 2 ===" << G4endl;
  for (auto const& pair : fDetector2Particles) {
    G4cout << pair.first << ": " << pair.second << G4endl;
  }
  G4cout << "==========================================" << G4endl;
  
  // Add a section for particles detected at detector 3
  G4cout << "\n=== Muons and Pions Detected at Detector 3 ===" << G4endl;
  for (auto const& pair : fDetector3Particles) {
    G4cout << pair.first << ": " << pair.second << G4endl;
  }
  G4cout << "==========================================" << G4endl;
  
  // Add a section for particles through RF cavity
  G4cout << "\n=== Muons and Pions Through RF Cavity ===" << G4endl;
  for (auto const& pair : fRFCavityParticles) {
    G4cout << pair.first << ": " << pair.second << G4endl;
  }
  G4cout << "==========================================" << G4endl;
  
  // Add a section for particles detected at detector 4
  G4cout << "\n=== Muons and Pions Detected at Detector 4 ===" << G4endl;
  for (auto const& pair : fDetector4Particles) {
    G4cout << pair.first << ": " << pair.second << G4endl;
  }
  G4cout << "==========================================" << G4endl;
}

// NEW: Initialize field logging file
void SteppingAction::InitializeFieldLogFile()
{
  if (fFileInitialized) return;
  
  fLogFile.open(fLogFileName, std::ios::out);
  if (fLogFile.is_open()) {
    // Write simplified CSV header - only position and field
    fLogFile << "X_mm,Y_mm,Z_mm,Bx_T,By_T,Bz_T,B_magnitude_T\n";
    fFileInitialized = true;
    G4cout << "Magnetic field logging initialized: " << fLogFileName << G4endl;
  } else {
    G4cout << "ERROR: Cannot open field log file: " << fLogFileName << G4endl;
    fLogField = false;
  }
}

// NEW: Log magnetic field at position - simplified version
void SteppingAction::LogFieldAtPosition(const G4ThreeVector& position, G4int trackID, G4int stepNumber, G4Track* track)
{
  if (!fLogFile.is_open()) return;
  
  // Get magnetic field at current position
  G4ThreeVector bField = GetMagneticFieldAtPosition(position, track);
  
  // Calculate field magnitude
  G4double bMagnitude = bField.mag();
  
  // Write only position and field to CSV file
  fLogFile << std::fixed << std::setprecision(4)
           << position.x()/mm << ","
           << position.y()/mm << ","
           << position.z()/mm << ","
           << std::setprecision(8)  // High precision for field values
           << bField.x()/tesla << ","
           << bField.y()/tesla << ","
           << bField.z()/tesla << ","
           << bMagnitude/tesla << "\n";
  
  // Flush frequently for continuous monitoring
  if (stepNumber % 100 == 0) {
    fLogFile.flush();
  }
}

// NEW: Get magnetic field at position
G4ThreeVector SteppingAction::GetMagneticFieldAtPosition(const G4ThreeVector& position, G4Track* track)
{
  // Get the field manager for current volume
  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()
                                  ->GetFieldManager();
  
  if (!fieldManager) {
    // Try to get local field manager
    if (track && track->GetVolume()) {
      G4LogicalVolume* logVol = track->GetVolume()->GetLogicalVolume();
      fieldManager = logVol->GetFieldManager();
    }
  }
  
  if (!fieldManager) {
    return G4ThreeVector(0, 0, 0);
  }
  
  const G4Field* field = fieldManager->GetDetectorField();
  if (!field) {
    return G4ThreeVector(0, 0, 0);
  }
  
  // Cast to magnetic field
  const G4MagneticField* magField = dynamic_cast<const G4MagneticField*>(field);
  if (!magField) {
    return G4ThreeVector(0, 0, 0);
  }
  
  // Get field value at position
  G4double point[4] = {position.x(), position.y(), position.z(), 0.0};
  G4double bField[3] = {0.0, 0.0, 0.0};
  
  magField->GetFieldValue(point, bField);
  
  return G4ThreeVector(bField[0], bField[1], bField[2]);
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // NEW: Initialize field logging on first call
  if (fLogField && !fFileInitialized) {
    InitializeFieldLogFile();
  }
  
  // NEW: Continuous magnetic field logging
  if (fLogField && fFileInitialized) {
    fFieldStepCounter++;
    
    // Log every step for continuous tracking
    if (fFieldStepCounter % fLogInterval == 0) {
      G4Track* track = step->GetTrack();
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
      LogFieldAtPosition(position, track->GetTrackID(), fFieldStepCounter, track);
    }
  }

  // ========================================================================
  // ALL YOUR EXISTING CODE BELOW REMAINS COMPLETELY UNCHANGED
  // ========================================================================

  // Initialize volumes if not already done
  if (!fScoringVolume || !fDetector1Volume || !fDetector2Volume || !fDetector3Volume || !fDetector4Volume || !fRFCavityVolume) { 
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();
    fDetector1Volume = detectorConstruction->GetDetector1Volume();
    fDetector2Volume = detectorConstruction->GetDetector2Volume();
    fDetector3Volume = detectorConstruction->GetDetector3Volume();
    fDetector4Volume = detectorConstruction->GetDetector4Volume();
    fRFCavityVolume = detectorConstruction->GetRFCavityVolume();
    
    G4cout << "Detector 1 position: " << detectorConstruction->GetDetector1Position()/cm << " cm" << G4endl;
    G4cout << "Detector 2 position: " << detectorConstruction->GetDetector2Position()/cm << " cm" << G4endl;
    G4cout << "Detector 3 position: " << detectorConstruction->GetDetector3Position()/cm << " cm" << G4endl;
    G4cout << "Detector 4 position: " << detectorConstruction->GetDetector4Position()/cm << " cm" << G4endl;
  }
    
  // Get the RunAction - using const_cast to handle the constness issue
  const G4UserRunAction* constRunAction = G4RunManager::GetRunManager()->GetUserRunAction();
  RunAction* runAction = const_cast<RunAction*>(dynamic_cast<const RunAction*>(constRunAction));

  // Get current track
  G4Track* track = step->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String particleName = particle->GetParticleName();
  G4double energy = track->GetKineticEnergy();

  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    
    // Check if a magnetic field is assigned
    if (fieldMgr) {
        const G4Field* field = fieldMgr->GetDetectorField();  // Get the field as const G4Field
        if (field) {
            // Cast the field to G4MagneticField
            const G4MagneticField* magField = dynamic_cast<const G4MagneticField*>(field);
            if (magField) {
                // Get the magnetic field at the current point (step position)
                G4double point[4] = {track->GetPosition().x(), track->GetPosition().y(), track->GetPosition().z(), 0.0};
                G4double fieldValue[3];  // Array to store the field values (fx, fy, fz)

                magField->GetFieldValue(point, fieldValue);  // Get the field values in the array

                // Create a stringstream to format the output as CSV
                std::stringstream ss;

                // Get the position and field values
                G4double x = track->GetPosition().x();
                G4double y = track->GetPosition().y();
                G4double z = track->GetPosition().z();
                G4double fx = fieldValue[0];  // x component of the field (always 0)
                G4double fy = fieldValue[1];  // y component of the field (always 0)
                G4double fz = fieldValue[2];  // z component of the field

                // Write the values to the stringstream in CSV format
                ss << x << "," << y << "," << z << "," << fx << "," << fy << "," << fz << "\n";

                // Open the CSV file in append mode
                std::ofstream outFile("magnetic_field.csv", std::ios_base::app);  // Open in append mode to add to the file

                // Check if the file was successfully opened
                if (!outFile.is_open()) {
                    G4cerr << "Error: Could not open CSV file for writing!" << G4endl;
                } else {
                    // Write the formatted values to the file
                    outFile << ss.str();
                    outFile.close();  // Close the file after writing
                    G4cout << "Data written to magnetic_field.csv successfully." << G4endl;
                }
            }
        }
    }
  
  // Check for pion decay specifically
  G4String processName = "Unknown";
  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  if (process) processName = process->GetProcessName();
  
  // If this is a decay process
  if (processName == "Decay") {
    // If the current particle is a pion
    if (particleName == "pi+" || particleName == "pi-") {
      // Get secondaries created in this step
      const std::vector<const G4Track*>* secondaries = step->GetSecondaryInCurrentStep();
      
      if (secondaries && secondaries->size() > 0) {
        for (const G4Track* secTrack : *secondaries) {
          G4String secName = secTrack->GetDefinition()->GetParticleName();
          
          // If a muon is created from pion decay
          if (secName == "mu+" || secName == "mu-") {
            G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
            G4double secEnergy = secTrack->GetKineticEnergy();
            
            G4cout << "\n!!! PION DECAY DETECTED !!!" << G4endl;
            G4cout << particleName << " → " << secName << G4endl;
            G4cout << "Position: " << position/mm << " mm" << G4endl;
            G4cout << "Parent Energy: " << energy/MeV << " MeV" << G4endl;
            G4cout << "Muon Energy: " << secEnergy/MeV << " MeV" << G4endl;
          }
        }
      }
    }
  }
  
  // Get current volume
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
  
  // Check for muons and charged pions in detector 1
  if (step->IsFirstStepInVolume() && volume == fDetector1Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum();
    G4double totalEnergy = step->GetPreStepPoint()->GetTotalEnergy();
    
    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons
      fDetector1Particles[particleName]++;
      runAction->RecordParticleToExcel(particleName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(1, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector1();
      }
      
      G4cout << "\n!!! MUON DETECTED IN DETECTOR 1 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions
      fDetector1Particles[particleName]++;
      runAction->RecordParticleToExcel(particleName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(1, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector1();
      }
      
      G4cout << "\n!!! PION DETECTED IN DETECTOR 1 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
  }

  // Check for muons and charged pions in detector 2
  if (step->IsFirstStepInVolume() && volume == fDetector2Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum();
    G4double totalEnergy = step->GetPreStepPoint()->GetTotalEnergy();
    
    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons at Detector 2
      fDetector2Particles[particleName]++;
      G4String recordName = "2" + particleName;
      runAction->RecordParticleToExcel(recordName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(2, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector2();
      }
      
      G4cout << "\n!!! MUON DETECTED IN DETECTOR 2 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions at Detector 2
      fDetector2Particles[particleName]++;
      G4String recordName = "2" + particleName;
      runAction->RecordParticleToExcel(recordName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(2, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector2();
      }
      
      G4cout << "\n!!! PION DETECTED IN DETECTOR 2 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
  }

  // Check for muons and charged pions in detector 3
  if (step->IsFirstStepInVolume() && volume == fDetector3Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum();
    G4double totalEnergy = step->GetPreStepPoint()->GetTotalEnergy();
    
    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons at Detector 3
      fDetector3Particles[particleName]++;
      G4String recordName = "3" + particleName;
      runAction->RecordParticleToExcel(recordName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(3, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector3();
      }
      
      G4cout << "\n!!! MUON DETECTED IN DETECTOR 3 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions at Detector 3
      fDetector3Particles[particleName]++;
      G4String recordName = "3" + particleName;
      runAction->RecordParticleToExcel(recordName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(3, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector3();
      }
      
      G4cout << "\n!!! PION DETECTED IN DETECTOR 3 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
  }
  
  // RF Cavity functionality - track particles entering the cavity
  if (step->IsFirstStepInVolume() && volume == fRFCavityVolume) {
    // Get position and momentum at cavity entrance
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum();
    G4double totalEnergy = step->GetPreStepPoint()->GetTotalEnergy();
    
    // Store initial properties for particles of interest
    if (particleName == "mu+" || particleName == "mu-" || 
        particleName == "pi+" || particleName == "pi-") {
        
      // Count particles
      fRFCavityParticles[particleName]++;
      
      // Store entry data for later comparison
      fRFCavityEntranceEnergy[track->GetTrackID()] = totalEnergy;
      fRFCavityEntranceMomentum[track->GetTrackID()] = momentum;
      
      // Record 6D vector for cavity entrance
      G4String recordName = "RF_in_" + particleName;
      runAction->Record6DVector(5, recordName, position, momentum, totalEnergy);
      
      G4cout << "\n!!! PARTICLE ENTERING RF CAVITY !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << totalEnergy/MeV << " MeV" << G4endl;
      G4cout << "Momentum: " << momentum.mag()/MeV << " MeV/c" << G4endl;
    }
  }
  
  // RF Cavity functionality - track particles exiting the cavity
  if (step->IsLastStepInVolume() && volume == fRFCavityVolume) {
    // Get position and momentum at cavity exit
    G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPostStepPoint()->GetMomentum();
    G4double totalEnergy = step->GetPostStepPoint()->GetTotalEnergy();
    G4int trackID = track->GetTrackID();
    
    // Check if this particle was tracked at entrance
    if (particleName == "mu+" || particleName == "mu-" || 
        particleName == "pi+" || particleName == "pi-") {
      
      // Record 6D vector for cavity exit
      G4String recordName = "RF_out_" + particleName;
      runAction->Record6DVector(6, recordName, position, momentum, totalEnergy);
      
      // Calculate energy gain if we tracked this particle at entrance
      if (fRFCavityEntranceEnergy.find(trackID) != fRFCavityEntranceEnergy.end()) {
        G4double initialEnergy = fRFCavityEntranceEnergy[trackID];
        G4ThreeVector initialMomentum = fRFCavityEntranceMomentum[trackID];
        G4double energyGain = totalEnergy - initialEnergy;
        
        G4cout << "\n!!! PARTICLE EXITING RF CAVITY !!!" << G4endl;
        G4cout << "Type: " << particleName << G4endl;
        G4cout << "Initial Energy: " << initialEnergy/MeV << " MeV" << G4endl;
        G4cout << "Final Energy: " << totalEnergy/MeV << " MeV" << G4endl;
        G4cout << "Energy Gain: " << energyGain/MeV << " MeV" << G4endl;
        G4cout << "Initial Momentum Z: " << initialMomentum.z()/MeV << " MeV/c" << G4endl;
        G4cout << "Final Momentum Z: " << momentum.z()/MeV << " MeV/c" << G4endl;
        
        // Clean up the tracking maps to avoid memory growth
        fRFCavityEntranceEnergy.erase(trackID);
        fRFCavityEntranceMomentum.erase(trackID);
      }
    }
  }

  // Check for muons and charged pions in detector 4
  if (step->IsFirstStepInVolume() && volume == fDetector4Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum();
    G4double totalEnergy = step->GetPreStepPoint()->GetTotalEnergy();
    
    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons at Detector 4
      fDetector4Particles[particleName]++;
      G4String recordName = "4" + particleName;
      runAction->RecordParticleToExcel(recordName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(4, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector4();
      }
      
      G4cout << "\n!!! MUON DETECTED IN DETECTOR 4 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions at Detector 4
      fDetector4Particles[particleName]++;
      G4String recordName = "4" + particleName;
      runAction->RecordParticleToExcel(recordName, energy);
      // Add 6D vector recording
      runAction->Record6DVector(4, particleName, position, momentum, totalEnergy);
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector4();
      }
      
      G4cout << "\n!!! PION DETECTED IN DETECTOR 4 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/MeV << " MeV" << G4endl;
    }
  }

  // Check if we are in scoring volume for energy deposition
  if (volume == fScoringVolume) {
    // Collect energy deposited in this step
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);
  }

}