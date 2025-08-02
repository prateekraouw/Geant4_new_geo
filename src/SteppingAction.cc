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
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh" 
#include <fstream>
#include <iomanip>
#include <memory>


SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(nullptr),
  fDetector1Volume(nullptr),
  fDetector2Volume(nullptr),
  fDetector3Volume(nullptr),
  fDetector4Volume(nullptr),
  fRFCavityVolume(nullptr)
{}

SteppingAction::~SteppingAction()
{
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


void SteppingAction::UserSteppingAction(const G4Step* step)
{
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

    /* 
    G4cout << "Detector 1 position: " << detectorConstruction->GetDetector1Position()/cm << " cm" << G4endl;
    G4cout << "Detector 2 position: " << detectorConstruction->GetDetector2Position()/cm << " cm" << G4endl;
    G4cout << "Detector 3 position: " << detectorConstruction->GetDetector3Position()/cm << " cm" << G4endl;
    G4cout << "Detector 4 position: " << detectorConstruction->GetDetector4Position()/cm << " cm" << G4endl;
    */
    
  }

  // Get the RunAction - using const_cast to handle the constness issue
  const G4UserRunAction* constRunAction = G4RunManager::GetRunManager()->GetUserRunAction();
  RunAction* runAction = const_cast<RunAction*>(dynamic_cast<const RunAction*>(constRunAction));

  // Get current track
  G4Track* track = step->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String particleName = particle->GetParticleName();
  G4double energy = track->GetKineticEnergy();


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
            G4cout << "Parent Energy: " << energy/GeV << " GeV" << G4endl;
            G4cout << "Muon Energy: " << secEnergy/GeV << " GeV" << G4endl;
          }
        }
      }
    }
  }

  // Get current volume
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
  G4String volumeName = volume->GetName();



  // Check for muons and charged pions in detector 1
  if (step->IsFirstStepInVolume() && volume == fDetector1Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum() /GeV ;
    G4double totalEnergy = step->GetTrack()->GetKineticEnergy();

    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons
      fDetector1Particles[particleName]++;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 1;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
        
        // Also record particle data to Excel
        runAct->RecordParticleToExcel(particleName, track->GetKineticEnergy());
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector1();
      }

      G4cout << "\n!!! MUON DETECTED IN DETECTOR 1 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions
      fDetector1Particles[particleName]++;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 1;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
        
        // Also record particle data to Excel
        runAct->RecordParticleToExcel(particleName, track->GetKineticEnergy());
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector1();
      }

      G4cout << "\n!!! PION DETECTED IN DETECTOR 1 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
  }

  // Check for muons and charged pions in detector 2
  if (step->IsFirstStepInVolume() && volume == fDetector2Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum() /GeV;
    G4double totalEnergy = step->GetTrack()->GetKineticEnergy();

    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons at Detector 2
      fDetector2Particles[particleName]++;
      G4String recordName = "2" + particleName;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 2;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
        
        // Also record particle data to Excel
        runAct->RecordParticleToExcel(particleName, track->GetKineticEnergy());
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector2();
      }

      G4cout << "\n!!! MUON DETECTED IN DETECTOR 2 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions at Detector 2
      fDetector2Particles[particleName]++;
      G4String recordName = "2" + particleName;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 2;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector2();
      }

      G4cout << "\n!!! PION DETECTED IN DETECTOR 2 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
  }

  // Check for muons and charged pions in detector 3
  if (step->IsFirstStepInVolume() && volume == fDetector3Volume) {
    // Get position and momentum for 6D vector
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum() /GeV;
    G4double totalEnergy = step->GetTrack()->GetKineticEnergy();

    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons at Detector 3
      fDetector3Particles[particleName]++;
      G4String recordName = "3" + particleName;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 3;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector3();
      }

      G4cout << "\n!!! MUON DETECTED IN DETECTOR 3 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions at Detector 3
      fDetector3Particles[particleName]++;
      G4String recordName = "3" + particleName;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 3;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector3();
      }

      G4cout << "\n!!! PION DETECTED IN DETECTOR 3 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
  }

  // RF Cavity functionality - track particles entering the cavity
  if (step->IsFirstStepInVolume() && volume == fRFCavityVolume) {
    // Get position and momentum at cavity entrance
    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum() /GeV;
    G4double totalEnergy = step->GetTrack()->GetKineticEnergy();

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
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 5;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }

      G4cout << "\n!!! PARTICLE ENTERING RF CAVITY !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << totalEnergy/GeV << " GeV" << G4endl;
      G4cout << "Momentum: " << momentum.mag()/GeV << " GeV/c" << G4endl;
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
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 6;
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }

      // Calculate energy gain if we tracked this particle at entrance
      if (fRFCavityEntranceEnergy.find(trackID) != fRFCavityEntranceEnergy.end()) {
        G4double initialEnergy = fRFCavityEntranceEnergy[trackID];
        G4ThreeVector initialMomentum = fRFCavityEntranceMomentum[trackID];
        G4double energyGain = totalEnergy - initialEnergy;

        G4cout << "\n!!! PARTICLE EXITING RF CAVITY !!!" << G4endl;
        G4cout << "Type: " << particleName << G4endl;
        G4cout << "Initial Energy: " << initialEnergy/GeV << " GeV" << G4endl;
        G4cout << "Final Energy: " << totalEnergy/GeV << " GeV" << G4endl;
        G4cout << "Energy Gain: " << energyGain/GeV << " GeV" << G4endl;
        G4cout << "Initial Momentum Z: " << initialMomentum.z()/GeV << " GeV/c" << G4endl;
        G4cout << "Final Momentum Z: " << momentum.z()/GeV << " GeV/c" << G4endl;

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
    G4ThreeVector momentum = step->GetPreStepPoint()->GetMomentum() /GeV;
    G4double totalEnergy = step->GetTrack()->GetKineticEnergy();

    if (particleName == "mu+" || particleName == "mu-") {
      // Count muons at Detector 4
      fDetector4Particles[particleName]++;
      G4String recordName = "4" + particleName;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 4;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddMuonAtDetector4();
      }

      G4cout << "\n!!! MUON DETECTED IN DETECTOR 4 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
    // Only count charged pions (pi+, pi-)
    else if (particleName == "pi+" || particleName == "pi-") {
      // Count charged pions at Detector 4
      fDetector4Particles[particleName]++;
      G4String recordName = "4" + particleName;
      //
      // Add 6D vector recording
      const auto* baseRun =
        G4RunManager::GetRunManager()->GetUserRunAction();
      
      // Try to down‐cast to *your* RunAction type
      const RunAction* constRun =
        dynamic_cast<const RunAction*>(baseRun);
      
      if ( constRun ) {
        // Cast away constness so we can call the non‐const member
        RunAction* runAct = const_cast<RunAction*>(constRun);
      
        // Now call with the correct detector ID, for example:
        //   1 for Detector 1, 2 for Detector 2, etc.
        constexpr G4int myDetID = 4;
      
        runAct->Record6DVector(
          myDetID,
          particleName,
          position,
          momentum,
          totalEnergy
        );
      }
      // Add to event counts
      if (fEventAction) {
        fEventAction->AddPionAtDetector4();
      }

      G4cout << "\n!!! PION DETECTED IN DETECTOR 4 !!!" << G4endl;
      G4cout << "Type: " << particleName << G4endl;
      G4cout << "Energy: " << track->GetKineticEnergy()/GeV << " GeV" << G4endl;
    }
  }
  // Check if we are in scoring volume for energy deposition
  if (volume == fScoringVolume) {
    // Collect energy deposited in this step
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);
  }
}

