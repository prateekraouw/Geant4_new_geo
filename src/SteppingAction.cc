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
#include "G4ElectroMagneticField.hh"
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
    if (!fDetector1Volume) {
      auto detCon = static_cast<const DetectorConstruction*>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      fDetector1Volume = detCon->GetDetector1Volume();
      fDetector2Volume = detCon->GetDetector2Volume();
      fDetector3Volume = detCon->GetDetector3Volume();
      fDetector4Volume = detCon->GetDetector4Volume();
    }
    
    // 2) record each crossing into detectors 1–4
    if (step->IsFirstStepInVolume()) {
      auto volume = step->GetPreStepPoint()
                         ->GetTouchableHandle()
                         ->GetVolume()
                         ->GetLogicalVolume();
      G4int detID = 0;
      if      (volume == fDetector1Volume) detID = 1;
      else if (volume == fDetector2Volume) detID = 2;
      else if (volume == fDetector3Volume) detID = 3;
      else if (volume == fDetector4Volume) detID = 4;
    
      if (detID > 0) {
        auto* track = step->GetTrack();
        const G4String name = track->GetDefinition()->GetParticleName();
      
        // only muons & pions (as you had)
              if (name=="mu+" || name=="mu-" || name=="pi+" || name=="pi-") {
                const auto* pre = step->GetPreStepPoint();
                const G4ThreeVector pos = pre->GetPosition();
                const G4ThreeVector mom = pre->GetMomentum();
                const G4double E = track->GetKineticEnergy();
                const G4double t = pre->GetGlobalTime(); // ns
      
                // log (B is in internal Geant4 units; divide by 'tesla' when writing CSV)
                if (auto* runAct = dynamic_cast<RunAction*>(
                      const_cast<G4UserRunAction*>(G4RunManager::GetRunManager()->GetUserRunAction())))
                {
                  runAct->Record6DVector(detID, name,
                                         pre->GetPosition(),
                                         pre->GetMomentum(),
                                         track->GetTrackID(),
                                         step->GetTrack()->GetCurrentStepNumber());
                }

        }
    }
}
}