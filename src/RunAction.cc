#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"

RunAction::RunAction()
: G4UserRunAction()
{
  // Open the 6D vector file
  Open6DVectorFile();
}

RunAction::~RunAction()
{
  if (fOutputFile.is_open()) {
    fOutputFile.close();
  }
  
  // Close the 6D vector file
  Close6DVectorFile();
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  fSecondaryParticles.clear();
  fParticleCounts.clear();
  
  // Open Excel file for particle data
  G4String fileName = "particle_data" + std::to_string(run->GetRunID()) + ".csv";
  fOutputFile.open(fileName);
  
  // Write CSV header with more information
  if (fOutputFile.is_open()) {
    fOutputFile << "ParticleType,Energy" << std::endl;
    G4cout << "Recording particle data to file: " << fileName << G4endl;
  } else {
    G4cerr << "ERROR: Could not open output file " << fileName << G4endl;
  }
  
  // If 6D vector file was closed, reopen it for this run
  if (!file6DVector.is_open()) {
    Open6DVectorFile();
  }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  // Print simple particle summary
  G4cout << "\n=== PARTICLE SUMMARY ===" << G4endl;
  for (const auto& pair : fParticleCounts) {
    G4cout << pair.first << ": " << pair.second << G4endl;
  }
  G4cout << "=========================" << G4endl;
  
  // Close Excel file
  if (fOutputFile.is_open()) {
    fOutputFile.close();
    G4cout << "Particle data saved to Excel file" << G4endl;
  }
  
  // Ensure 6D vector data is flushed to disk
  if (file6DVector.is_open()) {
    file6DVector.flush();
    G4cout << "6D vector data flushed to disk" << G4endl;
  }
}

void RunAction::RecordParticleToExcel(const G4String& name, 
                                     const G4double& kineticEnergy)
{
  if (fOutputFile.is_open()) {
    G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
    
    // Write to Excel with enhanced information
    fOutputFile << name << ","
                << kineticEnergy/MeV << std::endl;
  }
  
  // Count this particle type for the summary
  CountParticle(name);
}



// Function to open 6D vector file
void RunAction::Open6DVectorFile()
{
  // Create a new file for 6D phase space data
  std::string filename = "6D_vector.csv";
  
  // Open the file in write mode (overwrite if exists)
  file6DVector.open(filename, std::ios::out);
  
  // Add header to the CSV file
  if (file6DVector.is_open()) {
    file6DVector << "Detector,ParticleType,x[cm],px[MeV/c],y[cm],py[MeV/c],z[cm],pz[MeV/c],TotalEnergy[MeV]" << std::endl;
    G4cout << "Opened 6D vector file: " << filename << G4endl;
  } else {
    G4cout << "ERROR: Could not open 6D vector file: " << filename << G4endl;
  }
}

// Function to write 6D vector data for each particle
void RunAction::Record6DVector(G4int detectorID, 
                              const G4String& particleName, 
                              const G4ThreeVector& position, 
                              const G4ThreeVector& momentum,
                              G4double totalEnergy)
{
  if (file6DVector.is_open()) {
    file6DVector << detectorID << ","
                << particleName << ","
                << position.x()/cm << ","
                << momentum.x()/MeV << ","
                << position.y()/cm << ","
                << momentum.y()/MeV << ","
                << position.z()/cm << ","
                << momentum.z()/MeV << ","
                << totalEnergy/MeV << std::endl;
  }
}

// Function to close 6D vector file
void RunAction::Close6DVectorFile()
{
  if (file6DVector.is_open()) {
    file6DVector.close();
    G4cout << "Closed 6D vector file" << G4endl;
  }
}