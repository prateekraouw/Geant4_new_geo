#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UnitsTable.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"    // defines G4Mutex, G4MUTEX_INITIALIZER     // still needed for G4AutoLock
#include <fstream>
#include <string>


// Remove global mutex and header flag - each thread will have its own file
RunAction::RunAction()
: G4UserRunAction()
{
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
  G4int runID = run->GetRunID();
  G4int tid   = G4Threading::G4GetThreadId();    // 0 = master, >0 = worker

  // Debug output
  G4cout << "=== DEBUG: BeginOfRunAction ===" << G4endl;
  G4cout << "Run ID: " << runID << G4endl;
  G4cout << "Thread ID: " << tid << G4endl;
  G4cout << "Is Master Thread: " << (tid == 0 ? "YES" : "NO") << G4endl;
  G4cout << "===============================" << G4endl;

  //
  //——— particle data file ———
  //
  std::string pfile = "particle_data_run"
                    + std::to_string(runID)
                    + "_t" + std::to_string(tid)
                    + ".csv";
  fOutputFile.open(pfile);
  if (fOutputFile.is_open()) {
    fOutputFile << "ParticleType,Energy[MeV]\n";
    G4cout << "Thread " << tid
           << " writing particle data to " << pfile << G4endl;
  } else {
    G4cerr << "ERROR: Thread " << tid << " could not open file: " << pfile << G4endl;
  }

  //
  //——— 6D vector file ———
  std::ostringstream vfn;
   vfn << "6D_vector_run" << runID << "_t" << tid << ".csv";
   file6DVector.open(vfn.str(), std::ios::out);
   if (file6DVector.is_open()) {
     // write header for this thread's file
     file6DVector << "Detector,ParticleType,"
                  << "x[cm],px[MeV/c],"
                  << "y[cm],py[MeV/c],"
                  << "z[cm],pz[MeV/c],"
                  << "TotalEnergy[MeV]\n";
     G4cout << "[Thread " << tid << "] Opened 6D vector file: " 
            << vfn.str() << G4endl;
   } else {
     G4cerr << "[Thread " << tid << "] ERROR: Could not open 6D vector file: "
            << vfn.str() << G4endl;
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

  //SaveMagneticFieldAlongZ();
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



// Function to write 6D vector data for each particle
void RunAction::Record6DVector(G4int detectorID, 
                              const G4String& particleName, 
                              const G4ThreeVector& position, 
                              const G4ThreeVector& momentum,
                              G4double totalEnergy)
{
    G4int tid = G4Threading::G4GetThreadId();
    
    // Debug output
    G4cout << "[Thread " << tid << "] Record6DVector called for detector " 
           << detectorID << ", particle " << particleName << G4endl;
    G4cout << "[Thread " << tid << "] File open: " << (file6DVector.is_open() ? "YES" : "NO") << G4endl;
    
    // Make sure the file is open - no mutex needed since each thread has its own file
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
        G4cout << "[Thread " << tid << "] Data written successfully" << G4endl;
    } else {
        G4cerr << "[Thread " << tid << "] ERROR: File not open for writing!" << G4endl;
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

 // For file handling

// Function to save the magnetic field data along the Z-axis to a CSV file