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
#include <sstream>


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
  const G4int runID = run->GetRunID();

  if (G4Threading::IsMasterThread()) {
    G4cout << "[MASTER] BeginOfRunAction for run " << runID << G4endl;
    return; // Master doesn't generate events — don't open data files here.
  }

  const G4int tid = G4Threading::G4GetThreadId();
  // particle_data_...
  std::string pfile = "particle_data_run" + std::to_string(runID) +
                      "_t" + std::to_string(tid) + ".csv";
  fOutputFile.open(pfile, std::ios::out);
  if (fOutputFile) {
    fOutputFile << "ParticleType,Energy[MeV]\n";
    G4cout << "[T" << tid << "] writing particle data to " << pfile << G4endl;
  }

  // 6D_vector_...
  std::ostringstream vfn;
  vfn << "6D_vector_run" << runID << "_t" << tid << ".csv";
  file6DVector.open(vfn.str(), std::ios::out);
  if (file6DVector) {
    file6DVector << "Detector,ParticleType,"
                 << "x[mm],px[GeV/c],"
                 << "y[mm],py[GeV/c],"
                 << "z[mm],pz[GeV/c],"
                 << "Bx[T],By[T],Bz[T],"<<"track,step"<<"\n";
    G4cout << "[T" << tid << "] Opened 6D vector file: " << vfn.str() << G4endl;
  }
}

void RunAction::EndOfRunAction(const G4Run* run)
{
  if (G4Threading::IsMasterThread()) return;

  if (fOutputFile.is_open()) fOutputFile.close();
  if (file6DVector.is_open()) file6DVector.close(); // actually close here
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
void RunAction::Record6DVector(G4int detectorID, const G4String& particleName,
                             const G4ThreeVector& pos, const G4ThreeVector& mom,
                             G4int trackID, G4int stepNum)
{
    // Make sure the file is open - no mutex needed since each thread has its own file
    if (file6DVector.is_open()) {
        file6DVector << detectorID << ","
                  << particleName << ","
                  << pos.x()/cm << ","
                  << mom.x()/GeV << ","
                  << pos.y()/cm << ","
                  << mom.y()/GeV << ","
                  << pos.z()/cm << ","
                  << mom.z()/GeV << ","
                  << trackID << ","
                  << stepNum << std::endl;
    } else {
        G4int tid = G4Threading::G4GetThreadId();
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