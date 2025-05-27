#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"

int main(int argc, char** argv){
  // Remove this line to actually run the simulation
  // if(true)return 0;
  
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Construct the default run manager
  auto* runManager = G4RunManagerFactory::CreateRunManager();
  
  // Set mandatory initialization classes
  // Use DetectorConstruction instead of ChicaneConstruction
  runManager->SetUserInitialization(new DetectorConstruction());
  runManager->SetUserInitialization(new PhysicsList());
  runManager->SetUserInitialization(new ActionInitialization());
  
  // Initialize G4 kernel
  runManager->Initialize();
  
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();
  
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if (argc != 1) {
    // Batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command + fileName);
  }
  else {
    // Interactive mode
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }
  
  // Job termination
  delete visManager;
  delete runManager;
  
  return 0;
}