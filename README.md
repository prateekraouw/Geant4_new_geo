Geant4 Beamline Simulation

A multi-threaded Geant4 application for simulating a proton beamline with magnetic solenoids, RF cavities, and detector planes. Records 6D phase-space vectors for muons, pions, and protons as they cross detector volumes.

✨ Features

Gaussian Proton Beam: Configurable transverse sigma via /gun/sigma UI command.

Multi-Threaded Output: Per-thread CSV logging of 6D vectors (x, px, y, py, z, pz, E) for each particle crossing detectors.

Detector Planes: Four detector volumes with entry recording; RF cavity entry/exit monitoring.

Particle Filters: Record only mu+, mu-, pi+, pi-, and proton crossings.

Flexible Logging: Merge per-thread CSVs easily using provided shell commands.

Energy Scoring: Collect and summarize energy deposition in a scoring volume.

Decay Monitoring: Prints pion→muon decay events with kinematic details.

🚀 Prerequisites

C++17 compiler

Geant4 v10.7 or newer

CLHEP library (for RandGauss)

CMake (v3.10+)

Make or Ninja build system

🛠️ Build Instructions

mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=./install ../
make -j$(nproc)
make install

▶️ Running the Simulation

Create a macro file (run.mac):

# Set up beam energy and sigma
/gun/energy 10 GeV
/gun/sigma 2.5 mm
# Set number of threads and events
/run/initialize
/run/numberOfThreads 4
/run/beamOn 10000

Launch:

./bin/YourApp run.mac

Output CSVs appear as:

6D_vector_run<runID>_t<threadID>.csv

Merge with:

head -n1 6D_vector_run0_t0.csv > merged.csv
tail -q -n +2 6D_vector_run0_t*.csv >> merged.csv

📂 Project Structure

.
├── build
│   ├── 6D_vector.csv
│   ├── 6D_vector_run0_t0.csv
│   ├── 6D_vector_run0_t-1.csv
│   ├── 6D_vector_run0_t1.csv
│   ├── 6D_vector_run0_t2.csv
│   ├── 6D_vector_run0_t3.csv
│   ├── 6D_vector_run0_t4.csv
│   ├── CMakeCache.txt
│   ├── CMakeFiles/
│   ├── cmake_install.cmake
│   ├── init_vis.mac
│   ├── Makefile
│   ├── merged.csv
│   ├── particle_data_run0_t*.csv
│   ├── plot.ipynb
│   ├── run.mac
│   ├── tungsten_sim
│   └── vis.mac
├── CMakeLists.txt
├── include/
│   ├── ActionInitialization.hh
│   ├── DetectorConstruction.hh
│   ├── EventAction.hh
│   ├── MomentumChicane.hh
│   ├── PhysicsList.hh
│   ├── PrimaryGeneratorAction.hh
│   ├── RFCavityField.hh
│   ├── RunAction.hh
│   ├── SolenoidSystem.hh
│   └── SteppingAction.hh
├── src/
│   ├── ActionInitialization.cc
│   ├── DetectorConstruction.cc
│   ├── EventAction.cc
│   ├── MomentumChicane.cc
│   ├── PhysicsList.cc
│   ├── PrimaryGeneratorAction.cc
│   ├── RFCavityField.cc
│   ├── RunAction.cc
│   ├── SolenoidSystem.cc
│   ├── SteppingAction.cc
│   └── temp.txt
├── main
├── particle_data0.csv
├── particle_data_run1.csv
├── README.md
├── run.mac
├── tungsten_sim.cc
└── vis.mac
---

## 📖 Usage Tips

- Use `/gun/sigma` to tune beam spot size.
- Adjust detector positions in `DetectorConstruction`.
- Enable ntuple merging via Geant4 analysis manager as an alternative to manual CSV merge.

---

## 👤 Author & License

**Prateek Rao** – University of Wisconsin–Madison

Licensed under the MIT License. See [LICENSE](LICENSE) for details.
