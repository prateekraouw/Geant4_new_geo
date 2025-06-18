# Geant4 Solenoid Beamline Simulation

This project simulates a realistic particle beamline using **Geant4**, designed with a series of magnetic solenoids that guide and focus a charged particle beam. The simulation incorporates:
- Tapering solenoids
- Adjustable solenoid geometry
- Realistic magnetic fringe fields (using `tanh` transitions)
- Customizable gaps between solenoids

It is intended for use in high-energy physics simulations, especially muon or proton beam transport systems.

---

## 📁 Project Structure


---

## ⚙️ Prerequisites

- [Geant4 (v11+)](https://geant4.web.cern.ch/support/download) built with visualization and multithreading
- CMake (>= 3.16)
- C++17-compatible compiler (e.g., g++, clang++)
- Linux or macOS (tested on AlmaLinux and Ubuntu)

---

## 🔧 Building the Project

1. Source Geant4 environment:
   ```bash
    source /path/to/geant4-install/bin/geant4make.sh
   ```
   ```bash
   mkdir build
   ```
   ```bash
   cd build
   ```
   ```bash
   cmake ..
    ```
   ```bash
   make -j$(nproc)
   ```
## run in GUI mode
  ```bash
  ./tungstem_sim 
  ```
## run in command line mode
  ```bash
  ./tungstem_sim run.mac 
  ```
  
## Analysis
I. A Jupyter-Notebook plot.ipynb exists in the `build` directory.
- i. It has all the code build in for field and space-phase analysis
 
## Data for Analysis ( in `build` Directory)
1.6D_vector.csv 
- i. for space-phase Analysis
## 
2. all_23_solenoids.csv 
- ii. for solenoid field analysis
