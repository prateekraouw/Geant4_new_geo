# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/prateek/Geant4_new_geo

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/prateek/Geant4_new_geo/build

# Include any dependencies generated for this target.
include CMakeFiles/tungsten_sim.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/tungsten_sim.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/tungsten_sim.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tungsten_sim.dir/flags.make

CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o: /home/prateek/Geant4_new_geo/tungsten_sim.cc
CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o -MF CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o.d -o CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o -c /home/prateek/Geant4_new_geo/tungsten_sim.cc

CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/tungsten_sim.cc > CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.i

CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/tungsten_sim.cc -o CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.s

CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o: /home/prateek/Geant4_new_geo/src/DetectorConstruction.cc
CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o -MF CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o -c /home/prateek/Geant4_new_geo/src/DetectorConstruction.cc

CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/DetectorConstruction.cc > CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.i

CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/DetectorConstruction.cc -o CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.s

CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o: /home/prateek/Geant4_new_geo/src/PhysicsList.cc
CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o -MF CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o -c /home/prateek/Geant4_new_geo/src/PhysicsList.cc

CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/PhysicsList.cc > CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.i

CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/PhysicsList.cc -o CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.s

CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o: /home/prateek/Geant4_new_geo/src/ActionInitialization.cc
CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o -MF CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o -c /home/prateek/Geant4_new_geo/src/ActionInitialization.cc

CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/ActionInitialization.cc > CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.i

CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/ActionInitialization.cc -o CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.s

CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o: /home/prateek/Geant4_new_geo/src/PrimaryGeneratorAction.cc
CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o -MF CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o -c /home/prateek/Geant4_new_geo/src/PrimaryGeneratorAction.cc

CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/PrimaryGeneratorAction.cc > CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.i

CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/PrimaryGeneratorAction.cc -o CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.s

CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o: /home/prateek/Geant4_new_geo/src/RunAction.cc
CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o -MF CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o -c /home/prateek/Geant4_new_geo/src/RunAction.cc

CMakeFiles/tungsten_sim.dir/src/RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/RunAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/RunAction.cc > CMakeFiles/tungsten_sim.dir/src/RunAction.cc.i

CMakeFiles/tungsten_sim.dir/src/RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/RunAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/RunAction.cc -o CMakeFiles/tungsten_sim.dir/src/RunAction.cc.s

CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o: /home/prateek/Geant4_new_geo/src/EventAction.cc
CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o -MF CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o -c /home/prateek/Geant4_new_geo/src/EventAction.cc

CMakeFiles/tungsten_sim.dir/src/EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/EventAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/EventAction.cc > CMakeFiles/tungsten_sim.dir/src/EventAction.cc.i

CMakeFiles/tungsten_sim.dir/src/EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/EventAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/EventAction.cc -o CMakeFiles/tungsten_sim.dir/src/EventAction.cc.s

CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o: /home/prateek/Geant4_new_geo/src/SteppingAction.cc
CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o -MF CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o -c /home/prateek/Geant4_new_geo/src/SteppingAction.cc

CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/SteppingAction.cc > CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.i

CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/SteppingAction.cc -o CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.s

CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o: /home/prateek/Geant4_new_geo/src/RFCavityField.cc
CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o -MF CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o -c /home/prateek/Geant4_new_geo/src/RFCavityField.cc

CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/RFCavityField.cc > CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.i

CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/RFCavityField.cc -o CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.s

CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o: /home/prateek/Geant4_new_geo/src/SolenoidSystem.cc
CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o -MF CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o -c /home/prateek/Geant4_new_geo/src/SolenoidSystem.cc

CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/SolenoidSystem.cc > CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.i

CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/SolenoidSystem.cc -o CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.s

CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o: CMakeFiles/tungsten_sim.dir/flags.make
CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o: /home/prateek/Geant4_new_geo/src/MomentumChicane.cc
CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o: CMakeFiles/tungsten_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o -MF CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o.d -o CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o -c /home/prateek/Geant4_new_geo/src/MomentumChicane.cc

CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/prateek/Geant4_new_geo/src/MomentumChicane.cc > CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.i

CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/prateek/Geant4_new_geo/src/MomentumChicane.cc -o CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.s

# Object files for target tungsten_sim
tungsten_sim_OBJECTS = \
"CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o" \
"CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o"

# External object files for target tungsten_sim
tungsten_sim_EXTERNAL_OBJECTS =

tungsten_sim: CMakeFiles/tungsten_sim.dir/tungsten_sim.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/DetectorConstruction.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/PhysicsList.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/ActionInitialization.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/PrimaryGeneratorAction.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/RunAction.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/EventAction.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/SteppingAction.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/RFCavityField.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/SolenoidSystem.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/src/MomentumChicane.cc.o
tungsten_sim: CMakeFiles/tungsten_sim.dir/build.make
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4Tree.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4FR.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4GMocren.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4visHepRep.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4RayTracer.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4VRML.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4ToolsSG.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4OpenGL.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4vis_management.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4modeling.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4interfaces.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4mctruth.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4geomtext.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4error_propagation.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4readout.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4physicslists.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4run.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4event.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4tracking.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4parmodels.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4processes.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4digits_hits.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4track.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4particles.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4geometry.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4materials.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4graphics_reps.so
tungsten_sim: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.15.13
tungsten_sim: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.15.13
tungsten_sim: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.15.13
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4analysis.so
tungsten_sim: /usr/lib/x86_64-linux-gnu/libexpat.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4zlib.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4intercoms.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4global.so
tungsten_sim: /home/prateek/geant4-v11.3.2/geant4-install/lib/libG4ptl.so.3.0.0
tungsten_sim: /usr/local/lib/libCLHEP-2.4.7.1.so
tungsten_sim: CMakeFiles/tungsten_sim.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/home/prateek/Geant4_new_geo/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Linking CXX executable tungsten_sim"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tungsten_sim.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tungsten_sim.dir/build: tungsten_sim
.PHONY : CMakeFiles/tungsten_sim.dir/build

CMakeFiles/tungsten_sim.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tungsten_sim.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tungsten_sim.dir/clean

CMakeFiles/tungsten_sim.dir/depend:
	cd /home/prateek/Geant4_new_geo/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/prateek/Geant4_new_geo /home/prateek/Geant4_new_geo /home/prateek/Geant4_new_geo/build /home/prateek/Geant4_new_geo/build /home/prateek/Geant4_new_geo/build/CMakeFiles/tungsten_sim.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/tungsten_sim.dir/depend

