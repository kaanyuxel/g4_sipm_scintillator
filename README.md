# g4_sipm_scintillator


## Introduction
This source code is a [Geant4](https://github.com/Geant4/geant4) example for a high energy physics particle detector that simulates particle detection elements. The example defines a detector that consists of two 20cm x 20cm x 1cm scintillator planes and each plane has one silicon photomultiplier (SiPM) at the middle of the left edge of the scintillator. When the charged particle passes through the scintillator planes, it crates photons and those photons are detected by SiPM. The SiPMs give a signal at the end of the event proportional to the detected photons. Therefore, detector properties are important to detect particles in high energy physics. The example simulates different brands/models of SiPM and scintillator devices thanks to the G4Sipm and GoDDeSS libraries. All the properties of SiPM(thermal noise, crosstalk, etc.) and scintillators(photon wavelength, scintillation yield, etc.) can be simulated with G4Sipm and GoDDeSS libraries. Normally, simulating these kinds of properties is not easy in Geant4 and requires knowledge of coding and devices. This Geant4 example helps to make a realistic simulation of a fundamental particle detector using external of G4Sipm and GoDDeSS libraries.

### G4Sipm
The G4Sipm is simulation toolkit to simulate SiPMs. The G4Sipm has been developed from the laboratory measurements and existing publications. It can be used in Geant4 applications or tested externally. For further information, please check the following links;
* [Academic Article](http://dx.doi.org/10.1016/j.nima.2015.01.067)
* [Documantation Page](http://g4sipm.readthedocs.io/en/latest/index.html) 
* [Git Page](https://github.com/ntim/g4sipm)

Even though it can be installed externally, the source of G4Sipm has been added to the example under *externals/g4sipm* folder and it will be installed automatically.

### GODDeSS
The GODDeSS(Geant4 Objects for Detailed Detectors with Scintillators and SiPMs) is a Geant4 extension for easy modelling of optical detector components like scintillator, fibre, etc. One can define optical detector components inside Geant4. However, it is not easy for non-experienced users and this leads user to make mistakes in the research field. Moreover, the popular scintillator products in the high energy physics field has been defined in its database. User can define experiment setups using GODDeSS libraries in the Geant4 simulation source code like in this example. For further information, please check the following links;
* [Academic Article](https://iopscience.iop.org/article/10.1088/1748-0221/12/04/P04026)
* [Documantation Page](https://git.rwth-aachen.de/3pia/forge/goddess-package/-/wikis/Documentation) 
* [Git Page](https://git.rwth-aachen.de/3pia/forge/goddess-package)

Althoug it can be installed externally, the source of GODDeSS(v4.3) has been added to the example under the *externals/GODDeSS* folder and it will be installed automatically.

## Requirements
The programs below have to be installed before running the example. Otherwise, the code will not compile.
 
* [Geant4](https://github.com/Geant4/geant4) (4.10 or newer)
* [Boost](http://www.boost.org/) (1.50.0 or newer)
* [ROOT](https://root.cern.ch) (5.34 or newer)
* [CMake](https://cmake.org/) (3.3 or newer)
* [zlib](https://zlib.net/) (1.2.9 or newer)

## How to install and run g4_sipm_scintillator

If you already installed all the requirements, you can install and run the example. Let's start with cloning the code,
```
git clone https://github.com/kaanyuxel/g4_sipm_scintillator.git
```
Then, create folder with a name *build* inside of g4_sipm_scintillator and access that folder;
```
cd g4_sipm_scintillator
mkdir build
cd build
```
Finally you can install it with typing following commands,
```
cmake ..
make -j4
```
If it compiles without any error, it creates a executable called *sipm_scintillator_app* and you can run it;
```
./sipm_scintillator_app
```
It will open a Qt window and one can see two scintillator planes with a seperation. In order to start experiment, type the following command on *Qt Session* line;
```
run/beamOn 1
```
This command will send one 5 GeV muon to scintillator planes. If you want to send more, you just type *run/beamOn N* (N is number of particle that will be sent). At the end of run, you can see the ROOT file with a name *time.root* which contains signal related branches. 

## Check Results with Simple C++ Macro
To see the produced signals by SiPMs, you can run the macro which is called *DrawSignal.cpp* with a command;
```
root -l DrawSignal.cpp
```
It opens two canvases (one shows the signals with dots and the other with lines) and the signals can be seen on produced canvases. 
