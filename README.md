# g4_sipm_scintillator


## Introduction
It is a GEANT4 example of a high energy physics detector that simulates particle detection elements. The example defines a detector that consists of two 20cm x 20cm x 1cm scintillator planes and each plane has one silicon photomultiplier (SiPM) at the middle of the left edge of the scintillator. When the charged particle passes through the scintillator planes, it crates photons and those photons are detected by SiPM. The SiPMs give a signal at the end of the event proportional to the detected photons. Therefore, detector properties are important to detect particles in high energy physics. The example simulates different brands/models of SiPM and scintillator devices thanks to the G4Sipm and GoDDeSS libraries. All the properties of SiPM(thermal noise, crosstalk, etc.) and scintillators(photon wavelength, scintillation yield, etc.) can be simulated with G4Sipm and GoDDeSS libraries. Normally, simulating these kinds of properties is not easy in GEANT4 and requires knowledge of coding and devices. This GEANT4 example helps to make a realistic simulation of a fundamental particle detector using G4Sipm and GoDDeSS external libraries.

## Requirements

## How to run g4_sipm_scintillator

## Check Results with Simple C++ Macro
