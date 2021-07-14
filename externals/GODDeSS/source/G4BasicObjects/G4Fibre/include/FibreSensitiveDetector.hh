/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef FIBRESENSITIVEDETECTOR_HH_
#define FIBRESENSITIVEDETECTOR_HH_

#include <G4VSensitiveDetector.hh>

#include <GODDeSS_DataStorage.hh>



// class variables begin with capital letters, local variables with small letters



///  a class defining the action if a sentitive detector is hit
class FibreSensitiveDetector : public G4VSensitiveDetector
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 */
	FibreSensitiveDetector(G4String SDname, GODDeSS_DataStorage * dataStorage)
	// inheriting from other classes:
	: G4VSensitiveDetector(SDname)
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	, DataStorage(dataStorage)
	{
	}

	/**
	 *  Destructor (empty)
	 */
	~FibreSensitiveDetector()
	{
	}



	/**
	 *  Function to process the hits of the sensitive detector.\n
	 *  It will be executed by GEANT4, when a particle starts its step within the sensitive detector.
	 */
	G4bool ProcessHits(G4Step * theStep, G4TouchableHistory *);

	void Initialize(G4HCofThisEvent *)
	{
	}

	void EndOfEvent(G4HCofThisEvent *)
	{
	}

private:
	GODDeSS_DataStorage * DataStorage;

};

#endif /* FIBRESENSITIVEDETECTOR_HH_ */
