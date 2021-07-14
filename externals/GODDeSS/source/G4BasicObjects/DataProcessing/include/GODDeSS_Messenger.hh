/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef GODDeSS_Messenger_h
#define GODDeSS_Messenger_h 1

#include <GODDeSS_DataStorage.hh>
#include <PropertyToolsManager.hh>
#include <ScintillatorTileConstructor.hh>
#include <FibreConstructor.hh>
#include <PhotonDetectorConstructor.hh>



// class variables begin with capital letters, local variables with small letters



///  a class to store and provide variables which are needed in different parts of the GODDeSS
class GODDeSS_Messenger
{
public:

	/**
	 *  Constructor:
	 *  - sets class variables to default values
	 *  - creates a data storage object (GODDeSS_DataStorage)
	 *  - creates a property tools manager object (PropertyToolsManager)
	 */
	GODDeSS_Messenger(vector<G4double> energyRangeVector)
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: STConstructor(0)
	, FConstructor(0)
	, PDConstructor(0)
	{
		DataStorage = new GODDeSS_DataStorage();
		PropertyTools = new PropertyToolsManager(energyRangeVector);
	}

	/**
	 *  Destructor (empty)
	 */
	~GODDeSS_Messenger()
	{
	}



private:
	GODDeSS_DataStorage * DataStorage;
	PropertyToolsManager * PropertyTools;

	ScintillatorTileConstructor * STConstructor;
	FibreConstructor * FConstructor;
	PhotonDetectorConstructor * PDConstructor;

public:
	/**
	 *  @return pointer to the GODDeSS_DataStorage
	 */
	GODDeSS_DataStorage * GetDataStorage() const
	{ return DataStorage; }

	/**
	 *  @return pointer to the PropertyToolsManager
	 */
	PropertyToolsManager * GetPropertyToolsManager() const
	{ return PropertyTools; }

	/**
	 *  Set the pointer to the ScintillatorTileConstructor.
	 */
	void SetScintillatorTileConstructor(ScintillatorTileConstructor * scintillatorTileConstructor)
	{ STConstructor = scintillatorTileConstructor; }
	/**
	 *  @return pointer to the ScintillatorTileConstructor
	 */
	ScintillatorTileConstructor * GetScintillatorTileConstructor() const
	{ return STConstructor; }

	/**
	 *  Set the pointer to the FibreConstructor.
	 */
	void SetFibreConstructor(FibreConstructor * fibreConstructor)
	{ FConstructor = fibreConstructor; }
	/**
	 *  @return pointer to the FibreConstructor
	 */
	FibreConstructor * GetFibreConstructor() const
	{ return FConstructor; }

	/**
	 *  Set the pointer to the PhotonDetectorConstructor.
	 */
	void SetPhotonDetectorConstructor(PhotonDetectorConstructor * photonDetectorConstructor)
	{ PDConstructor = photonDetectorConstructor; }
	/**
	 *  @return pointer to the PhotonDetectorConstructor
	 */
	PhotonDetectorConstructor * GetPhotonDetectorConstructor() const
	{ return PDConstructor; }
};

#endif
