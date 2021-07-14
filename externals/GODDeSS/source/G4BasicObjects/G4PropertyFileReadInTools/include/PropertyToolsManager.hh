/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#ifndef TOOLS_HH_
#define TOOLS_HH_

#include <G4MaterialPropertyVector.hh>
#include <G4Material.hh>

#include "GoddessProperties.hh"

using std::vector;



// class variables begin with capital letters, local variables with small letters



///  a class making the processing of the properties etc. much easier.
class PropertyToolsManager
{
public:

	/**
	 *  Constructor:
	 *  - sets the energy range on which the property is to be defined
	 */
	PropertyToolsManager(std::vector<G4double> energyRangeVector)
	// initialising the variables (doing it with default values, "" or "0" is just to prevent errors from wrongly initialised variables), this has to be done in the order of their appearance in the hh-file:
	: EnergyRangeVector(energyRangeVector)
	{
	}

	/**
	 *  Destructor (empty)
	 */
	virtual ~PropertyToolsManager()
	{
	}



	void AddElementsFromTable( G4Material *& material,		/**< material to which the elements is to be added */
			GoddessProperties::tabular elementTable	/**< tabular (from GoddessProperties) with the chemical component data of the material */
				 );

	G4MaterialPropertyVector * GetPropertyDistribution( GoddessProperties properties,	/**< tabular (from GoddessProperties) of the desired property */
							     G4String propertyKey,	/**< key to identify the desired property */
							     G4bool invert = false,	/**< bool if the distribution is to be inverted (i.e. 1 - property value is used) */
							     G4double linearMatchingCoefficient = NAN,	/**< linear coefficient to adjust the property value (i.e. property value * (linearMatchingCoefficient + property value * quadraticMatchingCoefficient) is used)*/
							     G4double quadraticMatchingCoefficient = NAN	/**< quadratic coefficient to adjust the property value (i.e. property value * (linearMatchingCoefficient + property value * quadraticMatchingCoefficient) is used)*/
							   );

	G4MaterialPropertyVector * GetPropertyDistribution( G4double propertyValue   /**< constant vakue that is to be transformed into a distribution */
							   );

	G4MaterialPropertyVector * tabular_to_G4MaterialPropertyVector( G4String propertyKey,   /**< key of the desired property */
									 GoddessProperties::tabular tabular,   /**< tabular (from GoddessProperties) of the desired property */
									 G4bool invert = false,		/**< bool if the distribution is to be inverted (i.e. 1 - property value is used) */
									 G4double linearMatchingCoefficient = NAN,	/**< linear coefficient to adjust the property value (i.e. property value * (linearMatchingCoefficient + property value * quadraticMatchingCoefficient) is used)*/
									 G4double quadraticMatchingCoefficient = NAN	/**< quadratic coefficient to adjust the property value (i.e. property value * (linearMatchingCoefficient + property value * quadraticMatchingCoefficient) is used)*/
								       );

	void checkIfMaterialAlreadyExists( G4String materialName,	/**< material name that is to be checked */
					   G4Material *& material	/**< material that is to be checked */
					 );

	G4bool isConstantProperty( G4MaterialPropertyVector * propertyVector   /**< spectrum to be checked if it is constant */
				 );


private:
	G4double getValueFromLinearFit( G4double energy,
				        vector<double> propertyValues,
				        vector<double> energies,
				        G4int iter
				      );

	G4bool sameMaterials( G4Material * material1,	/**< 1. material to be compaired */
			      G4Material * material2	/**< 2. material to be compaired */
			    );

	G4bool properySpectrumEqual( G4MaterialPropertyVector * propertyVector1,	/**< 1. spectrum to be compaired */
				     G4MaterialPropertyVector * propertyVector2		/**< 2. spectrum to be compaired */
				   );


	std::vector<G4double> EnergyRangeVector;
};

#endif
