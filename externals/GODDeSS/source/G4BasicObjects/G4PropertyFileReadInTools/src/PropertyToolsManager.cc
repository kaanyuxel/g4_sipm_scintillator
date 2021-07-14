/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <VectorUtil.hh>

#include <boost/lexical_cast.hpp>

#include <PropertyToolsManager.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to add the elements (specified in the chemical component data in the material property file) to the material
 */
void PropertyToolsManager::AddElementsFromTable(G4Material*& material, GoddessProperties::tabular elementTable)
{
	// type cast vector<boost::any> to the needed vector
	vector<std::string> elements = VectorUtil::adapt<std::string>(elementTable["element"]);
	vector<double> massFractions = VectorUtil::adapt<double>(elementTable["mass_fraction"]);
	vector<double> atomsPerUnit = VectorUtil::adapt<double>(elementTable["atoms_per_unit"]);

	// get the elements and add them to the material
	if(massFractions.size() != 0)
	{
		G4double massFractionSum = 0.;
		G4String elementNames = "";

		for(unsigned int iter = 0; iter < elements.size(); ++iter)
		{
			G4Element* element = G4Element::GetElement(elements[iter].c_str(), false);   // without warning if element does not exist -> use your own warning message

			if(element)
			{
				massFractionSum += massFractions[iter];
				elementNames += element->GetName() + " ";
			}
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# PropertyToolsManager::AddElementsFromTable():" <<
				std::endl << "# The element \"" << elements[iter] << "\" has not been defined. You have to do that in your in the detector construction!" <<
				std::endl << "##########" <<
				std::endl << std::endl;

				// abort the whole programme, if a critical error occured
				exit(1);
			}
		}

		if(fabs(massFractionSum - 1.) < 1e-12)
		{
			for(unsigned int iter = 0; iter < elements.size(); ++iter)
			{
				G4Element* element = G4Element::GetElement(elements[iter].c_str(), false);   // without warning if element does not exist -> use your own warning message
				material->AddElement(element, massFractions[iter]);
			}
		}
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# PropertyToolsManager::AddElementsFromTable():" <<
			std::endl << "# The elements ( " << elementNames << ") add up to more than 100% mass fractions." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			// abort the whole programme, if a critical error occured
			exit(1);
		}
	}
	else if(atomsPerUnit.size() != 0)
	{
		G4bool atomsPerUnitBelowOne = false;
		G4String elementNames = "";

		for(unsigned int iter = 0; iter < elements.size(); ++iter)
		{
			G4Element* element = G4Element::GetElement(elements[iter].c_str(), false);   // without warning if element does not exist -> use your own warning message

			if(element)
			{
				if(atomsPerUnit[iter] < 1.) atomsPerUnitBelowOne = true;
				elementNames += element->GetName() + " ";
			}
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# PropertyToolsManager::AddElementsFromTable():" <<
				std::endl << "# The element \"" << elements[iter] << "\" has not been defined. You have to do that in your in the detector construction!" <<
				std::endl << "##########" <<
				std::endl << std::endl;

				// abort the whole programme, if a critical error occured
				exit(1);
			}
		}

		if( !atomsPerUnitBelowOne )
		{
			for(unsigned int iter = 0; iter < elements.size(); ++iter)
			{
				G4Element* element = G4Element::GetElement(elements[iter].c_str(), false);   // without warning if element does not exist -> use your own warning message
				material->AddElement(element, (G4int) atomsPerUnit[iter]);
			}
		}
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# PropertyToolsManager::AddElementsFromTable():" <<
			std::endl << "# At least one of the elements ( " << elementNames << ") has less than 1 atom per unit." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			// abort the whole programme, if a critical error occured
			exit(1);
		}
	}
}



/**
 *  Function to extract the spectrum from the Properties of the material property file:
 *  - if Properties contains a spectrum of the desired property, tabular_to_G4MaterialPropertyVector() is used
 *  - if Properties contains a constant of the desired property, it is changed to a spectrum
 *
 *  @return G4MaterialPropertyVector with the energy dependent spectrum of the desired property
 */
G4MaterialPropertyVector * PropertyToolsManager::GetPropertyDistribution(GoddessProperties properties, G4String propertyKey, G4bool invert, G4double linearMatchingCoefficient, G4double quadraticMatchingCoefficient)
{
	G4MaterialPropertyVector * propertyVector = 0;

	if(!isnan(linearMatchingCoefficient) && isnan(quadraticMatchingCoefficient)) quadraticMatchingCoefficient = 0.;

	if(properties.containsTabular(propertyKey))
	{
		propertyVector = new G4MaterialPropertyVector();
		propertyVector = tabular_to_G4MaterialPropertyVector(propertyKey, properties.getTabular(propertyKey), invert, linearMatchingCoefficient, quadraticMatchingCoefficient);
	}
	else if(properties.containsNumber(propertyKey))
	{
		G4double propertyValue = properties.getNumber(propertyKey);
		if(!isnan(linearMatchingCoefficient)) propertyValue = propertyValue * (linearMatchingCoefficient + propertyValue * quadraticMatchingCoefficient);

		propertyVector = new G4MaterialPropertyVector();

		for(unsigned int i = 0; i < EnergyRangeVector.size(); i++)
		{
			propertyVector->InsertValues(EnergyRangeVector[i], propertyValue);
		}
	}

	return propertyVector;   //returns a 0 pointer if the desired property could not be found
}



/**
 *  Function to transform a constant property into an energy dependent distribution.
 *
 *  @return G4MaterialPropertyVector with the energy dependent spectrum of the given property
 */
G4MaterialPropertyVector * PropertyToolsManager::GetPropertyDistribution(G4double propertyValue)
{
	G4MaterialPropertyVector * propertyVector = new G4MaterialPropertyVector();

	for(unsigned int i = 0; i < EnergyRangeVector.size(); i++)
	{
		propertyVector->InsertValues(EnergyRangeVector[i], propertyValue);
	}

	return propertyVector;
}



/**
 *  Function to extract the spectrum from a tabular (from Properties):
 *  - if necessary, calculates the energy from the wave length
 *  - fills the spectrum into a G4MaterialPropertyVector
 *
 *  @return G4MaterialPropertyVector with the energy dependent spectrum of the desired property
 */
G4MaterialPropertyVector * PropertyToolsManager::tabular_to_G4MaterialPropertyVector(G4String propertyKey, GoddessProperties::tabular tabular, G4bool invert, G4double linearMatchingCoefficient, G4double quadraticMatchingCoefficient)
{
	// type cast vector<boost::any> to the needed vector
	vector<double> wavelengths = VectorUtil::adapt<double>(tabular["wavelength"]);
	tabular.erase("wavelength");
	vector<double> energies = VectorUtil::adapt<double>(tabular["energy"]);
	tabular.erase("energy");
	vector<double> propertyValues = VectorUtil::adapt<double>((* tabular.begin()).second);

	// get the energy and and the property and add them to the G4MaterialPropertyVector
	if(energies.size() == 0 && wavelengths.size() != 0)
	{
		for(unsigned int iter = 0; iter < wavelengths.size(); ++iter)
		{
			energies.push_back(CLHEP::c_light * CLHEP::h_Planck / wavelengths[iter]);
		}
	}
	else if( !(energies.size() != 0 && wavelengths.size() == 0) )
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# PropertyToolsManager::tabular_to_G4MaterialPropertyVector():" <<
		std::endl << "# Something is wrong with the wavelengths/energies of the property tabular \"" << propertyKey << "\". Either both have entries or both have no entries." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		// abort the whole programme, if a critical error occured
		exit(1);
	}

	if( propertyValues.size() == 0 || !(propertyValues.size() == energies.size() || propertyValues.size() == wavelengths.size()) )
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# PropertyToolsManager::tabular_to_G4MaterialPropertyVector():" <<
		std::endl << "# There are no property entries in the property tabular \"" << propertyKey << "\" (or not exactly the same number as wavelength/energy entries)!" <<
		std::endl << "##########" <<
		std::endl << std::endl;

		// abort the whole programme, if a critical error occured
		exit(1);
	}

	G4MaterialPropertyVector * propertyVector = new G4MaterialPropertyVector();

	G4bool insideOpticalEnergyRange = false;
	G4bool insertLastEntry = false;
	G4bool alreadyPrintedCherenkovWarningForThisPropertyVector = false;	// only one warning per spectrum

	for(unsigned int iter = 0; iter < propertyValues.size(); ++iter)
	{
		G4double energy = NAN;
		G4double propertyValue = NAN;

		if(iter > 0)
		{
			if(energies[iter] <= energies[iter - 1])
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# PropertyToolsManager::tabular_to_G4MaterialPropertyVector():" <<
				std::endl << "# Property tabular \"" << propertyKey << "\": " << energies[iter]/CLHEP::eV << "eV is not larger than " << energies[iter - 1]/CLHEP::eV << "eV!" <<
				std::endl << "# The property entries have to be sorted by rising energy values!" <<
				std::endl << "##########" <<
				std::endl << std::endl;

				// abort the whole programme, if a critical error occured
				exit(1);
			}
			else if(  ! alreadyPrintedCherenkovWarningForThisPropertyVector
				 && propertyKey.find("n_ref") != std::string::npos	// for refractive indices only
				 && propertyKey.find("_real") == std::string::npos && propertyKey.find("_imag") == std::string::npos	// not for complex refractive indices
				 && propertyValues[iter] <= propertyValues[iter - 1]
				)
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# WARNING:" << propertyKey <<
				std::endl << "# PropertyToolsManager::tabular_to_G4MaterialPropertyVector():" <<
				std::endl << "# Property tabular \"" << propertyKey << "\": " << propertyValues[iter] << " (at " << energies[iter]/CLHEP::eV << "eV) is not larger than " << propertyValues[iter - 1] << " (at " << energies[iter - 1]/CLHEP::eV << "eV)!" <<
				std::endl << "# Be aware that Cherenkov radiation in Geant4 is designed for normal-dispersive materials only!" <<
				std::endl << "##########" <<
				std::endl << std::endl;

				alreadyPrintedCherenkovWarningForThisPropertyVector = true;
			}
		}

		if(energies[iter] > EnergyRangeVector[0] && energies[iter] < EnergyRangeVector[EnergyRangeVector.size() - 1])
		{
			if(!insideOpticalEnergyRange)
			{
				// estimating the property value for the minimal energy defined in EnergyRangeVector from rough linear fit (for all properties are defined in the same energy range)
				G4int iter_temp = iter;        // using the previous and the current data point for the linear fit
				if(iter_temp == 0) iter_temp++;   // if there is no previous: using the current and the next data point for the linear fit

				energy = EnergyRangeVector[0];
				propertyValue = getValueFromLinearFit(EnergyRangeVector[0], propertyValues, energies, iter_temp);

				if(invert) propertyValue = 1. - propertyValue;
				if(!isnan(linearMatchingCoefficient)) propertyValue = propertyValue * (linearMatchingCoefficient + propertyValue * quadraticMatchingCoefficient);

				propertyVector->InsertValues(energy, propertyValue);

				insideOpticalEnergyRange = true;
			}

			energy = energies[iter];
			propertyValue = propertyValues[iter];

			if(iter == propertyValues.size() - 1) insertLastEntry = true;
		}
		else if(insideOpticalEnergyRange && energies[iter] > EnergyRangeVector[EnergyRangeVector.size() - 1]) insertLastEntry = true;

		if(!isnan(energy) && !isnan(propertyValue))
		{
			if(invert) propertyValue = 1. - propertyValue;
			if(!isnan(linearMatchingCoefficient))
			{
				G4double propertyValue_temp = propertyValue * (linearMatchingCoefficient + propertyValue * quadraticMatchingCoefficient);
				if(isinf(propertyValue_temp)) propertyValue = copysign(std::numeric_limits<double>::max(), propertyValue);
				else propertyValue = propertyValue_temp;
			}

			propertyVector->InsertValues(energy, propertyValue);
		}

		if(insertLastEntry)
		{
			energy = EnergyRangeVector[EnergyRangeVector.size() - 1];
			propertyValue = getValueFromLinearFit(EnergyRangeVector[EnergyRangeVector.size() - 1], propertyValues, energies, iter);

			if(invert) propertyValue = 1. - propertyValue;
			if(!isnan(linearMatchingCoefficient))
			{
				G4double propertyValue_temp = propertyValue * (linearMatchingCoefficient + propertyValue * quadraticMatchingCoefficient);
				if(isinf(propertyValue_temp)) propertyValue = copysign(std::numeric_limits<double>::max(), propertyValue);
				else propertyValue = propertyValue_temp;
			}

			propertyVector->InsertValues(energy, propertyValue);

			break;
		}
	}

	return propertyVector;
}



G4double PropertyToolsManager::getValueFromLinearFit(G4double energy, vector<double> propertyValues, vector<double> energies, G4int iter)
{
	G4double value = 0.;

	if(fabs(propertyValues[iter - 1]) > 1e-12)
	{
		G4double deltaProperty = propertyValues[iter] - propertyValues[iter - 1];
		G4double deltaEnergy = energies[iter] - energies[iter - 1];
		G4double deltaEnergyRangeVector = energy - energies[iter];

		value = propertyValues[iter] + deltaProperty / deltaEnergy * deltaEnergyRangeVector;
	}

	return value;
}



/**
 *  Function to check if a <b> newly generated </b> material had already been registered to GEANT:
 *  - checks if a material with \em materialName already exists
 *    - <b> if not </b>, set the name of \em material to \em materialName
 *    - <b> else </b>, check if \em material and the existing material with \em materialName are the same (sameMaterials())
 *      - <b> if so </b>, remove the <b> > last < </b> material that was registered to GEANT
 *      - <b> else </b>, set the name of \em material to an extended \em materialName
 */
void PropertyToolsManager::checkIfMaterialAlreadyExists(G4String materialName, G4Material *& material)
{
	G4MaterialTable* matTable = (G4MaterialTable*)G4Material::GetMaterialTable();

	for(unsigned int iter = 0; iter < matTable->size(); iter++)
	{
		G4String tempMaterialName = materialName;
		if(iter != 0) tempMaterialName += ("_" + boost::lexical_cast<std::string>(iter)).c_str();

		G4Material* tempMaterial = 0;
		if(G4Material::GetMaterial(tempMaterialName, false))   // if a material with this name already exists, get it
		{
			tempMaterial = G4Material::GetMaterial(tempMaterialName);
		}
		else   // else rename the new one (and abort the function)
		{
			material->SetName(tempMaterialName);
			return;
		}

		if(sameMaterials(material, tempMaterial))   // if the new material is the same as the already existing one, delete it and keep the existing one (and abort the function)
		{
			matTable->pop_back();

// NOTE: Here I would like to delete the unused material to free disk space, but GEANT4 crashes (for some simulation setups) if I do so. Don't ask me why.....
// 			delete material->GetMaterialPropertiesTable();
// 			delete material;

			material = tempMaterial;
			return;
		}

		if( !(iter < matTable->size() - 1) )   // if the end of the material table is reached and the new material has not been found inside, rename it
		{
			material->SetName(tempMaterialName);
		}
	}
}



/**
 *  Function to check, if a spectrum is constant.
 */
G4bool PropertyToolsManager::isConstantProperty(G4MaterialPropertyVector * propertyVector)
{
	if(!propertyVector)
	{
		G4cout << "PropertyToolsManager::isConstantProperty(): propertyVector does not exist!" << G4endl;
		return false;
	}

	G4double propertyValue = propertyVector->Value(propertyVector->Energy(0));

	for(unsigned int i = 1; i < propertyVector->GetVectorLength(); i++)
	{
		if(fabs(propertyVector->Value(propertyVector->Energy(i)) - propertyValue) > 1e-12) return false;
	}

	return true;
}



/**
 *  Function to compair two material, if they are exactly the same.
 */
G4bool PropertyToolsManager::sameMaterials(G4Material * material1, G4Material * material2)
{
	if(   material1->GetNumberOfElements() != material2->GetNumberOfElements()
	   || fabs(material1->GetDensity() - material2->GetDensity()) > 1e-12   )
	{
		return false;
	}

	for(unsigned int iter = 0; iter < material1->GetNumberOfElements(); iter++)
	{
		if(material1->GetElement(iter) != material2->GetElement(iter)) return false;
	}

	G4MaterialPropertiesTable * propertiesTable1 = material1->GetMaterialPropertiesTable();
	G4MaterialPropertiesTable * propertiesTable2 = material2->GetMaterialPropertiesTable();


	// get the maps with the property spectra
	const std::map<G4String, G4MaterialPropertyVector*> * propertiesMap1 = propertiesTable1->GetPropertiesMap();
	const std::map<G4String, G4MaterialPropertyVector*> * propertiesMap2 = propertiesTable2->GetPropertiesMap();
	// get the maps with the constant properties
	const std::map<G4String, G4double> * constPropertiesMap1 = propertiesTable1->GetPropertiesCMap();
	const std::map<G4String, G4double> * constPropertiesMap2 = propertiesTable2->GetPropertiesCMap();

	// compair the number of properties
	if(propertiesMap1->size() != propertiesMap2->size() || constPropertiesMap1->size() != constPropertiesMap2->size()) return false;

	// compair the constant properties
	for(std::map<G4String, G4double>::const_iterator iterC1 = constPropertiesMap1->begin(); iterC1 != constPropertiesMap1->end(); iterC1++)
	{
		G4bool propertyFound = false;

		for(std::map<G4String, G4double>::const_iterator iterC2 = constPropertiesMap2->begin(); iterC2 != constPropertiesMap2->end(); iterC2++)
		{
			// compair the names of the constant properties
			if(iterC1->first == iterC2->first)
			{
				// compair the values of the constant properties
				if(iterC1->second != iterC2->second) return false;

				propertyFound = true;
				break;
			}
		}

		if(!propertyFound) return false;
	}

	// compair the property spectra
	for(std::map<G4String, G4MaterialPropertyVector*>::const_iterator iter1 = propertiesMap1->begin(); iter1 != propertiesMap1->end(); iter1++)
	{
		G4bool propertyFound = false;

		for(std::map<G4String, G4MaterialPropertyVector*>::const_iterator iter2 = propertiesMap2->begin(); iter2 != propertiesMap2->end(); iter2++)
		{
			// compair the names of the property spectra
			if(iter1->first == iter2->first)
			{
				// compair the values of the property spectra
				if(!properySpectrumEqual(iter1->second, iter2->second)) return false;

				propertyFound = true;
				break;
			}
		}

		if(!propertyFound) return false;
	}

	return true;
}



/**
 *  Function to compair two spectra, if they are exactly the same.
 */
G4bool PropertyToolsManager::properySpectrumEqual(G4MaterialPropertyVector * propertyVector1, G4MaterialPropertyVector * propertyVector2)
{
	if( (propertyVector1 && !propertyVector2) || (!propertyVector1 && propertyVector2) ) return false;

	if(propertyVector1 && propertyVector2)
	{
		if(fabs(propertyVector1->GetVectorLength() - propertyVector2->GetVectorLength()) > 1e-12) return false;

		for(unsigned int iter = 0; iter < propertyVector1->GetVectorLength(); iter++)
		{
			if(fabs(propertyVector1->Energy(iter) - propertyVector2->Energy(iter)) > 1e-12) return false;
			if(fabs(propertyVector1->Value(propertyVector1->Energy(iter)) - propertyVector2->Value(propertyVector2->Energy(iter))) > 1e-12) return false;
		}
	}

	return true;
}
