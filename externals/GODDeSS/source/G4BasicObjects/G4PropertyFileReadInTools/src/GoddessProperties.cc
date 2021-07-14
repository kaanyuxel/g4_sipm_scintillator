/*
 * GoddessProperties.cc
 *
 * @date Mar 16, 2012
 * @author Tim Niggemann, III Phys. Inst. A, RWTH Aachen University
 * @copyright Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Germany License
 */

#include "GoddessProperties.hh"

#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>

#include <G4UnitsTable.hh>

#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Units/PhysicalConstants.h>

#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>

GoddessProperties::GoddessProperties() {
	std::string spacePattern = "[[:space:]]*";
	std::string stringPattern = "([a-zA-Z0-9_-]+)";
	std::string numberPattern = "([+-]?[0-9]*[.]?[0-9]*[eE]?[+-]?[0-9]*)";
	std::string stringOrNumberPattern = "([a-zA-Z0-9_-]+|[+-]?[0-9]*[.]?[0-9]*[eE]?[+-]?[0-9]*)";
	std::string keyPattern = stringPattern + spacePattern + "[=:]{1}" + spacePattern;
	std::string unitPattern = spacePattern + "([1a-zA-Z%]{1}[a-zA-Z0-9/*]*)" + spacePattern;
	// Define matchers.
	emptyRx = boost::regex("[^[:w:]]*");
	commentRx = boost::regex("[^[:w:]]*#+.*");
	commentInLineRx = boost::regex("([^#]+)(#.*)");
	numberRx = boost::regex(spacePattern + numberPattern + spacePattern);
	keyNumberValuePairRx = boost::regex(keyPattern + numberPattern + spacePattern);
	keyNumberValueUnitTripletRx = boost::regex(keyPattern + numberPattern + spacePattern + "[*]{1}" + unitPattern);
	keyStringValuePairRx = boost::regex(keyPattern + stringPattern + spacePattern);
	tabularStartRx = boost::regex(keyPattern + "(tabular)" + spacePattern);
	descriptorRx = boost::regex(stringPattern + spacePattern + "/" + unitPattern);
}

GoddessProperties::~GoddessProperties() {
	//
}

bool GoddessProperties::matchesCelsius(std::string str) {
	return str == "Celsius";
}

bool GoddessProperties::matchesPerCent(std::string str) {
	return (str == "perCent" || str == "%");
}

bool GoddessProperties::matchesInfinity(std::string str) {
	boost::match_results<std::string::const_iterator> result;
	return regex_match(str, result, boost::regex("([+-]?[iI][nN][fF]([iI][nN][iI][tT][yY])?)"));
}

double GoddessProperties::convert(double value, std::string unitStr) {
	// Match empty unit string.
	if (unitStr == "") {
		return value;
	}
	// Try to match Celsius.
	if (matchesCelsius(unitStr)) {
		return value * G4UnitDefinition::GetValueOf("K") + CLHEP::STP_Temperature;
	}
	// Try to match %.
	if (matchesPerCent(unitStr)) {
		return value * CLHEP::perCent;
	}
	// Match with Geant4.
	std::vector<std::string> directMultiplicativeUnits;
	std::vector<std::string> inverseMultiplicativeUnits;

	while(unitStr.find("*") != std::string::npos || unitStr.find("/") != std::string::npos)
	{
		if(   (unitStr.find("*") != std::string::npos && unitStr.find("/") == std::string::npos)
		   || (unitStr.find("*") != std::string::npos && unitStr.find("/") != std::string::npos && unitStr.rfind("*") > unitStr.rfind("/"))   )
		{
			directMultiplicativeUnits.push_back(unitStr.substr(unitStr.rfind("*") + 1));
			unitStr = unitStr.substr(0, unitStr.rfind("*"));
		}
		else
		{
			inverseMultiplicativeUnits.push_back(unitStr.substr(unitStr.rfind("/") + 1));
			unitStr = unitStr.substr(0, unitStr.rfind("/"));
		}
	}

	directMultiplicativeUnits.push_back(unitStr);

	G4double unit = 1.;

	for (unsigned int iter = 0; iter < directMultiplicativeUnits.size(); iter++)
	{
		if(atof(directMultiplicativeUnits[iter].c_str())) unit *= atof(directMultiplicativeUnits[iter].c_str());
		else unit *= G4UnitDefinition::GetValueOf(G4String(directMultiplicativeUnits[iter]));
	}

	for (unsigned int iter = 0; iter < inverseMultiplicativeUnits.size(); iter++)
	{
		if(atof(inverseMultiplicativeUnits[iter].c_str())) unit /= atof(inverseMultiplicativeUnits[iter].c_str());
		else unit /= G4UnitDefinition::GetValueOf(G4String(inverseMultiplicativeUnits[iter]));
	}

	return value * unit;
}

void GoddessProperties::parseTabular(std::ifstream* in, std::string key) {
	tabular tab;
	// Temporary container.
	std::vector<std::string> columnNames;
	std::vector<std::string> columnUnits;
	std::vector<std::vector<boost::any> > columns;
	std::string line;
	boost::match_results<std::string::const_iterator> result;
	// Match header line.
	getline(*in, line);
	// Match comment at line end and cut it off.
	if (boost::regex_search(line, result, commentInLineRx, boost::match_all)) {
		line = result[1];
	}
	// Tokenizer definition.
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
	boost::char_separator<char> sep("\t");
	// Parse the header line.
	tokenizer headerTokens(line, sep);
	for (tokenizer::iterator it = headerTokens.begin(); it != headerTokens.end(); ++it) {
		std::string descriptor = (std::string) (*it);
		// Matched descriptor with unit.
		boost::match_results<std::string::const_iterator> regexResult;
		if (boost::regex_match(descriptor, regexResult, descriptorRx, boost::match_all)) {
			columnNames.push_back(regexResult[1]);
			columnUnits.push_back(regexResult[2]);
		} else {
			// Matched descriptor without.
			columnNames.push_back(descriptor);
			columnUnits.push_back("");
		}
		columns.push_back(std::vector<boost::any>());
	}
	// Exit if the header line could not be resolved.
	if (columnNames.size() == 0) {
		std::cerr << "GoddessProperties::parseTabular(in, key=\"" << key << "\"): could not parse tabular " << std::endl;
		return;
	}
	// Parse table content.
	while (getline(*in, line)) {
		// Match comment at line end and cut it off.
		if (boost::regex_search(line, result, commentInLineRx, boost::match_all)) {
			line = result[1];
		}

		tokenizer lineTokens(line, sep);
		// Determine number of tokens.
		unsigned int size = 0;
		for (tokenizer::iterator it = lineTokens.begin(); it != lineTokens.end(); ++it) {
			size++;
		}
		// Skip if the line was empty or the number of tokens mismatches the number of columns.
		if (line.size() == 0 || size != columns.size()) {
			// TODO(tim): reset file position.
			break;
		}
		// Insert values.
		int colId = 0;
		for (tokenizer::iterator it = lineTokens.begin(); it != lineTokens.end(); ++it) {
			std::string cell = (std::string) (*it);
			// Match number entry.
			if (regex_match(cell, result, numberRx)) {
				double value = atof(cell.c_str());
				columns[colId].push_back(convert(value, columnUnits[colId]));
			} else if(matchesInfinity(cell)) {
				double value = atof(cell.c_str());
				if(value < 0) value = - std::numeric_limits<double>::max();
				else value = std::numeric_limits<double>::max();
				columns[colId].push_back(convert(value, columnUnits[colId]));
			} else {
				// Match string entry.
				columns[colId].push_back(cell);
			}
			colId++;
		}
	}
	// Build up table.
	for (unsigned int i = 0; i < columnNames.size(); ++i) {
		tab[columnNames[i]] = columns[i];
	}
	entries[key] = tab;
}

bool GoddessProperties::load(std::string _filename) {
	entries.clear();
	// Open the file.
	std::ifstream in(_filename.c_str());
	if (!in.is_open()) {
		std::cerr << "GoddessProperties::load(filename = \"" << _filename << "\"): can not open file." << std::endl;
		return false;
	}
	// Iterate over all lines in the file.
	std::string line;
	boost::match_results<std::string::const_iterator> result;
	while (getline(in, line)) {
		// Match empty line.
		if (line.empty() || boost::regex_match(line, emptyRx)) {
			continue;
		}
		// Match comment line.
		if (regex_match(line, commentRx)) {
			continue;
		}
		// Match comment at line end and cut it off.
		if (boost::regex_search(line, result, commentInLineRx, boost::match_all)) {
			line = result[1];
		}
		// Match key-number-pair.
		if (boost::regex_search(line, result, keyNumberValuePairRx, boost::match_all)) {
			std::string key = result[1];
			std::string value = result[2];
			entries[key] = static_cast<double>(atof(value.c_str()));
			continue;
		}
		// Match key-number-unit-triplet.
		if (regex_search(line, result, keyNumberValueUnitTripletRx, boost::match_all)) {
			std::string key = result[1];
			std::string value = result[2];
			std::string unit = result[3];
			entries[key] = convert(static_cast<double>(atof(value.c_str())), unit);
			continue;
		}
		// Match tabular data.
		if (boost::regex_search(line, result, tabularStartRx, boost::match_all)) {
			parseTabular(&in, result[1]);
			continue;
		}
		// Match key-string-pair.
		if (boost::regex_search(line, result, keyStringValuePairRx, boost::match_all)) {
			std::string key = result[1];
			std::string value = result[2];
			entries[key] = value;
			continue;
		}
		// Line could not be matched, so break operation with an error.
		std::cerr << "GoddessProperties::load(filename = \"" << _filename << "\"): can not match \"" << line << "\"."
				<< std::endl;
		return false;
	}
	filename = _filename;
	return true;
}

double GoddessProperties::getNumber(std::string key) const {
	if (!containsNumber(key)) {
		std::cerr << "GoddessProperties::getNumber(key = \"" << key << "\"): does not exist." << std::endl;
		return NAN;
	}
	return boost::any_cast<double>(entries.at(key));
}

std::string GoddessProperties::getString(std::string key) const {
	if (!containsString(key)) {
		std::cerr << "GoddessProperties::getString(key = \"" << key << "\"): does not exist." << std::endl;
		return "";
	}
	return boost::any_cast<std::string>(entries.at(key));
}

bool GoddessProperties::containsNumber(std::string key) const {
	try {
		boost::any_cast<double>(entries.at(key));
		return true;
	} catch (...) {
		//
	}
	return false;
}

bool GoddessProperties::containsString(std::string key) const {
	try {
		boost::any_cast<std::string>(entries.at(key));
		return true;
	} catch (...) {
		//
	}
	return false;
}

GoddessProperties::tabular GoddessProperties::getTabular(std::string key) const {
	if (!containsTabular(key)) {
		std::cerr << "GoddessProperties::getTabular(key = \"" << key << "\"): does not exist." << std::endl;
		return tabular();
	}
	return boost::any_cast<tabular>(entries.at(key));
}

std::string GoddessProperties::toString() const {
	std::stringstream out;
	// Iterate over all properties.
	for (std::map<std::string, boost::any>::const_iterator it = entries.begin(); it != entries.end(); it++) {
		out << it->first << ":\t";
		const std::type_info& type = it->second.type();
		if (type == typeid(double)) {
			out << boost::any_cast<double>(it->second);
		} else if (type == typeid(std::string)) {
			out << boost::any_cast<std::string>(it->second);
		} else if (type == typeid(tabular)) {
			out << "tabular\n";
			tabular tab = boost::any_cast<tabular>(it->second);
			// Determine number of entries.
			unsigned int n = tab.begin()->second.size();
			// Print column names.
			out << "\t";
			for (tabular::const_iterator colIt = tab.begin(); colIt != tab.end(); colIt++) {
				out << colIt->first << "\t";
			}
			out << "\n";
			// Print values.
			for (unsigned int i = 0; i < n; i++) {
				out << "\t";
				for (tabular::const_iterator colIt = tab.begin(); colIt != tab.end(); colIt++) {
					boost::any entry = colIt->second.at(i);
					if (entry.type() == typeid(double)) {
						out << boost::any_cast<double>(entry);
					} else if (entry.type() == typeid(std::string)) {
						out << boost::any_cast<std::string>(entry);
					}
					out << "\t";
				}
				out << "\n";
			}
		}
		out << "\n";
	}
	return out.str();
}

void GoddessProperties::print() const {
	std::cout << toString() << std::endl;
}

bool GoddessProperties::containsTabular(std::string key) const {
	try {
		boost::any_cast<tabular>(entries.at(key));
		return true;
	} catch (...) {
		//
	}
	return false;
}

