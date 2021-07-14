/*
 * GoddessProperties.hh
 *
 * @date Mar 16, 2012
 * @author Tim Niggemann, III Phys. Inst. A, RWTH Aachen University
 * @copyright Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Germany License
 */

#ifndef GODDESSPROPERTIES_HH_
#define GODDESSPROPERTIES_HH_

#include <iostream>
#include <map>
#include <string>
#include <boost/any.hpp>
#include <boost/regex.hpp>

/**
 * This class can parse ".properties" files containing key-value-pairs and tables.
 * Uses boost::regex to match the contents.
 */
class GoddessProperties {
private:
	std::string filename;
	std::map<std::string, boost::any> entries;
	boost::regex emptyRx;
	boost::regex commentRx;
	boost::regex commentInLineRx;
	boost::regex numberRx;
	boost::regex keyNumberValuePairRx;
	boost::regex keyNumberValueUnitTripletRx;
	boost::regex keyStringValuePairRx;
	boost::regex tabularStartRx;
	boost::regex descriptorRx;

protected:
	/**
	 * @param str - the string to match.
	 * @return bool - true of matched "celsius" or "C".
	 */
	bool matchesCelsius(std::string str);

	/**
	 * @param str - the string to match.
	 * @return bool - true of matched "perCent" or "%".
	 */
	bool matchesPerCent(std::string str);

	/**
	 * @param str - the string to match.
	 * @return bool - true of matched "inf", "infinity", "Inf", "Infinity", "INF", or "INFINITY".
	 */
	bool matchesInfinity(std::string str);

	/**
	 * @param value - the numeric value.
	 * @param unitStr - the unit string.
	 * @return double - the converted value.
	 */
	double convert(double value, std::string unitStr);

	/**
	 * @param in - the input stream.
	 * @param key - the key of the tabular.
	 */
	void parseTabular(std::ifstream* in, std::string key);

public:
	/**
	 * Type represents a table in the property file.
	 */
	typedef std::map<std::string, std::vector<boost::any> > tabular;

	GoddessProperties();
	virtual ~GoddessProperties();

	/**
	 * @param filename - the properties file name.
	 * @return boolean - true if the file could be parsed completely.
	 */
	bool load(std::string filename);

	/**
	 * @param key - the key.
	 * @return double - the number.
	 */
	double getNumber(std::string key) const;

	/**
	 * @param key - the key.
	 * @return string - the string.
	 */
	std::string getString(std::string key) const;

	/**
	 * @param key - the key.
	 * @return map - the tabular data.
	 */
	tabular getTabular(std::string key) const;

	/**
	 * @param key - the key.
	 * @return boolean - true if the key-number-pair exists.
	 */
	bool containsNumber(std::string key) const;

	/**
	 * @param key - the key.
	 * @return boolean - true if the key-string-pair exists.
	 */
	bool containsString(std::string key) const;

	/**
	 * @param key - the key.
	 * @return boolean - true if the key-tabular-pair exists.
	 */
	bool containsTabular(std::string key) const;

	/**
	 * @return prints all properties to a string which is compatible with the file format.
	 */
	std::string toString() const;

	/**
	 * Print all read properties.
	 */
	void print() const;

	/**
	 * @return string - the file name.
	 */
	std::string getFilename() const {
		return filename;
	}
};

#endif /* GODDESSPROPERTIES_HH_ */
