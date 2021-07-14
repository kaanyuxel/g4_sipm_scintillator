/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <zlib.h>
#include <assert.h>

#include <PhotonListSourceMessenger.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to define the action of the commands.\n
 *  It will be executed by GEANT4, when a command (created by this class) is called in the interactive mode.
 */
void PhotonListSourceMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == InfilePathCmd)
	{
		if(IsEventAvailable()) ParticleInformations.clear();

		ReadInPhotonList(newValue);
	}
}



/**
 *  Function to read a the list with photon data from a compressed (.gz) text file.\n
 *  It will be executed, when a new list with photon data is provided (with the "/pls/photonListPath" command).\n
 *  The compressed (.gz) text file is read in, then the data is decompressed, processed, and saved within the PhotonListSourceMessenger class.
 */
void PhotonListSourceMessenger::ReadInPhotonList(G4String infilePath)
{
	FILE * infile;
	infile = fopen(infilePath.c_str(), "r");
	if(infile == 0) std::cerr << "Properties::load(filename = \"" << infilePath << "\"): can not open file." << std::endl;
	else
	{

		G4cout << "Properties::load(filename = \"" << infilePath << "\"): reading in file." << G4endl << G4endl;
// fist part of the function to read in the file and decompress it
	/* Decompress from file source to file dest until stream ends or EOF.
	   errors are:
	      Z_OK on success
	      Z_MEM_ERROR if memory could not be allocated for processing
	      Z_DATA_ERROR if the deflate data is invalid or incomplete
	      Z_VERSION_ERROR if the version of zlib.h and the version of the library linked do not match
	      Z_ERRNO if there is an error reading or writing the files. */

		int readinBufferSize = 524288;
		int processBufferSize = 512;

		int returnValue = readinBufferSize;
		z_stream zlibStream;

		// declare variable to read in the deflated source file:
		unsigned char * inChar = new unsigned char[readinBufferSize];

		// declare variable to temporarily save the inflated file:
		unsigned char * tempBufferChar = new unsigned char[processBufferSize+1];
		tempBufferChar[processBufferSize] = 0;
		G4String processedDataString = "";

		// allocate inflate state
		zlibStream.zalloc = Z_NULL;
		zlibStream.zfree = Z_NULL;
		zlibStream.opaque = Z_NULL;
		zlibStream.avail_in = 0;
		zlibStream.next_in = Z_NULL;
		returnValue = inflateInit2(&zlibStream, (16+MAX_WBITS));
		if(returnValue != Z_OK)
		{
			G4cout << G4endl << "#####" << G4endl << "# Error when trying to read in the file to memory! The data will not be considered!" << G4endl << "#####" << G4endl << G4endl;
		}
		else
		{
			setDefaultValues();


			// decompress until deflate stream ends or end of file
			do
			{
				// try to read in the next sourceSize characters from the source file and tell zlibStream how many characters could be read
				zlibStream.avail_in = fread(inChar, 1, readinBufferSize, infile);

				if(ferror(infile))
				{
					(void) inflateEnd(&zlibStream);
					G4cout << G4endl << "#####" << G4endl << "# Error when trying to read in the file to memory! The data will not be considered!" << G4endl << "#####" << G4endl << G4endl;
					break;
				}
				// if no character could be read from the source file...
				if(zlibStream.avail_in == 0) break;

				// pass the read in data to zlibStream
				zlibStream.next_in = inChar;

				do
				{
					// tell zlibStream the size of tempBufferChar
					zlibStream.avail_out = processBufferSize;
					// tell zlibStream where to save the inflated data
					zlibStream.next_out = tempBufferChar;
					// inflate data and save the output of inflate() for error treatment
					returnValue = inflate(&zlibStream, Z_NO_FLUSH);

					// error treatment
					assert(returnValue != Z_STREAM_ERROR);  // state not clobbered
					switch(returnValue)
					{
						case Z_NEED_DICT:
							returnValue = Z_DATA_ERROR;     // and fall through
						case Z_DATA_ERROR:
						case Z_MEM_ERROR:
							(void)inflateEnd(&zlibStream);
							G4cout << G4endl << "#####" << G4endl << "# Error when trying to decompress the file in to memory! The data will PARTLY (i.e. in a wrong way) be considered!" << G4endl << "#####" << G4endl << G4endl;
							break;
					}

					// if inflate() could not totally fill tempBufferChar again (there were not enough characters to be read), cut off everything behind the newly filled characters
					if(zlibStream.avail_out > 0) tempBufferChar[processBufferSize - zlibStream.avail_out] = 0;

					// save the date from an "unsigned char" to a "char *" which is easier to work with
					char * tempBufferCharArray = (char *) tempBufferChar;

					// get the position to the last line break
					char * pointerToLastLineEnd = strrchr(tempBufferCharArray,'\n');

					// replace the last line break with a '\0'
					if(pointerToLastLineEnd) pointerToLastLineEnd[0] = '\0';

					// save save everything until the first '\0' to the data string (the potential beginning of an incomplete line will be cut)
					processedDataString += tempBufferCharArray;
					if(!pointerToLastLineEnd) continue;
					else (* processedDataString.end()) = '\n';   /* a string does not like to end on '\0' */

// processing of read in data
					const char * dataCharArray = processedDataString.c_str();
					const char * dataCharArrayEnd = dataCharArray + processedDataString.size();

					while(dataCharArray <= dataCharArrayEnd)
					{

						if(strncmp(dataCharArray, " ", 1) == 0 || strncmp(dataCharArray, "\0", 1) == 0)
						{
							dataCharArray++;
						}
						else
						{
							if( getIntValue(dataCharArray, "EventID:", NumEvent) )
							{
								createNewEvent();

								continue;
							}

							if(!Start_time_ns_set)
							{
								if( getDoubleValue(dataCharArray, "t/ns:", Start_time_ns) )
								{
									Start_time_ns_set = true;

									savePhotonData();

									continue;
								}
							}

							if(!Start_position_mm_set)
							{
								if( getVectorValue(dataCharArray, "pos/mm:", Start_position_mm) )
								{
									Start_position_mm_set = true;

									savePhotonData();

									continue;
								}
							}

							if(!Start_momentum_MeV_set)
							{
								if( getVectorValue(dataCharArray, "momentum/MeV:", Start_momentum_MeV) )
								{
									Start_momentum_MeV_set = true;

									savePhotonData();

									continue;
								}
							}

							if(!Start_polarisation_set)
							{
								if( getVectorValue(dataCharArray, "pol:", Start_polarisation) )
								{
									Start_polarisation_set = true;

									savePhotonData();

									continue;
								}
							}

							relocateCharArrayPointer(dataCharArray, dataCharArrayEnd, "\n");
						}
					}

// second part of the function to read in the file and decompress it
					// save last (incomplete) line to data string (it will be processe in the next run)
					processedDataString = pointerToLastLineEnd + 1;

				// done when inflate() could not totally fill tempBufferChar again
				}while(zlibStream.avail_out == 0);

			// done when inflate() says it's done
			}while(returnValue != Z_STREAM_END);

			// clean up and return
			(void) inflateEnd(&zlibStream);

			if(returnValue != Z_STREAM_END) G4cout << G4endl << "#####" << G4endl << "# Error when trying to decompress the file in to memory! The data will PARTLY (i.e. in a wrong way) be considered!" << G4endl << "#####" << G4endl << G4endl;
		}
	}

	fclose(infile);
}



void PhotonListSourceMessenger::createNewEvent()
{
	std::vector< std::vector<G4ThreeVector> > tempVector;

	ParticleInformations.push_back(tempVector);
}



void PhotonListSourceMessenger::savePhotonData()
{
	if(Start_time_ns_set && Start_position_mm_set && Start_momentum_MeV_set && Start_polarisation_set && ParticleInformations.size())
	{
		std::vector<G4ThreeVector> tempVector;
		tempVector.push_back(G4ThreeVector(Start_time_ns, NAN, NAN));
		tempVector.push_back(Start_position_mm);
		tempVector.push_back(Start_momentum_MeV);
		tempVector.push_back(Start_polarisation);

		(ParticleInformations[ParticleInformations.size() - 1]).insert((ParticleInformations[ParticleInformations.size() - 1]).begin(), tempVector);

		setDefaultValues();
	}
}



G4bool PhotonListSourceMessenger::getIntValue(const char * & data, G4String keyword, G4int & output)
{
	if(strncmp(data, keyword.c_str(), strlen(keyword.c_str())) == 0)
	{
		G4double rfloat = NAN;
		G4int matchedStringlength = 0;

		sscanf(data, "%*s %lf%n", &rfloat, &matchedStringlength);
		output = (G4int) rfloat;

		data += matchedStringlength + 1;

		return true;
	}

	return false;
}



G4bool PhotonListSourceMessenger::getDoubleValue(const char * & data, G4String keyword, G4double & output)
{
	if(strncmp(data, keyword.c_str(), strlen(keyword.c_str())) == 0)
	{
		G4double rfloat = NAN;
		G4int matchedStringlength = 0;

		sscanf(data, "%*s %lf%n", &rfloat, &matchedStringlength);
		output = rfloat;

		data += matchedStringlength + 1;

		return true;
	}

	return false;
}



G4bool PhotonListSourceMessenger::getVectorValue(const char * & data, G4String keyword, G4ThreeVector & output)
{
	if(strncmp(data, keyword.c_str(), strlen(keyword.c_str())) == 0)
	{
		G4double vx = NAN;
		G4double vy = NAN;
		G4double vz = NAN;
		G4int matchedStringlength = 0;

		sscanf(data, "%*s (%lf,%lf,%lf)%n", &vx, &vy, &vz, &matchedStringlength);

		if(!isnan(vx) && !isnan(vy) && !isnan(vz))
		{
			G4ThreeVector vec = G4ThreeVector(vx, vy, vz);
			output = vec;
		}

		data += matchedStringlength + 1;

		return true;
	}

	return false;
}



void PhotonListSourceMessenger::relocateCharArrayPointer(const char * & charArray, const char * charArrayEnd, G4String relocationDestination)
{
	G4int relocationDestination_length = strlen(relocationDestination.c_str());
	G4bool stopSkipping = false;

	while(charArray <= charArrayEnd)
	{
		if(strncmp(charArray, relocationDestination.c_str(), relocationDestination_length) == 0) stopSkipping = true;
		charArray++;
		if(stopSkipping) break;
	}
}



void PhotonListSourceMessenger::setDefaultValues()
{
	Start_time_ns = NAN;
	Start_time_ns_set = false;
	Start_position_mm = G4ThreeVector(NAN, NAN, NAN);
	Start_position_mm_set = false;
	Start_momentum_MeV = G4ThreeVector(NAN, NAN, NAN);
	Start_momentum_MeV_set = false;
	Start_polarisation = G4ThreeVector(NAN, NAN, NAN);
	Start_polarisation_set = false;
}
