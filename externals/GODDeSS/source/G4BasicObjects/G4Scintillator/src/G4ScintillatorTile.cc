/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4PVPlacement.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4VisAttributes.hh>

#include <G4ScintillatorTile.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to construct all volumes and surfaces for scintillator tiles.
 */
void G4ScintillatorTile::ConstructScintiVolume()
{
// solids (dimensions):

	G4double twinTileAirGap = 0.;
	G4double twinTileReflectorWidth = 0.;
	if(TwinTile)
	{
		if(TwinTileProperties.containsNumber("thick_air")) twinTileAirGap = TwinTileProperties.getNumber("thick_air");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::ConstructScintiVolume():" <<
			std::endl << "# The scintillator twin tile's air gap thickness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(TwinTileProperties.containsNumber("thick_reflector")) twinTileReflectorWidth = TwinTileProperties.getNumber("thick_reflector");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::ConstructScintiVolume():" <<
			std::endl << "# The scintillator twin tile's reflector thickness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Dimensions += G4ThreeVector(0., Dimensions.y() + twinTileReflectorWidth + 2. * twinTileAirGap, 0.);
	}

	G4String scintillator_solid_name = ScintillatorName + "_solid";
	G4VSolid * scintillator_solid = new G4Box(scintillator_solid_name, Dimensions.x() / 2., Dimensions.y() / 2., Dimensions.z() / 2.);

	if(CutDimple)
	{
		G4VSolid * dimple_solid = new G4Sphere("dimple_solid", 0., DimpleRadius, 0. * CLHEP::deg, 360. * CLHEP::deg, 0. * CLHEP::deg,  180. * CLHEP::deg);
		scintillator_solid = new G4SubtractionSolid(scintillator_solid_name, scintillator_solid, dimple_solid, G4Transform3D(G4RotationMatrix(), DimpleCentrePosition));
	}

	G4Box * twinTileReflector_solid = 0;
	G4Box * twinTileAirGap_solid = 0;
	if(TwinTile)
	{
		twinTileReflector_solid = new G4Box("twinTileReflector_solid", Dimensions.getX() / 2., twinTileReflectorWidth / 2., Dimensions.getZ() / 2.);
		twinTileAirGap_solid = new G4Box("twinTileAirGap_solid", Dimensions.getX() / 2., twinTileReflectorWidth / 2. + twinTileAirGap, Dimensions.getZ() / 2.);
	}

// logical volumes (material):

	G4String scintillator_logical_name = ScintillatorName + "_logical";
	G4LogicalVolume * scintillator_logical = new G4LogicalVolume(scintillator_solid, Material_Scinti, scintillator_logical_name, 0, 0, 0);

	G4VisAttributes ScintillatorVisAtt(G4Colour::Black());
	ScintillatorVisAtt.SetForceWireframe(true);
	ScintillatorVisAtt.SetLineWidth(3.);
	scintillator_logical->SetVisAttributes(ScintillatorVisAtt);

	G4LogicalVolume * twinTileReflector_logical = 0;
	G4LogicalVolume * twinTileAirGap_logical = 0;
	if(TwinTile)
	{
		twinTileReflector_logical = new G4LogicalVolume(twinTileReflector_solid, Material_TwinTile, "twinTileReflector_logical", 0, 0, 0);
		twinTileReflector_logical->SetVisAttributes(G4Colour::Grey());

		twinTileAirGap_logical = new G4LogicalVolume(twinTileAirGap_solid, Material_Air, "twinTileAirGap_logical", 0, 0, 0);
		twinTileAirGap_logical->SetVisAttributes(G4Colour::Blue());
	}

	if(ConstructSensitiveDetector)
	{
		// try to find an already existing senitive detector
		ScintillatorSensitiveDetector * sensitiveDetector = (ScintillatorSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("ScintillatorSD", false);
		if(!sensitiveDetector)
		{
			// create a new senitive detector
			sensitiveDetector = new ScintillatorSensitiveDetector("ScintillatorSD", DataStorage);
			// register it to Geant4
			G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDetector);
		}

		// assign it to detector parts
		scintillator_logical->SetSensitiveDetector(sensitiveDetector);
	}

// physical volumes (placement):
//NOTE: G4PVPlacement puts the volume to be placed into EVERY physical volume emanating from the same logical volume (no matter whether the logical or physical volume is specified as mother volume)!
//      => If G4PVPlacement is to be able to distinguish physical volumes, for each physical volume a separate logical volume has to be created.
//      => rule of thumb: For each volume that might become a mother volume, a separate logical volume should to be created.

	Scintillator_physical = new G4PVPlacement(Transformation, ScintillatorName, scintillator_logical, MotherVolume_physical, false, 0, SearchOverlaps);

	if(TwinTile)
	{
		TwinTileAirGap_physical = new G4PVPlacement(G4Transform3D(), "twinTileAirGap", twinTileAirGap_logical, Scintillator_physical, false, 0, SearchOverlaps);
		TwinTileReflector_physical = new G4PVPlacement(G4Transform3D(), "twinTileReflector", twinTileReflector_logical, TwinTileAirGap_physical, false, 0, SearchOverlaps);
	}

// ######  surfaces:
//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.
	ConstructScintiSurface();
}



/**
 *  Function to define materials (compounds, alloys) and their properties:
 *  - creates new materials according to the specifications from the property file(s)
 *  - if the same materials already exist, the newly created ones are deleted
 */
void G4ScintillatorTile::DefineMaterials()
{
	// NOTE The materials which have already been define (e.g. in the DetectorConstruction) are used here

//---------- scintillator ----------//
	G4double density_Scinti;
	if(ScintiProperties.containsNumber("density")) density_Scinti = ScintiProperties.getNumber("density");
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineMaterials():" <<
		std::endl << "# The scintillator's density has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	GoddessProperties::tabular chemicalComponents_Scinti;
	if(ScintiProperties.containsTabular("chemical_components")) chemicalComponents_Scinti = ScintiProperties.getTabular("chemical_components");
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineMaterials():" <<
		std::endl << "# The scintillator's chemical components have not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4double BirksConstant_Scinti = NAN;
	if(ScintiProperties.containsNumber("birks_constant")) BirksConstant_Scinti = ScintiProperties.getNumber("birks_constant");

	Material_Scinti = new G4Material("tempName", density_Scinti, (int) chemicalComponents_Scinti["element"].size());

	PropertyTools->AddElementsFromTable(Material_Scinti, chemicalComponents_Scinti);
	if(!isnan(BirksConstant_Scinti)) Material_Scinti->GetIonisation()->SetBirksConstant(BirksConstant_Scinti);

	DefineScintillatorMaterialProperties();

	PropertyTools->checkIfMaterialAlreadyExists("ScintillatorMaterial", Material_Scinti);

//---------- internal reflector ----------//
	if(TwinTile)
	{
		G4double density_twinTile;
		if(TwinTileProperties.containsNumber("density")) density_twinTile = TwinTileProperties.getNumber("density");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineMaterials():" <<
			std::endl << "# The scintillator twin tile's reflector density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_twinTile;
		if(TwinTileProperties.containsTabular("chemical_components")) chemicalComponents_twinTile = TwinTileProperties.getTabular("chemical_components");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineMaterials():" <<
			std::endl << "# The scintillator twin tile's reflector chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_TwinTile = new G4Material("tempName", density_twinTile, (int) chemicalComponents_twinTile["element"].size());
		PropertyTools->AddElementsFromTable(Material_TwinTile, chemicalComponents_twinTile);

		DefineTwinTileMaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("TwinTileMaterial", Material_TwinTile);


		Material_Air = new G4Material("tempName", 1.293 * CLHEP::kg / CLHEP::m3, 2);
		Material_Air->AddElement(G4Element::GetElement("Nitrogen"), 70 * CLHEP::perCent);
		Material_Air->AddElement(G4Element::GetElement("Oxygen"), 30 * CLHEP::perCent);

		DefineAirMaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("Air", Material_Air);
	}
}



//NOTE: possible Properties:
//         "RINDEX":			(spectrum (in dependence of the photon energy))		(obligatory property!)
//		defines the refraction index of the material, used for boundary processes, Cerenkov radiation and Rayleigh scattering
//         "ABSLENGTH":			(spectrum (in dependence of the photon energy))
//		defines the absorption length (absorption spectrum) of the material, used for the "normal" absorption of optical photons (default is infinity, i.e. no absorption)
//		(the absorption length for the WLS process of WLS materials is specified by "WLSABSLENGTH", "ABSLENGTH" can be specified additionally to simulate a non-WLS absorption fraction)
//         "RAYLEIGH":			(spectrum (in dependence of the photon energy))
//		defines the absorption length of the material, used for the rayleigh scattering of optical photons (default is infinity, i.e. no scattering)
//
//         "SCINTILLATIONYIELD":	(constant value (energy independent))			(obligatory property for scintillator materials!)
//		defines the mean number of photons, emitted per MeV energy deposition in the scintillator material (the real number is Poisson/Gauss distributed)
//		(can also be specified separately for different particles by putting "ELECTRON...", "PROTON...", "DEUTERON...", "TRITON...", "ALPHA...", "ION..." infront of "SCINTILLATIONYIELD")
//		(default is 0, i.e. no scintillation process)
//         "RESOLUTIONSCALE":		(constant value (energy independent))
//		defines the intrinsic resolution of the scintillator material, used for the statistical distribution of the number of generated photons in the scintillation process
//		(values > 1 result in a wider distribution, values < 1 result in a narrower distribution -> 1 is to be chosen as default)
//		(default is 0)
//         "FASTCOMPONENT":		(spectrum (in dependence of the photon energy))		(at least one "...COMPONENT" is obligatory for scintillator materials!)
//		defines the emission spectrum of the material, used for the fast scintillation process	NOTE: emission spectra are NOT linearly extrapolated between two given points!
//         "SLOWCOMPONENT":		(spectrum (in dependence of the photon energy))		(at least one "...COMPONENT" is obligatory for scintillator materials!)
//		defines the emission spectrum of the material, used for the slow scintillation process	NOTE: emission spectra are NOT linearly extrapolated between two given points!
//         "FASTTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between energy deposition and photon emission), used for the fast scintillation process (default is 0)
//         "SLOWTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between energy deposition and photon emission), used for the slow scintillation process (default is 0)
//         "FASTSCINTILLATIONRISETIME":	(constant value (energy independent))
//		defines the rise time (time between the start of the emission and the emission peak), used for the fast scintillation process (default is 0)
//         "SLOWSCINTILLATIONRISETIME":	(constant value (energy independent))
//		defines the rise time (time between the start of the emission and the emission peak), used for the slow scintillation process (default is 0)
//         "YIELDRATIO":		(constant value (energy independent))			(obligatory property for scintillator materials, if both "...COMPONENT"s are specified!)
//		defines relative strength of the fast scintillation process as a fraction of total scintillation yield (default is 0)
//
//         "WLSABSLENGTH":		(spectrum (in dependence of the photon energy))		(obligatory property for WLS materials!)
//		defines the absorption length (absorption spectrum) of the material, used for the WLS process (default is infinity, i.e. no WLS process)
//         "WLSCOMPONENT":		(spectrum (in dependence of the photon energy))		(obligatory property for WLS materials!)
//		defines the emission spectrum of the material, used for the WLS process	NOTE: emission spectra are NOT linearly extrapolated between two given points!
//         "WLSTIMECONSTANT":		(constant value (energy independent))
//		defines the decay time (time between absorption and emission), used for the WLS process (default is 0)
//         "WLSMEANNUMBERPHOTONS":	(constant value (energy independent))
//		defines the mean number of photons, emitted for each photon that was absorbed by the WLS material
//		(if specified, the real number of emitted photons is Poisson distributed, else the real number of emitted photons is 1)

/**
 *  Function to define the following material properties of the scintillator (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 *  - emmission spectrum
 *  - decay time
 *  - rise time
 *  - scintillation yield
 *  - ratio between the fast and the slow scintillation process
 *  - intrinsic resolution factor
 */
void G4ScintillatorTile::DefineScintillatorMaterialProperties()
{
	G4double scintillationyield = 0.;
	if(ScintiProperties.containsNumber("lightOutput_rel") && ScintiProperties.containsNumber("lightOutput"))
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# WARNING:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# An absolut AND a relative scintillation yield are given. The absolute one will be used." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		scintillationyield = ScintiProperties.getNumber("lightOutput");
	}
	else if(ScintiProperties.containsNumber("lightOutput")) scintillationyield = ScintiProperties.getNumber("lightOutput");
	else if(ScintiProperties.containsNumber("lightOutput_rel")) scintillationyield = ScintiProperties.getNumber("lightOutput_rel") * 20000.;
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# Neither an absolut nor a relative scintillation yield is given." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4double intrinsicResolutionFactor = -1;
	if(ScintiProperties.containsNumber("res_intr")) intrinsicResolutionFactor = ScintiProperties.getNumber("res_intr");
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# WARNING:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# No intrinsic resolution factor of the scintillator is given. It will be set to 1." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		intrinsicResolutionFactor = 1.;
	}

	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Scinti = PropertyTools->GetPropertyDistribution(ScintiProperties, "n_ref");
	if(!refractiveIndex_Scinti)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# The scintillator's refraction index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_Scinti = PropertyTools->GetPropertyDistribution(ScintiProperties, "mu_att");
	if( !attenuationLength_Scinti || PropertyTools->isConstantProperty(attenuationLength_Scinti) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_Scinti = PropertyTools->GetPropertyDistribution(ScintiProperties, "n_ref_imag");

		if(refractiveIndex_imag_Scinti)
		{
			if( !attenuationLength_Scinti || (PropertyTools->isConstantProperty(attenuationLength_Scinti) && !PropertyTools->isConstantProperty(refractiveIndex_imag_Scinti)) )
			{
				if(attenuationLength_Scinti) delete attenuationLength_Scinti;
				attenuationLength_Scinti = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_Scinti->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_Scinti->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_Scinti->Value(energy) );

					attenuationLength_Scinti->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_Scinti)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# The scintillator's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}


	G4MaterialPropertyVector * fastLightOutput_Scinti = 0;
	if(ScintiProperties.containsTabular("epsilon_fast") && ScintiProperties.containsTabular("epsilon"))
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# An emission spectrum AND a fast emission spectrum of the scintillator are given. Which one is to be used?" <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}
	else
	{
		fastLightOutput_Scinti = PropertyTools->GetPropertyDistribution(ScintiProperties, "epsilon");
		if(!fastLightOutput_Scinti) fastLightOutput_Scinti = PropertyTools->GetPropertyDistribution(ScintiProperties, "epsilon_fast");
	}

	G4MaterialPropertyVector * slowLightOutput_Scinti = PropertyTools->GetPropertyDistribution(ScintiProperties, "epsilon_slow");

	if(!fastLightOutput_Scinti && !slowLightOutput_Scinti)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
		std::endl << "# The scintillator's emission spectrum has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4double fastScintiDecayTime;
	if(fastLightOutput_Scinti)
	{
		if(ScintiProperties.containsNumber("t_decay_fast") && ScintiProperties.containsNumber("t_decay"))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# A decay time AND a fast decay time of the scintillator are given. Which one is to be used?" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
		else if(ScintiProperties.containsNumber("t_decay")) fastScintiDecayTime = ScintiProperties.getNumber("t_decay");
		else if(ScintiProperties.containsNumber("t_decay_fast")) fastScintiDecayTime = ScintiProperties.getNumber("t_decay_fast");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# The scintillator's (fast) decay time has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	G4double slowScintiDecayTime;
	if(slowLightOutput_Scinti)
	{
		if(ScintiProperties.containsNumber("t_decay_slow")) slowScintiDecayTime = ScintiProperties.getNumber("t_decay_slow");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# The scintillator's slow decay time has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	G4double fastScintiRiseTime;
	if(fastLightOutput_Scinti)
	{
		if(ScintiProperties.containsNumber("t_rise_fast") && ScintiProperties.containsNumber("t_rise"))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# A rise time AND a fast rise time of the scintillator are given. Which one is to be used?" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
		else if(ScintiProperties.containsNumber("t_rise")) fastScintiRiseTime = ScintiProperties.getNumber("t_rise");
		else if(ScintiProperties.containsNumber("t_rise_fast")) fastScintiRiseTime = ScintiProperties.getNumber("t_rise_fast");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# The scintillator's (fast) rise time has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	G4double slowScintiRiseTime;
	if(slowLightOutput_Scinti)
	{
		if(ScintiProperties.containsNumber("t_rise_slow")) slowScintiRiseTime = ScintiProperties.getNumber("t_rise_slow");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# The scintillator's slow rise time has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	G4double fastScintiDecayFraction;
	if(fastLightOutput_Scinti && slowLightOutput_Scinti)
	{
		if(ScintiProperties.containsNumber("fraction_fast")) fastScintiDecayFraction = ScintiProperties.getNumber("fraction_fast");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::DefineScintillatorMaterialProperties():" <<
			std::endl << "# The scintillator's fraction of fast scintillation decay has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}


	G4MaterialPropertiesTable * MPT_Scinti = new G4MaterialPropertiesTable();
	MPT_Scinti->AddProperty("RINDEX", refractiveIndex_Scinti);
	MPT_Scinti->AddProperty("ABSLENGTH", attenuationLength_Scinti);
	if(fastLightOutput_Scinti) MPT_Scinti->AddProperty("FASTCOMPONENT", fastLightOutput_Scinti);
	if(slowLightOutput_Scinti) MPT_Scinti->AddProperty("SLOWCOMPONENT", slowLightOutput_Scinti);
	MPT_Scinti->AddConstProperty("SCINTILLATIONYIELD", scintillationyield);
	MPT_Scinti->AddConstProperty("RESOLUTIONSCALE", intrinsicResolutionFactor);
	if(fastLightOutput_Scinti && fastScintiDecayTime) MPT_Scinti->AddConstProperty("FASTTIMECONSTANT", fastScintiDecayTime);
	if(slowLightOutput_Scinti && slowScintiDecayTime) MPT_Scinti->AddConstProperty("SLOWTIMECONSTANT", slowScintiDecayTime);
	if(fastLightOutput_Scinti && fastScintiRiseTime) MPT_Scinti->AddConstProperty("FASTSCINTILLATIONRISETIME", fastScintiRiseTime);
	if(slowLightOutput_Scinti && slowScintiRiseTime) MPT_Scinti->AddConstProperty("SLOWSCINTILLATIONRISETIME", slowScintiRiseTime);
	if(fastLightOutput_Scinti && slowLightOutput_Scinti) MPT_Scinti->AddConstProperty("YIELDRATIO", fastScintiDecayFraction);
	Material_Scinti->SetMaterialPropertiesTable(MPT_Scinti);
}

/**
 *  Function to define the following material properties of the twin tile's reflective layer (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4ScintillatorTile::DefineTwinTileMaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "n_ref");
	if(!refractiveIndex_twinTile) refractiveIndex_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "n_ref_real");
	if(!refractiveIndex_twinTile)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineTwinTileMaterialProperties():" <<
		std::endl << "# The scintillator twin tile's reflector refraction index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "mu_att");   // getting the given attenuation length spectrum
	if( !attenuationLength_twinTile || PropertyTools->isConstantProperty(attenuationLength_twinTile) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "n_ref_imag");

		if(refractiveIndex_imag_twinTile)
		{
			if( !attenuationLength_twinTile || (PropertyTools->isConstantProperty(attenuationLength_twinTile) && !PropertyTools->isConstantProperty(refractiveIndex_imag_twinTile)) )
			{
				if(attenuationLength_twinTile) delete attenuationLength_twinTile;
				attenuationLength_twinTile = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_twinTile->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_twinTile->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_twinTile->Value(energy) );

					attenuationLength_twinTile->InsertValues(energy, attLength);
				}
			}
		}
	}
	if(!attenuationLength_twinTile)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::DefineTwinTileMaterialProperties():" <<
		std::endl << "# The scintillator twin tile's reflector attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * MPT_TwinTile = new G4MaterialPropertiesTable();
	MPT_TwinTile->AddProperty("RINDEX", refractiveIndex_twinTile);
	if(attenuationLength_twinTile) MPT_TwinTile->AddProperty("ABSLENGTH", attenuationLength_twinTile);
	Material_TwinTile->SetMaterialPropertiesTable(MPT_TwinTile);
}

/**
 *  Function to define the following material properties of the air gap:
 *  - refractive index spectrum
 */
void G4ScintillatorTile::DefineAirMaterialProperties()
{
	// by now, GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Air = PropertyTools->GetPropertyDistribution(1.);
	// the refractive index of air can approximated by 1 as it only varies from 1.000308 (230nm) to 1.00027417 (1000nm)
	// http://refractiveindex.info/?group=GASES&material=Air

	// absorption and rayleigh scattering can be neglected, as only relatively thin layers are used in this simulation

	G4MaterialPropertiesTable * mpt_Air = new G4MaterialPropertiesTable();
	mpt_Air->AddProperty("RINDEX", refractiveIndex_Air);
	Material_Air->SetMaterialPropertiesTable(mpt_Air);
}



/**
 *  Function to construct surfaces between different volumes and define their mechanical and optical properties (according to the specifications from the property file(s)). \n
 *  The following properties are defined:
 *  - surface type ("dielectric_metal" or "dielectric_dielectric")
 *  - surface roughness
 *  - reflectivity (DefineSurfacesProperties())
 */
void G4ScintillatorTile::ConstructScintiSurface()
{
//NOTE: !!!!!! In order to understand simulation of optical surface properties (and esspecially the relations between the different options and parameters),
//             you should definitely start with reading the diagram to be found at:
//             http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2012/05/23/23.20-70533-ces_in_geant4_revised.png
//             (in http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/445.html)
//             The following comments are only a excerpt and CANNOT compete with the precision of this diagram !!!!!!
//
//      possible reflection models:
//         glisur:  original and obsolete GEANT3 model
//         unified: unified model (provides a range of different reflection mechanisms, cf. "possible Properties" below)
//         LUT:     using the Look-Up-Table for the surface simulation
//
//      possible surface types:
//         dielectric_dielectric: if both materials are dielectric, i.e. the reflection and refraction probabilities are defined by the refractive indices via the Fresnel equations
//                                (the reflection probability CANNOT be defined by the material's reflectivity, cf. "REFLECTIVITY" below)
//         dielectric_metal:      if one material is dielectric and the other is metal-like (i.e. the photons can only be reflected or absorbed but not refracted)
//                                this surface type has to be used if the reflection probability is to be defined by the material's reflectivity
//                                (this does not work for dielectric_dielectric, cf. "REFLECTIVITY" below)
//         dielectric_LUT:        if Look-Up-Tables is to be used
//         firsov:                for Firsov Process (O.B. Firsov, “Reflection of fast ions from a dense medium at glancing angles”, Sov. Phys.-Docklady, vol. 11, No. 8, pp.732-733, 1967)
//         x_ray:                 for x-ray mirror process
//
//      possible surface finishs: (cf. http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2012/05/23/23.20-70533-ces_in_geant4_revised.png !!!)
//         for dielectric_metal:
//            polished:             perfectly polished surface
//                                  -> specular spike reflection
//                                     (requiring the reflectivity of the metal)
//            ground:               rough surface
//                                  -> specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the reflectivity of the metal and the angle alpha between a micro-facet normal and the average surface normal)
//
//         for dielectric_dielectric:
//            polished:             perfectly polished surface
//                                  -> specular spike reflection
//                                     (requiring the refraction indices of both materials)
//            ground:               rough surface
//                                  -> specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of both materials and the angle alpha between a micro-facet normal and the average surface normal)
//            polishedfrontpainted: the volume has a perfectly polished surface and an absorbing paint without air gap
//                                  -> specular spike reflection
//                                     (requiring the reflectivity of the paint)
//            groundfrontpainted:   the volume has a rough surface and an absorbing paint without air gap
//                                  -> (diffuse) Lambertian reflection
//                                     (requiring the reflectivity of the paint)
//            polishedbackpainted:  the volume has a rough surface and a perfectly polished, absorbing paint with air gap
//                                  -> volume-air surface: specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of the volume material and the surface as well as the angle alpha between a micro-facet normal and the average surface normal)
//                                  -> air-paint surface: specular spike reflection
//                                     (requiring the reflectivity of the paint)
//            groundbackpainted:    the volume has a rough surface and a rough, absorbing paint with air gap
//                                  -> volume-air surface: specular spike, specular lobe reflection, (diffuse) Lambertian reflection or back scattering
//                                     (requiring the refraction indices of the volume material and the surface as well as the angle alpha between a micro-facet normal and the average surface normal)
//                                  -> air-paint surface: (diffuse) Lambertian reflection
//                                     (requiring the reflectivity of the paint)
//
//         for the Look-Up-Table model:
//            polishedlumirrorair:
//            polishedlumirrorglue:
//            polishedair:
//            polishedteflonair:
//            polishedtioair:
//            polishedtyvekair:
//            polishedvm2000air:
//            polishedvm2000glue:
//            etchedlumirrorair:
//            etchedlumirrorglue:
//            etchedair:
//            etchedteflonair:
//            etchedtioair:
//            etchedtyvekair:
//            etchedvm2000air:
//            etchedvm2000glue:
//            groundlumirrorair:
//            groundlumirrorglue:
//            groundair:
//            groundteflonair:
//            groundtioair:
//            groundtyvekair:
//            groundvm2000air:
//            groundvm2000glue:
//
//      sigma alpha: defines the roughness of the surface (sigma alpha specifies the standard deviation of the distribution of the micro-facet normals in [rad])

	// create the surface and define its properties
	G4double scintiSigmaAlpha;
	if(ScintiProperties.containsNumber("roughness")) scintiSigmaAlpha = ScintiProperties.getNumber("roughness");
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4ScintillatorTile::ConstructScintiSurface():" <<
		std::endl << "# The scintillator's roughness has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	OptSurf_ScintiMother = new G4OpticalSurface("ScintiMotherSurface", SurfaceModel, SurfaceFinish, dielectric_dielectric, scintiSigmaAlpha);

	if(TwinTile)
	{
		if(TwinTileProperties.containsString("is_metal")) TwinTileReflectorIsMetal = (TwinTileProperties.getString("is_metal") == "true");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4ScintillatorTile::ConstructScintiSurface():" <<
			std::endl << "# It is not specified if the scintillator twin tile's reflector is metal-like. It will be set to false." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			TwinTileReflectorIsMetal = false;
		}

		G4double twinTileSigmaAlpha;
		if(TwinTileProperties.containsNumber("roughness")) twinTileSigmaAlpha = TwinTileProperties.getNumber("roughness");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4ScintillatorTile::ConstructScintiSurface():" <<
			std::endl << "# The scintillator twin tile's reflector roughness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(TwinTileReflectorIsMetal) OptSurf_reflective = new G4OpticalSurface("reflective surface", unified, ground, dielectric_metal, twinTileSigmaAlpha);
		else OptSurf_reflective = new G4OpticalSurface("reflective surface", unified, ground, dielectric_dielectric, twinTileSigmaAlpha);
	}


	// define the optical surface properties
	DefineSurfacesProperties();


//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.

	// specify between which volumes the previously defined surface is to be applied
	// the first volume is where the photon comes from, the second volume is where the photon heads for.
	new G4LogicalBorderSurface("Scintillator/Mother Surface", Scintillator_physical, MotherVolume_physical, OptSurf_ScintiMother);
	new G4LogicalBorderSurface("Mother/Scintillator Surface", MotherVolume_physical, Scintillator_physical, OptSurf_ScintiMother);

	if(TwinTile)
	{
		new G4LogicalBorderSurface("Scintillator/TileAirGap Surface", Scintillator_physical, TwinTileAirGap_physical, OptSurf_ScintiMother);
		new G4LogicalBorderSurface("TileAirGap/Scintillator Surface", TwinTileAirGap_physical, Scintillator_physical, OptSurf_ScintiMother);
		new G4LogicalSkinSurface("reflector Surface", TwinTileReflector_physical->GetLogicalVolume(), OptSurf_reflective);
	}
}

/**
 *  Function to define the optical properties of the reflective surfaces (according to the specifications from the property file(s)). \n
 *  The following properties are defined:
 *  - reflectivity
 *  - fraction of reflection at the average surface
 *  - fraction of reflection at a micro facet surface
 *  - fraction of back scattering
 *  - fraction of internal (diffuse) Lambertian reflection
 */
void G4ScintillatorTile::DefineSurfacesProperties()
{
//NOTE: possible Properties (for reflection model "unified"):
//         "RINDEX":				(spectrum (in dependence of the photon energy))			(obligatory property for surfaces with surface finish "...backpainted"!)
//		in case of surface finish "...backpainted", defines the refraction index of the thin material layer between the surface and the paint
//		(only used for boundary processes and only in case of surface finish "...backpainted")
//         "REFLECTIVITY":			(spectrum (in dependence of the photon energy))
//		does NOT IN ANY CASE define the reflectivity of the surface but defines (1 - absorption coefficient), which may make a difference:
//			for dielectric_metal:      it is really the reflectivity of the metal, i.e. the probability that the photon is reflected at all
//						   (every photon not absorbed by the metal will be reflected)
//			for dielectric_dielectric: it is NOT the reflectivity but defines the absorption coefficient (absorption coefficient = 1 - "reflectivity") of a dirty surface
//						   (every photon not absorbed by dirt will be reflected or refracted as normal for "dielectric_dielectric" surfaces and the reflectivity
//						   for this process can be specified via "TRANSMITTANCE" or is calculated via the Fresnel equations, cf. "TRANSMITTANCE")
//		(default is 1., i.e. nothing is absorbed)
//         "REALRINDEX":			(spectrum (in dependence of the photon energy))
//		defines the real part of the complex refraction index of the surface   //FIXME surface <-> material ???
//		(in case that "REFLECTIVITY" is not specified and if "REALRINDEX" and "IMAGINARYRINDEX" are both specified,
//		the refectivity is calculated from the complex refraction index via the Fresnel equations   //FIXME shouldn't it be two refraction indices???
//		-> therefore, the complex refraction index should not be specified for surface type "dielectric_dielectric", cf. "REFLECTIVITY"
//		   (if one wants to define the reflectivity of a "dielectric_dielectric" surface from the complex refraction index,
//		   one has to calculate the transmittance via the Fresnel equations and specify it with "TRANSMITTANCE"))
//         "IMAGINARYRINDEX":			(spectrum (in dependence of the photon energy))
//		defines the imaginary part of the complex refraction index of the surface   //FIXME surface <-> material ???
//		(in case that "REFLECTIVITY" is not specified and if "REALRINDEX" and "IMAGINARYRINDEX" are both specified,
//		the refectivity is calculated from the complex refraction index via the Fresnel equations   //FIXME shouldn't it be two refraction indices???
//		-> therefore, the complex refraction index should not be specified for surface type "dielectric_dielectric", cf. "REFLECTIVITY"
//		   (if one wants to define the reflectivity of a "dielectric_dielectric" surface from the complex refraction index,
//		   one has to calculate the transmittance via the Fresnel equations and specify it with "TRANSMITTANCE"))
//         "EFFICIENCY":			(spectrum (in dependence of the photon energy))
//		defines the detection efficiency of absorbed photons (default is 0)
//         "TRANSMITTANCE":			(spectrum (in dependence of the photon energy))
//		in case of "dielectric_dielectric" surfaces, defines the fraction of photons (reaching the surface despite of dirt (cf. "REFLECTIVITY"))
//		which are refracted by the surface instead of being reflected
//		(if "TRANSMITTANCE" is not specified, the transmittance is calculated from the (real) refraction indices of the two materials forming the surface via the Fresnel equations)
//		(only used for boundary processes and only in case of surface type "dielectric_dielectric")
//         "SPECULARSPIKECONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability for reflection at the average surface
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//         "SPECULARLOBECONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability for reflection at a micro facet surface, i.e. the direction is smeared around the direction of the specular spike reflection
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//         "BACKSCATTERCONSTANT":		(spectrum (in dependence of the photon energy))
//		defines the probability of back scattering, caused by several reflections within a deep grove
//		(only used in case of surface finish "ground" or "...backpainted", default is 0)
//	    !!! the probability of internal (diffuse) Lambertian reflection can not be specified directly but is defined via
//		    100% = "SPECULARSPIKECONSTANT" + "SPECULARLOBECONSTANT" + "BACKSCATTERCONSTANT" + Lambertian (-> default is 1) !!!

	// define the optical surface properties
	G4MaterialPropertyVector * specularspike_ScintiMother = PropertyTools->GetPropertyDistribution(ScintiProperties, "fraction_specularspike");
	if(!specularspike_ScintiMother)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The scintillator's fraction of specular-spike-reflection has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertyVector * specularlobe_ScintiMother = PropertyTools->GetPropertyDistribution(ScintiProperties, "fraction_specularlobe");
	if(!specularlobe_ScintiMother)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The scintillator's fraction of specular-lobe-reflection has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertyVector * backscattering_ScintiMother = PropertyTools->GetPropertyDistribution(ScintiProperties, "fraction_backscattering");
	if(!backscattering_ScintiMother)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The scintillator's fraction of backscattering-reflection has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * MPT_OptSurf_ScintiMother = new G4MaterialPropertiesTable();
// 	if(reflectivity_ScintiMother) MPT_OptSurf_ScintiMother->AddProperty("TRANSMITTANCE", reflectivity_ScintiMother); // surface reflectvity
	MPT_OptSurf_ScintiMother->AddProperty("SPECULARLOBECONSTANT", specularlobe_ScintiMother);
	MPT_OptSurf_ScintiMother->AddProperty("SPECULARSPIKECONSTANT", specularspike_ScintiMother);
	MPT_OptSurf_ScintiMother->AddProperty("BACKSCATTERCONSTANT", backscattering_ScintiMother);
	OptSurf_ScintiMother->SetMaterialPropertiesTable(MPT_OptSurf_ScintiMother);


	if(TwinTile)
	{
		G4MaterialPropertyVector * reflectivity_twinTile = 0;
		G4MaterialPropertyVector * refractiveIndex_real_twinTile = 0;
		G4MaterialPropertyVector * refractiveIndex_imag_twinTile = 0;

		if(TwinTileReflectorIsMetal)
		{
			// trying to get a given reflectivity spectrum
			reflectivity_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "reflec");
			if(!reflectivity_twinTile)
			{
				// using the complex reflective index to specify the reflectivity spectrum
				refractiveIndex_real_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "n_ref_real");
				refractiveIndex_imag_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "n_ref_imag");

				if(!refractiveIndex_real_twinTile || !refractiveIndex_imag_twinTile)
				{
					std::cerr <<
					std::endl << "##########" <<
					std::endl << "# ERROR:" <<
					std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
					std::endl << "# The scintillator twin tile's reflector reflectivity has not been specified." <<
					std::endl << "##########" <<
					std::endl << std::endl;

					CriticalErrorOccured = true;
				}
			}
		}

		G4MaterialPropertyVector * specularspike_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "fraction_specularspike");
		if(!specularspike_twinTile)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The scintillator twin tile's reflector fraction of specular-spike-reflection has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		G4MaterialPropertyVector * specularlobe_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "fraction_specularlobe");
		if(!specularlobe_twinTile)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The scintillator twin tile's reflector fraction of specular-lobe-reflection has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		G4MaterialPropertyVector * backscattering_twinTile = PropertyTools->GetPropertyDistribution(TwinTileProperties, "fraction_backscattering");
		if(!backscattering_twinTile)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The scintillator twin tile's reflector fraction of backscattering-reflection has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		G4MaterialPropertiesTable * MPT_OptSurf_reflective = new G4MaterialPropertiesTable();
		if(reflectivity_twinTile) MPT_OptSurf_reflective->AddProperty("REFLECTIVITY", reflectivity_twinTile);
		if(refractiveIndex_real_twinTile && refractiveIndex_imag_twinTile)
		{
			MPT_OptSurf_reflective->AddProperty("REALRINDEX", refractiveIndex_real_twinTile);
			MPT_OptSurf_reflective->AddProperty("IMAGINARYRINDEX", refractiveIndex_imag_twinTile);
		}
		MPT_OptSurf_reflective->AddProperty("SPECULARLOBECONSTANT", specularlobe_twinTile);
		MPT_OptSurf_reflective->AddProperty("SPECULARSPIKECONSTANT", specularspike_twinTile);
		MPT_OptSurf_reflective->AddProperty("BACKSCATTERCONSTANT", backscattering_twinTile);
		OptSurf_reflective->SetMaterialPropertiesTable(MPT_OptSurf_reflective);
	}
}



void G4ScintillatorTile::SetDefaults()
{
	CriticalErrorOccured = false;

	Material_Scinti = 0;
	Material_TwinTile = 0;
	Material_Air = 0;
	Transformation = G4Transform3D();
	Scintillator_physical = 0;
	TwinTileAirGap_physical = 0;
	TwinTileReflector_physical = 0;
	OptSurf_ScintiMother = 0;
	TwinTileReflectorIsMetal = false;
	OptSurf_reflective = 0;
}
