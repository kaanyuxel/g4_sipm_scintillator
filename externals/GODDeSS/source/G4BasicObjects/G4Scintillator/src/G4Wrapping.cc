/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4VisAttributes.hh>

#include <G4Wrapping.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to create a wrapping around a G4ScintillatorTile by using the surface properties and Look Up Tables.
 */
void G4Wrapping::ConstructLUTWrapping( G4ScintillatorTile* scinti_to_be_wrapped,
				       G4OpticalSurfaceModel surfaceModel,
				       G4OpticalSurfaceFinish surfaceFinish,
				       G4SurfaceType surfaceType,
				       G4double roughness,
				       G4VPhysicalVolume* surrounding_volume )
{
	G4OpticalSurface * OptSurf_MotherWrapping = new G4OpticalSurface("Wrapping1MotherSurface", surfaceModel, surfaceFinish, surfaceType, roughness);
	G4VPhysicalVolume * volume_to_be_wrapped = scinti_to_be_wrapped->GetScintillator_physicalVolume();

	if(surrounding_volume)
	{
		G4LogicalBorderSurface * OptBorderSurf_ScintiMother = G4LogicalBorderSurface::GetSurface(volume_to_be_wrapped, surrounding_volume);

		if(OptBorderSurf_ScintiMother)
		{
			OptBorderSurf_ScintiMother->SetSurfaceProperty(OptSurf_MotherWrapping);
			OptBorderSurf_ScintiMother->SetName("Scintillator/LUTWrapping/Mother");
		}
		else
		{
			new G4LogicalBorderSurface("Scintillator/LUTWrapping/Mother", volume_to_be_wrapped, surrounding_volume, OptSurf_MotherWrapping);
		}


		G4LogicalBorderSurface * OptBorderSurf_MotherScinti = G4LogicalBorderSurface::GetSurface(surrounding_volume, volume_to_be_wrapped);

		if(OptBorderSurf_MotherScinti)
		{
			OptBorderSurf_MotherScinti->SetSurfaceProperty(OptSurf_MotherWrapping);
			OptBorderSurf_MotherScinti->SetName("Mother/LUTWrapping/Scintillator");
		}
		else
		{
			new G4LogicalBorderSurface("Mother/LUTWrapping/Scintillator", surrounding_volume, volume_to_be_wrapped, OptSurf_MotherWrapping);
		}
	}
	else
	{
		new G4LogicalSkinSurface("LUTWrapping", volume_to_be_wrapped->GetLogicalVolume(), OptSurf_MotherWrapping);
	}
}



/**
 *  Function to construct all volumes and surfaces for wrapping volumes.
 */
void G4Wrapping::ConstructWrappingVolumes()
{
// solids (dimensions):

	G4String wrapping1_solid_name;
	G4String wrapping2_solid_name;

	// air gap between scintillator and wrapping
	G4ThreeVector dimensions_airGap = Dimensions + G4ThreeVector(2 * AirGapThickness, 2 * AirGapThickness, 2 * AirGapThickness);
	G4VSolid * airGap_solid = 0;
	if(AirGapThickness) airGap_solid = new G4Box("airGap_solid", dimensions_airGap.x() / 2., dimensions_airGap.y() / 2., dimensions_airGap.z() / 2.);

	// first wrapping layer
	G4ThreeVector dimensions_wrapping1 = dimensions_airGap + G4ThreeVector(2 * Wrapping1Thickness, 2 * Wrapping1Thickness, 2 * Wrapping1Thickness);

	if(Wrapping2Exists) wrapping1_solid_name = WrappingName + "_firstLayer_solid";
	else wrapping1_solid_name = WrappingName + "_solid";
	G4VSolid * wrapping1_solid = new G4Box(wrapping1_solid_name, dimensions_wrapping1.x() / 2., dimensions_wrapping1.y() / 2., dimensions_wrapping1.z() / 2.);

	// second wrapping layer
	G4VSolid * wrapping2_solid;
	if(Wrapping2Exists)
	{
		G4ThreeVector dimensions_wrapping2 = dimensions_wrapping1 + G4ThreeVector(2 * Wrapping2Thickness, 2 * Wrapping2Thickness, 2 * Wrapping2Thickness);

		wrapping2_solid_name = WrappingName + "_secondLayer_solid";
		wrapping2_solid = new G4Box(wrapping2_solid_name, dimensions_wrapping2.x() / 2., dimensions_wrapping2.y() / 2., dimensions_wrapping2.z() / 2.);
	}


	if(ScintiToBeWrapped_physical->HasDimple())
	{
		G4double dimpleRadius = ScintiToBeWrapped_physical->GetDimpleRadius();
		G4ThreeVector dimpleTranslation = ScintiToBeWrapped_physical->GetDimpleCentrePosition();

		if(airGap_solid)
		{
			G4VSolid * dimpleAirGap_solid = new G4Sphere("dimpleAirGap_solid", 0., dimpleRadius - AirGapThickness, 0. * CLHEP::deg, 360. * CLHEP::deg, 0. * CLHEP::deg,  180. * CLHEP::deg);
			airGap_solid = new G4SubtractionSolid("airGap_solid", airGap_solid, dimpleAirGap_solid, G4Transform3D(G4RotationMatrix(), dimpleTranslation));
		}

		G4VSolid * dimpleWrapping1_solid = new G4Sphere("dimpleWrapping1_solid", 0., dimpleRadius - AirGapThickness - Wrapping1Thickness, 0. * CLHEP::deg, 360. * CLHEP::deg, 0. * CLHEP::deg,  180. * CLHEP::deg);
		wrapping1_solid = new G4SubtractionSolid(wrapping1_solid_name, wrapping1_solid, dimpleWrapping1_solid, G4Transform3D(G4RotationMatrix(), dimpleTranslation));

		if(Wrapping2Exists)
		{
			G4VSolid * dimpleWrapping2_solid = new G4Sphere("dimpleWrapping2_solid", 0., dimpleRadius - AirGapThickness - Wrapping1Thickness - Wrapping2Thickness, 0. * CLHEP::deg, 360. * CLHEP::deg, 0. * CLHEP::deg,  180. * CLHEP::deg);
			wrapping2_solid = new G4SubtractionSolid(wrapping2_solid_name, wrapping2_solid, dimpleWrapping2_solid, G4Transform3D(G4RotationMatrix(), dimpleTranslation));
		}
	}

	// cut air gap and scintillator
	if(airGap_solid) wrapping1_solid = new G4SubtractionSolid(wrapping1_solid_name, wrapping1_solid, airGap_solid, G4Transform3D());
	else wrapping1_solid = new G4SubtractionSolid(wrapping1_solid_name, wrapping1_solid, ScintiToBeWrapped_physical->GetScintillator_solidVolume(), G4Transform3D());

	if(Wrapping2Exists)
	{
		if(airGap_solid) wrapping2_solid = new G4SubtractionSolid(wrapping2_solid_name, wrapping2_solid, airGap_solid, G4Transform3D());
		else wrapping2_solid = new G4SubtractionSolid(wrapping2_solid_name, wrapping2_solid, ScintiToBeWrapped_physical->GetScintillator_solidVolume(), G4Transform3D());
	}

	// cut volumes specified in CutVolumes
	for(size_t i=0; i < CutVolumes.size(); i++)
	{
		if(CutVolumes[i]->GetMotherLogical() && CutVolumes[i] != ScintiToBeWrapped_physical->GetScintillator_physicalVolume())   // CutVolume != World && CutVolume != Scinti
		{
			if(CutVolumes[i]->GetMotherLogical() == MotherVolume_physical->GetLogicalVolume())   //both volumes are placed within the same mother volume - not perfect due to problems if one of the volumes is parameterised!!!
			{
				G4VSolid * cutVolume_solid = CutVolumes[i]->GetLogicalVolume()->GetSolid();
				G4RotationMatrix relCutRotation = Transformation.getRotation().inverse() * CutVolumes[i]->GetObjectRotationValue();

				G4ThreeVector relCutTranslation = CutVolumes[i]->GetObjectTranslation() - Transformation.getTranslation();
				relCutTranslation.transform(Transformation.getRotation().inverse());

				G4Transform3D FibreTransformation_outsideMother = G4Transform3D(relCutRotation, relCutTranslation);

				if(			   G4IntersectionSolid("", wrapping1_solid, cutVolume_solid, FibreTransformation_outsideMother).GetCubicVolume()
				    || (Wrapping2Exists && G4IntersectionSolid("", wrapping2_solid, cutVolume_solid, FibreTransformation_outsideMother).GetCubicVolume()) )
				{

					while(cutVolume_solid->GetConstituentSolid(0))
					{
						cutVolume_solid = cutVolume_solid->GetConstituentSolid(0);
					}
					cutVolume_solid = cutVolume_solid->Clone();

					if(cutVolume_solid->GetEntityType() == "G4Box")
					{
						((G4Box*) cutVolume_solid)->SetXHalfLength( ((G4Box*) cutVolume_solid)->GetXHalfLength() + AirGapThickness * 0.9 );
						((G4Box*) cutVolume_solid)->SetYHalfLength( ((G4Box*) cutVolume_solid)->GetYHalfLength() + AirGapThickness * 0.9 );
						((G4Box*) cutVolume_solid)->SetZHalfLength( ((G4Box*) cutVolume_solid)->GetZHalfLength() + AirGapThickness * 0.9 );
					}
					if(cutVolume_solid->GetEntityType() == "G4Tubs")
					{
						((G4Tubs*) cutVolume_solid)->SetOuterRadius( ((G4Tubs*) cutVolume_solid)->GetOuterRadius() + AirGapThickness * 0.9 );
						((G4Tubs*) cutVolume_solid)->SetInnerRadius( 0 );
// 						((G4Tubs*) cutVolume_solid)->SetZHalfLength( ((G4Tubs*) cutVolume_solid)->GetZHalfLength() + AirGapThickness * 0.9 );
					}

					if(fabs(G4IntersectionSolid("", wrapping1_solid, cutVolume_solid, FibreTransformation_outsideMother).GetCubicVolume()) >= 1e-12)
					{
						wrapping1_solid = new G4SubtractionSolid(wrapping1_solid_name, wrapping1_solid, cutVolume_solid, FibreTransformation_outsideMother);
					}
					if(Wrapping2Exists && fabs(G4IntersectionSolid("", wrapping2_solid, cutVolume_solid, FibreTransformation_outsideMother).GetCubicVolume()) >= 1e-12)
					{
						wrapping2_solid = new G4SubtractionSolid(wrapping2_solid_name, wrapping2_solid, cutVolume_solid, FibreTransformation_outsideMother);
					}
				}
			}
		}
	}


// logical volumes (material):

	G4String wrapping1_logical_name;
	G4String wrapping2_logical_name;

	// first wrapping layer
	if(Wrapping2Exists) wrapping1_logical_name = WrappingName + "_firstLayer_logical";
	else wrapping1_logical_name = WrappingName + "_logical";
	G4LogicalVolume * wrapping1_logical = new G4LogicalVolume(wrapping1_solid, Material_Wrapping1, wrapping1_logical_name, 0, 0, 0);

	G4VisAttributes wrapping1VisAtt;
	if(Wrapping2Exists) wrapping1VisAtt.SetColour(G4Colour::Green());
	else wrapping1VisAtt.SetColour(G4Colour::Red());
	wrapping1VisAtt.SetForceWireframe(true);
	wrapping1VisAtt.SetLineWidth(3.);
	wrapping1_logical->SetVisAttributes(wrapping1VisAtt);

	// second wrapping layer
	G4LogicalVolume * wrapping2_logical;
	if(Wrapping2Exists)
	{
		wrapping2_logical_name = WrappingName + "_secondLayer_logical";
		wrapping2_logical = new G4LogicalVolume(wrapping2_solid, Material_Wrapping2, wrapping2_logical_name, 0, 0, 0);

		G4VisAttributes wrapping2VisAtt(G4Colour::Red());
		wrapping2VisAtt.SetForceWireframe(true);
		wrapping2VisAtt.SetLineWidth(3.);
		wrapping2_logical->SetVisAttributes(wrapping2VisAtt);
	}

	if(ConstructSensitiveDetector)
	{
		// try to find an already existing senitive detector
		WrappingSensitiveDetector * sensitiveDetector = (WrappingSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("WrappingSD", false);
		if(!sensitiveDetector)
		{
			// create a new senitive detector
			sensitiveDetector = new WrappingSensitiveDetector("WrappingSD", DataStorage);
			// register it to Geant4
			G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDetector);
		}

		// assign it to detector parts
		wrapping1_logical->SetSensitiveDetector(sensitiveDetector);
		if(Wrapping2Exists) wrapping2_logical->SetSensitiveDetector(sensitiveDetector);
	}

// physical volumes (placement):
//NOTE: G4PVPlacement puts the volume to be placed into EVERY physical volume emanating from the same logical volume (no matter whether the logical or physical volume is specified as mother volume)!
//      => If G4PVPlacement is to be able to distinguish physical volumes, for each physical volume a separate logical volume has to be created.
//      => rule of thumb: For each volume that might become a mother volume, a separate logical volume should to be created.

	G4String wrapping2_physical_name;
	G4String wrapping1_physical_name;
	if(Wrapping2Exists)
	{
	// second wrapping layer
		wrapping2_physical_name = WrappingName + "_secondLayer";
		Wrapping2_physical = new G4PVPlacement(Transformation, wrapping2_physical_name, wrapping2_logical, MotherVolume_physical, false, 0, SearchOverlaps);

	// first wrapping layer
		wrapping1_physical_name = WrappingName + "_firstLayer";
		Wrapping1_physical = new G4PVPlacement(G4Transform3D(), wrapping1_physical_name, wrapping1_logical, Wrapping2_physical, false, 0, SearchOverlaps);
	}
	else
	{
	// first wrapping layer
		Wrapping1_physical = new G4PVPlacement(Transformation, WrappingName, wrapping1_logical, MotherVolume_physical, false, 0, SearchOverlaps);
	}

// surfaces:
//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.
	ConstructWrappingSurface();
}



/**
 *  Function to define materials (compounds, alloys) and their properties:
 *  - creates new materials according to the specifications from the property file(s)
 *  - if the same materials already exist, the newly created ones are deleted
 */
void G4Wrapping::DefineWrappingMaterials()
{
	// NOTE The materials which have already been define (e.g. in the DetectorConstruction) are used here

//---------- external reflector / wrapping ----------//

	// for the first wrapping layer
	G4double density_Wrapping1;
	if(WrappingProperties.containsNumber("density_layer_1")) density_Wrapping1 = WrappingProperties.getNumber("density_layer_1");
	else if(WrappingProperties.containsNumber("density")) density_Wrapping1 = WrappingProperties.getNumber("density");
	else
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::DefineWrappingMaterials():" <<
		std::endl << "# The" << counterString << " wrapping layer's density has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	GoddessProperties::tabular chemicalComponents_Wrapping1;
	if(WrappingProperties.containsTabular("chemical_components_layer_1")) chemicalComponents_Wrapping1 = WrappingProperties.getTabular("chemical_components_layer_1");
	else if(WrappingProperties.containsTabular("chemical_components")) chemicalComponents_Wrapping1 = WrappingProperties.getTabular("chemical_components");
	else
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::DefineWrappingMaterials():" <<
		std::endl << "# The" << counterString << " wrapping layer's chemical components have not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	Material_Wrapping1 = new G4Material("tempName", density_Wrapping1, (int) chemicalComponents_Wrapping1["element"].size());
	PropertyTools->AddElementsFromTable(Material_Wrapping1, chemicalComponents_Wrapping1);

	DefineWrapping1MaterialProperties();

	PropertyTools->checkIfMaterialAlreadyExists("ScintillatorWrapping1Material", Material_Wrapping1);

	// for the second wrapping layer
	if(Wrapping2Exists)
	{
		G4double density_Wrapping2;
		if(WrappingProperties.containsNumber("density_layer_2")) density_Wrapping2 = WrappingProperties.getNumber("density_layer_2");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::DefineWrappingMaterials():" <<
			std::endl << "# The second wrapping layer's density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_Wrapping2;
		if(WrappingProperties.containsTabular("chemical_components_layer_2")) chemicalComponents_Wrapping2 = WrappingProperties.getTabular("chemical_components_layer_2");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::DefineWrappingMaterials():" <<
			std::endl << "# The second wrapping layer's chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_Wrapping2 = new G4Material("tempName", density_Wrapping2, (int) chemicalComponents_Wrapping2["element"].size());
		PropertyTools->AddElementsFromTable(Material_Wrapping2, chemicalComponents_Wrapping2);

		DefineWrapping2MaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("ScintillatorWrapping2Material", Material_Wrapping2);
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
 *  Function to define the following material properties of the inner wrapping layer (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Wrapping::DefineWrapping1MaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_layer_1");
	if(!refractiveIndex_Wrapping1) refractiveIndex_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref");
	if(!refractiveIndex_Wrapping1) refractiveIndex_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real_layer_1");
	if(!refractiveIndex_Wrapping1) refractiveIndex_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real");
	if(!refractiveIndex_Wrapping1)
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::DefineWrapping1MaterialProperties():" <<
		std::endl << "# The" << counterString << " wrapping layer's refraction index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "mu_att_layer_1");
	if(!attenuationLength_Wrapping1) attenuationLength_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "mu_att");
	if( !attenuationLength_Wrapping1 || PropertyTools->isConstantProperty(attenuationLength_Wrapping1) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag_layer_1");
		if(!refractiveIndex_imag_Wrapping1) refractiveIndex_imag_Wrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag");

		if(refractiveIndex_imag_Wrapping1)
		{
			if( !attenuationLength_Wrapping1 || (PropertyTools->isConstantProperty(attenuationLength_Wrapping1) && !PropertyTools->isConstantProperty(refractiveIndex_imag_Wrapping1)) )
			{
				if(attenuationLength_Wrapping1) delete attenuationLength_Wrapping1;
				attenuationLength_Wrapping1 = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_Wrapping1->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_Wrapping1->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_Wrapping1->Value(energy) );

					attenuationLength_Wrapping1->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_Wrapping1)
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::DefineWrapping1MaterialProperties():" <<
		std::endl << "# The" << counterString << " wrapping layer's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * MPT_Wrapping1 = new G4MaterialPropertiesTable();
	MPT_Wrapping1->AddProperty("RINDEX", refractiveIndex_Wrapping1);
	if(attenuationLength_Wrapping1) MPT_Wrapping1->AddProperty("ABSLENGTH", attenuationLength_Wrapping1);
	Material_Wrapping1->SetMaterialPropertiesTable(MPT_Wrapping1);
}

/**
 *  Function to define the following material properties of the outer wrapping layer (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Wrapping::DefineWrapping2MaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_layer_2");
	if(!refractiveIndex_Wrapping2) refractiveIndex_Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real_layer_2");
	if(!refractiveIndex_Wrapping2)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::DefineWrapping2MaterialProperties():" <<
		std::endl << "# The second wrapping layer's refraction index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "mu_att_layer_2");
	if( !attenuationLength_Wrapping2 || PropertyTools->isConstantProperty(attenuationLength_Wrapping2) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag_layer_2");

		if(refractiveIndex_imag_Wrapping2)
		{
			if( !attenuationLength_Wrapping2 || (PropertyTools->isConstantProperty(attenuationLength_Wrapping2) && !PropertyTools->isConstantProperty(refractiveIndex_imag_Wrapping2)) )
			{
				if(attenuationLength_Wrapping2) delete attenuationLength_Wrapping2;
				attenuationLength_Wrapping2 = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_Wrapping2->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_Wrapping2->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_Wrapping2->Value(energy) );

					attenuationLength_Wrapping2->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_Wrapping2)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::DefineWrapping2MaterialProperties():" <<
		std::endl << "# The second wrapping layer's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * MPT_Wrapping2 = new G4MaterialPropertiesTable();
	MPT_Wrapping2->AddProperty("RINDEX", refractiveIndex_Wrapping2);
	if(!refractiveIndex_Wrapping2) MPT_Wrapping2->AddProperty("ABSLENGTH", attenuationLength_Wrapping2);
	Material_Wrapping2->SetMaterialPropertiesTable(MPT_Wrapping2);
}



/**
 *  Function to construct surfaces between different volumes and define their mechanical and optical properties (according to the specifications from the property file(s)). \n
 *  The following properties are defined:
 *  - surface type ("dielectric_metal" or "dielectric_dielectric")
 *  - surface roughness
 *  - reflectivity (DefineWrappingSurfacesProperties())
 */
void G4Wrapping::ConstructWrappingSurface()
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

	// create the surface an define its properties
	G4double wrapping1SigmaAlpha = 0;
	if(WrappingProperties.containsNumber("roughness_layer_1")) wrapping1SigmaAlpha = WrappingProperties.getNumber("roughness_layer_1");
	else if(WrappingProperties.containsNumber("roughness")) wrapping1SigmaAlpha = WrappingProperties.getNumber("roughness");
	else
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The" << counterString << " wrapping layer's roughness has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4double wrapping2SigmaAlpha = 0;
	if(Wrapping2Exists)
	{
		if(WrappingProperties.containsNumber("roughness_layer_2")) wrapping2SigmaAlpha = WrappingProperties.getNumber("roughness_layer_2");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The second wrapping layer's roughness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	if(Wrapping1IsMetal) OptSurf_MotherWrapping1 = new G4OpticalSurface("Wrapping1MotherSurface", SurfaceModel, SurfaceFinish, dielectric_metal, wrapping1SigmaAlpha);
	else OptSurf_MotherWrapping1 = new G4OpticalSurface("Wrapping1MotherSurface", SurfaceModel, SurfaceFinish, dielectric_dielectric, wrapping1SigmaAlpha);

	if(Wrapping2Exists)
	{
		if(Wrapping2IsMetal)
		{
			OptSurf_MotherWrapping2 = new G4OpticalSurface("Wrapping2MotherSurface", SurfaceModel, SurfaceFinish, dielectric_metal, wrapping2SigmaAlpha);
			if(!Wrapping1IsMetal) OptSurf_Wrapping1Wrapping2 = new G4OpticalSurface("Wrapping1Wrapping2Surface", SurfaceModel, SurfaceFinish, dielectric_metal, wrapping2SigmaAlpha);
		}
		else
		{
			OptSurf_MotherWrapping2 = new G4OpticalSurface("Wrapping2MotherSurface", SurfaceModel, SurfaceFinish, dielectric_dielectric, wrapping2SigmaAlpha);
			if(Wrapping1IsMetal) OptSurf_Wrapping1Wrapping2 = new G4OpticalSurface("Wrapping1Wrapping2Surface", SurfaceModel, SurfaceFinish, dielectric_metal, wrapping2SigmaAlpha);
			else OptSurf_Wrapping1Wrapping2 = new G4OpticalSurface("Wrapping1Wrapping2Surface", SurfaceModel, SurfaceFinish, dielectric_dielectric, wrapping2SigmaAlpha);
		}
	}

	// define the optical surface properties
	DefineWrappingSurfacesProperties();





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
	new G4LogicalSkinSurface("Wrapping1 Surface", Wrapping1_physical->GetLogicalVolume(), OptSurf_MotherWrapping1);
	if(!AirGapThickness)
	{
		G4LogicalBorderSurface * OptBorderSurf_ScintiWrapping = G4LogicalBorderSurface::GetSurface(ScintiToBeWrapped_physical->GetScintillator_physicalVolume(), MotherVolume_physical);

		if(OptBorderSurf_ScintiWrapping)
		{
			OptBorderSurf_ScintiWrapping->SetSurfaceProperty(OptSurf_MotherWrapping1);
			OptBorderSurf_ScintiWrapping->SetPhysicalVolumes(ScintiToBeWrapped_physical->GetScintillator_physicalVolume(), Wrapping1_physical);
			OptBorderSurf_ScintiWrapping->SetName("Scintillator/Wrapping1 Surface");
		}
		else
		{
			new G4LogicalBorderSurface("Scintillator/Wrapping1 Surface", ScintiToBeWrapped_physical->GetScintillator_physicalVolume(), Wrapping1_physical, OptSurf_MotherWrapping1);
		}

		G4LogicalBorderSurface * OptBorderSurf_WrappingScinti = G4LogicalBorderSurface::GetSurface(MotherVolume_physical, ScintiToBeWrapped_physical->GetScintillator_physicalVolume());

		if(OptBorderSurf_WrappingScinti)
		{
			OptBorderSurf_WrappingScinti->SetSurfaceProperty(OptSurf_MotherWrapping1);
			OptBorderSurf_WrappingScinti->SetPhysicalVolumes(Wrapping1_physical, ScintiToBeWrapped_physical->GetScintillator_physicalVolume());
			OptBorderSurf_WrappingScinti->SetName("Scintillator/Wrapping1 Surface");
		}
		else
		{
			new G4LogicalBorderSurface("Wrapping1/Scintillator Surface", Wrapping1_physical, ScintiToBeWrapped_physical->GetScintillator_physicalVolume(), OptSurf_MotherWrapping1);
		}
	}

	if(Wrapping2Exists)
	{
		new G4LogicalSkinSurface("Wrapping2 Surface", Wrapping2_physical->GetLogicalVolume(), OptSurf_MotherWrapping2);
		new G4LogicalBorderSurface("Wrapping1/Wrapping2 Surface", Wrapping1_physical, Wrapping2_physical, OptSurf_Wrapping1Wrapping2);
		new G4LogicalBorderSurface("Wrapping2/Wrapping1 Surface", Wrapping2_physical, Wrapping1_physical, OptSurf_Wrapping1Wrapping2);
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
void G4Wrapping::DefineWrappingSurfacesProperties()
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
	G4MaterialPropertyVector * reflectivity_MotherWrapping1 = 0;
	G4MaterialPropertyVector * refractiveIndex_real_MotherWrapping1 = 0;
	G4MaterialPropertyVector * refractiveIndex_imag_MotherWrapping1 = 0;
	G4String reflectivityString_MotherWrapping1 = "";

	if(Wrapping1IsMetal)
	{
		reflectivity_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_1");
		if(!reflectivity_MotherWrapping1) reflectivity_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec");
	}
	else
	{
		reflectivity_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_1", true);
		if(!reflectivity_MotherWrapping1) reflectivity_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec", true);
	}

	if(reflectivity_MotherWrapping1)
	{
		if(Wrapping1IsMetal) reflectivityString_MotherWrapping1 = "REFLECTIVITY";
		else reflectivityString_MotherWrapping1 = "TRANSMITTANCE";
	}
	else if(Wrapping1IsMetal)
	{
		// using the complex reflective index to specify the reflectivity spectrum
		refractiveIndex_real_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real_layer_1");
		refractiveIndex_imag_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag_layer_1");
		if(!refractiveIndex_real_MotherWrapping1 || !refractiveIndex_imag_MotherWrapping1)
		{
			refractiveIndex_real_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real");
			refractiveIndex_imag_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag");
		}

		if(!refractiveIndex_real_MotherWrapping1 || !refractiveIndex_imag_MotherWrapping1)
		{
			G4String counterString = "";
			if(Wrapping2Exists) counterString = " first";

			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The" << counterString << " wrapping layer's reflectivity has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	G4MaterialPropertyVector * specularspike_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularspike_layer_1");
	if(!specularspike_MotherWrapping1) specularspike_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularspike");
	if(!specularspike_MotherWrapping1)
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The" << counterString << " wrapping layer's fraction of specular-spike-reflection has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertyVector * specularlobe_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularlobe_layer_1");
	if(!specularlobe_MotherWrapping1) specularlobe_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularlobe");
	if(!specularlobe_MotherWrapping1)
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The" << counterString << " wrapping layer's fraction of specular-lobe-reflection has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertyVector * backscattering_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_backscattering_layer_1");
	if(!backscattering_MotherWrapping1) backscattering_MotherWrapping1 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_backscattering");
	if(!backscattering_MotherWrapping1)
	{
		G4String counterString = "";
		if(Wrapping2Exists) counterString = " first";

		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
		std::endl << "# The" << counterString << " wrapping layer's fraction of backscattering-reflection has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * MPT_OptSurf_MotherWrapping1 = new G4MaterialPropertiesTable();
	if(reflectivity_MotherWrapping1) MPT_OptSurf_MotherWrapping1->AddProperty(reflectivityString_MotherWrapping1.c_str(), reflectivity_MotherWrapping1);
	if(refractiveIndex_real_MotherWrapping1 && refractiveIndex_imag_MotherWrapping1)
	{
		MPT_OptSurf_MotherWrapping1->AddProperty("REALRINDEX", refractiveIndex_real_MotherWrapping1);
		MPT_OptSurf_MotherWrapping1->AddProperty("IMAGINARYRINDEX", refractiveIndex_imag_MotherWrapping1);
	}
	MPT_OptSurf_MotherWrapping1->AddProperty("SPECULARSPIKECONSTANT", specularspike_MotherWrapping1);
	MPT_OptSurf_MotherWrapping1->AddProperty("SPECULARLOBECONSTANT", specularlobe_MotherWrapping1);
	MPT_OptSurf_MotherWrapping1->AddProperty("BACKSCATTERCONSTANT", backscattering_MotherWrapping1);
	OptSurf_MotherWrapping1->SetMaterialPropertiesTable(MPT_OptSurf_MotherWrapping1);


	if(Wrapping2Exists)
	{
		G4MaterialPropertyVector * reflectivity_Wrapping1Wrapping2 = 0;
		G4MaterialPropertyVector * refractiveIndex_real_Wrapping1Wrapping2 = 0;
		G4MaterialPropertyVector * refractiveIndex_imag_Wrapping1Wrapping2 = 0;
		G4MaterialPropertyVector * reflectivity_MotherWrapping2 = 0;
		G4MaterialPropertyVector * refractiveIndex_real_MotherWrapping2 = 0;
		G4MaterialPropertyVector * refractiveIndex_imag_MotherWrapping2 = 0;
		G4String reflectivityString_MotherWrapping2 = "";
		G4String reflectivityString_Wrapping1Wrapping2 = "";

		if(Wrapping2IsMetal) reflectivity_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_2");
		else reflectivity_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_2", true);

		if(reflectivity_MotherWrapping2)
		{
			if(Wrapping2IsMetal) reflectivityString_MotherWrapping2 = "REFLECTIVITY";
			else reflectivityString_MotherWrapping2 = "TRANSMITTANCE";
		}
		else if(Wrapping2IsMetal)
		{
			// using the complex reflective index to specify the reflectivity spectrum
			refractiveIndex_real_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real_layer_2");
			refractiveIndex_imag_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag_layer_2");

			if(!refractiveIndex_real_MotherWrapping2 || !refractiveIndex_imag_MotherWrapping2)
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
				std::endl << "# The second wrapping layer's reflectivity has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
		}

		if(Wrapping1IsMetal && !Wrapping2IsMetal)
		{
			reflectivity_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_1");
			if(!reflectivity_Wrapping1Wrapping2) reflectivity_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec");

			if(reflectivity_Wrapping1Wrapping2) reflectivityString_Wrapping1Wrapping2 = "REFLECTIVITY";
			else
			{
				// using the complex reflective index to specify the reflectivity spectrum
				refractiveIndex_real_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real_layer_1");
				refractiveIndex_imag_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag_layer_1");
				if(!refractiveIndex_real_Wrapping1Wrapping2 || !refractiveIndex_imag_Wrapping1Wrapping2)
				{
					refractiveIndex_real_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real");
					refractiveIndex_imag_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag");
				}

				if(!refractiveIndex_real_Wrapping1Wrapping2 || !refractiveIndex_imag_Wrapping1Wrapping2)
				{
					std::cerr <<
					std::endl << "##########" <<
					std::endl << "# ERROR:" <<
					std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
					std::endl << "# The reflectivity between the wrapping layers has not been specified." <<
					std::endl << "##########" <<
					std::endl << std::endl;

					CriticalErrorOccured = true;
				}
			}
		}
		else if(!Wrapping1IsMetal && Wrapping2IsMetal)
		{
			reflectivity_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_2");

			if(reflectivity_Wrapping1Wrapping2) reflectivityString_Wrapping1Wrapping2 = "REFLECTIVITY";
			else
			{
				// using the complex reflective index to specify the reflectivity spectrum
				refractiveIndex_real_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_real_layer_2");
				refractiveIndex_imag_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "n_ref_imag_layer_2");

				if(!refractiveIndex_real_Wrapping1Wrapping2 || !refractiveIndex_imag_Wrapping1Wrapping2)
				{
					std::cerr <<
					std::endl << "##########" <<
					std::endl << "# ERROR:" <<
					std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
					std::endl << "# The reflectivity between the wrapping layers has not been specified." <<
					std::endl << "##########" <<
					std::endl << std::endl;

					CriticalErrorOccured = true;
				}
			}
		}
		else if(!Wrapping1IsMetal && !Wrapping2IsMetal)
		{
			reflectivity_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_2", true);
			if(!reflectivity_Wrapping1Wrapping2) reflectivity_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec_layer_1", true);
			if(!reflectivity_Wrapping1Wrapping2) reflectivity_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "reflec", true);

			if(reflectivity_Wrapping1Wrapping2) reflectivityString_Wrapping1Wrapping2 = "TRANSMITTANCE";
		}

		if( !(Wrapping1IsMetal && Wrapping2IsMetal) )
		{
			G4MaterialPropertyVector * specularspike_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularspike_layer_2");
			if(!specularspike_Wrapping1Wrapping2)
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
				std::endl << "# The fraction of specular-spike-reflection between the wrapping layers has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}

			G4MaterialPropertyVector * specularlobe_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularlobe_layer_2");
			if(!specularlobe_Wrapping1Wrapping2)
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
				std::endl << "# The fraction of specular-lobe-reflection between the wrapping layers has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}

			G4MaterialPropertyVector * backscattering_Wrapping1Wrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_backscattering_layer_2");
			if(!backscattering_Wrapping1Wrapping2)
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
				std::endl << "# The fraction of backscattering-reflection between the wrapping layers has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}

			G4MaterialPropertiesTable * MPT_OptSurf_Wrapping1Wrapping2 = new G4MaterialPropertiesTable();
			if(reflectivity_Wrapping1Wrapping2) MPT_OptSurf_Wrapping1Wrapping2->AddProperty(reflectivityString_Wrapping1Wrapping2.c_str(), reflectivity_Wrapping1Wrapping2);
			if(refractiveIndex_real_Wrapping1Wrapping2 && refractiveIndex_imag_Wrapping1Wrapping2)
			{
				MPT_OptSurf_Wrapping1Wrapping2->AddProperty("REALRINDEX", refractiveIndex_real_Wrapping1Wrapping2);
				MPT_OptSurf_Wrapping1Wrapping2->AddProperty("IMAGINARYRINDEX", refractiveIndex_imag_Wrapping1Wrapping2);
			}
			MPT_OptSurf_Wrapping1Wrapping2->AddProperty("SPECULARSPIKECONSTANT", specularspike_Wrapping1Wrapping2);
			MPT_OptSurf_Wrapping1Wrapping2->AddProperty("SPECULARLOBECONSTANT", specularlobe_Wrapping1Wrapping2);
			MPT_OptSurf_Wrapping1Wrapping2->AddProperty("BACKSCATTERCONSTANT", backscattering_Wrapping1Wrapping2);
			OptSurf_Wrapping1Wrapping2->SetMaterialPropertiesTable(MPT_OptSurf_Wrapping1Wrapping2);
		}

		G4MaterialPropertyVector * specularspike_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularspike_layer_2");
		if(!specularspike_MotherWrapping2)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The second wrapping layer's fraction of specular-spike-reflection has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		G4MaterialPropertyVector * specularlobe_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_specularlobe_layer_2");
		if(!specularlobe_MotherWrapping2)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The second wrapping layer's fraction of specular-lobe-reflection has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		G4MaterialPropertyVector * backscattering_MotherWrapping2 = PropertyTools->GetPropertyDistribution(WrappingProperties, "fraction_backscattering_layer_2");
		if(!backscattering_MotherWrapping2)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Wrapping::ConstructWrappingSurface():" <<
			std::endl << "# The second wrapping layer's fraction of backscattering-reflection has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		G4MaterialPropertiesTable * MPT_OptSurf_MotherWrapping2 = new G4MaterialPropertiesTable();
		if(reflectivity_MotherWrapping2) MPT_OptSurf_MotherWrapping2->AddProperty(reflectivityString_MotherWrapping2.c_str(), reflectivity_MotherWrapping2);
		if(refractiveIndex_real_MotherWrapping2 && refractiveIndex_imag_MotherWrapping2)
		{
			MPT_OptSurf_MotherWrapping2->AddProperty("REALRINDEX", refractiveIndex_real_MotherWrapping2);
			MPT_OptSurf_MotherWrapping2->AddProperty("IMAGINARYRINDEX", refractiveIndex_imag_MotherWrapping2);
		}
		MPT_OptSurf_MotherWrapping2->AddProperty("SPECULARSPIKECONSTANT", specularspike_MotherWrapping2);
		MPT_OptSurf_MotherWrapping2->AddProperty("SPECULARLOBECONSTANT", specularlobe_MotherWrapping2);
		MPT_OptSurf_MotherWrapping2->AddProperty("BACKSCATTERCONSTANT", backscattering_MotherWrapping2);
		OptSurf_MotherWrapping2->SetMaterialPropertiesTable(MPT_OptSurf_MotherWrapping2);
	}
}



void G4Wrapping::SetDefaults()
{
	CriticalErrorOccured = false;

	Material_Wrapping1 = 0;
	Material_Wrapping2 = 0;
	Wrapping2Exists = false;
	Wrapping1IsMetal = false;
	Wrapping2IsMetal = false;
	AirGapThickness = 0.;
	Wrapping1Thickness = NAN;
	Wrapping2Thickness = NAN;
	Dimensions = G4ThreeVector(NAN, NAN, NAN);
	Transformation = G4Transform3D();
	Wrapping1_physical = 0;
	Wrapping2_physical = 0;
	OptSurf_MotherWrapping1 = 0;
	OptSurf_Wrapping1Wrapping2 = 0;
	OptSurf_MotherWrapping2 = 0;
}
