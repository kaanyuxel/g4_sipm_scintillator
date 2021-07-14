/*
 * author:      Erik Dietz-Laursonn
 * institution: Physics Institute 3A, RWTH Aachen University, Aachen, Germany
 * copyright:   Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
 */


#include <boost/any.hpp>
#include <boost/lexical_cast.hpp>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Torus.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4VisAttributes.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4PhysicalVolumeStore.hh>

#include <G4Fibre.hh>



// class variables begin with capital letters, local variables with small letters



/**
 *  Function to construct all volumes and surfaces for straight fibres.
 */
void G4Fibre::ConstructVolumes(G4String fibreType)
{
	if(!OnlyInsideMother) GrandMotherAndAuntVolumes = findGrandMotherAndAuntVolumes(CutAuntVolumesWithDaughters);


// embedment
	if(Glued)
	{
		FibreEmbedment_physical = ConstructEmbedmentPhysical(FibreNamePrefix + "_embedment", fibreType);
	}

// absorbing coating
	if(CoatingExists)
	{
		G4Colour colour_insideMother = G4Colour::Grey();
		G4Colour colour_outsideMother = G4Colour::Black();

		Coating_physical = ConstructFibreLayerPhysical(Coating_mothers_physical, FibreTransformation_insideMother, FibreTransformation_outsideMother, FibreNamePrefix + "_absorbingCoating", RelativeFibreRadius_CoatingMin, RelativeFibreRadius_CoatingMax, Material_Coating, colour_insideMother, colour_outsideMother, fibreType);
	}

// third (outer) fibre-cladding
	if(Cladding3Exists)
	{
		G4Colour colour_insideMother = G4Colour(1.0, 0.5, 0.5);   // a different red for fibre parts inside the scintillator
		G4Colour colour_outsideMother = G4Colour::Red();

		Cladding3_physical = ConstructFibreLayerPhysical(Cladding3_mothers_physical, FibreTransformation_insideMother, FibreTransformation_outsideMother, FibreNamePrefix + "_cladding3", RelativeFibreRadius_Cladding2, RelativeFibreRadius_Cladding3, Material_Cladding3, colour_insideMother, colour_outsideMother, fibreType);

		// for rough/smooth FibreStartPoint
		if(!CreateReflectiveStartPointVolumes)
		{
			RoughenedStartPoint_cladding3_physical = ConstructFibreLayerPhysical(RoughenedStartPoint_cladding3_mothers_physical, RoughenedStartPointTransformation_insideMother, RoughenedStartPointTransformation_outsideMother, FibreNamePrefix + "_roughenedStartPoint_cladding3", RelativeFibreRadius_Cladding2, RelativeFibreRadius_Cladding3, Material_Cladding3, colour_insideMother, colour_outsideMother, fibreType, "roughenedStartPoint");
		}

		// for rough/smooth FibreEndPoint
		if(!CreateReflectiveEndPointVolumes)
		{
			RoughenedEndPoint_cladding3_physical = ConstructFibreLayerPhysical(RoughenedEndPoint_cladding3_mothers_physical, RoughenedEndPointTransformation_insideMother, RoughenedEndPointTransformation_outsideMother, FibreNamePrefix + "_roughenedEndPoint_cladding3", RelativeFibreRadius_Cladding2, RelativeFibreRadius_Cladding3, Material_Cladding3, colour_insideMother, colour_outsideMother, fibreType, "roughenedEndPoint");
		}
	}


// second fibre-cladding
	if(Cladding2Exists)
	{
		G4Colour colour_insideMother = G4Colour(0.5, 0.5, 1.0);   // a different blue for fibre parts inside the scintillator
		G4Colour colour_outsideMother = G4Colour::Blue();

		Cladding2_physical = ConstructFibreLayerPhysical(Cladding2_mothers_physical, FibreTransformation_insideMother, FibreTransformation_outsideMother, FibreNamePrefix + "_cladding2", RelativeFibreRadius_Cladding1, RelativeFibreRadius_Cladding2, Material_Cladding2, colour_insideMother, colour_outsideMother, fibreType);

		// for rough/smooth FibreStartPoint
		if(!CreateReflectiveStartPointVolumes)
		{
			RoughenedStartPoint_cladding2_physical = ConstructFibreLayerPhysical(RoughenedStartPoint_cladding2_mothers_physical, RoughenedStartPointTransformation_insideMother, RoughenedStartPointTransformation_outsideMother, FibreNamePrefix + "_roughenedStartPoint_cladding2", RelativeFibreRadius_Cladding1, RelativeFibreRadius_Cladding2, Material_Cladding2, colour_insideMother, colour_outsideMother, fibreType, "roughenedStartPoint");
		}

		// for rough/smooth FibreEndPoint
		if(!CreateReflectiveEndPointVolumes)
		{
			RoughenedEndPoint_cladding2_physical = ConstructFibreLayerPhysical(RoughenedEndPoint_cladding2_mothers_physical, RoughenedEndPointTransformation_insideMother, RoughenedEndPointTransformation_outsideMother, FibreNamePrefix + "_roughenedEndPoint_cladding2", RelativeFibreRadius_Cladding1, RelativeFibreRadius_Cladding2, Material_Cladding2, colour_insideMother, colour_outsideMother, fibreType, "roughenedEndPoint");
		}
	}


// first (inner) fibre-cladding
	if(Cladding1Exists)
	{
		G4Colour colour_insideMother = G4Colour(0.5, 1.0, 0.5);   // a different green for fibre parts inside the scintillator
		G4Colour colour_outsideMother = G4Colour::Green();

		Cladding1_physical = ConstructFibreLayerPhysical(Cladding1_mothers_physical, FibreTransformation_insideMother, FibreTransformation_outsideMother, FibreNamePrefix + "_cladding1", RelativeFibreRadius_Core, RelativeFibreRadius_Cladding1, Material_Cladding1, colour_insideMother, colour_outsideMother, fibreType);

		// for rough/smooth FibreStartPoint
		if(!CreateReflectiveStartPointVolumes)
		{
			RoughenedStartPoint_cladding1_physical = ConstructFibreLayerPhysical(RoughenedStartPoint_cladding1_mothers_physical, RoughenedStartPointTransformation_insideMother, RoughenedStartPointTransformation_outsideMother, FibreNamePrefix + "_roughenedStartPoint_cladding1", RelativeFibreRadius_Core, RelativeFibreRadius_Cladding1, Material_Cladding1, colour_insideMother, colour_outsideMother, fibreType, "roughenedStartPoint");
		}

		// for rough/smooth FibreEndPoint
		if(!CreateReflectiveEndPointVolumes)
		{
			RoughenedEndPoint_cladding1_physical = ConstructFibreLayerPhysical(RoughenedEndPoint_cladding1_mothers_physical, RoughenedEndPointTransformation_insideMother, RoughenedEndPointTransformation_outsideMother, FibreNamePrefix + "_roughenedEndPoint_cladding1", RelativeFibreRadius_Core, RelativeFibreRadius_Cladding1, Material_Cladding1, colour_insideMother, colour_outsideMother, fibreType, "roughenedEndPoint");
		}
	}


// fibre core
	G4Colour colour_insideMother = G4Colour(1.0, 0.5, 0.5);   // a different red for fibre parts inside the scintillator
	G4Colour colour_outsideMother = G4Colour::Red();

	FibreCore_physical = ConstructFibreLayerPhysical(FibreCore_mothers_physical, FibreTransformation_insideMother, FibreTransformation_outsideMother, FibreNamePrefix + "_fibreCore", 0., RelativeFibreRadius_Core, Material_FibreCore, colour_insideMother, colour_outsideMother, fibreType);

	// for rough/smooth FibreStartPoint
	if(!CreateReflectiveStartPointVolumes)
	{
		RoughenedStartPoint_fibreCore_physical = ConstructFibreLayerPhysical(RoughenedStartPoint_fibreCore_mothers_physical, RoughenedStartPointTransformation_insideMother, RoughenedStartPointTransformation_outsideMother, FibreNamePrefix + "_roughenedStartPoint_fibreCore", 0., RelativeFibreRadius_Core, Material_FibreCore, colour_insideMother, colour_outsideMother, fibreType, "roughenedStartPoint");
	}

	// for rough/smooth FibreEndPoint
	if(!CreateReflectiveEndPointVolumes)
	{
		RoughenedEndPoint_fibreCore_physical = ConstructFibreLayerPhysical(RoughenedEndPoint_fibreCore_mothers_physical, RoughenedEndPointTransformation_insideMother, RoughenedEndPointTransformation_outsideMother, FibreNamePrefix + "_roughenedEndPoint_fibreCore", 0., RelativeFibreRadius_Core, Material_FibreCore, colour_insideMother, colour_outsideMother, fibreType, "roughenedEndPoint");
	}


// reflective volume at FibreStartPoint
	if(CreateReflectiveStartPointVolumes)
	{
		ReflectiveStartPoint_physical = ConstructFibreLayerPhysical(ReflectiveStartPoint_mothers_physical, ReflectiveStartPointTransformation_insideMother, ReflectiveStartPointTransformation_outsideMother, FibreNamePrefix + "_reflectiveStartPoint", 0., 1., Material_FibreCore, G4Colour::Black(), G4Colour::Black(), fibreType, "reflectiveStartPoint");  //FIXME Material?
	}
// reflective volume at FibreEndPoint
	if(CreateReflectiveEndPointVolumes)
	{
		ReflectiveEndPoint_physical = ConstructFibreLayerPhysical(ReflectiveEndPoint_mothers_physical, ReflectiveEndPointTransformation_insideMother, ReflectiveEndPointTransformation_outsideMother, FibreNamePrefix + "_reflectiveEndPoint", 0., 1., Material_FibreCore, G4Colour::Black(), G4Colour::Black(), fibreType, "reflectiveEndPoint");  //FIXME Material?
	}


// surfaces

//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.
	ConstructSurface();
}



/**
 *  Function that places the logical volume created by ConstructFibreLayerLogical(G4Transform3D transformationInsideMother, G4Transform3D transformationOutsideMother, G4String nameBase, G4double fibreRadiusMin_rel, G4double fibreRadiusMax_rel, G4bool roughenedPart, G4Material * material, G4Colour colour_in, G4Colour colour_out).
 *
 *  @return <b> if the desired part of the layer exists </b> the pointer to the created physical volume
 *  @return <b> else </b> 0
 */
std::vector<G4VPhysicalVolume *> G4Fibre::ConstructFibreLayerPhysical( std::vector<G4VPhysicalVolume *> &mothers,
								       G4Transform3D cutTransformation,
								       G4Transform3D transformation,
								       G4String nameBase,
								       G4double fibreRadiusMin_rel,
								       G4double fibreRadiusMax_rel,
								       G4Material * material,
								       G4Colour colour_in,
								       G4Colour colour_out,
								       G4String fibreType,
								       G4String fibrePart
								     )
{
	std::vector<boost::any> logicals = ConstructFibreLayerLogical(cutTransformation, transformation, nameBase, fibreRadiusMin_rel, fibreRadiusMax_rel, material, colour_in, colour_out, fibreType, fibrePart);
	std::vector<G4VPhysicalVolume *> physicals;


	for(unsigned int iter = 0; iter < logicals.size(); iter += 4)
	{
		if(!boost::any_cast<G4LogicalVolume *>(logicals[iter])) continue;

		G4String fibreLayer_physical_name = boost::any_cast<G4String>(logicals[iter + 1]);

		transformation = boost::any_cast<G4Transform3D>(logicals[iter + 3]);

		physicals.push_back(new G4PVPlacement(transformation, fibreLayer_physical_name, boost::any_cast<G4LogicalVolume *>(logicals[iter]), boost::any_cast<G4VPhysicalVolume *>(logicals[iter + 2]), false, 0, SearchOverlaps));
		mothers.push_back(boost::any_cast<G4VPhysicalVolume *>(logicals[iter + 2]));
	}

	return physicals;
}



/**
 *  Function that creates the logical volume for one layer of the G4Fibre% (either inside or outside the mother volume).
 *
 *  @return <b> if the desired part of the layer exists </b> the pointer to the created logical volume
 *  @return <b> else </b> 0
 */
std::vector<boost::any> G4Fibre::ConstructFibreLayerLogical( G4Transform3D transformationInsideMother,
							     G4Transform3D transformationOutsideMother,
							     G4String nameBase,
							     G4double fibreRadiusMin_rel,
							     G4double fibreRadiusMax_rel,
							     G4Material * material,
							     G4Colour colour_in,
							     G4Colour colour_out,
							     G4String fibreType,
							     G4String fibrePart )
{

// solids (dimensions):

	G4VSolid * fibreLayer_solid = 0;
	std::vector<boost::any> solids;
	G4String volumeName = "";

	if(fibreType == "bent")
	{
		G4double startAngle = BendingStartAngle;
		G4double deltaAngle = BendingDeltaAngle;

		if(fibrePart == "roughenedEndPoint")
		{
			startAngle += deltaAngle - Rough_bendingDeltaAngle;
			deltaAngle = Rough_bendingDeltaAngle;
		}
		else if(fibrePart == "roughenedStartPoint")
		{
			deltaAngle = Rough_bendingDeltaAngle;
		}
		else if(fibrePart == "reflectiveEndPoint")
		{
			startAngle += deltaAngle - Reflective_bendingDeltaAngle;
			deltaAngle = Reflective_bendingDeltaAngle;
		}
		else if(fibrePart == "reflectiveStartPoint")
		{
			deltaAngle = Reflective_bendingDeltaAngle;
		}
		else
		{
			if(CreateReflectiveStartPointVolumes)
			{
				startAngle += Reflective_bendingDeltaAngle;
				deltaAngle -= Reflective_bendingDeltaAngle;
			}
			else
			{
				startAngle += Rough_bendingDeltaAngle;
				deltaAngle -= Rough_bendingDeltaAngle;
			}

			if(CreateReflectiveEndPointVolumes)
			{
				deltaAngle -= Reflective_bendingDeltaAngle;
			}
			else
			{
				deltaAngle -= Rough_bendingDeltaAngle;
			}
		}

		if(FibreRound)
		{
			// volume of round fibre
			G4double rMin = fibreRadiusMin_rel * FibreRadius;
			G4double rMax = fibreRadiusMax_rel * FibreRadius;

			fibreLayer_solid = new G4Torus(nameBase + "_round_solid", rMin, rMax, BendingRadius, startAngle, deltaAngle);
		}
		else
		{
			G4double edgeLengthMin = fibreRadiusMin_rel * FibreEdgeLength;
			G4double edgeLengthMax = fibreRadiusMax_rel * FibreEdgeLength;

			// volume of square fibre
			G4VSolid * solid = new G4Tubs("solid", BendingRadius - edgeLengthMax / 2., BendingRadius + edgeLengthMax / 2., edgeLengthMax / 2., startAngle, deltaAngle);
			if(fabs(edgeLengthMin) >= 1e-12)
			{
				G4VSolid * cutSolid = new G4Tubs("cutSolid", BendingRadius - edgeLengthMin / 2., BendingRadius + edgeLengthMin / 2., edgeLengthMin / 2., 0. * CLHEP::deg, 360. * CLHEP::deg);
				fibreLayer_solid = new G4SubtractionSolid(nameBase + "_square_solid", solid, cutSolid, G4Transform3D());
			}
			else
			{
				fibreLayer_solid = solid;
				fibreLayer_solid->SetName((nameBase + "_square_solid").c_str());
			}
		}
	}
	else if(fibreType == "straight")
	{
		G4double length = Length;

		if(fibrePart == "roughenedStartPoint" || fibrePart == "roughenedEndPoint")
		{
			length = Rough_thickness;
		}
		else if(fibrePart == "reflectiveStartPoint" || fibrePart == "reflectiveEndPoint")
		{
			length = Reflective_thickness;
		}
		else
		{
			if(CreateReflectiveStartPointVolumes) length -= Reflective_thickness;
			else length -= Rough_thickness;
			if(CreateReflectiveEndPointVolumes) length -= Reflective_thickness;
			else length -= Rough_thickness;
		}

		if(FibreRound)
		{
			G4double phi_start = 0.0 * CLHEP::rad;
			G4double delta_Phi = 2.0 * CLHEP::pi * CLHEP::rad;

			// volume of round fibre
			G4double rMin = fibreRadiusMin_rel * FibreRadius;
			G4double rMax = fibreRadiusMax_rel * FibreRadius;
			fibreLayer_solid = new G4Tubs(nameBase + "_round_solid", rMin, rMax, length / 2., phi_start, delta_Phi);
		}
		else
		{
			G4double edgeLengthMin = fibreRadiusMin_rel * FibreEdgeLength;
			G4double edgeLengthMax = fibreRadiusMax_rel * FibreEdgeLength;

			// volume of square fibre
			G4VSolid * solid = new G4Box("solid", edgeLengthMax / 2., edgeLengthMax / 2., length / 2.);
			if(fabs(edgeLengthMin) >= 1e-12)
			{
				G4VSolid * cutSolid = new G4Box("cutSolid", edgeLengthMin / 2., edgeLengthMin / 2., length / 2. + 1. * CLHEP::mm);
				fibreLayer_solid = new G4SubtractionSolid(nameBase + "_square_solid", solid, cutSolid, G4Transform3D());
			}
			else
			{
				fibreLayer_solid = solid;
				fibreLayer_solid->SetName((nameBase + "_square_solid").c_str());
			}
		}
	}


	// divide fibre into parts outside and inside the mother volume (scintillator):
	volumeName = nameBase + "_(inside " + MotherVolume_physical->GetName() + ")";
	G4VSolid * fibreLayerIntersection_solid = new G4IntersectionSolid((volumeName + "_solid").c_str(), fibreLayer_solid, MotherVolume_solid, transformationInsideMother.inverse());

	if(fabs(fibreLayerIntersection_solid->GetCubicVolume()) < 1e-12) fibreLayerIntersection_solid = 0;

	solids.push_back(fibreLayerIntersection_solid);
	solids.push_back(volumeName);
	solids.push_back(MotherVolume_physical);
	solids.push_back(transformationInsideMother);


	if(!OnlyInsideMother && GrandMotherAndAuntVolumes[0])
	{
		fibreLayer_solid = new G4SubtractionSolid((nameBase + "_(outside " + MotherVolume_physical->GetName() + ")_solid").c_str(), fibreLayer_solid, MotherVolume_solid, transformationInsideMother.inverse());

		if(fabs(fibreLayer_solid->GetCubicVolume()) >= 1e-12)
		{
			// cut away everything outside the grand mother volume
			fibreLayerIntersection_solid = new G4IntersectionSolid((nameBase + "_(outside " + MotherVolume_physical->GetName() + " inside " + GrandMotherAndAuntVolumes[0]->GetName() + ")_solid").c_str(), fibreLayer_solid, GrandMotherAndAuntVolumes[0]->GetLogicalVolume()->GetSolid(), transformationOutsideMother.inverse());

			if(fabs(fibreLayerIntersection_solid->GetCubicVolume() - fibreLayer_solid->GetCubicVolume()) >= 1e-12) fibreLayer_solid = fibreLayerIntersection_solid;

			if(GrandMotherAndAuntVolumes.size() > 1)
			{
				G4bool onlyCut = false;

				for(unsigned int iter = 1; iter < GrandMotherAndAuntVolumes.size(); iter++)
				{
					if(! GrandMotherAndAuntVolumes[iter])
					{
						onlyCut = true;
						continue;
					}

					G4RotationMatrix rotation = GrandMotherAndAuntVolumes[iter]->GetObjectRotationValue().inverse() * transformationOutsideMother.getRotation();

					// WARNING-NOTE: as transform() changes the vector it is applied to, the following 3 commands are needed instead of a 1 line command
					G4ThreeVector translation = transformationOutsideMother.getTranslation();
					translation -= GrandMotherAndAuntVolumes[iter]->GetObjectTranslation();
					translation.transform(GrandMotherAndAuntVolumes[iter]->GetObjectRotationValue().inverse());

					volumeName = nameBase + "_(inside " + GrandMotherAndAuntVolumes[iter]->GetName() + ")";
					fibreLayerIntersection_solid = new G4IntersectionSolid((volumeName + "_solid").c_str(), fibreLayer_solid, GrandMotherAndAuntVolumes[iter]->GetLogicalVolume()->GetSolid(), G4Transform3D(rotation, translation).inverse());


					if(fabs(fibreLayerIntersection_solid->GetCubicVolume()) >= 1e-12)
					{
						if(!onlyCut)
						{
							solids.push_back(fibreLayerIntersection_solid);
							solids.push_back(volumeName);
							solids.push_back(GrandMotherAndAuntVolumes[iter]);
							solids.push_back(G4Transform3D(rotation, translation));
						}

						fibreLayer_solid = new G4SubtractionSolid(("tempSubtraction_" + boost::lexical_cast<G4String>(iter)).c_str(), fibreLayer_solid, GrandMotherAndAuntVolumes[iter]->GetLogicalVolume()->GetSolid(), G4Transform3D(rotation, translation).inverse());
					}
// 					else
// 					{
// 						GrandMotherAndAuntVolumes.erase(GrandMotherAndAuntVolumes.begin() + iter);
// 						iter--;
// 					}
				}

				if(fabs(fibreLayer_solid->GetCubicVolume()) >= 1e-12)
				{
					volumeName = nameBase + "_(inside " + GrandMotherAndAuntVolumes[0]->GetName() + ")";
					fibreLayer_solid->SetName((volumeName + "_solid").c_str());

					solids.push_back(fibreLayer_solid);
					solids.push_back(volumeName);
					solids.push_back(GrandMotherAndAuntVolumes[0]);
					solids.push_back(transformationOutsideMother);
				}
			}
			else
			{
				volumeName = nameBase + "_(inside " + GrandMotherAndAuntVolumes[0]->GetName() + ")";
				fibreLayer_solid->SetName((volumeName + "_solid").c_str());

				solids.push_back(fibreLayer_solid);
				solids.push_back(volumeName);
				solids.push_back(GrandMotherAndAuntVolumes[0]);
				solids.push_back(transformationOutsideMother);
			}
		}
// 		else
// 		{
// 			GrandMotherAndAuntVolumes.clear();
// 		}
	}
	else
	{
		GrandMotherAndAuntVolumes.clear();
	}

// logical volumes (material):

	std::vector<boost::any> logicals;

	for(unsigned int iter = 0; iter < solids.size(); iter += 4)
	{
		if(!boost::any_cast<G4VSolid *>(solids[iter])) continue;

		G4Colour colour = colour_out;
		if(iter == 0) colour = colour_in;

		G4String nameOfVolume = boost::any_cast<G4String>(solids[iter + 1]);


		G4LogicalVolume * fibreLayer_logical = new G4LogicalVolume(boost::any_cast<G4VSolid *>(solids[iter]), material, (nameOfVolume + "_logical").c_str(), 0, 0, 0);

		G4VisAttributes fibreLayer_vis = G4VisAttributes(colour);
		fibreLayer_vis.SetForceSolid(true);
		fibreLayer_logical->SetVisAttributes(fibreLayer_vis);

		logicals.push_back(fibreLayer_logical);
		logicals.push_back(solids[iter + 1]);
		logicals.push_back(solids[iter + 2]);
		logicals.push_back(solids[iter + 3]);


		if(ConstructSensitiveDetector)
		{
			// try to find an already existing senitive detector
			FibreSensitiveDetector * sensitiveDetector = (FibreSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("FibreSD", false);
			if(!sensitiveDetector)
			{
				// create a new senitive detector
				sensitiveDetector = new FibreSensitiveDetector("FibreSD", DataStorage);
				// register it to Geant4
				G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDetector);
			}

			// assign it to detector parts
			fibreLayer_logical->SetSensitiveDetector(sensitiveDetector);
		}
	}

	return logicals;
}



/**
 *  Function that places the logical volume created by ConstructFibreLayerLogical(G4Transform3D transformationInsideMother, G4Transform3D transformationOutsideMother, G4String nameBase, G4double fibreRadiusMin_rel, G4double fibreRadiusMax_rel, G4bool roughenedPart, G4Material * material, G4Colour colour_in, G4Colour colour_out).
 *
 *  @return <b> if the desired part of the layer exists </b> the pointer to the created physical volume
 *  @return <b> else </b> 0
 */
std::vector<G4VPhysicalVolume *> G4Fibre::ConstructEmbedmentPhysical( G4String nameBase,
								      G4String fibreType
								     )
{
	std::vector<boost::any> logicals = ConstructEmbedmentLogical(nameBase, fibreType);
	std::vector<G4VPhysicalVolume *> physicals;


	for(unsigned int iter = 0; iter < logicals.size(); iter += 4)
	{
		if(!boost::any_cast<G4LogicalVolume *>(logicals[iter])) continue;

		G4String fibreLayer_physical_name = boost::any_cast<G4String>(logicals[iter + 1]);

		G4Transform3D transformation = boost::any_cast<G4Transform3D>(logicals[iter + 3]);

		physicals.push_back(new G4PVPlacement(transformation, fibreLayer_physical_name + " (optical cement)", boost::any_cast<G4LogicalVolume *>(logicals[iter]), boost::any_cast<G4VPhysicalVolume *>(logicals[iter + 2]), false, 0, SearchOverlaps));
	}

	return physicals;
}



/**
 *  Function that creates the logical volume for one layer of the G4Fibre% (either inside or outside the mother volume).
 *
 *  @return <b> if the desired part of the layer exists </b> the pointer to the created logical volume
 *  @return <b> else </b> 0
 */
std::vector<boost::any> G4Fibre::ConstructEmbedmentLogical( G4String nameBase,
							    G4String fibreType )
{

// solids (dimensions):

	G4VSolid * fibreLayer_solid = 0;
	G4Transform3D cutTransformation = G4Transform3D();
	std::vector<boost::any> solids;
	G4String volumeName = "";

	if(fibreType == "bent")
	{
		G4double embedment_bendingDeltaAngle = 2. * asin(EmbedmentThickness / 2. / BendingRadius);

		G4double embedmentStartAngle = BendingStartAngle;
		G4double embedmentDeltaAngle = BendingDeltaAngle;
		G4double cutStartAngle = BendingStartAngle;
		G4double cutDeltaAngle = BendingDeltaAngle;
		if(WhereIsEmbedmentToBeFlushWithFibre == "")
		{
			embedmentStartAngle -= embedment_bendingDeltaAngle;
			embedmentDeltaAngle += 2 * embedment_bendingDeltaAngle;
		}
		else if(WhereIsEmbedmentToBeFlushWithFibre == "S")
		{
			embedmentDeltaAngle += embedment_bendingDeltaAngle;
			cutStartAngle -= embedment_bendingDeltaAngle;
			cutDeltaAngle += embedment_bendingDeltaAngle;
		}
		else if(WhereIsEmbedmentToBeFlushWithFibre == "E")
		{
			embedmentStartAngle -= embedment_bendingDeltaAngle;
			embedmentDeltaAngle += embedment_bendingDeltaAngle;
			cutDeltaAngle += embedment_bendingDeltaAngle;
		}
		else if(WhereIsEmbedmentToBeFlushWithFibre == "SE")
		{
			cutStartAngle -= embedment_bendingDeltaAngle;
			cutDeltaAngle += 2 * embedment_bendingDeltaAngle;
		}


		G4VSolid * solid = 0;
		G4VSolid * cutSolid = 0;


		if(EmbedmentProfile == "round")
		{
			G4double rMax = 0.;

			if(FibreRound) rMax = FibreRadius + EmbedmentThickness;
			else rMax = sqrt(0.5) * FibreEdgeLength + EmbedmentThickness;

			solid = new G4Torus("solid", 0, rMax, BendingRadius, embedmentStartAngle, embedmentDeltaAngle);
		}
		else if(EmbedmentProfile == "quadratic")
		{
			G4double edgeLengthMax = 0.;

			if(FibreRound) edgeLengthMax = 2 * FibreRadius + 2 * EmbedmentThickness;
			else edgeLengthMax = FibreEdgeLength + 2 * EmbedmentThickness;

			solid = new G4Tubs("solid", BendingRadius - edgeLengthMax / 2., BendingRadius + edgeLengthMax / 2., edgeLengthMax / 2., embedmentStartAngle, embedmentDeltaAngle);
		}

		if(FibreRound)
		{
			cutSolid = new G4Torus("cutSolid", 0, FibreRadius, BendingRadius, cutStartAngle, cutDeltaAngle);
			fibreLayer_solid = new G4SubtractionSolid(nameBase + "_round_solid", solid, cutSolid, G4Transform3D());
		}
		else
		{
			cutSolid = new G4Tubs("cutSolid", BendingRadius - FibreEdgeLength / 2., BendingRadius + FibreEdgeLength / 2., FibreEdgeLength / 2., cutStartAngle, cutDeltaAngle);
			fibreLayer_solid = new G4SubtractionSolid(nameBase + "_square_solid", solid, cutSolid, G4Transform3D());
		}
	}
	else if(fibreType == "straight")
	{
		G4double phi_start = 0.0 * CLHEP::rad;
		G4double delta_Phi = 2.0 * CLHEP::pi * CLHEP::rad;

		G4double embedmentLength = Length;
		G4double cutLength = Length;
		if(WhereIsEmbedmentToBeFlushWithFibre == "")
		{
			embedmentLength += 2. * EmbedmentThickness;
		}
		else if(WhereIsEmbedmentToBeFlushWithFibre == "S")
		{
			embedmentLength += EmbedmentThickness;
			cutLength += EmbedmentThickness;
			cutTransformation = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, EmbedmentThickness));
		}
		else if(WhereIsEmbedmentToBeFlushWithFibre == "E")
		{
			embedmentLength += EmbedmentThickness;
			cutLength += EmbedmentThickness;
			cutTransformation = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, 0, -EmbedmentThickness));
		}
		else if(WhereIsEmbedmentToBeFlushWithFibre == "SE")
		{
			cutLength += 2. * EmbedmentThickness;
		}


		G4VSolid * solid = 0;
		G4VSolid * cutSolid = 0;


		if(EmbedmentProfile == "round")
		{
			G4double rMax = 0.;

			if(FibreRound) rMax = FibreRadius + EmbedmentThickness;
			else rMax = sqrt(0.5) * FibreEdgeLength + EmbedmentThickness;

			solid = new G4Tubs("solid", 0, rMax, embedmentLength / 2., phi_start, delta_Phi);
		}
		else if(EmbedmentProfile == "quadratic")
		{
			G4double edgeLengthMax = 0.;

			if(FibreRound) edgeLengthMax = 2 * FibreRadius + 2 * EmbedmentThickness;
			else edgeLengthMax = FibreEdgeLength + 2 * EmbedmentThickness;

			solid = new G4Box("solid", edgeLengthMax / 2., edgeLengthMax / 2., embedmentLength / 2.);
		}

		if(FibreRound)
		{
			cutSolid = new G4Tubs("cutSolid", 0, FibreRadius, cutLength / 2., phi_start, delta_Phi);
			fibreLayer_solid = new G4SubtractionSolid(nameBase + "_round_solid", solid, cutSolid, cutTransformation.inverse());
		}
		else
		{
			cutSolid = new G4Box("cutSolid", FibreEdgeLength / 2., FibreEdgeLength / 2., cutLength / 2.);
			fibreLayer_solid = new G4SubtractionSolid(nameBase + "_square_solid", solid, cutSolid, cutTransformation.inverse());
		}
	}


	if(VolumesWithEmbedment.size() == 0)
	{
		volumeName = nameBase + "_(inside " + MotherVolume_physical->GetName() + ")";
		G4VSolid * fibreLayerIntersection_solid = new G4IntersectionSolid((volumeName + "_solid").c_str(), fibreLayer_solid, MotherVolume_solid, FibreTransformation_insideMother.inverse());

		if(fabs(fibreLayerIntersection_solid->GetCubicVolume()) < 1e-12) fibreLayerIntersection_solid = 0;

		solids.push_back(fibreLayerIntersection_solid);
		solids.push_back(volumeName);
		solids.push_back(MotherVolume_physical);

		// WARNING-NOTE: as transform() changes the vector it is applied to, the following 3 commands are needed instead of a 1 line command
		G4ThreeVector translation = cutTransformation.getTranslation() / 2.;
		translation.transform(FibreTransformation_insideMother.getRotation());
		translation += FibreTransformation_insideMother.getTranslation();
		solids.push_back(G4Transform3D(FibreTransformation_insideMother.getRotation(), translation));
	}
	else
	{
		G4bool embedmentInMother = false;
		G4bool embedmentInGrandMother = false;
		for(unsigned int iterE = 0; iterE < VolumesWithEmbedment.size(); iterE++)
		{
			if(VolumesWithEmbedment[iterE] == MotherVolume_physical)
			{
				embedmentInMother = true;

				VolumesWithEmbedment.erase(VolumesWithEmbedment.begin() + iterE);
				iterE--;

				if(embedmentInGrandMother) break;
				else continue;
			}

			if(GrandMotherAndAuntVolumes.size() && VolumesWithEmbedment[iterE] == GrandMotherAndAuntVolumes[0])
			{
				embedmentInGrandMother = true;

				VolumesWithEmbedment.erase(VolumesWithEmbedment.begin() + iterE);
				iterE--;

				if(embedmentInMother) break;
				else continue;
			}
		}


		if(embedmentInMother || embedmentInGrandMother)
		{
			G4VSolid * fibreLayerIntersection_solid = new G4IntersectionSolid((volumeName + "_solid").c_str(), fibreLayer_solid, MotherVolume_solid, FibreTransformation_insideMother.inverse());
			volumeName = nameBase + "_(inside " + MotherVolume_physical->GetName() + ")";

			if(fabs(fibreLayerIntersection_solid->GetCubicVolume()) >= 1e-12)
			{
				if(embedmentInGrandMother) fibreLayer_solid = new G4SubtractionSolid("tempSubtraction", fibreLayer_solid, MotherVolume_solid, FibreTransformation_insideMother.inverse());
				if(! embedmentInMother) fibreLayerIntersection_solid = 0;
			}
			else
			{
				fibreLayerIntersection_solid = 0;
			}

			solids.push_back(fibreLayerIntersection_solid);
			solids.push_back(volumeName);
			solids.push_back(MotherVolume_physical);

			// WARNING-NOTE: as transform() changes the vector it is applied to, the following 3 commands are needed instead of a 1 line command
			G4ThreeVector translation = cutTransformation.getTranslation() / 2.;
			translation.transform(FibreTransformation_insideMother.getRotation());
			translation += FibreTransformation_insideMother.getTranslation();
			solids.push_back(G4Transform3D(FibreTransformation_insideMother.getRotation(), translation));
		}


		if(VolumesWithEmbedment.size() > 0 || embedmentInGrandMother)
		{
			G4bool onlyCut = false;

			for(unsigned int iter = 1; iter < GrandMotherAndAuntVolumes.size(); iter++)
			{
				if(! GrandMotherAndAuntVolumes[iter])
				{
					onlyCut = true;
					continue;
				}

				G4RotationMatrix rotation = GrandMotherAndAuntVolumes[iter]->GetObjectRotationValue().inverse() * FibreTransformation_outsideMother.getRotation();

				// WARNING-NOTE: as transform() changes the vector it is applied to, the following 5 commands are needed instead of a 1 line command
				G4ThreeVector translation = cutTransformation.getTranslation() / 2.;
				translation.transform(FibreTransformation_insideMother.getRotation());
				translation += FibreTransformation_outsideMother.getTranslation();
				translation -= GrandMotherAndAuntVolumes[iter]->GetObjectTranslation();
				translation.transform(GrandMotherAndAuntVolumes[iter]->GetObjectRotationValue().inverse());

				volumeName = nameBase + "_(inside " + GrandMotherAndAuntVolumes[iter]->GetName() + ")";
				G4VSolid * fibreLayerIntersection_solid = new G4IntersectionSolid((volumeName + "_solid").c_str(), fibreLayer_solid, GrandMotherAndAuntVolumes[iter]->GetLogicalVolume()->GetSolid(), G4Transform3D(rotation, translation).inverse());


				if(fabs(fibreLayerIntersection_solid->GetCubicVolume()) >= 1e-12)
				{

					G4bool saveVolume = false;
					if(!onlyCut)
					{
						for(unsigned int iterE = 0; iterE < VolumesWithEmbedment.size(); iterE++)
						{
							if(VolumesWithEmbedment[iterE] == GrandMotherAndAuntVolumes[iter])
							{
								VolumesWithEmbedment.erase(VolumesWithEmbedment.begin() + iterE);

								saveVolume = true;

								break;
							}
						}
					}

					if(saveVolume)
					{
						solids.push_back(fibreLayerIntersection_solid);
						solids.push_back(volumeName);
						solids.push_back(GrandMotherAndAuntVolumes[iter]);
						solids.push_back(G4Transform3D(rotation, translation));
					}

					if(embedmentInGrandMother) fibreLayer_solid = new G4SubtractionSolid(("tempSubtraction_" + boost::lexical_cast<G4String>(iter)).c_str(), fibreLayer_solid, GrandMotherAndAuntVolumes[iter]->GetLogicalVolume()->GetSolid(), G4Transform3D(rotation, translation).inverse());
				}
			}
		}


		if(embedmentInGrandMother && fabs(fibreLayer_solid->GetCubicVolume()) >= 1e-12)
		{
			volumeName = nameBase + "_(inside " + GrandMotherAndAuntVolumes[0]->GetName() + ")";
			fibreLayer_solid->SetName((volumeName + "_solid").c_str());

			solids.push_back(fibreLayer_solid);
			solids.push_back(volumeName);
			solids.push_back(GrandMotherAndAuntVolumes[0]);

			// WARNING-NOTE: as transform() changes the vector it is applied to, the following 3 commands are needed instead of a 1 line command
			G4ThreeVector translation = cutTransformation.getTranslation() / 2.;
			translation.transform(FibreTransformation_outsideMother.getRotation());
			translation += FibreTransformation_outsideMother.getTranslation();
			solids.push_back(G4Transform3D(FibreTransformation_outsideMother.getRotation(), translation));
		}
	}



// logical volumes (material):

	std::vector<boost::any> logicals;

	for(unsigned int iter = 0; iter < solids.size(); iter += 4)
	{
		if(!boost::any_cast<G4VSolid *>(solids[iter])) continue;

		G4String nameOfVolume = boost::any_cast<G4String>(solids[iter + 1]);


		G4LogicalVolume * fibreLayer_logical = new G4LogicalVolume(boost::any_cast<G4VSolid *>(solids[iter]), Material_OpticalCement, (nameOfVolume + "_logical").c_str(), 0, 0, 0);

		G4VisAttributes fibreLayer_vis = G4VisAttributes(G4Colour::Cyan());
		fibreLayer_vis.SetForceWireframe(true);
		fibreLayer_logical->SetVisAttributes(fibreLayer_vis);

		logicals.push_back(fibreLayer_logical);
		logicals.push_back(solids[iter + 1]);
		logicals.push_back(solids[iter + 2]);
		logicals.push_back(solids[iter + 3]);


		if(ConstructSensitiveDetector)
		{
			// try to find an already existing senitive detector
			FibreSensitiveDetector * sensitiveDetector = (FibreSensitiveDetector *) G4SDManager::GetSDMpointer()->FindSensitiveDetector("FibreSD", false);
			if(!sensitiveDetector)
			{
				// create a new senitive detector
				sensitiveDetector = new FibreSensitiveDetector("FibreSD", DataStorage);
				// register it to Geant4
				G4SDManager::GetSDMpointer()->AddNewDetector(sensitiveDetector);
			}

			// assign it to detector parts
			fibreLayer_logical->SetSensitiveDetector(sensitiveDetector);
		}
	}

	return logicals;
}



/**
 *  Function to initialise variables:
 *  - which claddings is to be created?
 *  - which profile does the fibre have?
 *  - determines the relative radii of the core/claddings
 */
void G4Fibre::InitialiseVariables()
{
	MotherVolume_solid = MotherVolume_physical->GetLogicalVolume()->GetSolid();

	if(Glued)
	{
		if(OpticalCementProperties.containsNumber("d_layer")) EmbedmentThickness = OpticalCementProperties.getNumber("d_layer");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# The embedment layer's thickness has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	// which claddings is to be created?
	if(FibreProperties.containsNumber("R_rel_cladding_1") || FibreProperties.containsNumber("R_rel_cladding")) Cladding1Exists = true;   // if R_rel_cladding_XYZ is given, i.e. a the cladding exists
	else Cladding1Exists = false;

	if(FibreProperties.containsNumber("R_rel_cladding_2"))
	{
		if(Cladding1Exists) Cladding2Exists = true;
		else
		{
			Cladding2Exists = false;

			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is no first but a second cladding???" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	if(FibreProperties.containsNumber("R_rel_cladding_3"))
	{
		if(Cladding2Exists) Cladding3Exists = true;
		else
		{
			Cladding3Exists = false;

			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is no second but a third cladding???" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	if(FibreProperties.containsNumber("R_relMax_coating")) CoatingExists = true;

	// which profile is to be used?
	G4String fibre_cross_section = FibreProperties.getString("profile");
	if(fibre_cross_section == "round") FibreRound = true;
	else if(fibre_cross_section == "quadratic") FibreRound = false;
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# WARNING:" <<
		std::endl << "# G4Fibre::InitialiseVariables():" <<
		std::endl << "# The fibre's cross section has not been specified properly. It will be set to \"round\"." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		FibreRound = true;
	}

	if(FibreRound)
	{
		if(FibreProperties.containsNumber("radius")) FibreRadius = FibreProperties.getNumber("radius");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# The fibre's radius has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(FibreBent && BendingRadius - FibreRadius < 1e-12)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - The fibre is bent and the fibre's radius is bigger than (or equal to) its bending radius." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else
	{
		if(FibreProperties.containsNumber("edge_length")) FibreEdgeLength = FibreProperties.getNumber("edge_length");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# The fibre's edge length has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(FibreBent && BendingRadius - FibreEdgeLength / 2. < 1e-12)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - The fibre is bent and the fibre's edge length is bigger than (or equal to) its bending diameter." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	// get the relative radii of the core/claddings
	if(Cladding1Exists)
	{
		if(FibreProperties.containsNumber("R_rel_cladding_1")) RelativeFibreRadius_Cladding1 = FibreProperties.getNumber("R_rel_cladding_1");
		else RelativeFibreRadius_Cladding1 = FibreProperties.getNumber("R_rel_cladding");

		if(! RelativeFibreRadius_Cladding1)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is a first cladding but with thickness 0???" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else RelativeFibreRadius_Cladding1 = 0.;

	if(Cladding2Exists)
	{
		RelativeFibreRadius_Cladding2 = FibreProperties.getNumber("R_rel_cladding_2");

		if(! RelativeFibreRadius_Cladding2)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is a second cladding but with thickness 0???" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else RelativeFibreRadius_Cladding2 = 0.;

	if(Cladding3Exists)
	{
		RelativeFibreRadius_Cladding3 = FibreProperties.getNumber("R_rel_cladding_3");

		if(! RelativeFibreRadius_Cladding3)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is a third cladding but with thickness 0???" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else RelativeFibreRadius_Cladding3 = 0.;

	RelativeFibreRadius_Core = 1. - (RelativeFibreRadius_Cladding1 + RelativeFibreRadius_Cladding2 + RelativeFibreRadius_Cladding3);
	if(Cladding1Exists) RelativeFibreRadius_Cladding1 = 1. - (RelativeFibreRadius_Cladding2 + RelativeFibreRadius_Cladding3);
	if(Cladding2Exists) RelativeFibreRadius_Cladding2 = 1. - RelativeFibreRadius_Cladding3;
	if(Cladding3Exists) RelativeFibreRadius_Cladding3 = 1.;


	if(CoatingExists)
	{
		RelativeFibreRadius_CoatingMin = FibreProperties.getNumber("R_relMin_coating");
		RelativeFibreRadius_CoatingMax = FibreProperties.getNumber("R_relMax_coating");

		if(isnan(RelativeFibreRadius_CoatingMax))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is a coating but its outer radius has not been specified correctly!!!" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(isnan(RelativeFibreRadius_CoatingMin))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is a coating but its inner radius has not been specified correctly. It will be set to the value of the fibre radius." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			RelativeFibreRadius_CoatingMin = 1.;
		}

		if(RelativeFibreRadius_CoatingMin < 1.)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# The inner radius of the fibre coating is smaller than the fibre radius. It will be set to the value of the fibre radius." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			RelativeFibreRadius_CoatingMin = 1.;
		}

		if(RelativeFibreRadius_CoatingMax < RelativeFibreRadius_CoatingMin)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::InitialiseVariables():" <<
			std::endl << "# Something went wrong - There is a coating but its inner radius is smaller than its outer radius!!!" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
}



std::vector<G4VPhysicalVolume *> G4Fibre::findGrandMotherAndAuntVolumes(G4bool cutAuntVolumesWithDaughters)
{
	// the default cut_volumes are all physical volumes
	std::vector<G4VPhysicalVolume *> allVolumes = *G4PhysicalVolumeStore::GetInstance();

	std::vector<G4VPhysicalVolume *> grandMotherAndAuntVolumes;
	std::vector<G4VPhysicalVolume *> cutonlyAuntVolumes;


	// first the grandMother (in the order as they were created)
	if(MotherVolume_physical->GetMotherLogical())
	{
		for(unsigned int iter = 0; iter < allVolumes.size(); iter++)
		{
			if(allVolumes[iter]->GetLogicalVolume() == MotherVolume_physical->GetMotherLogical())
			{
				grandMotherAndAuntVolumes.push_back(allVolumes[iter]);
				break;
			}
		}
	}
	else
	{
		grandMotherAndAuntVolumes.push_back(0);
		return grandMotherAndAuntVolumes;
	}


	// then the aunts
	for(unsigned int iter = 0; iter < allVolumes.size(); iter++)
	{
		if(allVolumes[iter]->GetMotherLogical())
		{
			if(allVolumes[iter]->GetMotherLogical() == grandMotherAndAuntVolumes[0]->GetLogicalVolume() && allVolumes[iter]->GetLogicalVolume() != MotherVolume_physical->GetLogicalVolume())
			{
				if(! (cutAuntVolumesWithDaughters && allVolumes[iter]->GetLogicalVolume()->GetNoDaughters()) )
				{
					grandMotherAndAuntVolumes.push_back(allVolumes[iter]);
				}
				else
				{
					cutonlyAuntVolumes.push_back(allVolumes[iter]);
				}
			}
		}
	}

	// and the aunts that should not be used
	if(cutonlyAuntVolumes.size())
	{
		grandMotherAndAuntVolumes.push_back(0);

		for(unsigned int iter = 0; iter < cutonlyAuntVolumes.size(); iter++)
		{
			grandMotherAndAuntVolumes.push_back(cutonlyAuntVolumes[iter]);
		}
	}

	// in the order they were created
	return grandMotherAndAuntVolumes;
}




/**
 *  Function to generate the transformations needed:
 *  - the transformation for cutting the fibre into parts inside and outside the mother volume
 *  - the transformation for placing the fibre parts inside the mother volume
 *  - the transformation for placing the fibre parts outside the mother volume
 *
 *  For this purpose, the function considers:
 *  - the transformation of the fibre relative to the reference volume
 *  - the transformation of the reference volume relative to the mother volume
 */
void G4Fibre::GenerateTransformation(G4String fibreType)
{
	// get the fibre's transformation relative to the reference volume
	G4RotationMatrix fibreRotation_rel = G4RotationMatrix();
	G4ThreeVector fibreTranslation_rel = G4ThreeVector(0., 0., 0.);

	if(!isnan(Length))   //i.e. a straight fibre is defined by Transformation and length
	{
		fibreTranslation_rel = FibreTransformation_rel.getTranslation();
		fibreRotation_rel = FibreTransformation_rel.getRotation();
	}
	else   //i.e. the fibre is bent or defined by start and end point -> FibreTransformation_rel has to be calculated
	{
		G4ThreeVector fibreAxis;
		if(fibreType == "straight")
		{
			fibreAxis = EndPoint - StartPoint;   // vector from the start to the end point
			Length = fibreAxis.mag();
			fibreTranslation_rel = (StartPoint + EndPoint) / 2.;
		}
		else if(fibreType == "bent")
		{
			fibreAxis = BendingAxis;
			fibreTranslation_rel = BendingCircularCentre;
		}

		// get rotation parameters for positioning the fibre
		G4ThreeVector originAxis(0, 0, 1);   // the fibre is originally generated parallel to (straight fibre) or circularly around (bent fibre) the z-axis
		G4ThreeVector rotationAxis = originAxis.cross(fibreAxis);

		G4double alpha = 0.;
		G4ThreeVector referenceAxis(0, 1, 0);
		if(fabs(rotationAxis.mag()) > 1e-12)
		{
			alpha = fibreAxis.angle(originAxis);
		}
		else if(fibreAxis / fibreAxis.mag() == - originAxis)
		{
			alpha = 180. * CLHEP::deg;
			rotationAxis = referenceAxis;
		}
		fibreRotation_rel.rotate(alpha, rotationAxis);

		if(fibreType == "straight")
		{
			// get get rotation parameters for the right orientation of the fibre
			if(fabs(referenceAxis.cross(rotationAxis).mag()) > 1e-12)
			{
				G4ThreeVector surfaceNormal_requested = fibreAxis.cross(referenceAxis);
				surfaceNormal_requested = surfaceNormal_requested.cross(fibreAxis);

				G4ThreeVector surfaceNormal_real = referenceAxis.rotate(alpha, rotationAxis);

				G4double beta = surfaceNormal_real.angle(surfaceNormal_requested);
				fibreRotation_rel.rotate(beta, fibreAxis);
			}
		}
		else if(fibreType == "bent")
		{
			// calculate the the angle at which the bent fibre starts (if it is not given)
			if(isnan(BendingStartAngle))
			{
				// NOTE: the fibre is originally generated circularly around the z-axis, 0° is parallel to the x-axis and the angle is counted mathematically right-handed around the z-axis

				// WARNING-NOTE: as transform() changes the vector it is applied to, the following 3 commands are needed instead of a 1 line command
				G4ThreeVector startPoint_rel = StartPoint - BendingCircularCentre;
				startPoint_rel /= startPoint_rel.mag();
				startPoint_rel.transform(fibreRotation_rel.inverse());

				if(startPoint_rel.y() < 0) BendingStartAngle = 180. * CLHEP::deg + acos(- startPoint_rel.x());
				else BendingStartAngle = acos(startPoint_rel.x());
			}
		}
	}


	// considering corrections due to the circular shape of bent fibres
	G4ThreeVector centreStartEndCorretionVector = G4ThreeVector(0, 0, 0);
	if(fibreType == "bent")
	{
		G4ThreeVector centreStartEndPoint_rel = BendingRadius * fabs(cos(BendingDeltaAngle / 2.)) * G4ThreeVector(cos(BendingStartAngle + BendingDeltaAngle / 2.), sin(BendingStartAngle + BendingDeltaAngle / 2.) ,0);
		if(fabs(centreStartEndPoint_rel.mag()) > 1e-12)
		{
			// WARNING-NOTE: as transform() changes the vector it is applied to, the following 6 commands are needed instead of a 1 line command
			centreStartEndCorretionVector += centreStartEndPoint_rel;
			centreStartEndCorretionVector.rotate(-90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
			centreStartEndCorretionVector.transform(ReplicationTransformation.getRotation());
			centreStartEndCorretionVector.rotate(90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
			centreStartEndCorretionVector *= -1;
			centreStartEndCorretionVector += centreStartEndPoint_rel;
		}
	}


	// get the fibre's transformation inside the mother volume, which is also needed to cut the fibre into pieces fitting into the different volumina of the mother volume (scintillator, wrapping)
	G4RotationMatrix fibreRotation_insideMother = G4RotationMatrix();
	G4ThreeVector fibreTranslation_insideMother = G4ThreeVector(0., 0., 0.);

	// considering the transformation of the reference volume
	if(fibreType == "straight")
	{
		fibreRotation_insideMother = ReplicationTransformation.getRotation();
	}
	else if(fibreType == "bent")
	{
		fibreRotation_insideMother.rotate(-90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
		fibreRotation_insideMother = ReplicationTransformation.getRotation() * fibreRotation_insideMother;
		fibreRotation_insideMother.rotate(90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
	}
	fibreRotation_insideMother = MotherVolume_physical->GetObjectRotationValue().inverse() * ReferenceVolume_physical->GetObjectRotationValue() * ReferenceVolume_transformation.getRotation() * fibreRotation_rel * fibreRotation_insideMother;

// 	// WARNING-NOTE: as transform() changes the vector it is applied to, the following 11 commands are needed instead of a 1 line command
	fibreTranslation_insideMother = ReplicationTransformation.getTranslation();
	if(fibreType == "bent")
	{
		fibreTranslation_insideMother.rotate(90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
		if(fabs(centreStartEndCorretionVector.mag()) > 1e-12) fibreTranslation_insideMother += centreStartEndCorretionVector;
	}
	fibreTranslation_insideMother.transform(fibreRotation_rel);
	fibreTranslation_insideMother += fibreTranslation_rel;
	fibreTranslation_insideMother.transform(ReferenceVolume_physical->GetObjectRotationValue());
	fibreTranslation_insideMother += ReferenceVolume_physical->GetObjectTranslation();
	fibreTranslation_insideMother.transform(ReferenceVolume_transformation.getRotation());
	fibreTranslation_insideMother += ReferenceVolume_transformation.getTranslation();
	fibreTranslation_insideMother -= MotherVolume_physical->GetObjectTranslation();
	fibreTranslation_insideMother.transform(MotherVolume_physical->GetObjectRotationValue().inverse());

	if(fibreType == "straight")
	{
		if(CreateReflectiveStartPointVolumes) fibreTranslation_insideMother += G4ThreeVector(0., 0., Reflective_thickness / 2.).transform(fibreRotation_insideMother);
		else fibreTranslation_insideMother += G4ThreeVector(0., 0., Rough_thickness / 2.).transform(fibreRotation_insideMother);
		if(CreateReflectiveEndPointVolumes) fibreTranslation_insideMother -= G4ThreeVector(0., 0., Reflective_thickness / 2.).transform(fibreRotation_insideMother);
		else fibreTranslation_insideMother -= G4ThreeVector(0., 0., Rough_thickness / 2.).transform(fibreRotation_insideMother);
	}

	FibreTransformation_insideMother = G4Transform3D(fibreRotation_insideMother, fibreTranslation_insideMother);


	// get the fibre's transformation outside of the mother volume
	G4RotationMatrix fibreRotation_outsideMother = G4RotationMatrix();
	G4ThreeVector fibreTranslation_outsideMother = G4ThreeVector(0., 0., 0.);

	// considering the transformation of the reference volume
	if(fibreType == "straight")
	{
		fibreRotation_outsideMother = ReplicationTransformation.getRotation();
	}
	else if(fibreType == "bent")
	{
		fibreRotation_outsideMother.rotate(-90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
		fibreRotation_outsideMother = ReplicationTransformation.getRotation() * fibreRotation_outsideMother;
		fibreRotation_outsideMother.rotate(90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
	}
	fibreRotation_outsideMother = ReferenceVolume_physical->GetObjectRotationValue() * ReferenceVolume_transformation.getRotation() * fibreRotation_rel * fibreRotation_outsideMother;

	// WARNING-NOTE: as transform() changes the vector it is applied to, the following 9 commands are needed instead of a 1 line command
	fibreTranslation_outsideMother = ReplicationTransformation.getTranslation();
	if(fibreType == "bent")
	{
		fibreTranslation_outsideMother.rotate(90 * CLHEP::deg, G4ThreeVector(1, 0, 0));
		if(fabs(centreStartEndCorretionVector.mag()) > 1e-12) fibreTranslation_outsideMother += centreStartEndCorretionVector;
	}
	fibreTranslation_outsideMother.transform(fibreRotation_rel);
	fibreTranslation_outsideMother += fibreTranslation_rel;
	fibreTranslation_outsideMother.transform(ReferenceVolume_physical->GetObjectRotationValue());
	fibreTranslation_outsideMother += ReferenceVolume_physical->GetObjectTranslation();
	fibreTranslation_outsideMother.transform(ReferenceVolume_transformation.getRotation());
	fibreTranslation_outsideMother += ReferenceVolume_transformation.getTranslation();

	if(fibreType == "straight")
	{
		if(CreateReflectiveStartPointVolumes) fibreTranslation_outsideMother += G4ThreeVector(0., 0., Reflective_thickness / 2.).transform(fibreRotation_outsideMother);
		else fibreTranslation_outsideMother += G4ThreeVector(0., 0., Rough_thickness / 2.).transform(fibreRotation_outsideMother);
		if(CreateReflectiveEndPointVolumes) fibreTranslation_outsideMother -= G4ThreeVector(0., 0., Reflective_thickness / 2.).transform(fibreRotation_outsideMother);
		else fibreTranslation_outsideMother -= G4ThreeVector(0., 0., Rough_thickness / 2.).transform(fibreRotation_outsideMother);
	}

	FibreTransformation_outsideMother = G4Transform3D(fibreRotation_outsideMother, fibreTranslation_outsideMother);


	if(fibreType == "straight")
	{
		G4double length = Length;
		if(CreateReflectiveStartPointVolumes) length -= Reflective_thickness;
		else length -= Rough_thickness;
		if(CreateReflectiveEndPointVolumes) length -= Reflective_thickness;
		else length -= Rough_thickness;

		if(CreateReflectiveStartPointVolumes)
		{
			ReflectiveStartPointTransformation_outsideMother = G4Transform3D(fibreRotation_outsideMother, fibreTranslation_outsideMother + G4ThreeVector(0., 0., (- length - Reflective_thickness) / 2.).transform(fibreRotation_outsideMother));
			ReflectiveStartPointTransformation_insideMother = G4Transform3D(FibreTransformation_insideMother.getRotation(), FibreTransformation_insideMother.getTranslation() + G4ThreeVector(0., 0., (- length - Reflective_thickness) / 2.).transform(fibreRotation_insideMother));
		}
		else
		{
			RoughenedStartPointTransformation_outsideMother = G4Transform3D(fibreRotation_outsideMother, fibreTranslation_outsideMother + G4ThreeVector(0., 0., (- length - Rough_thickness) / 2.).transform(fibreRotation_outsideMother));
			RoughenedStartPointTransformation_insideMother = G4Transform3D(FibreTransformation_insideMother.getRotation(), FibreTransformation_insideMother.getTranslation() + G4ThreeVector(0., 0., (- length - Rough_thickness) / 2.).transform(fibreRotation_insideMother));
		}

		if(CreateReflectiveEndPointVolumes)
		{
			ReflectiveEndPointTransformation_outsideMother = G4Transform3D(fibreRotation_outsideMother, fibreTranslation_outsideMother + G4ThreeVector(0., 0., (length + Reflective_thickness) / 2.).transform(fibreRotation_outsideMother));
			ReflectiveEndPointTransformation_insideMother = G4Transform3D(FibreTransformation_insideMother.getRotation(), FibreTransformation_insideMother.getTranslation() + G4ThreeVector(0., 0., (length + Reflective_thickness) / 2.).transform(fibreRotation_insideMother));
		}
		else
		{
			RoughenedEndPointTransformation_outsideMother = G4Transform3D(fibreRotation_outsideMother, fibreTranslation_outsideMother + G4ThreeVector(0., 0., (length + Rough_thickness) / 2.).transform(fibreRotation_outsideMother));
			RoughenedEndPointTransformation_insideMother = G4Transform3D(FibreTransformation_insideMother.getRotation(), FibreTransformation_insideMother.getTranslation() + G4ThreeVector(0., 0., (length + Rough_thickness) / 2.).transform(fibreRotation_insideMother));
		}
	}
	else if(fibreType == "bent")
	{
		if(CreateReflectiveStartPointVolumes)
		{
			ReflectiveStartPointTransformation_outsideMother = FibreTransformation_outsideMother;
			ReflectiveStartPointTransformation_insideMother = FibreTransformation_insideMother;
		}
		else
		{
			RoughenedStartPointTransformation_outsideMother = FibreTransformation_outsideMother;
			RoughenedStartPointTransformation_insideMother = FibreTransformation_insideMother;
		}

		if(CreateReflectiveEndPointVolumes)
		{
			ReflectiveEndPointTransformation_outsideMother = FibreTransformation_outsideMother;
			ReflectiveEndPointTransformation_insideMother = FibreTransformation_insideMother;
		}
		else
		{
			RoughenedEndPointTransformation_outsideMother = FibreTransformation_outsideMother;
			RoughenedEndPointTransformation_insideMother = FibreTransformation_insideMother;
		}
	}
}



/**
 *  Function to define materials (compounds, alloys) and their properties:
 *  - creates new materials according to the specifications from the property file(s)
 *  - if the same materials already exist, the newly created ones are deleted
 */
void G4Fibre::DefineMaterials()
{
	// NOTE The materials which have already been define (e.g. in the DetectorConstruction) are used here

//---------- optical cement ----------//
	if(Glued)
	{
		G4double density_OpticalCement;
		if(OpticalCementProperties.containsNumber("density")) density_OpticalCement = OpticalCementProperties.getNumber("density");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The optical cement's density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_OpticalCement;
		if(OpticalCementProperties.containsTabular("chemical_components")) chemicalComponents_OpticalCement = OpticalCementProperties.getTabular("chemical_components");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The optical cement's chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_OpticalCement = new G4Material("tempName", density_OpticalCement, (int) chemicalComponents_OpticalCement["element"].size());
		PropertyTools->AddElementsFromTable(Material_OpticalCement, chemicalComponents_OpticalCement);

		DefineOpticalCementMaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("OpticalCementMaterial", Material_OpticalCement);
	}

//---------- fibre-core ----------//
	G4double density_FibreCore;
	if(FibreProperties.containsNumber("density_core")) density_FibreCore = FibreProperties.getNumber("density_core");
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineMaterials():" <<
		std::endl << "# The fibre core's density has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	GoddessProperties::tabular chemicalComponents_FibreCore;
	if(FibreProperties.containsTabular("chemical_components_core")) chemicalComponents_FibreCore = FibreProperties.getTabular("chemical_components_core");
	else
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineMaterials():" <<
		std::endl << "# The fibre core's chemical components have not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	Material_FibreCore = new G4Material("tempName", density_FibreCore, (int) chemicalComponents_FibreCore["element"].size());
	PropertyTools->AddElementsFromTable(Material_FibreCore, chemicalComponents_FibreCore);

	DefineCoreMaterialProperties();

	PropertyTools->checkIfMaterialAlreadyExists("FibreCoreMaterial", Material_FibreCore);

//---------- first (inner) fibre-cladding ----------//
	if(Cladding1Exists)
	{
		G4double density_Cladding1;
		if(FibreProperties.containsNumber("density_cladding_1")) density_Cladding1 = FibreProperties.getNumber("density_cladding_1");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre first cladding's density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_Cladding1;
		if(FibreProperties.containsTabular("chemical_components_cladding_1")) chemicalComponents_Cladding1 = FibreProperties.getTabular("chemical_components_cladding_1");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre first cladding's chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_Cladding1 = new G4Material("tempName", density_Cladding1, (int) chemicalComponents_Cladding1["element"].size());
		PropertyTools->AddElementsFromTable(Material_Cladding1, chemicalComponents_Cladding1);

		DefineCladding1MaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("FibreCladding1Material", Material_Cladding1);
	}

//---------- second fibre-cladding ----------//
	if(Cladding2Exists)
	{
		G4double density_Cladding2;
		if(FibreProperties.containsNumber("density_cladding_2")) density_Cladding2 = FibreProperties.getNumber("density_cladding_2");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre second cladding's density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_Cladding2;
		if(FibreProperties.containsTabular("chemical_components_cladding_2")) chemicalComponents_Cladding2 = FibreProperties.getTabular("chemical_components_cladding_2");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre second cladding's chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_Cladding2 = new G4Material("tempName", density_Cladding2, (int) chemicalComponents_Cladding2["element"].size());
		PropertyTools->AddElementsFromTable(Material_Cladding2, chemicalComponents_Cladding2);

		DefineCladding2MaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("FibreCladding2Material", Material_Cladding2);
	}

//---------- third (outer) fibre-cladding ----------//
	if(Cladding3Exists)
	{
		G4double density_Cladding3;
		if(FibreProperties.containsNumber("density_cladding_3")) density_Cladding3 = FibreProperties.getNumber("density_cladding_3");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre third cladding's density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_Cladding3;
		if(FibreProperties.containsTabular("chemical_components_cladding_3")) chemicalComponents_Cladding3 = FibreProperties.getTabular("chemical_components_cladding_3");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre third cladding's chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_Cladding3 = new G4Material("tempName", density_Cladding3, (int) chemicalComponents_Cladding3["element"].size());
		PropertyTools->AddElementsFromTable(Material_Cladding3, chemicalComponents_Cladding3);

		DefineCladding3MaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("FibreCladding3Material", Material_Cladding3);
	}

//---------- absorbing coating ----------//
	if(CoatingExists)
	{
		G4double density_Coating;
		if(FibreProperties.containsNumber("density_coating")) density_Coating = FibreProperties.getNumber("density_coating");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre coating's density has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		GoddessProperties::tabular chemicalComponents_Coating;
		if(FibreProperties.containsTabular("chemical_components_coating")) chemicalComponents_Coating = FibreProperties.getTabular("chemical_components_coating");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineMaterials():" <<
			std::endl << "# The fibre coating's chemical components have not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		Material_Coating = new G4Material("tempName", density_Coating, (int) chemicalComponents_Coating["element"].size());
		PropertyTools->AddElementsFromTable(Material_Coating, chemicalComponents_Coating);

		DefineCoatingMaterialProperties();

		PropertyTools->checkIfMaterialAlreadyExists("FibreCoatingMaterial", Material_Coating);
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
 *  Function to define the following material properties of the optical cement (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Fibre::DefineOpticalCementMaterialProperties()
{
	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_opticalCement = PropertyTools->GetPropertyDistribution(OpticalCementProperties, "mu_att");
	if( !attenuationLength_opticalCement || PropertyTools->isConstantProperty(attenuationLength_opticalCement) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_opticalCement = PropertyTools->GetPropertyDistribution(OpticalCementProperties, "n_ref_imag");

		if(refractiveIndex_imag_opticalCement)
		{
			if( !attenuationLength_opticalCement || (PropertyTools->isConstantProperty(attenuationLength_opticalCement) && !PropertyTools->isConstantProperty(refractiveIndex_imag_opticalCement)) )
			{
				if(attenuationLength_opticalCement) delete attenuationLength_opticalCement;
				attenuationLength_opticalCement = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_opticalCement->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_opticalCement->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_opticalCement->Value(energy) );

					attenuationLength_opticalCement->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_opticalCement)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineOpticalCementMaterialProperties():" <<
		std::endl << "# The optical cement's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_OpticalCement = PropertyTools->GetPropertyDistribution(OpticalCementProperties, "n_ref");
	if(!refractiveIndex_OpticalCement)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineOpticalCementMaterialProperties():" <<
		std::endl << "# The optical cement's refractive index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * mpt_OpticalCement = new G4MaterialPropertiesTable();
	mpt_OpticalCement->AddProperty("RINDEX", refractiveIndex_OpticalCement);
	mpt_OpticalCement->AddProperty("ABSLENGTH", attenuationLength_opticalCement);
	Material_OpticalCement->SetMaterialPropertiesTable(mpt_OpticalCement);
}

/**
 *  Function to define the following material properties of the fibre's core (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 *  - emmission spectrum (for scintillating and wave-length-shifting fibres)
 *  - decay time (for scintillating and wave-length-shifting fibres)
 *  - rise time (for scintillating fibres)
 *  - scintillation yield (for scintillating fibres)
 *  - ratio between the fast and the slow scintillation process
 *  - intrinsic resolution factor
 */
void G4Fibre::DefineCoreMaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_FibreCore = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_core");
	if(!refractiveIndex_FibreCore)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
		std::endl << "# The fibre core's refractive index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_FibreCore = PropertyTools->GetPropertyDistribution(FibreProperties, "mu_att_core", false, FibreProperties.getNumber("mu_att_linear"), FibreProperties.getNumber("mu_att_quadratic"));
	if( !attenuationLength_FibreCore || PropertyTools->isConstantProperty(attenuationLength_FibreCore) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_FibreCore = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_imag_core");

		if(refractiveIndex_imag_FibreCore)
		{
			if( !attenuationLength_FibreCore || (PropertyTools->isConstantProperty(attenuationLength_FibreCore) && !PropertyTools->isConstantProperty(refractiveIndex_imag_FibreCore)) )
			{
				if(attenuationLength_FibreCore) delete attenuationLength_FibreCore;
				attenuationLength_FibreCore = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_FibreCore->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_FibreCore->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_FibreCore->Value(energy) );

					attenuationLength_FibreCore->InsertValues(energy, attLength);
				}
			}
		}
	}

	G4MaterialPropertyVector * wlsAttenuationLength_FibreCore = 0;
	if(IsWLS)
	{
		wlsAttenuationLength_FibreCore = PropertyTools->GetPropertyDistribution(FibreProperties, "mu_wls_att_core");

		if(!wlsAttenuationLength_FibreCore)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# The fibre core's WLS attenuation length has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else if(!attenuationLength_FibreCore)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
		std::endl << "# The fibre core's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}


	// for scintillating fibres
	G4double scintillationYield = 0.;
	G4double intrinsicResolutionFactor = -1;
	G4MaterialPropertyVector * fastLightOutput_Scinti = 0;
	G4MaterialPropertyVector * slowLightOutput_Scinti = 0;
	G4double fastDecayTime_Scinti = -1;
	G4double slowDecayTime_Scinti = -1;
	G4double fastRiseTime_Scinti = -1;
	G4double slowRiseTime_Scinti = -1;
	G4double fastDecayFraction_Scinti = -1.;
	// for WLS fibres
	G4MaterialPropertyVector * lightOutput_WLS = 0;
	G4double decayTime_WLS = 0.;
	if(IsScinti)
	{
		if(FibreProperties.containsNumber("lightOutput_rel") && FibreProperties.containsNumber("lightOutput"))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# An absolut AND a relative scintillation yield are given. The absolute one will be used." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			scintillationYield = FibreProperties.getNumber("lightOutput");
		}
		else if(FibreProperties.containsNumber("lightOutput")) scintillationYield = FibreProperties.getNumber("lightOutput");
		else if(FibreProperties.containsNumber("lightOutput_rel")) scintillationYield = FibreProperties.getNumber("lightOutput_rel") * 20000.;
		else std::cerr << std::endl << "Neither an absolut nor a relative scintillation yield is given." << std::endl << std::endl;
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# Neither an absolut nor a relative scintillation yield has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(FibreProperties.containsNumber("res_intr")) intrinsicResolutionFactor = FibreProperties.getNumber("res_intr");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# No intrinsic resolution factor of the scintillator is given. It will be set to 1." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			intrinsicResolutionFactor = 1.;
		}

		if(FibreProperties.containsTabular("epsilon_fast") && FibreProperties.containsTabular("epsilon"))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# An emission spectrum AND a fast emission spectrum of the scintillator are given. Which one is to be used?" <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
		else
		{
			fastLightOutput_Scinti = PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon");
			if(!fastLightOutput_Scinti) fastLightOutput_Scinti = PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon_fast");
		}

		slowLightOutput_Scinti = PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon_slow");

		if(!fastLightOutput_Scinti && !slowLightOutput_Scinti)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# The fibre core's scintillation emission spectrum has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(fastLightOutput_Scinti)
		{
			if(FibreProperties.containsNumber("t_decay_fast") && FibreProperties.containsNumber("t_decay"))
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# A decay time AND a fast decay time of the scintillator are given. Which one is to be used?" <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
			else if(FibreProperties.containsNumber("t_decay")) fastDecayTime_Scinti = FibreProperties.getNumber("t_decay");
			else if(FibreProperties.containsNumber("t_decay_fast")) fastDecayTime_Scinti = FibreProperties.getNumber("t_decay_fast");
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# The fibre core's (fast) scintillation decay time has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
		}

		if(slowLightOutput_Scinti)
		{
			if(FibreProperties.containsNumber("t_decay_slow")) slowDecayTime_Scinti = FibreProperties.getNumber("t_decay_slow");
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# The fibre core's slow scintillation decay time has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
		}

		if(fastLightOutput_Scinti)
		{
			if(FibreProperties.containsNumber("t_rise_fast") && FibreProperties.containsNumber("t_rise"))
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# A rise time AND a fast rise time of the scintillator are given. Which one is to be used?" <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
			else if(FibreProperties.containsNumber("t_rise")) fastRiseTime_Scinti = FibreProperties.getNumber("t_rise");
			else if(FibreProperties.containsNumber("t_rise_fast")) fastRiseTime_Scinti = FibreProperties.getNumber("t_rise_fast");
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# The fibre core's (fast) scintillation rise time has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
		}

		if(slowLightOutput_Scinti)
		{
			if(FibreProperties.containsNumber("t_rise_slow")) slowRiseTime_Scinti = FibreProperties.getNumber("t_rise_slow");
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# The fibre core's slow scintillation rise time has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
		}

		if(fastLightOutput_Scinti && slowLightOutput_Scinti)
		{
			if(FibreProperties.containsNumber("fraction_fast")) fastDecayFraction_Scinti = FibreProperties.getNumber("fraction_fast");
			else
			{
				std::cerr <<
				std::endl << "##########" <<
				std::endl << "# ERROR:" <<
				std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
				std::endl << "# The fibre core's fraction of fast scintillation decay has not been specified." <<
				std::endl << "##########" <<
				std::endl << std::endl;

				CriticalErrorOccured = true;
			}
		}
	}
	else if(IsWLS)
	{
		lightOutput_WLS = PropertyTools->GetPropertyDistribution(FibreProperties, "epsilon");
		if(!lightOutput_WLS)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# The fibre core's WLS emission spectrum has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}

		if(FibreProperties.containsNumber("t_decay")) decayTime_WLS = FibreProperties.getNumber("t_decay");
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineCoreMaterialProperties():" <<
			std::endl << "# The fibre core's WLS decay time has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}

	G4MaterialPropertiesTable * mpt_FibreCore = new G4MaterialPropertiesTable();
	mpt_FibreCore->AddProperty("RINDEX", refractiveIndex_FibreCore);
	if(IsScinti)
	{
		mpt_FibreCore->AddProperty("ABSLENGTH", attenuationLength_FibreCore);
		if(fastLightOutput_Scinti) mpt_FibreCore->AddProperty("FASTCOMPONENT", fastLightOutput_Scinti);
		if(slowLightOutput_Scinti) mpt_FibreCore->AddProperty("SLOWCOMPONENT", slowLightOutput_Scinti);
		mpt_FibreCore->AddConstProperty("SCINTILLATIONYIELD", scintillationYield);
		mpt_FibreCore->AddConstProperty("RESOLUTIONSCALE", intrinsicResolutionFactor);
		if(fastLightOutput_Scinti && fastDecayTime_Scinti) mpt_FibreCore->AddConstProperty("FASTTIMECONSTANT", fastDecayTime_Scinti);
		if(slowLightOutput_Scinti && slowDecayTime_Scinti) mpt_FibreCore->AddConstProperty("SLOWTIMECONSTANT", slowDecayTime_Scinti);
		if(fastLightOutput_Scinti && fastRiseTime_Scinti) mpt_FibreCore->AddConstProperty("FASTSCINTILLATIONRISETIME", fastRiseTime_Scinti);
		if(slowLightOutput_Scinti && slowRiseTime_Scinti) mpt_FibreCore->AddConstProperty("SLOWSCINTILLATIONRISETIME", slowRiseTime_Scinti);
		if(fastLightOutput_Scinti && slowLightOutput_Scinti) mpt_FibreCore->AddConstProperty("YIELDRATIO", fastDecayFraction_Scinti);
	}
	else if(IsWLS)
	{
		if(attenuationLength_FibreCore) mpt_FibreCore->AddProperty("ABSLENGTH", attenuationLength_FibreCore);
		mpt_FibreCore->AddProperty("WLSABSLENGTH", wlsAttenuationLength_FibreCore);
		mpt_FibreCore->AddProperty("WLSCOMPONENT", lightOutput_WLS);
		mpt_FibreCore->AddConstProperty("WLSTIMECONSTANT", decayTime_WLS);
	}
	else mpt_FibreCore->AddProperty("ABSLENGTH", attenuationLength_FibreCore);
	Material_FibreCore->SetMaterialPropertiesTable(mpt_FibreCore);
}

/**
 *  Function to define the following material properties of the fibre's 1. cladding (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Fibre::DefineCladding1MaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Cladding1 = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_cladding_1");
	if(!refractiveIndex_Cladding1)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCladding1MaterialProperties():" <<
		std::endl << "# The fibre first cladding's refractive index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_Cladding1 = PropertyTools->GetPropertyDistribution(FibreProperties, "mu_att_cladding_1", false, FibreProperties.getNumber("mu_att_linear"), FibreProperties.getNumber("mu_att_quadratic"));
	if( !attenuationLength_Cladding1 || PropertyTools->isConstantProperty(attenuationLength_Cladding1) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_Cladding1 = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_imag_cladding_1");

		if(refractiveIndex_imag_Cladding1)
		{
			if( !attenuationLength_Cladding1 || (PropertyTools->isConstantProperty(attenuationLength_Cladding1) && !PropertyTools->isConstantProperty(refractiveIndex_imag_Cladding1)) )
			{
				if(attenuationLength_Cladding1) delete attenuationLength_Cladding1;
				attenuationLength_Cladding1 = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_Cladding1->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_Cladding1->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_Cladding1->Value(energy) );

					attenuationLength_Cladding1->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_Cladding1)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCladding1MaterialProperties():" <<
		std::endl << "# The fibre first cladding's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * mpt_Cladding1 = new G4MaterialPropertiesTable();
	mpt_Cladding1->AddProperty("RINDEX", refractiveIndex_Cladding1);
	mpt_Cladding1->AddProperty("ABSLENGTH", attenuationLength_Cladding1);
	Material_Cladding1->SetMaterialPropertiesTable(mpt_Cladding1);
}

/**
 *  Function to define the following material properties of the fibre's 2. cladding (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Fibre::DefineCladding2MaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Cladding2 = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_cladding_2");
	if(!refractiveIndex_Cladding2)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCladding2MaterialProperties():" <<
		std::endl << "# The fibre second cladding's refractive index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_Cladding2 = PropertyTools->GetPropertyDistribution(FibreProperties, "mu_att_cladding_2", false, FibreProperties.getNumber("mu_att_linear"), FibreProperties.getNumber("mu_att_quadratic"));
	if( !attenuationLength_Cladding2 || PropertyTools->isConstantProperty(attenuationLength_Cladding2) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_Cladding2 = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_imag_cladding_2");

		if(refractiveIndex_imag_Cladding2)
		{
			if( !attenuationLength_Cladding2 || (PropertyTools->isConstantProperty(attenuationLength_Cladding2) && !PropertyTools->isConstantProperty(refractiveIndex_imag_Cladding2)) )
			{
				if(attenuationLength_Cladding2) delete attenuationLength_Cladding2;
				attenuationLength_Cladding2 = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_Cladding2->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_Cladding2->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_Cladding2->Value(energy) );

					attenuationLength_Cladding2->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_Cladding2)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCladding2MaterialProperties():" <<
		std::endl << "# The fibre second cladding's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * mpt_Cladding2 = new G4MaterialPropertiesTable();
	mpt_Cladding2->AddProperty("RINDEX", refractiveIndex_Cladding2);
	mpt_Cladding2->AddProperty("ABSLENGTH", attenuationLength_Cladding2);
	Material_Cladding2->SetMaterialPropertiesTable(mpt_Cladding2);
}

/**
 *  Function to define the following material properties of the fibre's 3. cladding (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Fibre::DefineCladding3MaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Cladding3 = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_cladding_3");
	if(!refractiveIndex_Cladding3)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCladding3MaterialProperties():" <<
		std::endl << "# The fibre third cladding's refractive index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	// by now (GEANT 4.9.5), AddConstProperty("ABSLENGTH", G4double) DOES NOT WORK (the program runs, but the attenuation length does not seem to have any impact)!
	G4MaterialPropertyVector * attenuationLength_Cladding3 = PropertyTools->GetPropertyDistribution(FibreProperties, "mu_att_cladding_3", false, FibreProperties.getNumber("mu_att_linear"), FibreProperties.getNumber("mu_att_quadratic"));
	if( !attenuationLength_Cladding3 || PropertyTools->isConstantProperty(attenuationLength_Cladding3) )
	{
		// calculating the attenuation length spectrum from the given complex refractive index spectrum
		G4MaterialPropertyVector * refractiveIndex_imag_Cladding3 = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_imag_cladding_3");

		if(refractiveIndex_imag_Cladding3)
		{
			if( !attenuationLength_Cladding3 || (PropertyTools->isConstantProperty(attenuationLength_Cladding3) && !PropertyTools->isConstantProperty(refractiveIndex_imag_Cladding3)) )
			{
				if(attenuationLength_Cladding3) delete attenuationLength_Cladding3;
				attenuationLength_Cladding3 = new G4MaterialPropertyVector();
				for(unsigned int i = 0; i < refractiveIndex_imag_Cladding3->GetVectorLength(); i++)
				{
					G4double energy = refractiveIndex_imag_Cladding3->Energy(i);
					G4double attLength = CLHEP::c_light * CLHEP::h_Planck / (2. * CLHEP::pi * energy * refractiveIndex_imag_Cladding3->Value(energy) );

					attenuationLength_Cladding3->InsertValues(energy, attLength);
				}
			}
		}
	}

	if(!attenuationLength_Cladding3)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCladding3MaterialProperties():" <<
		std::endl << "# The fibre third cladding's attenuation length has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * mpt_Cladding3 = new G4MaterialPropertiesTable();
	mpt_Cladding3->AddProperty("RINDEX", refractiveIndex_Cladding3);
	mpt_Cladding3->AddProperty("ABSLENGTH", attenuationLength_Cladding3);
	Material_Cladding3->SetMaterialPropertiesTable(mpt_Cladding3);
}

/**
 *  Function to define the following material properties of the fibre's coating (according to the specifications from the property file):
 *  - refractive index spectrum
 *  - attenuation length / absorption spectrum
 */
void G4Fibre::DefineCoatingMaterialProperties()
{
	// by now (GEANT 4.9.5), GEANT expects that a material property with the name identifier RINDEX is wavelength/energy dependent (http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/379.html)!
	G4MaterialPropertyVector * refractiveIndex_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_coating");
	if(!refractiveIndex_Coating)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# ERROR:" <<
		std::endl << "# G4Fibre::DefineCoatingMaterialProperties():" <<
		std::endl << "# The fibre coating's refractive index has not been specified." <<
		std::endl << "##########" <<
		std::endl << std::endl;

		CriticalErrorOccured = true;
	}

	G4MaterialPropertiesTable * mpt_Coating = new G4MaterialPropertiesTable();
	mpt_Coating->AddProperty("RINDEX", refractiveIndex_Coating);
	Material_Coating->SetMaterialPropertiesTable(mpt_Coating);
}



/**
 *  Function to construct surfaces between different volumes and define their mechanical and optical properties (according to the specifications from the property file(s)). \n
 *  The following properties are defined:
 *  - surface type ("dielectric_metal" or "dielectric_dielectric")
 *  - surface roughness
 *  - reflectivity
 */
void G4Fibre::ConstructSurface()
{
//NOTE: A surface has to be defined for the borders between every two volumes.
//      It is needed to simulate optical boundary processes.
//      Exclusively in case of a perfectly smooth surface between two dielectic materials
//      (and only refractive indices needed to describe), no surface has to be defined.
//
//      In order to define different surface properties for different borders of the same volumes,
//      one has to define adjacent volumes in the simulation that butt up with the various sides
//      and are made of the same physical material as the mother.


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


// fibre core
	G4double fibreSigmaAlpha_Core = 0.;

	if(FibreProperties.containsNumber("roughness_core"))
	{
		fibreSigmaAlpha_Core = FibreProperties.getNumber("roughness_core");
	}
	else if(!Cladding1Exists)
	{
		std::cerr <<
		std::endl << "##########" <<
		std::endl << "# Warning:" <<
		std::endl << "# G4Fibre::ConstructSurface():" <<
// 		std::endl << "# The fibre core's roughness has not been specified. Roughness set to 0." <<
		std::endl << "# The roughness of the fibre's outer surface has not been specified. Roughness set to 0." <<
		std::endl << "##########" <<
		std::endl << std::endl;
//		CriticalErrorOccured = true;
	}

	G4OpticalSurface * OptSurf_Fibre_Core = 0;
	if(fibreSigmaAlpha_Core >= 1e-12)
	{
		OptSurf_Fibre_Core = new G4OpticalSurface("FibreCoreSurface", unified, ground, dielectric_dielectric, fibreSigmaAlpha_Core);
	}
	else
	{
		OptSurf_Fibre_Core = new G4OpticalSurface("FibreCoreSurface", unified, polished, dielectric_dielectric, 0.);
	}

	G4MaterialPropertiesTable *MPT_OptSurf_Fibre_Core = new G4MaterialPropertiesTable();
	MPT_OptSurf_Fibre_Core->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(1.));
	MPT_OptSurf_Fibre_Core->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(0.));
	MPT_OptSurf_Fibre_Core->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
	OptSurf_Fibre_Core->SetMaterialPropertiesTable(MPT_OptSurf_Fibre_Core);

	if(Cladding1Exists)
	{
		for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
		{
			// specify between which volumes the previously defined surface is to be applied	TODO
			// the first volume is where the photon comes from, the second volume is where the photon heads for.
			new G4LogicalBorderSurface("roughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", FibreCore_physical[counter], Cladding1_physical[counter], OptSurf_Fibre_Core);
			new G4LogicalBorderSurface("roughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", Cladding1_physical[counter], FibreCore_physical[counter], OptSurf_Fibre_Core);
		}
	}
	else
	{
		for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
		{
			new G4LogicalSkinSurface("roughened Surface (around " + FibreCore_physical[counter]->GetName() + ")", FibreCore_physical[counter]->GetLogicalVolume(), OptSurf_Fibre_Core);
		}
	}

// first (inner) cladding
	if(Cladding1Exists)
	{
		G4double fibreSigmaAlpha_Cladding1 = 0.;

		if(FibreProperties.containsNumber("roughness_cladding_1"))
		{
			fibreSigmaAlpha_Cladding1 = FibreProperties.getNumber("roughness_cladding_1");
		}
		else if(!Cladding2Exists)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# Warning:" <<
			std::endl << "# G4Fibre::ConstructSurface():" <<
// 			std::endl << "# The roughness of the first cladding has not been specified. Roughness set to 0." <<
			std::endl << "# The roughness of the fibre's outer surface has not been specified. Roughness set to 0." <<
			std::endl << "##########" <<
			std::endl << std::endl;
	//		CriticalErrorOccured = true;
		}

		G4OpticalSurface * OptSurf_Fibre_Cladding_1 = 0;
		if(fibreSigmaAlpha_Cladding1 >= 1e-12)
		{
			OptSurf_Fibre_Cladding_1 = new G4OpticalSurface("FibreCladding1Surface", unified, ground, dielectric_dielectric, fibreSigmaAlpha_Cladding1);
		}
		else
		{
			OptSurf_Fibre_Cladding_1 = new G4OpticalSurface("FibreCladding1Surface", unified, polished, dielectric_dielectric, 0.);
		}

		G4MaterialPropertiesTable *MPT_OptSurf_Fibre_Cladding_1 = new G4MaterialPropertiesTable();
		MPT_OptSurf_Fibre_Cladding_1->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(1.));
		MPT_OptSurf_Fibre_Cladding_1->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(0.));
		MPT_OptSurf_Fibre_Cladding_1->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
		OptSurf_Fibre_Cladding_1->SetMaterialPropertiesTable(MPT_OptSurf_Fibre_Cladding_1);

		if(Cladding2Exists)
		{
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				// specify between which volumes the previously defined surface is to be applied
				// the first volume is where the photon comes from, the second volume is where the photon heads for.
				new G4LogicalBorderSurface("roughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", Cladding1_physical[counter], Cladding2_physical[counter], OptSurf_Fibre_Cladding_1);
				new G4LogicalBorderSurface("roughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", Cladding2_physical[counter], Cladding1_physical[counter], OptSurf_Fibre_Cladding_1);
			}
		}
		else
		{
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalSkinSurface("roughened Surface (around " + Cladding1_physical[counter]->GetName() + ")", Cladding1_physical[counter]->GetLogicalVolume(), OptSurf_Fibre_Cladding_1);
			}
		}
	}

// second cladding
	if(Cladding2Exists)
	{
		G4double fibreSigmaAlpha_Cladding2 = 0.;

		if(FibreProperties.containsNumber("roughness_cladding_2"))
		{
			fibreSigmaAlpha_Cladding2 = FibreProperties.getNumber("roughness_cladding_2");
		}
		else if(!Cladding3Exists)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# Warning:" <<
			std::endl << "# G4Fibre::ConstructSurface():" <<
// 			std::endl << "# The roughness of the second cladding has not been specified. Roughness set to 0." <<
			std::endl << "# The roughness of the fibre's outer surface has not been specified. Roughness set to 0." <<
			std::endl << "##########" <<
			std::endl << std::endl;
	//		CriticalErrorOccured = true;
		}

		G4OpticalSurface *OptSurf_Fibre_Cladding_2 = 0;
		if(fibreSigmaAlpha_Cladding2 >= 1e-12)
		{
			OptSurf_Fibre_Cladding_2 = new G4OpticalSurface("FibreCladding2Surface", unified, ground, dielectric_dielectric, fibreSigmaAlpha_Cladding2);
		}
		else
		{
			OptSurf_Fibre_Cladding_2 = new G4OpticalSurface("FibreCladding2Surface", unified, polished, dielectric_dielectric, 0.);
		}

		G4MaterialPropertiesTable *MPT_OptSurf_Fibre_Cladding_2 = new G4MaterialPropertiesTable();
		MPT_OptSurf_Fibre_Cladding_2->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(1.));
		MPT_OptSurf_Fibre_Cladding_2->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(0.));
		MPT_OptSurf_Fibre_Cladding_2->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
		OptSurf_Fibre_Cladding_2->SetMaterialPropertiesTable(MPT_OptSurf_Fibre_Cladding_2);

		if(Cladding3Exists)
		{
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				// specify between which volumes the previously defined surface is to be applied
				// the first volume is where the photon comes from, the second volume is where the photon heads for.
				new G4LogicalBorderSurface("roughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", Cladding2_physical[counter], Cladding3_physical[counter], OptSurf_Fibre_Cladding_2);
				new G4LogicalBorderSurface("roughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", Cladding3_physical[counter], Cladding2_physical[counter], OptSurf_Fibre_Cladding_2);
			}
		}
		else
		{
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalSkinSurface("roughened Surface (around " + Cladding2_physical[counter]->GetName() + ")", Cladding2_physical[counter]->GetLogicalVolume(), OptSurf_Fibre_Cladding_2);
			}
		}
	}

// third (outer) cladding
	if(Cladding3Exists)
	{
		G4double fibreSigmaAlpha_Cladding3 = 0.;

		if(FibreProperties.containsNumber("roughness_cladding_3"))
		{
			fibreSigmaAlpha_Cladding3 = FibreProperties.getNumber("roughness_cladding_3");
		}
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# Warning:" <<
			std::endl << "# G4Fibre::ConstructSurface():" <<
// 			std::endl << "# The roughness of the third cladding has not been specified. Roughness set to 0." <<
			std::endl << "# The roughness of the fibre's outer surface has not been specified. Roughness set to 0." <<
			std::endl << "##########" <<
			std::endl << std::endl;
	//		CriticalErrorOccured = true;
		}
	}

	if(CoatingExists)
	{
		G4MaterialPropertyVector * reflectivity_Coating = 0;
		G4MaterialPropertyVector * refractiveIndex_real_Coating = 0;
		G4MaterialPropertyVector * refractiveIndex_imag_Coating = 0;

		if(FibreProperties.containsNumber("reflec_coating"))
		{
			reflectivity_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "reflec_coating");
		}
		else if(FibreProperties.containsNumber("n_ref_real_coating") && FibreProperties.containsNumber("n_ref_imag_coating"))
		{
			refractiveIndex_real_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_real_coating");
			refractiveIndex_imag_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "n_ref_imag_coating");
		}
		else
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::ConstructSurface():" <<
			std::endl << "# The reflectivity of the fibre's coating has not been specified. Reflectivity set to 0 (100% absorption)." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			reflectivity_Coating = PropertyTools->GetPropertyDistribution(0.);
		}

		G4MaterialPropertyVector * specularspike_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "fraction_specularspike_coating");
		if(!specularspike_Coating)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::ConstructWrappingSurface():" <<
			std::endl << "# The fibre coating's fraction of specular-spike-reflection has not been specified. Is set to 0." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			specularspike_Coating = PropertyTools->GetPropertyDistribution(0.);
		}

		G4MaterialPropertyVector * specularlobe_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "fraction_specularlobe_coating");
		if(!specularlobe_Coating)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::ConstructWrappingSurface():" <<
			std::endl << "# The fibre coating's wrapping layer's fraction of specular-lobe-reflection has not been specified. Is set to 1." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			specularlobe_Coating = PropertyTools->GetPropertyDistribution(1.);
		}

		G4MaterialPropertyVector * backscattering_Coating = PropertyTools->GetPropertyDistribution(FibreProperties, "fraction_backscattering_coating");
		if(!backscattering_Coating)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::ConstructWrappingSurface():" <<
			std::endl << "# The fibre coating's wrapping layer's fraction of backscattering-reflection has not been specified. Is set to 0." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			backscattering_Coating = PropertyTools->GetPropertyDistribution(0.);
		}

		G4double fibreSigmaAlpha_Coating = FibreProperties.getNumber("roughness_coating");
		if(isnan(fibreSigmaAlpha_Coating))
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# WARNING:" <<
			std::endl << "# G4Fibre::ConstructSurface():" <<
			std::endl << "# The roughness of the fibre's coating has not been specified. Roughness set to 0." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			fibreSigmaAlpha_Coating = 0.;
		}


		G4OpticalSurface *OptSurf_Fibre_Coating = new G4OpticalSurface("FibreCoatingSurface", unified, ground, dielectric_metal, fibreSigmaAlpha_Coating);

		G4MaterialPropertiesTable *MPT_OptSurf_Fibre_Coating = new G4MaterialPropertiesTable();
		if(reflectivity_Coating) MPT_OptSurf_Fibre_Coating->AddProperty("REFLECTIVITY", reflectivity_Coating);
		if(refractiveIndex_real_Coating && refractiveIndex_imag_Coating)
		{
			MPT_OptSurf_Fibre_Coating->AddProperty("REALRINDEX", refractiveIndex_real_Coating);
			MPT_OptSurf_Fibre_Coating->AddProperty("IMAGINARYRINDEX", refractiveIndex_imag_Coating);
		}
		MPT_OptSurf_Fibre_Coating->AddProperty("SPECULARLOBECONSTANT", specularspike_Coating);
		MPT_OptSurf_Fibre_Coating->AddProperty("SPECULARSPIKECONSTANT", specularlobe_Coating);
		MPT_OptSurf_Fibre_Coating->AddProperty("BACKSCATTERCONSTANT", backscattering_Coating);
		OptSurf_Fibre_Coating->SetMaterialPropertiesTable(MPT_OptSurf_Fibre_Coating);


		//NOTE: A surface has to be defined for the borders between every two volumes.
		//      It is needed to simulate optical boundary processes.
		//      Exclusively in case of a perfectly smooth surface between two dielectic materials
		//      (and only refractive indices needed to describe), no surface has to be defined.
		//
		//      In order to define different surface properties for different borders of the same volumes,
		//      one has to define adjacent volumes in the simulation that butt up with the various sides
		//      and are made of the same physical material as the mother.

		// specify between which volumes the previously defined surface should be applied
		// the first volume is where the photon comes from, the second volume is where the photon heads for.

		for(unsigned int counter = 0; counter < Coating_physical.size(); counter++)
		{
			new G4LogicalBorderSurface("Coating Surface (between " + Coating_physical[counter]->GetName() + " and " + Coating_mothers_physical[counter]->GetName() + ")", Coating_physical[counter], Coating_mothers_physical[counter], OptSurf_Fibre_Coating);
			new G4LogicalBorderSurface("Coating Surface (between " + Coating_mothers_physical[counter]->GetName() + " and " + Coating_physical[counter]->GetName() + ")", Coating_mothers_physical[counter], Coating_physical[counter], OptSurf_Fibre_Coating);

			if(Cladding3Exists)
			{
				new G4LogicalBorderSurface("Coating Surface (between " + Coating_physical[counter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", Coating_physical[counter], Cladding3_physical[counter], OptSurf_Fibre_Coating);
				new G4LogicalBorderSurface("Coating Surface (between " + Cladding3_physical[counter]->GetName() + " and " + Coating_physical[counter]->GetName() + ")", Cladding3_physical[counter], Coating_physical[counter], OptSurf_Fibre_Coating);
			}
			else if(Cladding2Exists)
			{
				new G4LogicalBorderSurface("Coating Surface (between " + Coating_physical[counter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", Coating_physical[counter], Cladding2_physical[counter], OptSurf_Fibre_Coating);
				new G4LogicalBorderSurface("Coating Surface (between " + Cladding2_physical[counter]->GetName() + " and " + Coating_physical[counter]->GetName() + ")", Cladding2_physical[counter], Coating_physical[counter], OptSurf_Fibre_Coating);
			}
			else if(Cladding1Exists)
			{
				new G4LogicalBorderSurface("Coating Surface (between " + Coating_physical[counter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", Coating_physical[counter], Cladding1_physical[counter], OptSurf_Fibre_Coating);
				new G4LogicalBorderSurface("Coating Surface (between " + Cladding1_physical[counter]->GetName() + " and " + Coating_physical[counter]->GetName() + ")", Cladding1_physical[counter], Coating_physical[counter], OptSurf_Fibre_Coating);
			}
			else
			{
				new G4LogicalBorderSurface("Coating Surface (between " + Coating_physical[counter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", Coating_physical[counter], FibreCore_physical[counter], OptSurf_Fibre_Coating);
				new G4LogicalBorderSurface("Coating Surface (between " + FibreCore_physical[counter]->GetName() + " and " + Coating_physical[counter]->GetName() + ")", FibreCore_physical[counter], Coating_physical[counter], OptSurf_Fibre_Coating);
			}
		}
	}

// fibre start & end points
	G4MaterialPropertyVector * reflectivity = 0;

	// perfectly smooth surface for in between fibre parts
	OptSurf_perfectlySmooth = new G4OpticalSurface("perfectly smooth surface", unified, polished, dielectric_dielectric, 0.);

	G4MaterialPropertiesTable * MPT_OptSurf_perfectlySmooth = new G4MaterialPropertiesTable();
	MPT_OptSurf_perfectlySmooth->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(0.));
	MPT_OptSurf_perfectlySmooth->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(1.));
	MPT_OptSurf_perfectlySmooth->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
	OptSurf_perfectlySmooth->SetMaterialPropertiesTable(MPT_OptSurf_perfectlySmooth);


// fibre start point
	if(CreateReflectiveStartPointVolumes)
	{
		OptSurf_startPoint = new G4OpticalSurface("reflective start point surface", unified, polished, dielectric_metal, 0.);

		reflectivity = PropertyTools->GetPropertyDistribution(Reflectivity_startPoint);
		if(!reflectivity)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineReflectiveSurfacesProperties():" <<
			std::endl << "# The fibre's start point reflectivity has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else
	{
		if(Roughness_startPoint >= 1e-12)
		{
			OptSurf_startPoint = new G4OpticalSurface("roughened start point surface", unified, ground, dielectric_dielectric, Roughness_startPoint);
		}
		else
		{
			OptSurf_startPoint = new G4OpticalSurface("smooth start point surface", unified, polished, dielectric_dielectric, 0.);
		}
	}

	G4MaterialPropertiesTable * MPT_OptSurf_startPoint = new G4MaterialPropertiesTable();
	if(CreateReflectiveStartPointVolumes) MPT_OptSurf_startPoint->AddProperty("REFLECTIVITY", reflectivity);
	if(!CreateReflectiveStartPointVolumes && Roughness_startPoint >= 1e-12)
	{
		MPT_OptSurf_startPoint->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(1.));
		MPT_OptSurf_startPoint->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(0.));
	}
	else
	{
		MPT_OptSurf_startPoint->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(0.));
		MPT_OptSurf_startPoint->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(1.));
	}
	MPT_OptSurf_startPoint->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
	OptSurf_startPoint->SetMaterialPropertiesTable(MPT_OptSurf_startPoint);


	if(CreateReflectiveStartPointVolumes)
	{
		for(unsigned int iter = 0; iter < ReflectiveStartPoint_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("reflective start point Surface (around " + ReflectiveStartPoint_physical[iter]->GetName() + ")", ReflectiveStartPoint_physical[iter]->GetLogicalVolume(), OptSurf_startPoint);

			for(unsigned int counter = 0; counter < Cladding3_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveStartPoint_physical[iter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", ReflectiveStartPoint_physical[iter], Cladding3_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + ReflectiveStartPoint_physical[iter]->GetName() + ")", Cladding3_physical[counter], ReflectiveStartPoint_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveStartPoint_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", ReflectiveStartPoint_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + ReflectiveStartPoint_physical[iter]->GetName() + ")", Cladding2_physical[counter], ReflectiveStartPoint_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveStartPoint_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", ReflectiveStartPoint_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + ReflectiveStartPoint_physical[iter]->GetName() + ")", Cladding1_physical[counter], ReflectiveStartPoint_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveStartPoint_physical[iter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", ReflectiveStartPoint_physical[iter], FibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + ReflectiveStartPoint_physical[iter]->GetName() + ")", FibreCore_physical[counter], ReflectiveStartPoint_physical[iter], OptSurf_perfectlySmooth);
			}
		}
	}
	else
	{
		// roughened SkinSurface for all parts of the roughened fibre part
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding3_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened start point Surface (around " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding3_physical[iter]->GetLogicalVolume(), OptSurf_startPoint);
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding2_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened start point Surface (around " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter]->GetLogicalVolume(), OptSurf_startPoint);
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding1_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened start point Surface (around " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter]->GetLogicalVolume(), OptSurf_startPoint);
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_fibreCore_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened start point Surface (around " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + ")", RoughenedStartPoint_fibreCore_physical[iter]->GetLogicalVolume(), OptSurf_startPoint);
		}

		// perfectly smooth BorderSurface between the different parts of the roughened fibre layer
		if(RoughenedStartPoint_cladding3_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding3_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedStartPoint_cladding3_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + " and " + RoughenedStartPoint_cladding3_physical[iter + counter]->GetName() + ")", RoughenedStartPoint_cladding3_physical[iter], RoughenedStartPoint_cladding3_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding3_physical[iter + counter]->GetName() + " and " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding3_physical[iter + counter], RoughenedStartPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}
		if(RoughenedStartPoint_cladding2_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding2_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedStartPoint_cladding2_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[iter + counter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter], RoughenedStartPoint_cladding2_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[iter + counter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter + counter], RoughenedStartPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}
		if(RoughenedStartPoint_cladding1_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding1_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedStartPoint_cladding1_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[iter + counter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter], RoughenedStartPoint_cladding1_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[iter + counter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter + counter], RoughenedStartPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}
		if(RoughenedStartPoint_fibreCore_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedStartPoint_fibreCore_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedStartPoint_fibreCore_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + " and " + RoughenedStartPoint_fibreCore_physical[iter + counter]->GetName() + ")", RoughenedStartPoint_fibreCore_physical[iter], RoughenedStartPoint_fibreCore_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_fibreCore_physical[iter + counter]->GetName() + " and " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + ")", RoughenedStartPoint_fibreCore_physical[iter + counter], RoughenedStartPoint_fibreCore_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}

		// perfectly smooth BorderSurface between the different layers of the roughened fibre part
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding3_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < RoughenedStartPoint_cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding3_physical[iter], RoughenedStartPoint_cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[counter], RoughenedStartPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding2_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < RoughenedStartPoint_cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter], RoughenedStartPoint_cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[counter], RoughenedStartPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding1_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < RoughenedStartPoint_fibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + " and " + RoughenedStartPoint_fibreCore_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter], RoughenedStartPoint_fibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_fibreCore_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + ")", RoughenedStartPoint_fibreCore_physical[counter], RoughenedStartPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
		}

		// perfectly smooth BorderSurface to all other fibre parts
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding3_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding3_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding3_physical[iter], Cladding3_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + ")", Cladding3_physical[counter], RoughenedStartPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding3_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding3_physical[iter]->GetName() + ")", Cladding2_physical[counter], RoughenedStartPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding2_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding3_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter], Cladding3_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + ")", Cladding3_physical[counter], RoughenedStartPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + ")", Cladding2_physical[counter], RoughenedStartPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding2_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding2_physical[iter]->GetName() + ")", Cladding1_physical[counter], RoughenedStartPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_cladding1_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + ")", Cladding2_physical[counter], RoughenedStartPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + ")", Cladding1_physical[counter], RoughenedStartPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", RoughenedStartPoint_cladding1_physical[iter], FibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + RoughenedStartPoint_cladding1_physical[iter]->GetName() + ")", FibreCore_physical[counter], RoughenedStartPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedStartPoint_fibreCore_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", RoughenedStartPoint_fibreCore_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + ")", Cladding1_physical[counter], RoughenedStartPoint_fibreCore_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", RoughenedStartPoint_fibreCore_physical[iter], FibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + RoughenedStartPoint_fibreCore_physical[iter]->GetName() + ")", FibreCore_physical[counter], RoughenedStartPoint_fibreCore_physical[iter], OptSurf_perfectlySmooth);
			}
		}
	}


// fibre end point
	if(CreateReflectiveEndPointVolumes)
	{
		OptSurf_endPoint = new G4OpticalSurface("reflective end point surface", unified, polished, dielectric_metal, 0.);

		reflectivity = PropertyTools->GetPropertyDistribution(Reflectivity_endPoint);
		if(!reflectivity)
		{
			std::cerr <<
			std::endl << "##########" <<
			std::endl << "# ERROR:" <<
			std::endl << "# G4Fibre::DefineReflectiveSurfacesProperties():" <<
			std::endl << "# The fibre's end point reflectivity has not been specified." <<
			std::endl << "##########" <<
			std::endl << std::endl;

			CriticalErrorOccured = true;
		}
	}
	else
	{
		if(Roughness_endPoint >= 1e-12)
		{
			OptSurf_endPoint = new G4OpticalSurface("roughened end point surface", unified, ground, dielectric_dielectric, Roughness_endPoint);
		}
		else
		{
			OptSurf_endPoint = new G4OpticalSurface("smooth end point surface", unified, polished, dielectric_dielectric, 0.);
		}
	}

	G4MaterialPropertiesTable * MPT_OptSurf_endPoint = new G4MaterialPropertiesTable();
	if(CreateReflectiveEndPointVolumes) MPT_OptSurf_endPoint->AddProperty("REFLECTIVITY", reflectivity);
	if(!CreateReflectiveEndPointVolumes && Roughness_endPoint >= 1e-12)
	{
		MPT_OptSurf_endPoint->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(1.));
		MPT_OptSurf_endPoint->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(0.));
	}
	else
	{
		MPT_OptSurf_endPoint->AddProperty("SPECULARLOBECONSTANT", PropertyTools->GetPropertyDistribution(0.));
		MPT_OptSurf_endPoint->AddProperty("SPECULARSPIKECONSTANT", PropertyTools->GetPropertyDistribution(1.));
	}
	MPT_OptSurf_endPoint->AddProperty("BACKSCATTERCONSTANT", PropertyTools->GetPropertyDistribution(0.));
	OptSurf_endPoint->SetMaterialPropertiesTable(MPT_OptSurf_endPoint);


	if(CreateReflectiveEndPointVolumes)
	{
		for(unsigned int iter = 0; iter < ReflectiveEndPoint_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("reflective start point Surface (around " + ReflectiveEndPoint_physical[iter]->GetName() + ")", ReflectiveEndPoint_physical[iter]->GetLogicalVolume(), OptSurf_endPoint);

			if(Cladding3Exists)
			{
				for(unsigned int counter = 0; counter < Cladding3_physical.size(); counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveEndPoint_physical[iter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", ReflectiveEndPoint_physical[iter], Cladding3_physical[counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + ReflectiveEndPoint_physical[iter]->GetName() + ")", Cladding3_physical[counter], ReflectiveEndPoint_physical[iter], OptSurf_perfectlySmooth);
				}
			}
			if(Cladding2Exists)
			{
				for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveEndPoint_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", ReflectiveEndPoint_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + ReflectiveEndPoint_physical[iter]->GetName() + ")", Cladding2_physical[counter], ReflectiveEndPoint_physical[iter], OptSurf_perfectlySmooth);
				}
			}
			if(Cladding1Exists)
			{
				for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveEndPoint_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", ReflectiveEndPoint_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + ReflectiveEndPoint_physical[iter]->GetName() + ")", Cladding1_physical[counter], ReflectiveEndPoint_physical[iter], OptSurf_perfectlySmooth);
				}
			}
			for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + ReflectiveEndPoint_physical[iter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", ReflectiveEndPoint_physical[iter], FibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + ReflectiveEndPoint_physical[iter]->GetName() + ")", FibreCore_physical[counter], ReflectiveEndPoint_physical[iter], OptSurf_perfectlySmooth);
			}
		}
	}
	else
	{
		// roughened SkinSurface for all parts of the roughened fibre part
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding3_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened Surface (around " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding3_physical[iter]->GetLogicalVolume(), OptSurf_endPoint);
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding2_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened Surface (around " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter]->GetLogicalVolume(), OptSurf_endPoint);
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding1_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened Surface (around " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter]->GetLogicalVolume(), OptSurf_endPoint);
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_fibreCore_physical.size(); iter++)
		{
			new G4LogicalSkinSurface("roughened Surface (around " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + ")", RoughenedEndPoint_fibreCore_physical[iter]->GetLogicalVolume(), OptSurf_endPoint);
		}

		// perfectly smooth BorderSurface between the different parts of the roughened fibre layer
		if(RoughenedEndPoint_cladding3_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding3_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedEndPoint_cladding3_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + " and " + RoughenedEndPoint_cladding3_physical[iter + counter]->GetName() + ")", RoughenedEndPoint_cladding3_physical[iter], RoughenedEndPoint_cladding3_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding3_physical[iter + counter]->GetName() + " and " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding3_physical[iter + counter], RoughenedEndPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}
		if(RoughenedEndPoint_cladding2_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding2_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedEndPoint_cladding2_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[iter + counter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter], RoughenedEndPoint_cladding2_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[iter + counter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter + counter], RoughenedEndPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}
		if(RoughenedEndPoint_cladding1_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding1_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedEndPoint_cladding1_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[iter + counter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter], RoughenedEndPoint_cladding1_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[iter + counter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter + counter], RoughenedEndPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}
		if(RoughenedEndPoint_fibreCore_physical.size())
		{
			for(unsigned int iter = 0; iter < RoughenedEndPoint_fibreCore_physical.size() - 1; iter++)
			{
				for(unsigned int counter = 1; counter < RoughenedEndPoint_fibreCore_physical.size() - iter; counter++)
				{
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + " and " + RoughenedEndPoint_fibreCore_physical[iter + counter]->GetName() + ")", RoughenedEndPoint_fibreCore_physical[iter], RoughenedEndPoint_fibreCore_physical[iter + counter], OptSurf_perfectlySmooth);
					new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_fibreCore_physical[iter + counter]->GetName() + " and " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + ")", RoughenedEndPoint_fibreCore_physical[iter + counter], RoughenedEndPoint_fibreCore_physical[iter], OptSurf_perfectlySmooth);
				}
			}
		}

		// perfectly smooth BorderSurface between the different layers of the roughened fibre part
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding3_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < RoughenedEndPoint_cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding3_physical[iter], RoughenedEndPoint_cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[counter], RoughenedEndPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding2_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < RoughenedEndPoint_cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter], RoughenedEndPoint_cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[counter], RoughenedEndPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding1_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < RoughenedEndPoint_fibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + " and " + RoughenedEndPoint_fibreCore_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter], RoughenedEndPoint_fibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_fibreCore_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + ")", RoughenedEndPoint_fibreCore_physical[counter], RoughenedEndPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
		}

	// perfectly smooth BorderSurface to all other fibre parts
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding3_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding3_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding3_physical[iter], Cladding3_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + ")", Cladding3_physical[counter], RoughenedEndPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding3_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding3_physical[iter]->GetName() + ")", Cladding2_physical[counter], RoughenedEndPoint_cladding3_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding2_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding3_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + " and " + Cladding3_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter], Cladding3_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding3_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + ")", Cladding3_physical[counter], RoughenedEndPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + ")", Cladding2_physical[counter], RoughenedEndPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding2_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding2_physical[iter]->GetName() + ")", Cladding1_physical[counter], RoughenedEndPoint_cladding2_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_cladding1_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding2_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + " and " + Cladding2_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter], Cladding2_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding2_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + ")", Cladding2_physical[counter], RoughenedEndPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + ")", Cladding1_physical[counter], RoughenedEndPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", RoughenedEndPoint_cladding1_physical[iter], FibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + RoughenedEndPoint_cladding1_physical[iter]->GetName() + ")", FibreCore_physical[counter], RoughenedEndPoint_cladding1_physical[iter], OptSurf_perfectlySmooth);
			}
		}
		for(unsigned int iter = 0; iter < RoughenedEndPoint_fibreCore_physical.size(); iter++)
		{
			for(unsigned int counter = 0; counter < Cladding1_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + " and " + Cladding1_physical[counter]->GetName() + ")", RoughenedEndPoint_fibreCore_physical[iter], Cladding1_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + Cladding1_physical[counter]->GetName() + " and " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + ")", Cladding1_physical[counter], RoughenedEndPoint_fibreCore_physical[iter], OptSurf_perfectlySmooth);
			}
			for(unsigned int counter = 0; counter < FibreCore_physical.size(); counter++)
			{
				new G4LogicalBorderSurface("nonroughened Surface (between " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + " and " + FibreCore_physical[counter]->GetName() + ")", RoughenedEndPoint_fibreCore_physical[iter], FibreCore_physical[counter], OptSurf_perfectlySmooth);
				new G4LogicalBorderSurface("nonroughened Surface (between " + FibreCore_physical[counter]->GetName() + " and " + RoughenedEndPoint_fibreCore_physical[iter]->GetName() + ")", FibreCore_physical[counter], RoughenedEndPoint_fibreCore_physical[iter], OptSurf_perfectlySmooth);
			}
		}
	}
}



void G4Fibre::SetDefaults()
{
	CriticalErrorOccured = false;

	MotherVolume_solid = 0;

// Materials & Elements
	Material_OpticalCement = 0;
	Material_FibreCore = 0;
	Material_Cladding1 = 0;
	Material_Cladding2 = 0;
	Material_Cladding3 = 0;
	Material_Coating = 0;

// Fibre
	FibreRound = false;
	IsScinti = false;
	IsWLS = false;
	Cladding1Exists = false;
	Cladding2Exists = false;
	Cladding3Exists = false;
	CoatingExists = false;
	FibreRadius = NAN;
	RelativeFibreRadius_Core = NAN;
	RelativeFibreRadius_Cladding1 = NAN;
	RelativeFibreRadius_Cladding2 = NAN;
	RelativeFibreRadius_Cladding3 = NAN;
	RelativeFibreRadius_CoatingMin = NAN;
	RelativeFibreRadius_CoatingMax = NAN;
	FibreEdgeLength = NAN;
	FibreTransformation_insideMother = G4Transform3D();
	FibreTransformation_outsideMother = G4Transform3D();

// Embedment
	EmbedmentThickness = NAN;


	OptSurf_startPoint = 0;
	OptSurf_endPoint = 0;
	OptSurf_perfectlySmooth = 0;

// Reflective fibre end
	Reflective_thickness = 1 * CLHEP::nm;

	ReflectiveStartPointTransformation_insideMother = G4Transform3D();
	ReflectiveStartPointTransformation_outsideMother = G4Transform3D();
	if(isnan(Reflectivity_startPoint)) CreateReflectiveStartPointVolumes = false;
	else CreateReflectiveStartPointVolumes = true;

	ReflectiveEndPointTransformation_insideMother = G4Transform3D();
	ReflectiveEndPointTransformation_outsideMother = G4Transform3D();
	if(isnan(Reflectivity_endPoint)) CreateReflectiveEndPointVolumes = false;
	else CreateReflectiveEndPointVolumes = true;

// Roughened fibre end
	Rough_thickness = 1 * CLHEP::nm;

	RoughenedStartPointTransformation_insideMother = G4Transform3D();
	RoughenedStartPointTransformation_outsideMother = G4Transform3D();

	RoughenedEndPointTransformation_insideMother = G4Transform3D();
	RoughenedEndPointTransformation_outsideMother = G4Transform3D();
}
