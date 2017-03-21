//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: WGR16DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file WGR16DetectorConstruction.cc
/// \brief Implementation of the WGR16DetectorConstruction class
#include "WGR16PMTSD.hh"
#include "WGR16DetectorConstruction.hh"
#include "WGR16MagneticField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"
#include "G4ThreeVector.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4VisExtent.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"

#include <cmath>
#include <stdio.h>
#include <float.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal WGR16MagneticField* WGR16DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WGR16DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16DetectorConstruction::WGR16DetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fVisAttributes(),fScoringVolume(0)
{
	DefineCommands();
}


WGR16DetectorConstruction::~WGR16DetectorConstruction()
{
	delete fMessenger;
	
	for (G4int i=0; i<G4int(fVisAttributes.size()); ++i)
	{
		delete fVisAttributes[i];
	}
}

G4VPhysicalVolume* WGR16DetectorConstruction::Construct()
{
	G4bool checkOverlaps = false;

	// -- meterials -- //
	ConstructMaterials();
	G4Material* vac = G4Material::GetMaterial("G4_Galactic");
	G4Material* cu  = G4Material::GetMaterial("G4_Cu");

	// const G4double cu_radlen = cu->GetRadlen(); // -- radiation length -- //

	// -- world -- //
	G4VSolid* worldSolid 
	= new G4Box("worldBox",10.*m,10.*m,10.*m);
	G4LogicalVolume* worldLogical
	= new G4LogicalVolume(worldSolid,vac,"worldLogical");
	G4VPhysicalVolume* worldPhysical
	= new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0, false,0,checkOverlaps);

	//////////////////
	// -- Cu box -- //
	//////////////////
	G4double pi = 3.14159265359;
	G4double radius = 1.8; // -- unit: m -- //
	G4double Circumference = 2 * pi * radius;
	// G4double nTower_PhiDir = 283;
	G4double nTower_PhiDir = 4;
	G4double CuLen_PhiDir = (Circumference / nTower_PhiDir)*m;

	G4double CuLen_EtaDir = CuLen_PhiDir*2.0;
	G4double CuLen_H = 2.5*m;

	G4cout << "[Cu] (PhiDir, EtaDir, Height (unit:m)) = (" << CuLen_PhiDir << ", " << CuLen_EtaDir << ", " << CuLen_H << ", )" << G4endl;

	G4Box* CuBox 
	= new G4Box("CuBox", CuLen_H/2.0, CuLen_PhiDir/2., CuLen_EtaDir/2.0);
	G4LogicalVolume *CuLogical
	= new G4LogicalVolume(CuBox, cu, "CuLogical");

	// -- default unit for the rotation -- //
	G4double dPhi = (2*pi) / nTower_PhiDir;
	G4double half_dPhi = 0.5*dPhi;
	for(G4int i_cu=0; i_cu<nTower_PhiDir; i_cu++)
	{
		G4double phi = i_cu*dPhi;
		G4RotationMatrix rotM  = G4RotationMatrix();
		rotM.rotateY(90*deg);
		rotM.rotateZ(phi);

		G4ThreeVector Unit_Z = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
		G4ThreeVector position = (radius + 0.5*CuLen_H)*Unit_Z; // -- radius = size of the vector -- //
		G4Transform3D transform = G4Transform3D(rotM,position);

		new G4PVPlacement(transform, CuLogical, "CuPhysical", worldLogical, false, i_cu, checkOverlaps );

		// G4ThreeVector origin(x,y,z);
		// G4RotationMatrix* RotMatrix = new G4RotationMatrix();

		// RotMatrix->rotateZ(90*deg);
		// RotMatrix->rotateZ(-i*theta_unit_ZRot);
		// RotMatrix->rotateX(90*deg);
		// RotMatrix->rotateX(-theta_unit*(copyNo+0.5));

		// -- place it -- //
		// new G4PVPlacement( RotMatrix, G4ThreeVector(), CuLogical, "CuPhysical", worldLogical, false, 0, checkOverlaps );
	}

	// -- visualization -- //
	G4VisAttributes* visAttr;
	visAttr = new G4VisAttributes( G4Colour(0.0,1.0,1.0) );
	visAttr->SetVisibility(true);
	CuLogical->SetVisAttributes(visAttr);
	fVisAttributes.push_back(visAttr);

	return worldPhysical;
}

void WGR16DetectorConstruction::ConstructSDandField()
{

}

void WGR16DetectorConstruction::ConstructMaterials()
{
	G4NistManager* nistManager = G4NistManager::Instance();
	
	// -- Vacuum: "Galactic" -- //
	nistManager->FindOrBuildMaterial("G4_Galactic");
	nistManager->FindOrBuildMaterial("G4_Cu");

	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void WGR16DetectorConstruction::DefineCommands()
{
	// Define /WGR16/detector command directory using generic messenger class
	// fMessenger = new G4GenericMessenger(this, 
	//                                     "/WGR16/detector/", 
	//                                     "Detector control");

	// G4GenericMessenger::Command& B2ECollAngleCmd
	//   = fMessenger->DeclareMethodWithUnit("Block2Angle","deg",
	//       &WGR16DetectorConstruction::SetBlock2Angle,
	//       "Set rotation angle of exit collimator and corresponding detector");
	// B2ECollAngleCmd.SetParameterName("angle",true);
	// B2ECollAngleCmd.SetRange("angle>=-90. && angle<=90");
	// B2ECollAngleCmd.SetDefaultValue("60.");	
}





















