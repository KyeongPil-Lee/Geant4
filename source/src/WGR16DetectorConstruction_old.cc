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
: G4VUserDetectorConstruction(), 
  fMessenger(0),
  fVisAttributes(),fScoringVolume(0)
{
    // define commands for this class
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16DetectorConstruction::~WGR16DetectorConstruction()
{
    delete fMessenger;
    
    for (G4int i=0; i<G4int(fVisAttributes.size()); ++i) 
    {
      delete fVisAttributes[i];
    }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WGR16DetectorConstruction::Construct()
{

    G4bool checkOverlaps = false;

    ConstructMaterials();
    G4Material* vac = G4Material::GetMaterial("G4_Galactic");

    G4Material* vac = G4Material::GetMaterial("G4_Galactic");
    G4Material* cu  = G4Material::GetMaterial("G4_Cu");

    const G4double cu_radlen = cu->GetRadlen();

    G4VSolid* worldSolid 
      = new G4Box("worldBox",10.*m,10.*m,10.*m);
    G4LogicalVolume* worldLogical
      = new G4LogicalVolume(worldSolid,vac,"worldLogical");
    G4VPhysicalVolume* worldPhysical
      = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
          false,0,checkOverlaps);

    return worldPhysical;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16DetectorConstruction::ConstructSDandField()
{

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16DetectorConstruction::ConstructMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();
  
  // Vacuum "Galactic"
  nistManager->FindOrBuildMaterial("G4_Galactic");
  nistManager->FindOrBuildMaterial("G4_Cu");

  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16DetectorConstruction::DefineCommands()
{
    // Define /WGR16/detector command directory using generic messenger class
/*    fMessenger = new G4GenericMessenger(this, 
                                        "/WGR16/detector/", 
                                        "Detector control");

    G4GenericMessenger::Command& B2ECollAngleCmd
      = fMessenger->DeclareMethodWithUnit("Block2Angle","deg",
          &WGR16DetectorConstruction::SetBlock2Angle,
          "Set rotation angle of exit collimator and corresponding detector");
    B2ECollAngleCmd.SetParameterName("angle",true);
    B2ECollAngleCmd.SetRange("angle>=-90. && angle<=90");
    B2ECollAngleCmd.SetDefaultValue("60.");
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......