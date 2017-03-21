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
// $Id$
//
/// \file B2bChamberParameterisation.cc
/// \brief Implementation of the B2bChamberParameterisation class

#include "WGR16TowerParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <stdio.h>
#include <string.h>

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16TowerParameterisation::WGR16TowerParameterisation(  
        G4int    noZRotation, 
		  G4int	  noTowers,
        G4double translationlength,
		  string type,
		  G4int no
         )
 : G4VPVParameterisation()
{
   fNoZRotation = noZRotation;
	fNoTowers = noTowers;
   fTranslationlength = translationlength;
	fType	= type;
	fNoZ = no;// this is the variable for barrel i.e tower number in x,y plane
   G4double pi = 3.141592653589793238462643383279;	

	theta_unit = 2*pi/(G4double)fNoZRotation;
/*
	//calculate series of ThreeVectors	
	if(fType == "Barrel"){
		origin_barrel = new G4ThreeVector(fTranslationlength*cos(fNoZ*theta_unit),fTranslationlength*sin(copyNo*theta_unit),0);

		for(int i=0;i<fNoTowers;i++){
		x0.push_back(new G4ThreeVector(,,));

		}
	}
*/


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16TowerParameterisation::~WGR16TowerParameterisation()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16TowerParameterisation::ComputeTransformation
(const G4int copyNo, G4VPhysicalVolume* physVol) const
{

		G4double pi =3.141592653589793238462643383279;

	if(fType == "center"){
		
		G4ThreeVector origin(fTranslationlength*cos(copyNo*theta_unit),fTranslationlength*sin(copyNo*theta_unit),0);
		// Note: copyNo will start with zero!

		//G4double Zposition = fStartZ + copyNo * fSpacing;
		//G4ThreeVector origin(0,0,0);
		physVol->SetTranslation(origin);
		G4RotationMatrix * RotMatrix = new G4RotationMatrix();
		RotMatrix->rotateY(-90*deg);
		//RotMatrix->rotateX((-10*copyNo)*deg);
		RotMatrix->rotateX(theta_unit*copyNo);
		physVol->SetRotation(RotMatrix);
	}
	
	if(fType == "barrelL"){

	}

	if(fType == "barrelR"){

	}

	if(fType == "endcapL"){

	}

	if(fType == "endcapR"){

	}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16TowerParameterisation::ComputeDimensions
(G4Trd& unitTower, const G4int copyNo, const G4VPhysicalVolume*) const
{
  // Note: copyNo will start with zero!
 /* G4double rmax = fRmaxFirst + copyNo * fRmaxIncr;
  unitTower.SetInnerRadius(0);
  unitTower.SetOuterRadius(rmax);
  unitTower.SetZHalfLength(fHalfWidth);
  unitTower.SetStartPhiAngle(0.*deg);
  unitTower.SetDeltaPhiAngle(360.*deg);*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
