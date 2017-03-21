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
// $Id: WGR16DetectorConstruction.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16DetectorConstruction.hh
/// \brief Definition of the WGR16DetectorConstruction class

#ifndef WGR16DetectorConstruction_h
#define WGR16DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RotationMatrix.hh"
#include "G4FieldManager.hh"
#include <string.h>
#include <vector>

class WGR16MagneticField;

class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;
class G4VPhysicalVolume;

/// Detector construction

class WGR16DetectorConstruction : public G4VUserDetectorConstruction
{
public:
	WGR16DetectorConstruction();
	virtual ~WGR16DetectorConstruction();
	
	virtual G4VPhysicalVolume* Construct();
	virtual void ConstructSDandField();

	void ConstructMaterials();
	G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

	G4MaterialPropertiesTable* MaterialPropertyTable_PMMA();
	G4MaterialPropertiesTable* MaterialPropertyTable_FS();
	G4MaterialPropertiesTable* MaterialPropertyTable_Glass();
	G4MaterialPropertiesTable* MaterialPropertyTable_Air();

	bool IsFiberC(G4int i_EtaDir, G4int i_PhiDir);

private:
	void DefineCommands();

	G4GenericMessenger* fMessenger;
	
	static G4ThreadLocal WGR16MagneticField* fMagneticField;
	static G4ThreadLocal G4FieldManager* fFieldMgr;
	std::vector<G4VisAttributes*> fVisAttributes;

	// const G4int nEntries = 50;
	// G4double PhotonEnergy[nEntries] = {
	// 2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	// 2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	// 2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	// 2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	// 2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	// 2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	// 2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	// 3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	// 3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	// 3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

protected:
	G4LogicalVolume* fScoringVolume;
};

#endif