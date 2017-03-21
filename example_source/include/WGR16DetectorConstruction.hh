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

private:
    void DefineCommands();

    G4GenericMessenger* fMessenger;
    
    static G4ThreadLocal WGR16MagneticField* fMagneticField;
	 static G4ThreadLocal G4FieldManager* fFieldMgr;
	
	 G4LogicalVolume* cuPlateLogical;
	 G4LogicalVolume* fMagneticfieldLogical;
	
	 std::vector<G4VisAttributes*> fVisAttributes;
    //G4VPhysicalVolume* fSecondArmPhys;

protected:
	 	G4LogicalVolume* fScoringVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
