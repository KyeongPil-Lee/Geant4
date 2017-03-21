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
// $Id: WGR16HodoscopeHit.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16HodoscopeHit.cc
/// \brief Implementation of the WGR16HodoscopeHit class

#include "WGR16PMTHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <iomanip>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4Allocator<WGR16PMTHit>* WGR16PMTHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//WGR16PMTHit::WGR16PMTHit(G4int i,G4double t)
WGR16PMTHit::WGR16PMTHit()
: G4VHit(), fTrackID(-1),fEdep(0.),fPos(G4ThreeVector())
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16PMTHit::~WGR16PMTHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16PMTHit::WGR16PMTHit(const WGR16PMTHit &right)
: G4VHit() {
	fTrackID   = right.fTrackID;
	//fChamberNb = right.fChamberNb;
	fEdep      = right.fEdep;
	fPos       = right.fPos;
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const WGR16PMTHit& WGR16PMTHit::operator=(const WGR16PMTHit &right)
{
	fTrackID   = right.fTrackID;
	//fChamberNb = right.fChamberNb;
	fEdep      = right.fEdep;
	fPos       = right.fPos;
	return *this;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int WGR16PMTHit::operator==(const WGR16PMTHit &/*right*/) const
{
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16PMTHit::Draw()
{
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	 if(pVVisManager)
	 {
		 G4Circle circle(fPos);
		 circle.SetScreenSize(4.);
		 circle.SetFillStyle(G4Circle::filled);
		 G4Colour colour(0.,1.,1.);
		 G4VisAttributes attribs(colour);
		 circle.SetVisAttributes(attribs);
		 pVVisManager->Draw(circle);
	 }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*const std::map<G4String,G4AttDef>* WGR16HodoscopeHit::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("WGR16HodoscopeHit",isNew);

    if (isNew) {
        (*store)["HitType"] 
          = G4AttDef("HitType","Hit Type","Physics","","G4String");
        
        (*store)["ID"] 
          = G4AttDef("ID","ID","Physics","","G4int");
        
        (*store)["Time"] 
          = G4AttDef("Time","Time","Physics","G4BestUnit","G4double");
        
        (*store)["Pos"] 
          = G4AttDef("Pos","Position","Physics","G4BestUnit","G4ThreeVector");
        
        (*store)["LVol"] 
          = G4AttDef("LVol","Logical Volume","Physics","","G4String");
    }
    return store;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*std::vector<G4AttValue>* WGR16HodoscopeHit::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;
    
    values
      ->push_back(G4AttValue("HitType","HodoscopeHit",""));
    values
      ->push_back(G4AttValue("ID",G4UIcommand::ConvertToString(fId),""));
    values
      ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
    values
      ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));
    
    if (fPLogV)
        values->push_back(G4AttValue("LVol",fPLogV->GetName(),""));
    else
        values->push_back(G4AttValue("LVol"," ",""));
    
    return values;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16PMTHit::Print()
{
	G4cout
		<< "  trackID: " << fTrackID //<< " chamberNb: " << fChamberNb
		<< "Edep: "
		<< std::setw(7) << G4BestUnit(fEdep,"Energy")
		<< " Position: "
		<< std::setw(7) << G4BestUnit( fPos,"Length")
		<< G4endl;


	 
	 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
