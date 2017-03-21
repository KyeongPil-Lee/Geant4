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
// $Id: WGR16EventAction.cc 94486 2015-11-19 08:33:37Z gcosmo $
//
/// \file WGR16EventAction.cc
/// \brief Implementation of the WGR16EventAction class

#include "WGR16DetectorConstruction.hh"
#include "WGR16EventAction.hh"
#include "WGR16Analysis.hh"
#include "WGR16PMTHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include <string.h>
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16EventAction::WGR16EventAction()
: G4UserEventAction()
  ,fHitcnt(0),fPMTNum(0),
	fEdep(0),PMTHCID(0)
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16EventAction::~WGR16EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16EventAction::BeginOfEventAction(const G4Event*)
{
	edep = 0;  
   ClearVectors();
	G4SDManager* sdManager = G4SDManager::GetSDMpointer();

	PMTHCID = sdManager->GetCollectionID("PMTColl");// get the ID of HitCollection
}     

void WGR16EventAction::ClearVectors()
{
  fHitcnt.clear();
  fPMTNum.clear();
  fEdep.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16EventAction::EndOfEventAction(const G4Event* event)
{
	G4HCofThisEvent* hce = event->GetHCofThisEvent();
	if (!hce) 
	{
		G4ExceptionDescription msg;
		msg << "No hits collection of this event found." << G4endl; 
		G4Exception("WGR16EventAction::EndOfEventAction()",
				"WGR16Code001", JustWarning, msg);
		return;
	}   
	G4AnalysisManager* aM = G4AnalysisManager::Instance();
	
	WGR16PMTHitsCollection* PMTHC
		//= static_cast<WGR16PMTHitsCollection*>(hce->GetHC(0));
		= static_cast<WGR16PMTHitsCollection*>(hce->GetHC(PMTHCID));
	fHitcnt.push_back(PMTHC->entries());
   
	for(G4int i=0;i<PMTHC->entries();i++){
		WGR16PMTHit* hit = (*PMTHC)[i];
		fEdep.push_back(hit->GetEdep());
	}  // hit information
	
	G4double a=0;
	for(int i =0;i<fEdep.size();i++){
		a=a+ fEdep.at(i);
	}
	std::cout<<a<<" " <<edep<<std::endl;
	
	aM->FillNtupleDColumn(1,edep);// from Stepping Action
	
	aM->AddNtupleRow(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
