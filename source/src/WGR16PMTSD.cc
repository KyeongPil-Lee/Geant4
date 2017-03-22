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
// $Id: WGR16HodoscopeSD.cc 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16HodoscopeSD.cc
/// \brief Implementation of the WGR16HodoscopeSD class

#include "WGR16PMTSD.hh"
#include "WGR16PMTHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

WGR16PMTSD::WGR16PMTSD(const G4String& name, const G4String& hitsCollectionName):
G4VSensitiveDetector(name), 
fHitsCollection(0), fHCID(-1)
{
	collectionName.insert(hitsCollectionName);
}

WGR16PMTSD::~WGR16PMTSD()
{

}

void WGR16PMTSD::Initialize(G4HCofThisEvent* HCOfThisEvent)
{
	fHitsCollection = new WGR16PMTHitsCollection(SensitiveDetectorName,collectionName[0]);
	if( fHCID < 0 )
		fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
	HCOfThisEvent->AddHitsCollection(fHCID, fHitsCollection);
}

G4bool WGR16PMTSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
	WGR16PMTHit* newHit = new WGR16PMTHit();
	fHitsCollection->insert( newHit );
	return true;
}

void WGR16PMTSD::EndOfEvent(G4HCofThisEvent*)
{
	if( verboseLevel > 1 )
	{
		G4int nofHits = fHitsCollection->entries();
		G4cout << G4endl
		<< "-------->Hits Collection: in this event they are " << nofHits
		<< " hits in the tracker chambers: " << G4endl;

		for( G4int i=0; i<nofHits; i++ ) 
			(*fHitsCollection)[i]->Print();
	}
}


