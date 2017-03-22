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
// $Id: WGR16HodoscopeSD.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16HodoscopeSD.hh
/// \brief Definition of the WGR16HodoscopeSD class

#ifndef WGR16PMTSD_h
#define WGR16PMTSD_h 1

#include "G4VSensitiveDetector.hh"
#include "WGR16PMTHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

/// Hodoscope sensitive detector

class WGR16PMTSD : public G4VSensitiveDetector
{
public:
	//WGR16PMTSD(G4String name);
	WGR16PMTSD(const G4String& name, const G4String& hitsCollectionName);
	virtual ~WGR16PMTSD();

	virtual void Initialize(G4HCofThisEvent* HCOfThisEvent);
	virtual G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
	virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

private:
	WGR16PMTHitsCollection* fHitsCollection;
	G4int fHCID;
};

#endif