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
// $Id: WGR16EventAction.hh 94486 2015-11-19 08:33:37Z gcosmo $
//
/// \file WGR16EventAction.hh
/// \brief Definition of the WGR16EventAction class

#ifndef WGR16EventAction_h
#define WGR16EventAction_h 1


#include "G4UserEventAction.hh"
#include "globals.hh"

#include <vector>

class WGR16EventAction : public G4UserEventAction
{
public:
	WGR16EventAction();
	virtual ~WGR16EventAction();

	virtual void BeginOfEventAction(const G4Event*);
	virtual void EndOfEventAction(const G4Event*);

	std::vector<G4int>& GetRofHitcntVec() { return fHitcnt; }
	std::vector<G4int>& GetRofPMTNumVec() { return fPMTNum; }
	std::vector<G4double>& GetRofEdepVec() { return fEdep; }

private:
	void ClearVectors();

	std::vector<G4int> fPMTNum;
	std::vector<G4int> fHitcnt;
	std::vector<G4double> fEdep;

	G4double edep;
	G4int PMTHCID;
};

#endif