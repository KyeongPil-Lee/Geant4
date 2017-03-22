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
// $Id: WGR16HodoscopeHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16HodoscopeHit.hh
/// \brief Definition of the WGR16HodoscopeHit class

#ifndef WGR16PMTHit_h
#define WGR16PMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Hodoscope hit
///
/// It records:
/// - the strip ID
/// - the particle time
/// - the strip logical volume, its position and rotation

class WGR16PMTHit : public G4VHit
{
public:
	WGR16PMTHit();
	//WGR16PMTHit(G4int i,G4double t);
	WGR16PMTHit(const WGR16PMTHit &right);
	virtual ~WGR16PMTHit();

	const WGR16PMTHit& operator=(const WGR16PMTHit &right);
	int operator==(const WGR16PMTHit &right) const;
	
	inline void *operator new(size_t);
	inline void operator delete(void*aHit);
	
	void Draw();
	//virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
	//virtual std::vector<G4AttValue>* CreateAttValues() const;
	void Print();

	///////////////////////////////////
	// -- methods from base class -- //
	///////////////////////////////////
	// Set methods
	void SetTrackID  (G4int track)      { fTrackID = track; };
	//void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
	void SetEdep     (G4double de)      { fEdep = de; };
	void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

	// Get methods
	G4int GetTrackID() const     { return fTrackID; };
	//G4int GetChamberNb() const   { return fChamberNb; };
	G4double GetEdep() const     { return fEdep; };
	G4ThreeVector GetPos() const { return fPos; };

private:
	G4int         fTrackID;
	//G4int         fChamberNb;
	G4double      fEdep;
	G4ThreeVector fPos;
};

typedef G4THitsCollection<WGR16PMTHit> WGR16PMTHitsCollection;

extern G4ThreadLocal G4Allocator<WGR16PMTHit>* WGR16PMTHitAllocator;

inline void* WGR16PMTHit::operator new(size_t)
{
	if(!WGR16PMTHitAllocator)
		WGR16PMTHitAllocator = new G4Allocator<WGR16PMTHit>;

	return (void*)WGR16PMTHitAllocator->MallocSingle();
}

inline void WGR16PMTHit::operator delete(void*aHit)
{
	WGR16PMTHitAllocator->FreeSingle((WGR16PMTHit*) aHit);
}

#endif