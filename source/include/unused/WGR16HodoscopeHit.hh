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

#ifndef WGR16HodoscopeHit_h
#define WGR16HodoscopeHit_h 1

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

class WGR16HodoscopeHit : public G4VHit
{
public:
    WGR16HodoscopeHit(G4int i,G4double t);
    WGR16HodoscopeHit(const WGR16HodoscopeHit &right);
    virtual ~WGR16HodoscopeHit();

    const WGR16HodoscopeHit& operator=(const WGR16HodoscopeHit &right);
    int operator==(const WGR16HodoscopeHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void*aHit);
    
    void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    void Print();
    
    G4int GetID() const { return fId; }

    void SetTime(G4double val) { fTime = val; }
    G4double GetTime() const { return fTime; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }

    void SetLogV(G4LogicalVolume* val) { fPLogV = val; }
    const G4LogicalVolume* GetLogV() const { return fPLogV; }
    
private:
    G4int fId;
    G4double fTime;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
    const G4LogicalVolume* fPLogV;
};

typedef G4THitsCollection<WGR16HodoscopeHit> WGR16HodoscopeHitsCollection;

extern G4ThreadLocal G4Allocator<WGR16HodoscopeHit>* WGR16HodoscopeHitAllocator;

inline void* WGR16HodoscopeHit::operator new(size_t)
{
    if (!WGR16HodoscopeHitAllocator)
        WGR16HodoscopeHitAllocator = new G4Allocator<WGR16HodoscopeHit>;
    return (void*)WGR16HodoscopeHitAllocator->MallocSingle();
}

inline void WGR16HodoscopeHit::operator delete(void*aHit)
{
    WGR16HodoscopeHitAllocator->FreeSingle((WGR16HodoscopeHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
