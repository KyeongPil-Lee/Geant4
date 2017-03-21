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
// $Id: WGR16HadCalorimeterHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16HadCalorimeterHit.hh
/// \brief Definition of the WGR16HadCalorimeterHit class

#ifndef WGR16HadCalorimeterHit_h
#define WGR16HadCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Hadron Calorimeter hit
///
/// It records:
/// - the cell column ID and row ID
/// - the energy deposit 
/// - the cell position and rotation

class WGR16HadCalorimeterHit : public G4VHit
{
public:
    WGR16HadCalorimeterHit();
    WGR16HadCalorimeterHit(G4int iCol,G4int iRow);
    WGR16HadCalorimeterHit(const WGR16HadCalorimeterHit &right);
    virtual ~WGR16HadCalorimeterHit();

    const WGR16HadCalorimeterHit& operator=(const WGR16HadCalorimeterHit &right);
    int operator==(const WGR16HadCalorimeterHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();
    
    void SetColumnID(G4int z) { fColumnID = z; }
    G4int GetColumnID() const { return fColumnID; }

    void SetRowID(G4int z) { fRowID = z; }
    G4int GetRowID() const { return fRowID; }

    void SetEdep(G4double de) { fEdep = de; }
    void AddEdep(G4double de) { fEdep += de; }
    G4double GetEdep() const { return fEdep; }

    void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    G4ThreeVector GetPos() const { return fPos; }

    void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
    G4RotationMatrix GetRot() const { return fRot; }
    
private:
    G4int fColumnID;
    G4int fRowID;
    G4double fEdep;
    G4ThreeVector fPos;
    G4RotationMatrix fRot;
};

typedef G4THitsCollection<WGR16HadCalorimeterHit> WGR16HadCalorimeterHitsCollection;

extern G4ThreadLocal G4Allocator<WGR16HadCalorimeterHit>* WGR16HadCalorimeterHitAllocator;

inline void* WGR16HadCalorimeterHit::operator new(size_t)
{
    if (!WGR16HadCalorimeterHitAllocator)
        WGR16HadCalorimeterHitAllocator = new G4Allocator<WGR16HadCalorimeterHit>;
    return (void*)WGR16HadCalorimeterHitAllocator->MallocSingle();
}

inline void WGR16HadCalorimeterHit::operator delete(void* aHit)
{
    WGR16HadCalorimeterHitAllocator->FreeSingle((WGR16HadCalorimeterHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
