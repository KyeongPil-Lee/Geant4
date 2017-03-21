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
// $Id: WGR163DHit.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR163DHit.hh
/// \brief Definition of the WGR163DHit class

#ifndef WGR163DHit_h
#define WGR163DHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class WGR163DHit : public G4VHit
{
public:
    enum kParticleKind{kGamma,kPositron,kElectron};

    WGR163DHit();
    WGR163DHit(G4int z);
    WGR163DHit(const WGR163DHit &right);
    virtual ~WGR163DHit();

    const WGR163DHit& operator=(const WGR163DHit &right);
    int operator==(const WGR163DHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    inline void SetWorldPos(G4ThreeVector xyz) { fWorldPos = xyz; }
    inline G4ThreeVector GetWorldPos() const { return fWorldPos; }

    inline void SetLocalPos(G4ThreeVector xyz) { fLocalPos = xyz; }
    inline G4ThreeVector GetLocalPos() const { return fLocalPos; }

    inline void SetKinEnergy(G4double in) { kinEnergy = in; }
    inline G4double GetKinEnergy() const { return kinEnergy; }

    inline void SetParticleKind(kParticleKind a){fPK = a;}
    inline kParticleKind GetParticleKind() const {return fPK;}

private:
    G4ThreeVector fLocalPos;
    G4ThreeVector fWorldPos;
    G4double kinEnergy;
    kParticleKind fPK;
};

typedef G4THitsCollection<WGR163DHit> WGR163DHitsCollection;

extern G4ThreadLocal G4Allocator<WGR163DHit>* WGR163DHitAllocator;

inline void* WGR163DHit::operator new(size_t)
{
    if (!WGR163DHitAllocator)
        WGR163DHitAllocator = new G4Allocator<WGR163DHit>;
    return (void*)WGR163DHitAllocator->MallocSingle();
}

inline void WGR163DHit::operator delete(void* aHit)
{
    WGR163DHitAllocator->FreeSingle((WGR163DHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
