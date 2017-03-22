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
// $Id: WGR16PrimaryGeneratorAction.hh 76474 2013-11-11 10:36:34Z gcosmo $
//
/// \file WGR16PrimaryGeneratorAction.hh
/// \brief Definition of the WGR16PrimaryGeneratorAction class

#ifndef WGR16PrimaryGeneratorAction_h
#define WGR16PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

class WGR16PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	WGR16PrimaryGeneratorAction();
	virtual ~WGR16PrimaryGeneratorAction();
	
	virtual void GeneratePrimaries(G4Event*);

	void SetEta(G4double eta) { feta = eta;}
	G4double GetEta() const { return feta;}
	void SetPhi(G4double phi) { fphi = phi;}
	G4double GetPhi() const { return fphi;}

private:
	void DefineCommands();

	G4ParticleGun* fParticleGun;
	G4GenericMessenger* fMessenger;
	G4ParticleDefinition* fElectron;
	G4ParticleDefinition* fPositron;
	G4ParticleDefinition* fMuon;
	G4ParticleDefinition* fPion;
	G4ParticleDefinition* fKaon;
	G4ParticleDefinition* fProton;
	G4ParticleDefinition* fOptGamma;

	G4double feta;
	G4double fphi;
};

#endif