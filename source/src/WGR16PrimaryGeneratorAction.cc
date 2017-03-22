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
// $Id: WGR16PrimaryGeneratorAction.cc 77781 2013-11-28 07:54:07Z gcosmo $
//
/// \file WGR16PrimaryGeneratorAction.cc
/// \brief Implementation of the WGR16PrimaryGeneratorAction class

#include "WGR16PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

WGR16PrimaryGeneratorAction::WGR16PrimaryGeneratorAction():
G4VUserPrimaryGeneratorAction(),
fParticleGun(0), fMessenger(0), 
fElectron(0), fPositron(0), fMuon(0), fPion(0), fKaon(0), fProton(0), fOptGamma(0),
feta(0),fphi(0)
{
	G4int n_particle = 1;
	fParticleGun  = new G4ParticleGun(n_particle);
	
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName;
	fElectron = particleTable->FindParticle(particleName="e-");
	fPositron = particleTable->FindParticle(particleName="e+");
	fMuon = particleTable->FindParticle(particleName="mu+");
	fPion = particleTable->FindParticle(particleName="pi+");
	fKaon = particleTable->FindParticle(particleName="kaon+");
	fProton = particleTable->FindParticle(particleName="proton");
	fOptGamma = particleTable->FindParticle(particleName="opticalphoton");
	
	// default particle kinematics
	G4double height = (1.8 + 2.5)*m + 1*mm;
	fParticleGun->SetParticlePosition(G4ThreeVector(height, 0.*m, 0.*m));
	
	// define commands for this class
	this->DefineCommands();
}

WGR16PrimaryGeneratorAction::~WGR16PrimaryGeneratorAction()
{
	delete fParticleGun;
	delete fMessenger;
}

void WGR16PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
	fParticleGun->SetParticleDefinition(fOptGamma);
	fParticleGun->SetParticleEnergy(2*eV);
	G4ThreeVector Direction;
	Direction.setREtaPhi(1,feta,fphi);
	fParticleGun->SetParticleMomentumDirection(Direction);
	fParticleGun->GeneratePrimaryVertex(event);
}

void WGR16PrimaryGeneratorAction::DefineCommands()
{
	fMessenger = new G4GenericMessenger(this, "/WGR16/generator/", "Primary generator control");

	 G4GenericMessenger::Command& etaCmd =
	 fMessenger->DeclareMethodWithUnit("eta","rad",&WGR16PrimaryGeneratorAction::SetEta,"eta of beam");
	 etaCmd.SetParameterName("eta",true);
	 etaCmd.SetDefaultValue("0.");

	 G4GenericMessenger::Command& phiCmd =
	 fMessenger->DeclareMethodWithUnit("phi","rad",&WGR16PrimaryGeneratorAction::SetPhi,"phi of beam");
	 phiCmd.SetParameterName("phi",true);
	 phiCmd.SetDefaultValue("0.");
}