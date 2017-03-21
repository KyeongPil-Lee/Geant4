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
// $Id: WGR16ActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file WGR16ActionInitialization.cc
/// \brief Implementation of the WGR16ActionInitialization class

#include "WGR16ActionInitialization.hh"
#include "WGR16PrimaryGeneratorAction.hh"
#include "WGR16RunAction.hh"
#include "WGR16EventAction.hh"
#include "WGR16SteppingAction.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16ActionInitialization::WGR16ActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16ActionInitialization::~WGR16ActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16ActionInitialization::BuildForMaster() const
{
  WGR16EventAction* eventAction = 0;
  SetUserAction(new WGR16RunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16ActionInitialization::Build() const
{
  SetUserAction(new WGR16PrimaryGeneratorAction);

  WGR16EventAction* eventAction = new WGR16EventAction;
  SetUserAction(eventAction);
  
  SetUserAction(new WGR16RunAction(eventAction));
 
  SetUserAction(new WGR16SteppingAction(eventAction));//
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
