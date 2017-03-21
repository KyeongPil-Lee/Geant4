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
// $Id: WGR16RunAction.cc 74204 2013-10-01 07:04:43Z ihrivnac $
//
/// \file WGR16RunAction.cc
/// \brief Implementation of the WGR16RunAction class

#include "WGR16DetectorConstruction.hh"
#include "WGR16RunAction.hh"
#include "WGR16EventAction.hh"
#include "WGR16Analysis.hh"
#include "WGR16PrimaryGeneratorAction.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16RunAction::WGR16RunAction(WGR16EventAction* eventAction)
 : G4UserRunAction(),
   fEventAction(eventAction)
{ 
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in WGR16Analysis.hh
    
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetVerboseLevel(1);

  G4RunManager::GetRunManager()->SetPrintProgress(10000);
  // Book histograms, ntuple
  //
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16RunAction::~WGR16RunAction()
{
  //delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  // Creating ntuple
  //
  if ( fEventAction ) {
    analysisManager->CreateNtuple("WGR16", "Hit Data");  // Branch Id = 0
    analysisManager->CreateNtupleIColumn("Hitcnt",fEventAction->GetRofHitcntVec());
	 analysisManager->CreateNtupleDColumn("EdepfromSteppingAction");
	 analysisManager->CreateNtupleDColumn("Edep",fEventAction->GetRofEdepVec());
	 analysisManager->FinishNtuple();
  }
  // Open an output file 
  // The default file name is set in WGR16RunAction::WGR16RunAction(),
  // it can be overwritten in a macro
//  for(unsigned int i = 2; i < std::strlen(tempS1);i++) if(tempS1[i]=='.') tempS1[i]='-';
  //analysisManager->OpenFile(outputname);
  analysisManager->OpenFile("example");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WGR16RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
  delete analysisManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
