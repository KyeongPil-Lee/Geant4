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
// $Id: WGR16DetectorConstruction.cc 77656 2013-11-27 08:52:57Z gcosmo $
//
/// \file WGR16DetectorConstruction.cc
/// \brief Implementation of the WGR16DetectorConstruction class
#include "WGR16PMTSD.hh"
#include "WGR16DetectorConstruction.hh"
#include "WGR16MagneticField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
// #include "G4MaterialPropertiesTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4PVParameterised.hh"
#include "G4ThreeVector.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"
#include "G4VisExtent.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"

#include <cmath>
#include <stdio.h>
#include <float.h>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal WGR16MagneticField* WGR16DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WGR16DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WGR16DetectorConstruction::WGR16DetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fVisAttributes(),fScoringVolume(0)
{
	DefineCommands();
}


WGR16DetectorConstruction::~WGR16DetectorConstruction()
{
	delete fMessenger;
	
	for (G4int i=0; i<G4int(fVisAttributes.size()); ++i)
	{
		delete fVisAttributes[i];
	}
}

// // -- calculate maximum theta value in barrel -- //
// // -- solve: radius*tan(theta) + CuLen_EtaDir / cos(theta) = BarrelLen_Half -- //
// G4double Calc_ThetaMax_Barrel( G4double	BarrelLen_Half, G4double radius, G4double CuLen_EtaDir )
// {
// 	G4double FirstTerm = std::asin( CuLen_EtaDir / sqrt(BarrelLen_Half*BarrelLen_Half + radius*radius) )
// 	G4double SecondTerm = std::atan( -BarrelLen_Half / radius );

// 	return FirstTerm - SecondTerm;
// }

// Theta equation using proper geometry
G4double Theta_Eq(G4double Pre_Theta, G4double Curr_Theta, G4double radius, G4double CuLen_EtaDir){
	return radius*(std::tan(Curr_Theta)-std::tan(Pre_Theta))+(CuLen_EtaDir/2)*(1/std::cos(Curr_Theta)-1/std::cos(Pre_Theta)-2*std::cos(Curr_Theta)-2*std::sin(Curr_Theta)*std::tan(Pre_Theta));
}

// solving increasing equation by binary search
G4double Solve_Eq(G4double Pre_Theta,G4double Low, G4double Max, G4double radius, G4double CuLen_EtaDir){
	G4double theta=(Low+Max)/2;
	G4double temp=Theta_Eq(Pre_Theta,theta,radius,CuLen_EtaDir);
	if (abs(temp)<0.000001) return theta;
	else if (temp>0.000001) return Solve_Eq(Pre_Theta,Low,theta,radius,CuLen_EtaDir);
	else return Solve_Eq(Pre_Theta,theta,Max,radius,CuLen_EtaDir);
}

G4VPhysicalVolume* WGR16DetectorConstruction::Construct()
{
	G4bool checkOverlaps = false;

	//////////////////////
	// -- meterials -- //
	/////////////////////
	ConstructMaterials();
	G4Material* vac = G4Material::GetMaterial("G4_Galactic");
	G4Material* vac_PMTHouse = G4Material::GetMaterial("G4_Galactic");
	G4Material* cu  = G4Material::GetMaterial("G4_Cu");
	// const G4double cu_radlen = cu->GetRadlen(); // -- radiation length -- //

	G4double z, a, density, fractionmass;
	G4int ncomponents, natoms;
	G4String symbol;
	G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
	G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
	G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
	G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z= 8., a= 16.00*g/mole);
	G4Element* F  = new G4Element("Fluorine",symbol="F" , z= 9., a= 18.9984*g/mole);
	// G4Element* Si = new G4Element("Silicon" ,symbol="Si", z= 14., a= 28.09*g/mole);/

	// -- for PMT Cathod -- //	
	G4Material* Al 
	= new G4Material("Aluminium", z=13., a=26.98*g/mole, density=2.700*g/cm3);

	// G4double PMTT = 1*mm;

	// -- Photocathod property -- //
	G4MaterialPropertiesTable* mpPMTPC = this->MaterialPropertyTable_PMTPC();
	Al->SetMaterialPropertiesTable( mpPMTPC );
	G4OpticalSurface* OpSurf_PMTPC = new G4OpticalSurface("OpSurf_PMTPC",glisur,polished,dielectric_metal);
	OpSurf_PMTPC->SetMaterialPropertiesTable(mpPMTPC);

	G4MaterialPropertiesTable* mpPMTHouse = this->MaterialPropertyTable_PMTHouse();
	vac_PMTHouse->SetMaterialPropertiesTable( mpPMTHouse );
	G4OpticalSurface* OpSurf_PMTHouse = new G4OpticalSurface("OpSurf_PMTHouse",unified,polished,dielectric_metal);
	OpSurf_PMTHouse->SetMaterialPropertiesTable(mpPMTHouse);




	// new G4Material("Lead"     , z=82., a=207.19*g/mole, density=11.35*g/cm3);
	// new G4Material("Copper"   , z=29., a=63.546*g/mole, density=8.96*g/cm3);

	// -- for PMT Glass -- //
	G4Material* Glass = new G4Material("Glass", density=1.032*g/cm3,2);
	Glass->AddElement(C,91.533*perCent);
	Glass->AddElement(H,8.467*perCent);

	G4MaterialPropertiesTable* mpGlass = this->MaterialPropertyTable_Glass();
	Glass->SetMaterialPropertiesTable(mpGlass);


	// -- for scintillation fiber core -- //
	G4Material* polystyrene
	= new G4Material("Polystyrene",density= 1.05*g/cm3, ncomponents=2);
	polystyrene->AddElement(C, natoms=8);
	polystyrene->AddElement(H, natoms=8);


	// -- for cladding (scintillation fibers) -- //
	G4Material* pmma_clad
	= new G4Material("PMMA_Clad",density= 1.19*g/cm3, ncomponents=3);
	pmma_clad->AddElement(C, natoms=5);
	pmma_clad->AddElement(H, natoms=8);
	pmma_clad->AddElement(O, natoms=2);


	// -- for Cerenkov fiber core -- //
	G4Material* pmma 
	= new G4Material("PMMA",density= 1.19*g/cm3, ncomponents=3);
	pmma->AddElement(C, natoms=5);
	pmma->AddElement(H, natoms=8);
	pmma->AddElement(O, natoms=2);

	G4MaterialPropertiesTable* mpPMMA = this->MaterialPropertyTable_PMMA();
	pmma->SetMaterialPropertiesTable(mpPMMA);


	// -- for cladding (Cerenkov fibers) -- //
	G4Material* fluorinatedPolymer 
	= new G4Material("Fluorinated_Polymer", density= 1.43*g/cm3, ncomponents=2);
	fluorinatedPolymer->AddElement(C,2);
	fluorinatedPolymer->AddElement(F,2);

	G4MaterialPropertiesTable* mpFS = this->MaterialPropertyTable_FS();
	fluorinatedPolymer->SetMaterialPropertiesTable(mpFS);


	// -- Air -- //
	G4Material* Air 
	= new G4Material("Air", density= 1.290*mg/cm3, ncomponents=2);
	Air->AddElement(N, fractionmass=0.7);
	Air->AddElement(O, fractionmass=0.3);

	G4MaterialPropertiesTable* mpAir = this->MaterialPropertyTable_Air();
	Air->SetMaterialPropertiesTable(mpAir);

	// G4MaterialPropertiesTable* mpPS
	// G4MaterialPropertiesTable* mpPMTPC;

	/////////////////
	// -- world -- //
	/////////////////
	G4VSolid* worldSolid 
	= new G4Box("worldBox",10.*m,10.*m,10.*m);
	G4LogicalVolume* worldLogical
	= new G4LogicalVolume(worldSolid,vac,"worldLogical");
	G4VPhysicalVolume* worldPhysical
	= new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0, false,0,checkOverlaps);

	//////////////////
	// -- Cu box -- //
	//////////////////
	G4double radius = 1.8*m; // -- size of inner tracker: it should be empty space to avoid any overlap with tracker system -- //
	G4double nTower_PhiDir = 283;
	// G4double nTower_PhiDir = 10;

	G4double pi = 3.14159265359;
	G4double dPhi = (2*pi) / nTower_PhiDir;
	G4double half_dPhi = 0.5*dPhi;
	G4double CuLen_PhiDir = 2*radius*std::tan(half_dPhi);
	G4double CuLen_EtaDir = CuLen_PhiDir;
	G4double CuLen_H = 2.5*m;

	bool DrawOneUnitTower = false;
	if( DrawOneUnitTower )
	{
		if( nTower_PhiDir != 283 ) cout << "nTower_PhiDir = " << nTower_PhiDir << " should be 283 for correct geometry" << endl;
		nTower_PhiDir = 1;
		CuLen_EtaDir = 40*mm;
		CuLen_PhiDir = CuLen_EtaDir;
		CuLen_H = 2.5*m;
	}

	cout << "[Cu] (PhiDir, EtaDir, Height) = (" << CuLen_PhiDir << ", " << CuLen_EtaDir << ", " << CuLen_H << ")" << endl;

	G4Box* CuBox 
	= new G4Box("CuBox", CuLen_EtaDir/2.0, CuLen_PhiDir/2., CuLen_H/2.0);
	G4LogicalVolume *CuLogical
	= new G4LogicalVolume(CuBox, cu, "CuLogical");

	///////////////////////////////////
	// -- Cu trapezoid (triangle) -- //
	///////////////////////////////////
	G4double CuTrdLen_PhiDir = 2*CuLen_H*std::sin(half_dPhi);
	G4double CuTrdLen_EtaDir = CuLen_EtaDir;
	G4double CuTrdLen_H = CuLen_H*std::cos(half_dPhi);
	G4Trd *CuTrd
	= new G4Trd("CuTrd", CuTrdLen_EtaDir/2.0, CuTrdLen_EtaDir/2.0, 0, CuTrdLen_PhiDir/2.0, CuTrdLen_H/2.0);
	G4LogicalVolume *CuTrdLogical
	= new G4LogicalVolume(CuTrd, cu, "CuTrdLogical");


	//////////////////
	// -- Fibers -- //
	//////////////////
	// -- Materials for Cerenkov fiber -- //
	G4Material *clad_C_Material = fluorinatedPolymer;
	G4Material *core_C_Material = pmma;

	// -- Materials for Scintillation fiber -- //
	G4Material *clad_S_Material = pmma_clad;
	G4Material *core_S_Material = polystyrene;

	// // -- Material for PMT glass -- //
	// G4Material *Glass_Material = Glass;

	// // -- Material for PMT Photocathod -- //
	// G4Material *PMTPC_Material = Al;

	// -- fibre parameters -- //
	G4double clad_C_rMin = 0*mm;
	G4double clad_C_rMax = 0.50*mm;
	G4double clad_C_Dz   = 2.5*m;
	G4double clad_C_Sphi = 0.;
	G4double clad_C_Dphi = 2.*M_PI;

	G4double core_C_rMin = 0.*mm;
	G4double core_C_rMax = 0.49*mm;
	G4double core_C_Dz   = 2.5*m;
	G4double core_C_Sphi = 0.;
	G4double core_C_Dphi = 2.*M_PI;

	// G4double clad_S_rMin = 0.485*mm;
	// G4double clad_S_rMax = 0.50*mm;
	// G4double clad_S_Dz   = 2.5*m;
	// G4double clad_S_Sphi = 0.;
	// G4double clad_S_Dphi = 2.*M_PI;

	G4double core_S_rMin = 0.*mm;
	G4double core_S_rMax = 0.485*mm;
	G4double core_S_Dz   = 2.5*m;
	G4double core_S_Sphi = 0.;
	G4double core_S_Dphi = 2.*M_PI;

	G4double dist_btwCore = 1.5*mm;

	bool Do_test = false;
	if( nTower_PhiDir == 10 ) Do_test = true;
	if( Do_test )
	{
		if( nTower_PhiDir != 10 )
			cout << "This test setting may not work for this # tower = " << nTower_PhiDir << endl;
		dist_btwCore = 250*mm; // -- when # tower = 10 -- //
		clad_C_rMax = dist_btwCore / 3.0;
		core_C_rMax = clad_C_rMax * 0.98;
		core_S_rMax = clad_C_rMax * 0.98;
	}

	const G4int nFiber_PhiDir = floor( CuLen_PhiDir / dist_btwCore ) - 1;
	const G4int nFiber_EtaDir = floor( CuLen_EtaDir / dist_btwCore ) - 1;
	G4double dist_edge_PhiDir = ( CuLen_PhiDir - (nFiber_PhiDir-1)*dist_btwCore ) / 2.0;
	G4double dist_edge_EtaDir = ( CuLen_EtaDir - (nFiber_EtaDir-1)*dist_btwCore ) / 2.0;
	//counting fiber for each direaction of triangle structure
	const G4int tri_nFiber_PhiDir = floor( CuTrdLen_PhiDir*std::cos(half_dPhi) / dist_btwCore ) - 1;
	const G4int tri_nFiber_EtaDir = floor( CuTrdLen_EtaDir / dist_btwCore ) - 1;
	G4double tri_dist_edge_PhiDir = ( CuTrdLen_PhiDir - (tri_nFiber_PhiDir-1)*dist_btwCore ) / 2.0;
	G4double tri_dist_edge_EtaDir = ( CuTrdLen_EtaDir*std::cos(half_dPhi) - (tri_nFiber_EtaDir-1)*dist_btwCore ) / 2.0;

	cout << "nFiber_PhiDir: " << nFiber_PhiDir << endl;
	cout << "nFiber_EtaDir: " << nFiber_EtaDir << endl;
	cout << "dist_edge_PhiDir: " << dist_edge_PhiDir << ", dist_edge_EtaDir: " << dist_edge_EtaDir << endl;
	cout << "tri_nFiber_PhiDir: " << tri_nFiber_PhiDir << endl;
	cout << "tri_nFiber_EtaDir: " << tri_nFiber_EtaDir << endl;
	cout << "tri_dist_edge_PhiDir: " << tri_dist_edge_PhiDir << ", dist_edge_EtaDir: " << tri_dist_edge_EtaDir << endl;

	// -- Solid -- //
	G4VSolid* fiberClad = new G4Tubs("fiberClad", 0.*mm, clad_C_rMax, CuLen_H/2., clad_C_Sphi, clad_C_Dphi); // -- S is the same -- //
	G4VSolid* fiberCoreC = new G4Tubs("fiberCoreC", 0.*mm, core_C_rMax, CuLen_H/2., core_C_Sphi, core_C_Dphi);
	G4VSolid* fiberCoreS = new G4Tubs("fiberCoreS", 0.*mm, core_S_rMax, CuLen_H/2., core_S_Sphi, core_S_Dphi);

	// -- PMT House, glass and PhotoCathod -- //
	G4double PMTHouseLen_H = 3*mm;
	G4double PMTHouseLen_EtaDir = CuLen_EtaDir;
	G4double PMTHouseLen_PhiDir = CuLen_PhiDir;

	G4double PMTGlassLen_H = 2*mm;
	G4double PMTGlassLen_EtaDir	= PMTHouseLen_EtaDir;
	G4double PMTGlassLen_PhiDir = PMTHouseLen_PhiDir;

	G4double PMTPCLen_H = 1*mm;
	G4double PMTPCLen_EtaDir = PMTHouseLen_EtaDir;
	G4double PMTPCLen_PhiDir = PMTHouseLen_PhiDir;

	if( Do_test )
	{
		PMTHouseLen_H *= 10;
		PMTGlassLen_H *= 10;
		PMTPCLen_H *= 10;
	}

	G4Box* PMTHouseBox 
	= new G4Box("PMTHouseBox", PMTHouseLen_EtaDir/2.0, PMTHouseLen_PhiDir/2., PMTHouseLen_H/2.0);
	G4LogicalVolume *PMTHouseBox_Logic
	= new G4LogicalVolume(PMTHouseBox, vac_PMTHouse, "PMTHouseBox_Logic");

	G4Box* PMTGlassBox 
	= new G4Box("PMTGlassBox", PMTGlassLen_EtaDir/2.0, PMTGlassLen_PhiDir/2., PMTGlassLen_H/2.0);
	G4LogicalVolume *PMTGlassBox_Logic
	= new G4LogicalVolume(PMTGlassBox, Glass, "PMTGlassBox_Logic");

	G4Box* PMTPCBox 
	= new G4Box("PMTPCBox", PMTPCLen_EtaDir/2.0, PMTPCLen_PhiDir/2., PMTPCLen_H/2.0);
	this->PMTPCBox_Logic
	= new G4LogicalVolume(PMTPCBox, Al, "PMTPCBox_Logic");

	//finding set of etas and saving it Eta
	// G4double Eta_Max=std::atan(2.5*m/radius);
	// G4double Eta[52];
	// Eta[0]=0;
	// const G4int nTower_EtaDir=51;
	// for (G4int i=1;i<=nTower_EtaDir;++i){
	// 	Eta[i]=Solve_Eq(Eta[i-1],Eta[i-1],Eta_Max,radius,CuLen_EtaDir);
	// }

	// G4double Theta_Max = Calc_ThetaMax_Barrel( 2.5*m, radius, CuLen_EtaDir );
	G4double Theta_Limit = M_PI / 2.0; // -- less than 90 degree -- //
	vector< G4double > vec_Theta;
	vec_Theta.push_back( 0.0 ); // -- first value: theta = eta = 0 -- //
	for(G4int i=0; i<100; i++) // -- up to 100: arbitrary number -- //
	{
		cout << "Inclined angle(degree) for " << i << "th tower in r-z plane: " << vec_Theta[i]*(180.0/M_PI) << endl; 
		G4double Theta_next = Solve_Eq(vec_Theta[i], vec_Theta[i], Theta_Limit, radius, CuLen_EtaDir);
		G4double length = radius*std::tan( Theta_next ) + CuLen_EtaDir / std::cos(Theta_next);

		if( length > 2.5*m )
			break;
		else
			vec_Theta.push_back( Theta_next );
	}

	G4int nTower_EtaDir = vec_Theta.size(); // -- including theta=0 tower -- //
	G4double Theta_Max = vec_Theta[nTower_EtaDir-1];
	G4double Eta_Max = (-1)*std::log( std::tan(Theta_Max / 2.0) );
	G4double length_Max = radius*std::tan( Theta_Max ) + CuLen_EtaDir / std::cos(Theta_Max);

	cout << "nTower_EtaDir: " << nTower_EtaDir << endl;
	cout << "Maximum inclined angle(degree) in barrel region: " << Theta_Max*(180.0/M_PI) << " (eta = " << Eta_Max << ")" << endl;
	cout << "\tCorresponding maximum length: " << length_Max << endl;

	// new G4LogicalSkinSurface("SkinSurf_PMTHouse", PMTHouseBox_Logic, OpSurf_PMTHouse);
	new G4LogicalSkinSurface("SkinSurf_PMTPC", PMTPCBox_Logic, OpSurf_PMTPC);

	// -- iteration for eta direction -- //
	// for(G4int i_barrel=-nTower_EtaDir-1; i_barrel<=nTower_EtaDir-1; i_barrel++) // -- nTower_EtaDir-1: # towers without the tower @ theta=0 -- //
	// {
	// 	G4int i_cu=0;
	// 	// -- iteration for phi direction -- //
	// 	// for(G4int i_cu=0; i_cu<nTower_PhiDir; i_cu++)
	// 	{
	// 		//////////////////
	// 		// -- Cu box -- //
	// 		//////////////////
	// 		G4double dTheta;
	// 		G4int I_factor; //correct error in sign of sin(dTheta)*tan(dTheta)
	// 		if (i_barrel<0) {
	// 			dTheta=-vec_Theta[-i_barrel];
	// 			I_factor=-1;
	// 		}
	// 		else {
	// 			dTheta=vec_Theta[i_barrel];
	// 			I_factor=1;
	// 		}

	// 		G4double phi = i_cu*dPhi;
	// 		//rotating matrix for Cu box
	// 		G4RotationMatrix rotM  = G4RotationMatrix();
	// 		rotM.rotateY(90*deg-dTheta);
	// 		rotM.rotateZ(phi);
	// 		//translational matrix for Cu box
	// 		G4double Phi_move=radius+0.5*CuLen_H*std::cos(dTheta)+I_factor*0.5*CuLen_EtaDir*std::sin(dTheta);
	// 		G4double Eta_move=radius*std::tan(dTheta)+0.5*CuLen_H*std::sin(dTheta)+I_factor*0.5*CuLen_EtaDir*std::sin(dTheta)*tan(dTheta);
	// 		G4ThreeVector Trans_Vector = G4ThreeVector(Phi_move*std::cos(phi), Phi_move*std::sin(phi),Eta_move);

	// 		G4Transform3D transform = G4Transform3D(rotM,Trans_Vector);
	// 		new G4PVPlacement(transform, CuLogical, "CuPhysical", worldLogical, false, i_cu, checkOverlaps );
	// 	}
	// }

	// -- iteration for phi direction -- //
	// i_cu = 0 // -- to test without phi rotation -- //
	for(G4int i_cu=0; i_cu<nTower_PhiDir; i_cu++)
	// for(G4int i_cu=0; i_cu<0; i_cu++)
	{
		G4int i_barrel = 0; // for a test -- //
		//////////////////
		// -- Cu box -- //
		//////////////////
		G4double dTheta;
		G4int I_factor; //correct error in sign of sin(dTheta)*tan(dTheta)
		if (i_barrel<0) {
			dTheta=-vec_Theta[-i_barrel];
			I_factor=-1;
		}
		else {
			dTheta=vec_Theta[i_barrel];
			I_factor=1;
		}

		G4double phi = i_cu*dPhi;
		G4RotationMatrix rotM  = G4RotationMatrix();
		rotM.rotateY(90*deg);
		rotM.rotateZ(phi);

		G4ThreeVector Unit_Z = G4ThreeVector(std::cos(phi),  std::sin(phi),0.);
		G4ThreeVector position = (radius + 0.5*CuLen_H)*Unit_Z; // -- multiply the size of the vector -- //
		G4Transform3D transform = G4Transform3D(rotM,position);

		new G4PVPlacement(transform, CuLogical, "CuPhysical", worldLogical, false, i_cu, checkOverlaps ); 

		// -- PMTs -- //
		G4ThreeVector position_PMTHouse = (radius + CuLen_H + 0.5*PMTHouseLen_H)*Unit_Z;
		G4Transform3D transform_PMTHouse = G4Transform3D(rotM,position_PMTHouse);
		new G4PVPlacement(transform_PMTHouse, PMTHouseBox_Logic, "PMTHouseBox_Phys", worldLogical, false, i_cu, checkOverlaps );

		new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5*PMTHouseLen_H + 0.5*PMTGlassLen_H), PMTGlassBox_Logic, "PMTGlassBox_Phys", PMTHouseBox_Logic, false, i_cu, checkOverlaps );
		new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5*PMTHouseLen_H - 0.5*PMTPCLen_H), this->PMTPCBox_Logic, "PMTPCBox_Phys", PMTHouseBox_Logic, false, i_cu, checkOverlaps );

		// -- fibers -- //
		G4int i_total = 0;
		for(G4int i_EtaDir=0; i_EtaDir<nFiber_EtaDir; i_EtaDir++)
		{
			G4double x_EtaDir = ((-1)*CuLen_EtaDir / 2.0) + dist_edge_EtaDir + dist_btwCore*i_EtaDir;
			for(G4int i_PhiDir=0; i_PhiDir<nFiber_PhiDir; i_PhiDir++)
			{
				i_total++;
				G4double x_PhiDir = ((-1)*CuLen_PhiDir / 2.0) + dist_edge_PhiDir + dist_btwCore*i_PhiDir;

				// -- cladding: same shape for both C and S fiber -- //
				G4VSolid* FiberClad_ith 
				= new G4IntersectionSolid("fiberClad", CuBox, fiberClad, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));

				G4LogicalVolume *FiberClad_Logic_ith
				= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");

				new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", CuLogical, false, i_total, checkOverlaps);

				// -- Cores -- //
				G4VSolid* FiberCore_ith;
				G4LogicalVolume *FiberCore_Logic_ith;

				// -- s, c, s, c, ... -- //
				// -- c, s, c, s, ... -- //
				bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);
				if( isFiberC )
				{
					FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuBox, fiberCoreC, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
					FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
				}
				else
				{
					FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuBox, fiberCoreS, 0, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
					FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
				}

				new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, i_total, checkOverlaps);

				G4VisAttributes* visAttr = new G4VisAttributes();
				if( isFiberC )
					visAttr->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
				else
					visAttr->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
				visAttr->SetForceSolid(true);
				visAttr->SetVisibility(true);
				FiberCore_Logic_ith->SetVisAttributes(visAttr);
			}
		}

		// ////////////////////////
		// // -- Cu trapezoid -- //
		// ////////////////////////
		// G4double Trd_phi = half_dPhi + i_cu*dPhi;
		// // G4double Trd_phi = half_dPhi + i_cu*dPhi + dPhi;
		// G4RotationMatrix Trd_rotM = G4RotationMatrix();
		// Trd_rotM.rotateY(90*deg-dTheta);
		// Trd_rotM.rotateZ(Trd_phi);
		
		// //this codes making overlap... need to be corrected
		// G4double move_tri=radius/std::cos(half_dPhi) + 0.5*CuLen_H*std::cos(dTheta)+0.5*CuLen_EtaDir*std::sin(dTheta);
		// G4double tri_z_move=radius*std::tan(dTheta)/std::cos(half_dPhi)+0.5*CuLen_H*std::sin(dTheta)+0.5*CuLen_EtaDir*std::sin(dTheta)*tan(dTheta);
		// G4ThreeVector Trd_Unit_Z = G4ThreeVector(move_tri*std::cos(Trd_phi), move_tri*std::sin(Trd_phi),tri_z_move);
		// //G4ThreeVector Trd_position = (radius/std::cos(half_dPhi) + 0.5*CuTrdLen_H)*Trd_Unit_Z; // -- multiply the size of the vector -- //
		// G4Transform3D Trd_transform = G4Transform3D(Trd_rotM,Trd_Unit_Z);
		
		// new G4PVPlacement(Trd_transform, CuTrdLogical, "CuTrdPhysical", worldLogical, false, i_cu, checkOverlaps );

		// // -- tri_fibers -- //
		// G4int tri_i_total = 0;
		// G4RotationMatrix *Trdfi_rotM = new G4RotationMatrix();
		// Trdfi_rotM->rotateX(-half_dPhi);
		// for(G4int i_EtaDir=0; i_EtaDir<nFiber_EtaDir; i_EtaDir++)
		// {
		// 	G4double x_EtaDir = ((-1)*CuTrdLen_EtaDir / 2.0) + tri_dist_edge_EtaDir + dist_btwCore*i_EtaDir;
		// 	for(G4int i_PhiDir=0; i_PhiDir<tri_nFiber_PhiDir; i_PhiDir++)
		// 	{
		// 		tri_i_total++;
		// 		G4double x_PhiDir =((-1)*CuTrdLen_PhiDir / 4.0) + tri_dist_edge_PhiDir/std::cos(half_dPhi) + dist_btwCore*i_PhiDir/std::cos(half_dPhi);
				
		// 		// -- cladding: same shape for both C and S fiber -- //
		// 		G4VSolid* FiberClad_ith
		// 		= new G4IntersectionSolid("fiberClad", CuTrd, fiberClad, Trdfi_rotM, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
				
		// 		G4LogicalVolume *FiberClad_Logic_ith
		// 		= new G4LogicalVolume(FiberClad_ith, clad_C_Material, "FiberClad_Logic");
				
		// 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberClad_Logic_ith, "FiberClad_Phys", CuTrdLogical, false, tri_i_total, checkOverlaps);
				
		// 		// -- Cores -- //
		// 		G4VSolid* FiberCore_ith;
		// 		G4LogicalVolume *FiberCore_Logic_ith;
				
		// 		// -- s, c, s, c, ... -- //
		// 		// -- c, s, c, s, ... -- //
		// 		bool isFiberC = this->IsFiberC(i_EtaDir, i_PhiDir);

		// 		if( isFiberC )
		// 		{
		// 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrd, fiberCoreC, Trdfi_rotM, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
		// 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_C_Material, "FiberCore_Logic");
		// 		}
		// 		else
		// 		{
		// 			FiberCore_ith = new G4IntersectionSolid( "fiberCore", CuTrd, fiberCoreS, Trdfi_rotM, G4ThreeVector(x_EtaDir, x_PhiDir, 0));
		// 			FiberCore_Logic_ith = new G4LogicalVolume(FiberCore_ith, core_S_Material, "FiberCore_Logic");
		// 		}
				
		// 		new G4PVPlacement(0, G4ThreeVector(0,0,0), FiberCore_Logic_ith, "FiberCore_Phys", FiberClad_Logic_ith, false, tri_i_total, checkOverlaps);
				
		// 		G4VisAttributes* visAttr1 = new G4VisAttributes();
		// 		if( isFiberC )
		// 			visAttr1->SetColour( G4Colour(0.0,0.0,1.0) ); // -- blue -- //
		// 		else
		// 			visAttr1->SetColour( G4Colour(1.0,1.0,0.0) );  // -- yellow -- //
		// 		visAttr1->SetForceSolid(true);
		// 		visAttr1->SetVisibility(true);
		// 		FiberCore_Logic_ith->SetVisAttributes(visAttr1);
		// 	}

		// }

		// ////////////////////////
		// // -- Cu trapezoid -- //
		// ////////////////////////
		// G4double Trd_phi = half_dPhi + i_cu*dPhi;
		// // G4double Trd_phi = half_dPhi + i_cu*dPhi + dPhi;
		// G4RotationMatrix Trd_rotM = G4RotationMatrix();
		// Trd_rotM.rotateY(90*deg);
		// Trd_rotM.rotateZ(Trd_phi);

		// G4ThreeVector Trd_Unit_Z = G4ThreeVector(std::cos(Trd_phi),  std::sin(Trd_phi),0.);
		// G4ThreeVector Trd_position = (radius/std::cos(half_dPhi) + 0.5*CuTrdLen_H)*Trd_Unit_Z; // -- multiply the size of the vector -- //
		// G4Transform3D Trd_transform = G4Transform3D(Trd_rotM,Trd_position);

		// new G4PVPlacement(Trd_transform, CuTrdLogical, "CuTrdPhysical", worldLogical, false, i_cu, checkOverlaps );




		// G4ThreeVector origin(x,y,z);
		// G4RotationMatrix* RotMatrix = new G4RotationMatrix();

		// RotMatrix->rotateZ(90*deg);
		// RotMatrix->rotateZ(-i*theta_unit_ZRot);
		// RotMatrix->rotateX(90*deg);
		// RotMatrix->rotateX(-theta_unit*(copyNo+0.5));

		// -- place it -- //
		// new G4PVPlacement( RotMatrix, G4ThreeVector(), CuLogical, "CuPhysical", worldLogical, false, 0, checkOverlaps );
	}

	// -- visualization -- //
	G4VisAttributes* visAttr;
	visAttr = new G4VisAttributes( G4Colour(0.0,1.0,1.0) ); // -- cyan -- //
	visAttr->SetVisibility(true);
	CuLogical->SetVisAttributes(visAttr);
	fVisAttributes.push_back(visAttr);

	G4VisAttributes* visAttr2 = new G4VisAttributes();
	visAttr2->SetColour( G4Colour(1.0,0.0,0.0) ); // -- red -- //
	visAttr2->SetForceSolid(true);
	visAttr2->SetVisibility(true);
	this->PMTPCBox_Logic->SetVisAttributes(visAttr2);
	fVisAttributes.push_back(visAttr2);

	G4VisAttributes* visAttr3 = new G4VisAttributes();
	visAttr3->SetColour( G4Colour(0.0,1.0,0.0) ); // -- green -- //
	// visAttr3->SetForceSolid(true);
	visAttr3->SetVisibility(true);
	PMTGlassBox_Logic->SetVisAttributes(visAttr3);
	fVisAttributes.push_back(visAttr3);

	return worldPhysical;
}

void WGR16DetectorConstruction::ConstructSDandField()
{
	G4SDManager* SDManager = G4SDManager::GetSDMpointer();

	G4String PMTName = "WGR16/PMTSD";
	G4String HitCollectionName = "PMTColl";
	G4VSensitiveDetector* SD_PMTPC = new WGR16PMTSD(PMTName, HitCollectionName);

	SDManager->AddNewDetector(SD_PMTPC);
	this->PMTPCBox_Logic->SetSensitiveDetector(SD_PMTPC);
}

void WGR16DetectorConstruction::ConstructMaterials()
{
	G4NistManager* nistManager = G4NistManager::Instance();
	
	// -- Vacuum: "Galactic" -- //
	nistManager->FindOrBuildMaterial("G4_Galactic");
	nistManager->FindOrBuildMaterial("G4_Cu");

	G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
	G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void WGR16DetectorConstruction::DefineCommands()
{
	// Define /WGR16/detector command directory using generic messenger class
	// fMessenger = new G4GenericMessenger(this, 
	//                                     "/WGR16/detector/", 
	//                                     "Detector control");

	// G4GenericMessenger::Command& B2ECollAngleCmd
	//   = fMessenger->DeclareMethodWithUnit("Block2Angle","deg",
	//       &WGR16DetectorConstruction::SetBlock2Angle,
	//       "Set rotation angle of exit collimator and corresponding detector");
	// B2ECollAngleCmd.SetParameterName("angle",true);
	// B2ECollAngleCmd.SetRange("angle>=-90. && angle<=90");
	// B2ECollAngleCmd.SetDefaultValue("60.");	
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_PMMA()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	// -- PMMA -- //
	G4double RefractiveIndex_PMMA[nEntries] =
	{
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
	};

	// G4double ABSLength_PMMA[nEntries] =
	// {
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m
	// };

	// G4double Reflectivity_PMMA[nEntries] =
	// {
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
	// };

	G4MaterialPropertiesTable* mpPMMA = new G4MaterialPropertiesTable();
	mpPMMA->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_PMMA,nEntries);
	//mpPMMA->AddProperty("ABSLENGTH",PhotonEnergy,ABSLength_PMMA,nEntries);
	//mpPMMA->AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity_PMMA,nEntries);

	return mpPMMA;
}


G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_FS()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//--- Fluorinated Polymer (FS) ---
	G4double RefractiveIndex_FluorinatedPolymer[nEntries] =
	{
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42,
	    1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42, 1.42
	};

	// G4double ABSLength_FS[nEntries] =
	// {
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m,
	//     8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m, 8.9*m
	// };

	// G4double Reflectivity_FS[nEntries] =
	// {
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
	//     1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
	// };

	G4MaterialPropertiesTable* mpFS = new G4MaterialPropertiesTable();
	mpFS->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_FluorinatedPolymer,nEntries);

	//mpFS->AddProperty("ABSLENGTH",PhotonEnergy,ABSLength_FS,nEntries);
	//mpFS->AddProperty("REFLECTIVITY",PhotonEnergy,Reflectivity_FS,nEntries);

	return mpFS;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_Glass()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//
	//Glass
	//
	G4double RefractiveIndex_Glass[nEntries] =
	{   1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49,
	    1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49, 1.49
	};
	
	// -- For test -- //
	// G4double RefractiveIndex_Glass[nEntries] =
	// {   1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	//     1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
	// };

	//Whenever I add some more detail optical properties, the GEANT4 will take lots of time to calculate
	//additional optical properties. I use only refractive indices to save time.
	// G4double Glass_AbsLength[nEntries] =
	// {
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	//     420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,420.*cm,
	// };

	G4MaterialPropertiesTable* mpGlass = new G4MaterialPropertiesTable();
	mpGlass->AddProperty("RINDEX",PhotonEnergy,RefractiveIndex_Glass,nEntries);
	//mpGlass->AddProperty("ABSLENGTH",PhotonEnergy,Glass_AbsLength,nEntries);

	return mpGlass;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_Air()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//
	//Air
	//
	G4double RefractiveIndex_Air[nEntries] =
	{
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
	};

	G4MaterialPropertiesTable* mpAir = new G4MaterialPropertiesTable();
	mpAir->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);

	return mpAir;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_PMTPC()
{
	G4double PhotonEnergy[] = {
	2.00*eV,2.03*eV,2.06*eV,2.09*eV,2.12*eV,
	2.15*eV,2.18*eV,2.21*eV,2.24*eV,2.27*eV,
	2.30*eV,2.33*eV,2.36*eV,2.39*eV,2.42*eV,
	2.45*eV,2.48*eV,2.51*eV,2.54*eV,2.57*eV,
	2.60*eV,2.63*eV,2.66*eV,2.69*eV,2.72*eV,
	2.75*eV,2.78*eV,2.81*eV,2.84*eV,2.87*eV,
	2.90*eV,2.93*eV,2.96*eV,2.99*eV,3.02*eV,
	3.05*eV,3.08*eV,3.11*eV,3.14*eV,3.17*eV,
	3.20*eV,3.23*eV,3.26*eV,3.29*eV,3.32*eV,
	3.35*eV,3.38*eV,3.41*eV,3.44*eV,3.47*eV};

	const G4int nEntries = sizeof(PhotonEnergy) / sizeof(G4double);

	//
	//Air
	//
	G4double RefractiveIndex_Air[nEntries] =
	{
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00
	};

	G4double p_mppc[2] = {2.00*eV, 3.47*eV};
	G4double refl_mppc[2] = {0.0, 0.0};
	G4double effi_mppc[2] = {0.11, 0.11}; // mimic Quantum Efficiency
	// G4double photocath_ReR[2] = {1.92, 1.92};
	// G4double photocath_ImR[2] = {1.69, 1.69};

	G4MaterialPropertiesTable* mpPMTPC = new G4MaterialPropertiesTable();
	mpPMTPC->AddProperty("REFLECTIVITY",p_mppc,refl_mppc,2);
	mpPMTPC->AddProperty("EFFICIENCY",p_mppc,effi_mppc,2);
	// mpPMTPC->AddProperty("REALINDEX",p_mppc,photocath_ReR,2);
	// mpPMTPC->AddProperty("IMAGINARYINDEX",p_mppc,photocath_ImR,2);
	mpPMTPC->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex_Air, nEntries);

	return mpPMTPC;
}

G4MaterialPropertiesTable* WGR16DetectorConstruction::MaterialPropertyTable_PMTHouse()
{
	G4double ephoton[] = {2.00*eV, 3.47*eV};
	const G4int num = sizeof(ephoton) / sizeof(G4double);

	G4double reficy[] = {0.0,0.0};
	G4double efficy[] = {1.0,1.0};

	G4MaterialPropertiesTable* mpPMTHouse = new G4MaterialPropertiesTable();
	mpPMTHouse->AddProperty("REFLECTIVITY",ephoton,reficy,num);
	mpPMTHouse->AddProperty("EFFICIENCY",ephoton,efficy,num);

	return mpPMTHouse;
}

bool WGR16DetectorConstruction::IsFiberC(G4int i_EtaDir, G4int i_PhiDir)
{
	bool Flag = false;
	if( i_EtaDir % 2 == 0 ) // start with c fiber -- //
	{
		if( i_PhiDir % 2 == 0 ) // -- c fiber -- //
			Flag = true;
		else
			Flag = false;
	}
	else // -- start with s fiber -- //
	{
		if( i_PhiDir % 2 == 0 ) // -- c fiber -- //
			Flag = false;
		else
			Flag = true;

	}

	return Flag;
}




