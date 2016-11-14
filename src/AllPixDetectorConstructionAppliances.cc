/**
 * Author: John Idarraga <idarraga@cern.ch> , 2010
 *
 */

#include "AllPixDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4UnionSolid.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"

#include"G4NistManager.hh"

void AllPixDetectorConstruction::BuildAppliances(int){

	switch(m_Appliances_type){
	case 0:
	{
		// Through the comand
		// --> /allpix/extras/setAppliancePosition
		// you can fill the vector "m_posVectorAppliances" available in this scope.
		// This vector holds the positions of the appliances volumes which can be placed with
		// respect to 'm_wrapper_log[detId]'.  This way your appliance properly rotates
		// with the detector.

		// Through the comand
		// --> /allpix/extras/setWrapperEnhancement
		// you can enhance the size of the wrapper so daughter volumens of 'm_wrapper_log[index]'
		// fit in.

		// Example
		G4cout << "Build Appliance " << G4endl;

		G4NistManager * nistman = G4NistManager::Instance();
		//G4Material * mat = nistman->FindOrBuildMaterial("G4_C");
		G4Material * mat = nistman->FindOrBuildMaterial("G4_Al");

		//G4_C
		//G4_Al
		//G4_LITHIUM_FLUORIDE
		//G4_POLYETHYLENE

		//std::vector<G4String> mats = nistman->GetNistMaterialNames();
		//std::vector<G4String>::iterator itr = mats.begin();
		//for( ; itr != mats.end() ; itr++) G4cout << *itr << G4endl;

//		G4Box * appliance_box = new G4Box("Appliance_box",
//				(100/2.)*mm,
//				(100/2.)*mm,
//				(20/2.)*mm);

		G4VisAttributes * applianceAtt = new G4VisAttributes(G4Color(0,0,0,0.6));
		applianceAtt->SetLineWidth(2);
		applianceAtt->SetForceSolid(true);
		applianceAtt->SetVisibility(true);

		// take into account wrapper for placemente
		map<int, G4ThreeVector>::iterator aplItr = m_posVectorAppliances.begin();
		int detId = 0;
		G4String appLogS = "Appliance_";
		G4String appPhysS = "Appliance_";
		char temp[128];

		G4Box *boxSup=new G4Box("boxSup",87*mm/2,79*mm/2,5*mm);
		G4Box *boxSupn=new G4Box("boxSupn",72*mm/2,54*mm/2,8*mm);
		G4Box *boxSupn2=new G4Box("boxSupn",52*mm/2,54*mm/2,5*mm);

		G4SubtractionSolid *supporttmp = new G4SubtractionSolid("BoxSup-BoxSupn",boxSup,boxSupn);
		G4SubtractionSolid *support = new G4SubtractionSolid("BoxSup-BoxSupn",supporttmp,boxSupn2,0,G4ThreeVector(0,44.5*mm,4*mm));

		for( ; aplItr != m_posVectorAppliances.end() ; aplItr++) {

			detId = (*aplItr).first;
			appLogS =  "Appliance_";
			appPhysS =  "Appliance_";
			sprintf(temp, "%d", detId);
			appLogS  += temp;
			appPhysS += temp;
			appLogS  += "_log";
			appPhysS += "_phys";

			G4LogicalVolume * appliance_log = new G4LogicalVolume(
					support,
					mat,
					appLogS);
			appliance_log->SetVisAttributes(applianceAtt);

			// If you mother vol is m_wrapper_log[detId],
			//  it will rotate with the medipix
			new G4PVPlacement(0,
					m_posVectorAppliances[detId]-G4ThreeVector(0,-10.25*mm,0*mm),
					appliance_log,
					appPhysS,
					m_wrapper_log[detId], // mother volume
					false,
					0,
					true);

			// Make this volume a sensitive device
			// I get automatically a hits ROOT file for each device
			//		sprintf(temp, "NeutronSD_%d", detId);
			//		G4String sdname = temp;
			//		AllPixTrackerSD * aTrackerSD = new AllPixTrackerSD( sdname,
			//				m_posVectorAppliances[detId],
			//				0);
			//
			//		SDman->AddNewDetector( aTrackerSD );
			//		appliance_log->SetSensitiveDetector( aTrackerSD );

		}
		break;
	}
	case 1:
	{
		// Aluminium box

		G4NistManager * nistman = G4NistManager::Instance();
		G4Material * alu = nistman->FindOrBuildMaterial("G4_Al");

		G4VisAttributes * applianceAtt = new G4VisAttributes(G4Color(255,0,0,0.7));
		applianceAtt->SetLineWidth(2);
		applianceAtt->SetForceSolid(true);
		applianceAtt->SetVisibility(true);

		G4Box * boxOut = new G4Box("outerBox", 54*mm/2,   94.25*mm/2, 12*mm/2);
		G4Box * boxIn  = new G4Box("innerBox", 52.5*mm/2, 92.5*mm/2,  12*mm/2);
		G4Box * window = new G4Box("window",   10*mm,     10*mm,      1.5*mm);

		G4SubtractionSolid * emptybox = new G4SubtractionSolid("emptybox",boxOut,boxIn,0,G4ThreeVector(0,0,-1.5*mm));
		G4SubtractionSolid * box = new G4SubtractionSolid("box",emptybox,window,0,G4ThreeVector(0,-22.25*mm,6*mm));

		// take into account wrapper for placement
		map<int, G4ThreeVector>::iterator aplItr = m_posVectorAppliances.begin();
		int detId = 0;
		G4String appLogS = "Appliance_";
		G4String appPhysS = "Appliance_";
		char temp[128];

		for( ; aplItr != m_posVectorAppliances.end() ; aplItr++) {

			detId = (*aplItr).first;
			appLogS =  "Appliance_";
			appPhysS =  "Appliance_";
			sprintf(temp, "%d", detId);
			appLogS  += temp;
			appPhysS += temp;
			appLogS  += "_log";
			appPhysS += "_phys";

			G4LogicalVolume * box_log = new G4LogicalVolume(box,alu,appLogS);
			box_log->SetVisAttributes(applianceAtt);

			new G4PVPlacement(0,
					m_posVectorAppliances[detId]+G4ThreeVector(0,0*mm,11.12*mm),
					box_log,
					appPhysS,
					m_wrapper_log[detId], // mother volume
					false,
					0,
					true);
		}
		break;
	}
	
	case 99 :
	{
		int detId = 350;
	
		G4cout << "Building Alibava box" << G4endl ;

		// materials
		G4NistManager* nistman = G4NistManager::Instance();
		G4Material * kaptonmaterial = nistman->FindOrBuildMaterial("G4_KAPTON");

		// create the armaflex material
		G4double z, a, fractionmass, density;
		G4String name, symbol;
		G4int ncomponents;

		a = 28.09*g/mole;
		G4Element* elSi  = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

		a = 16.00*g/mole;
		G4Element* elO  = new G4Element(name="Oxygen"  ,symbol="O" , z= 8., a);

		a = 1.01*g/mole;
		G4Element* elH  = new G4Element(name="Hydrogen"  ,symbol="H" , z= 1., a);

		a = 12.01*g/mole;
		G4Element* elC  = new G4Element(name="Carbon"  ,symbol="C" , z= 6., a);

		density = 85.629*mg/cm3;
		G4Material* armaflexmaterial = new G4Material(name="Armaflex",density,ncomponents=4);
		armaflexmaterial->AddElement(elSi, fractionmass=37.9*perCent);
		armaflexmaterial->AddElement(elO, fractionmass=21.6*perCent);
		armaflexmaterial->AddElement(elH, fractionmass=8.2*perCent);
		armaflexmaterial->AddElement(elC, fractionmass=32.4*perCent);

		// sizes in mm
		G4double kapthick = 0.1;
		G4double armathick = 5.0;
		G4double innerdim = 25.0;

		// first create the outer kapton box
		G4Box *kaptonbox=new G4Box("KaptonBox",(2.0*innerdim+armathick+kapthick)*mm,(3.0*innerdim+armathick+kapthick)*mm,(innerdim+armathick+kapthick)*mm);

		// then the armaflex box
		G4Box *armabox=new G4Box("ArmaBox",(2.0*innerdim+armathick)*mm,(3.0*innerdim+armathick)*mm,(innerdim+armathick)*mm);

		// remove the armaflex box volume from the kaptonbox volume
		G4SubtractionSolid *outerbox = new G4SubtractionSolid("KaptonBox-ArmaBox",kaptonbox,armabox);

		// the empty inner box to be removed from the armaflex box
		G4Box *inside=new G4Box("InnerBox",(2.0*innerdim)*mm,(3.0*innerdim)*mm,(innerdim)*mm);

		G4SubtractionSolid *innerbox = new G4SubtractionSolid("OuterBox-InnerBox",outerbox,inside);

		m_TestStructure_log = new G4LogicalVolume(outerbox,kaptonmaterial,"outerbox",0,0,0);
		G4VisAttributes * kapvis = new G4VisAttributes(G4Color(1,1,1,0.5));
		kapvis->SetLineWidth(1);
		kapvis->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(kapvis);
		G4VisAttributes * armavis = new G4VisAttributes(G4Color(1,1,1,0.5));
		armavis->SetLineWidth(1);
		armavis->SetForceSolid(true);
		m_TestStructure_log->SetVisAttributes(armavis);

		m_TestStructure_log = new G4LogicalVolume(innerbox,armaflexmaterial,"innerbox",0,0,0);

		m_TestStructure_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0),m_TestStructure_log,"ali_phys",m_wrapper_log[detId],false,0,true);

		G4cout << "Done Building Alibava box! Materials used:" << G4endl ;

		G4cout << armaflexmaterial << G4endl;

		G4cout << kaptonmaterial << G4endl;

	}

	default:
	{
		G4cout << "Unknown Appliance Type" << G4endl;
		break;
	}
	}
	

}
