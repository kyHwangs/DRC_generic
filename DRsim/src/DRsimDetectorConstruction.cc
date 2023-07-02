#include "DRsimDetectorConstruction.hh"
#include "DRsimCellParameterisation.hh"
#include "DRsimFilterParameterisation.hh"
#include "DRsimMirrorParameterisation.hh"
#include "DRsimSiPMSD.hh"

#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4IntersectionSolid.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4GeometryManager.hh"

#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

using namespace std;

G4ThreadLocal DRsimMagneticField* DRsimDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* DRsimDetectorConstruction::fFieldMgr = 0;

int DRsimDetectorConstruction::fNofRow = 1;
int DRsimDetectorConstruction::fNofModules = fNofRow * fNofRow;

DRsimDetectorConstruction::DRsimDetectorConstruction()
: G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  clad_C_rMin = 0.49*mm;
  clad_C_rMax = 0.50*mm;
  clad_C_Dz   = 2.5*m;
  clad_C_Sphi = 0.;
  clad_C_Dphi = 2.*M_PI;

  core_C_rMin = 0.*mm;
  core_C_rMax = 0.49*mm;
  core_C_Dz   = 2.5*m;
  core_C_Sphi = 0.;
  core_C_Dphi = 2.*M_PI;

  clad_S_rMin = 0.485*mm;
  clad_S_rMax = 0.50*mm;
  clad_S_Dz   = 2.5*m;
  clad_S_Sphi = 0.;
  clad_S_Dphi = 2.*M_PI;

  core_S_rMin = 0.*mm;
  core_S_rMax = 0.485*mm;
  core_S_Dz   = 2.5*m;
  core_S_Sphi = 0.;
  core_S_Dphi = 2.*M_PI;

  PMTT = 0.3*mm;
  filterT = 0.01*mm;
  reflectorT = 0.03*mm;

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0,0.5,0.,1.0));
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0.,0.,1.0,1.0));
  fVisAttrBlue->SetVisibility(true);
  fVisAttrGray = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.3));
  fVisAttrGray->SetVisibility(true);
  fVisAttrGreen = new G4VisAttributes(G4Colour(0.3,0.7,0.3));
  fVisAttrGreen->SetVisibility(true);
}

DRsimDetectorConstruction::~DRsimDetectorConstruction() {
  delete fMessenger;
  delete fMaterials;

  delete fVisAttrOrange;
  delete fVisAttrBlue;
  delete fVisAttrGray;
  delete fVisAttrGreen;
}

void DRsimDetectorConstruction::DefineMaterials() {
  fMaterials = DRsimMaterials::GetInstance();
}

G4VPhysicalVolume* DRsimDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  checkOverlaps = false;

  G4VSolid* worldSolid             = new G4Box("worldBox",10.*m,10.*m,10.*m);
  worldLogical                     = new G4LogicalVolume(worldSolid,FindMaterial("G4_Galactic"),"worldLogical");
  G4VPhysicalVolume* worldPhysical = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,false,0,checkOverlaps);

  float moduleUnitDimension = 122.;

  fFrontL     = 0.;
  fTowerDepth = 100.; 
  fModuleH    = moduleUnitDimension;
  fModuleW    = moduleUnitDimension;
  fFiberUnitH = 1.;


  // fRandomSeed = 1;

  doFiber     = true;
  doReflector = false;
  doPMT       = false;
  doPlace     = true;

  fiberUnit   = new G4Box("fiber_SQ", (fFiberUnitH/2) *mm, (1./2) *mm, (fTowerDepth/2) *mm);
  fiberClad   = new G4Tubs("fiber",  0, clad_C_rMax, fTowerDepth/2., 0 *deg, 360. *deg);   // S is the same
  fiberCoreC  = new G4Tubs("fiberC", 0, core_C_rMax, fTowerDepth/2., 0 *deg, 360. *deg);
  fiberCoreS  = new G4Tubs("fiberS", 0, core_S_rMax, fTowerDepth/2., 0 *deg, 360. *deg);

  dimCalc = new dimensionCalc();
  dimCalc->SetFrontL(fFrontL);
  dimCalc->SetTower_height(fTowerDepth);
  dimCalc->SetPMTT(PMTT+filterT);
  dimCalc->SetReflectorT(reflectorT);
  dimCalc->SetNofModules(fNofModules);
  dimCalc->SetNofRow(fNofRow);
  dimCalc->SetModuleHeight(fModuleH);
  dimCalc->SetModuleWidth(fModuleW);

  ModuleBuild(ModuleLogical,PMTGLogical,PMTfilterLogical,PMTcellLogical,PMTcathLogical,ReflectorMirrorLogical,fiberUnitIntersection,fiberCladIntersection,fiberCoreIntersection,fModuleProp);

  delete dimCalc;
  return worldPhysical;
}

void DRsimDetectorConstruction::ConstructSDandField() {
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4String SiPMName = "SiPMSD";

  // ! Not a memory leak - SDs are deleted by G4SDManager. Deleting them manually will cause double delete!
  if ( doPMT ) {
    for (int i = 0; i < fNofModules; i++) {
      DRsimSiPMSD* SiPMSDmodule = new DRsimSiPMSD("Module"+std::to_string(i), "ModuleC"+std::to_string(i), fModuleProp.at(i));
      SDman->AddNewDetector(SiPMSDmodule);
      PMTcathLogical[i]->SetSensitiveDetector(SiPMSDmodule);
    }
  }
}

void DRsimDetectorConstruction::ModuleBuild(G4LogicalVolume* ModuleLogical_[], 
                                            G4LogicalVolume* PMTGLogical_[], 
                                            G4LogicalVolume* PMTfilterLogical_[], 
                                            G4LogicalVolume* PMTcellLogical_[], 
                                            G4LogicalVolume* PMTcathLogical_[], 
                                            G4LogicalVolume* ReflectorMirrorLogical_[],
                                            std::vector<G4LogicalVolume*> fiberUnitIntersection_[], 
                                            std::vector<G4LogicalVolume*> fiberCladIntersection_[], 
                                            std::vector<G4LogicalVolume*> fiberCoreIntersection_[], 
                                            std::vector<DRsimInterface::DRsimModuleProperty>& ModuleProp_)
{

  std::vector<float> fModuleWidth;
  std::vector<float> fModuleHeight;
  std::vector<float> fModuleDepth;

  // %%%%% Module 0 : Center of 3D printing
  G4Box* solid_Module0 = new G4Box("Module0", 30.*mm / 2.,30.*mm / 2.,500.*mm / 2.);
  fModuleSolid[0] = static_cast<G4VSolid*>(solid_Module0);
  ModuleLogical_[0] = new G4LogicalVolume(fModuleSolid[0], FindMaterial("Copper"), "Module0");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(0), dimCalc->GetOrigin(0), ModuleLogical_[0], "Module0", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(30.);
  fModuleWidth.push_back(30.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 1 : Down SFHS
  G4Box* solid_Module1 = new G4Box("Module1", 36.*mm / 2.,45.*mm / 2.,500.*mm / 2.);
  fModuleSolid[1] = static_cast<G4VSolid*>(solid_Module1);
  ModuleLogical_[1] = new G4LogicalVolume(fModuleSolid[1], FindMaterial("Copper"), "Module1");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(1), dimCalc->GetOrigin(1), ModuleLogical_[1], "Module1", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(36.);
  fModuleWidth.push_back(45.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 2 : Left SFHS
  G4Box* solid_Module2 = new G4Box("Module2", 45.*mm / 2.,36.*mm / 2.,500.*mm / 2.);
  fModuleSolid[2] = static_cast<G4VSolid*>(solid_Module2);
  ModuleLogical_[2] = new G4LogicalVolume(fModuleSolid[2], FindMaterial("Copper"), "Module2");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(2), dimCalc->GetOrigin(2), ModuleLogical_[2], "Module2", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(45.);
  fModuleWidth.push_back(36.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 3 : Right SFHS
  G4Box* solid_Module3 = new G4Box("Module3", 45.*mm / 2.,36.*mm / 2.,500.*mm / 2.);
  fModuleSolid[3] = static_cast<G4VSolid*>(solid_Module3);
  ModuleLogical_[3] = new G4LogicalVolume(fModuleSolid[3], FindMaterial("Copper"), "Module3");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(3), dimCalc->GetOrigin(3), ModuleLogical_[3], "Module3", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(45.);
  fModuleWidth.push_back(36.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 4 : Up SFHS Wedge
  G4ThreeVector pt_Module4[8] = {G4ThreeVector()};
  pt_Module4[0] = G4ThreeVector(-22.5, -15., -250.1125255);
  pt_Module4[1] = G4ThreeVector(22.5, -15., -250.1125255);
  pt_Module4[2] = G4ThreeVector(-22.5, 15., -250.1125255);
  pt_Module4[3] = G4ThreeVector(22.5, 15., -250.1125255);
  pt_Module4[4] = G4ThreeVector(-22.5, -22.5, 250.1125255);
  pt_Module4[5] = G4ThreeVector(22.5, -22.5, 250.1125255);
  pt_Module4[6] = G4ThreeVector(-22.5, 22.5, 250.1125255);
  pt_Module4[7] = G4ThreeVector(22.5, 22.5, 250.1125255);

  G4Trap* solid_Module4 = new G4Trap("Module4", pt_Module4);
  fModuleSolid[4] = static_cast<G4VSolid*>(solid_Module4);
  ModuleLogical_[4] = new G4LogicalVolume(fModuleSolid[4], FindMaterial("Copper"), "Module4");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(4), dimCalc->GetOrigin(4), ModuleLogical_[4], "Module4", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(45.);
  fModuleWidth.push_back(45.);
  fModuleDepth.push_back(250.1125255 * 2);

  // %%%%% Module 5 : Right Down LEGO
  G4Box* solid_Module5 = new G4Box("Module5", 24.*mm / 2.,24.*mm / 2.,500.*mm / 2.);
  fModuleSolid[5] = static_cast<G4VSolid*>(solid_Module5);
  ModuleLogical_[5] = new G4LogicalVolume(fModuleSolid[5], FindMaterial("Copper"), "Module5");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(5), dimCalc->GetOrigin(5), ModuleLogical_[5], "Module5", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(24.);
  fModuleWidth.push_back(24.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 6 : Left Down LEGO
  G4Box* solid_Module6 = new G4Box("Module6", 24.*mm / 2.,24.*mm / 2.,500.*mm / 2.);
  fModuleSolid[6] = static_cast<G4VSolid*>(solid_Module6);
  ModuleLogical_[6] = new G4LogicalVolume(fModuleSolid[6], FindMaterial("Copper"), "Module6");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(6), dimCalc->GetOrigin(6), ModuleLogical_[6], "Module6", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(24.);
  fModuleWidth.push_back(24.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 7 : Right Up LEGO
  G4Box* solid_Module7 = new G4Box("Module7", 24.*mm / 2.,24.*mm / 2.,500.*mm / 2.);
  fModuleSolid[7] = static_cast<G4VSolid*>(solid_Module7);
  ModuleLogical_[7] = new G4LogicalVolume(fModuleSolid[7], FindMaterial("Copper"), "Module7");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(7), dimCalc->GetOrigin(7), ModuleLogical_[7], "Module7", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(24.);
  fModuleWidth.push_back(24.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 8 : Left Up LEGO
  G4Box* solid_Module8 = new G4Box("Module8", 24.*mm / 2.,24.*mm / 2.,500.*mm / 2.);
  fModuleSolid[8] = static_cast<G4VSolid*>(solid_Module8);
  ModuleLogical_[8] = new G4LogicalVolume(fModuleSolid[8], FindMaterial("Copper"), "Module8");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(8), dimCalc->GetOrigin(8), ModuleLogical_[8], "Module8", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(24.);
  fModuleWidth.push_back(24.);
  fModuleDepth.push_back(500.);

  // %%%%% Module 9 : Wing of 3D printing

  G4ThreeVector pt_Module9[8] = {G4ThreeVector()};
  pt_Module9[0] = G4ThreeVector(-15., -15., -250.);
  pt_Module9[1] = G4ThreeVector(15., -15., -250.);
  pt_Module9[2] = G4ThreeVector(-15., 15., -250.);
  pt_Module9[3] = G4ThreeVector(15., 15., -250.);
  pt_Module9[4] = G4ThreeVector(-30., -22.5, 250.);
  pt_Module9[5] = G4ThreeVector(30., -22.5, 250.);
  pt_Module9[6] = G4ThreeVector(-30., 22.5, 250.);
  pt_Module9[7] = G4ThreeVector(30., 22.5, 250.);

  G4Trap* aTrapTmp = new G4Trap("tmpTrap", pt_Module9);
  G4Box* aBoxTmp = new G4Box("tmpBox", 100.*mm / 2.,100.*mm / 2.,600.*mm / 2.);
  G4SubtractionSolid* subSolidTmp = new G4SubtractionSolid("tmpSubtrack", aTrapTmp, aBoxTmp, new G4RotationMatrix(), G4ThreeVector(-65., 0, 0));
  G4Box* centreBoxTmp = new G4Box("tmpCenterBox", 30.*mm / 2.,30.*mm / 2.,600.*mm / 2.);

  G4Box* removeBox = new G4Box("removeBox", 3.*mm / 2.,3.*mm / 2.,600.*mm / 2.);
  G4SubtractionSolid* solid_Wing_origin = new G4SubtractionSolid("Wing_origin", subSolidTmp, centreBoxTmp, new G4RotationMatrix(), G4ThreeVector());

  G4Box* removeBigBox = new G4Box("removeBox", 70.*mm / 2.,70.*mm / 2.,600.*mm / 2.);
  G4Box* removeBigBoxWidth = new G4Box("removeBox", 30.*mm / 2.,70.*mm / 2.,600.*mm / 2.);
  G4Box* removeSmallBox = new G4Box("removeBox", 15.*mm / 2.,3.*mm / 2.,600.*mm / 2.);

  G4SubtractionSolid* solid_Wing_top_1 = new G4SubtractionSolid("Wing_top_1", solid_Wing_origin, removeBigBoxWidth, new G4RotationMatrix(), G4ThreeVector(0., 0., 0.));
  G4SubtractionSolid* solid_Wing_top_2  = new G4SubtractionSolid("Wing_top_2", solid_Wing_top_1, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 6., 1.5 * -11., 0.));
  G4SubtractionSolid* solid_Wing_top_3  = new G4SubtractionSolid("Wing_top_3", solid_Wing_top_2, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 8., 1.5 * -12., 0.));
  G4SubtractionSolid* solid_Wing_top_4  = new G4SubtractionSolid("Wing_top_4", solid_Wing_top_3, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 10., 1.5 * -13., 0.));
  G4SubtractionSolid* solid_Wing_top_5  = new G4SubtractionSolid("Wing_top_5", solid_Wing_top_4, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 12., 1.5 * -14., 0.));
  G4SubtractionSolid* solid_Wing_top_6  = new G4SubtractionSolid("Wing_top_6", solid_Wing_top_5, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 14., 1.5 * -15., 0.));
  G4SubtractionSolid* solid_Wing_top_7  = new G4SubtractionSolid("Wing_top_7", solid_Wing_top_6, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 6., 1.5 * 11., 0.));
  G4SubtractionSolid* solid_Wing_top_8  = new G4SubtractionSolid("Wing_top_8", solid_Wing_top_7, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 8., 1.5 * 12., 0.));
  G4SubtractionSolid* solid_Wing_top_9  = new G4SubtractionSolid("Wing_top_9", solid_Wing_top_8, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 10., 1.5 * 13., 0.));
  G4SubtractionSolid* solid_Wing_top_10 = new G4SubtractionSolid("Wing_top_10", solid_Wing_top_9, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 12., 1.5 * 14., 0.));
  
  G4SubtractionSolid* solid_Module9 = new G4SubtractionSolid("Module9", solid_Wing_top_10, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 14., 1.5 * 15., 0.));
  fModuleSolid[9] = static_cast<G4VSolid*>(solid_Module9);
  ModuleLogical_[9] = new G4LogicalVolume(fModuleSolid[9], FindMaterial("Copper"), "Module9");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(9), dimCalc->GetOrigin(9), ModuleLogical_[9], "Module9", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(60.);
  fModuleWidth.push_back(45.);
  fModuleDepth.push_back(500.);

  G4SubtractionSolid* solid_Wing_left_1 = new G4SubtractionSolid("Wing_left_1", solid_Wing_origin, removeBigBox, new G4RotationMatrix(), G4ThreeVector(0., 20., 0.));
  G4SubtractionSolid* solid_Wing_left_2 = new G4SubtractionSolid("Wing_left_2", solid_Wing_left_1, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 16., 1.5 * -10., 0.));
  G4SubtractionSolid* solid_Wing_left_3 = new G4SubtractionSolid("Wing_left_3", solid_Wing_left_2, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 18., 1.5 * -11., 0.));
  G4SubtractionSolid* solid_Wing_left_4 = new G4SubtractionSolid("Wing_left_4", solid_Wing_left_3, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 20., 1.5 * -12., 0.));
  G4SubtractionSolid* solid_Wing_left_5 = new G4SubtractionSolid("Wing_left_5", solid_Wing_left_4, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 22., 1.5 * -13., 0.));
  G4SubtractionSolid* solid_Wing_left_6 = new G4SubtractionSolid("Wing_left_6", solid_Wing_left_5, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 24., 1.5 * -14., 0.));
  
  G4SubtractionSolid* solid_Module10 = new G4SubtractionSolid("Module10", solid_Wing_left_6, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 24., 1.5 * -15., 0.));
  fModuleSolid[10] = static_cast<G4VSolid*>(solid_Module10);
  ModuleLogical_[10] = new G4LogicalVolume(fModuleSolid[10], FindMaterial("Copper"), "Module10");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(9), dimCalc->GetOrigin(9), ModuleLogical_[10], "Module10", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(60.);
  fModuleWidth.push_back(45.);
  fModuleDepth.push_back(500.);

  G4SubtractionSolid* solid_Wing_right_1 = new G4SubtractionSolid("Wing_right_1", solid_Wing_origin, removeBigBox, new G4RotationMatrix(), G4ThreeVector(0., -20., 0.));
  G4SubtractionSolid* solid_Wing_right_2 = new G4SubtractionSolid("Wing_right_2", solid_Wing_right_1, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 16., 1.5 * 10., 0.));
  G4SubtractionSolid* solid_Wing_right_3 = new G4SubtractionSolid("Wing_right_3", solid_Wing_right_2, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 18., 1.5 * 11., 0.));
  G4SubtractionSolid* solid_Wing_right_4 = new G4SubtractionSolid("Wing_right_4", solid_Wing_right_3, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 20., 1.5 * 12., 0.));
  G4SubtractionSolid* solid_Wing_right_5 = new G4SubtractionSolid("Wing_right_5", solid_Wing_right_4, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 22., 1.5 * 13., 0.));
  G4SubtractionSolid* solid_Wing_right_6 = new G4SubtractionSolid("Wing_right_6", solid_Wing_right_5, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 24., 1.5 * 14., 0.));
  
  G4SubtractionSolid* solid_Module11 = new G4SubtractionSolid("Module11", solid_Wing_right_6, removeSmallBox, new G4RotationMatrix(), G4ThreeVector(1.5 * 24., 1.5 * 15., 0.));
  fModuleSolid[11] = static_cast<G4VSolid*>(solid_Module11);
  ModuleLogical_[11] = new G4LogicalVolume(fModuleSolid[11], FindMaterial("Copper"), "Module11");
  if (doPlace) new G4PVPlacement(dimCalc->GetRM(9), dimCalc->GetOrigin(9), ModuleLogical_[11], "Module11", worldLogical, false, 0, checkOverlaps);
  fModuleHeight.push_back(60.);
  fModuleWidth.push_back(45.);
  fModuleDepth.push_back(500.);

  FiberImplementSingle(ModuleLogical_, fModuleSolid, fModuleWidth, fModuleHeight, fModuleDepth, fiberUnitIntersection_, fiberCladIntersection_, fiberCoreIntersection_);

  //   dimCalc->SetisModule(true);
  //   module = new G4Box("Mudule", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (fTowerDepth/2.) *mm );
  //   ModuleLogical_[i] = new G4LogicalVolume(module,FindMaterial("Copper"),moduleName);
  //   // G4VPhysicalVolume* modulePhysical = new G4PVPlacement(0,dimCalc->GetOrigin(i),ModuleLogical_[i],moduleName,worldLogical,false,0,checkOverlaps);
  //   new G4PVPlacement(0,dimCalc->GetOrigin(i),ModuleLogical_[i],moduleName,worldLogical,false,0,checkOverlaps);

  //   if ( doPMT ) {
  //     dimCalc->SetisModule(false);  
  //     pmtg = new G4Box("PMTG", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (PMTT+filterT)/2. *mm );
  //     PMTGLogical_[i]  = new G4LogicalVolume(pmtg,FindMaterial("G4_AIR"),moduleName);
  //     new G4PVPlacement(0,dimCalc->GetOrigin_PMTG(i),PMTGLogical_[i],moduleName,worldLogical,false,0,checkOverlaps);
  //   }

  //   FiberImplement(i,ModuleLogical_,fiberUnitIntersection_,fiberCladIntersection_,fiberCoreIntersection_);

  //   DRsimInterface::DRsimModuleProperty ModulePropSingle;
  //   ModulePropSingle.towerXY   = fTowerXY;
  //   ModulePropSingle.ModuleNum = i;
  //   ModuleProp_.push_back(ModulePropSingle);

    // if ( doPMT ) {
    //   G4VSolid* SiPMlayerSolid = new G4Box("SiPMlayerSolid", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (PMTT/2.) *mm );
    //   G4LogicalVolume* SiPMlayerLogical = new G4LogicalVolume(SiPMlayerSolid,FindMaterial("G4_AIR"),"SiPMlayerLogical");
    //   new G4PVPlacement(0,G4ThreeVector(0.,0.,filterT/2.),SiPMlayerLogical,"SiPMlayerPhysical",PMTGLogical_[i],false,0,checkOverlaps);

    //   G4VSolid* filterlayerSolid = new G4Box("filterlayerSolid", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (filterT/2.) *mm );
    //   G4LogicalVolume* filterlayerLogical = new G4LogicalVolume(filterlayerSolid,FindMaterial("Glass"),"filterlayerLogical");
    //   new G4PVPlacement(0,G4ThreeVector(0.,0.,-PMTT/2.),filterlayerLogical,"filterlayerPhysical",PMTGLogical_[i],false,0,checkOverlaps);

    //   G4VSolid* PMTcellSolid = new G4Box("PMTcellSolid", 1.2/2. *mm, 1.2/2. *mm, PMTT/2. *mm );
    //   PMTcellLogical_[i] = new G4LogicalVolume(PMTcellSolid,FindMaterial("Glass"),"PMTcellLogical_");

    //   // DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fTowerXY.first,fTowerXY.second);
    //   DRsimCellParameterisation* PMTcellParam = new DRsimCellParameterisation(fFiberX, fFiberY, fFiberWhich);
    //   G4PVParameterised* PMTcellPhysical = new G4PVParameterised("PMTcellPhysical",PMTcellLogical_[i],SiPMlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second,PMTcellParam);

    //   G4VSolid* PMTcathSolid = new G4Box("PMTcathSolid", 1.2/2. *mm, 1.2/2. *mm, filterT/2. *mm );
    //   PMTcathLogical_[i] = new G4LogicalVolume(PMTcathSolid,FindMaterial("Silicon"),"PMTcathLogical_");
    //   new G4PVPlacement(0,G4ThreeVector(0.,0.,(PMTT-filterT)/2.*mm),PMTcathLogical_[i],"PMTcathPhysical",PMTcellLogical_[i],false,0,checkOverlaps);
    //   new G4LogicalSkinSurface("Photocath_surf",PMTcathLogical_[i],FindSurface("SiPMSurf"));

    //   G4VSolid* filterSolid = new G4Box("filterSolid", 1.2/2. *mm, 1.2/2. *mm, filterT/2. *mm );
    //   PMTfilterLogical_[i] = new G4LogicalVolume(filterSolid,FindMaterial("Gelatin"),"PMTfilterLogical_");

    //   int filterNo = (int)(fTowerXY.first * fTowerXY.second) / 2;
    //   if ( fTowerXY.first % 2 == 1)
    //     filterNo++;

    //   // DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fTowerXY.first,fTowerXY.second);
    //   DRsimFilterParameterisation* filterParam = new DRsimFilterParameterisation(fFiberX, fFiberY, fFiberWhich);
    //   G4PVParameterised* filterPhysical = new G4PVParameterised("filterPhysical",PMTfilterLogical_[i],filterlayerLogical,kXAxis,filterNo,filterParam);
    //   new G4LogicalBorderSurface("filterSurf",filterPhysical,PMTcellPhysical,FindSurface("FilterSurf"));
          
    //   PMTcathLogical_[i]->SetVisAttributes(fVisAttrGreen);
    //   PMTfilterLogical_[i]->SetVisAttributes(fVisAttrOrange);
    // }

    // if ( doReflector ) {
    //   G4VSolid* ReflectorlayerSolid = new G4Box("ReflectorlayerSolid", (fModuleH/2.) *mm, (fModuleW/2.) *mm, (reflectorT/2.) *mm );
    //   G4LogicalVolume* ReflectorlayerLogical = new G4LogicalVolume(ReflectorlayerSolid,FindMaterial("G4_Galactic"),"ReflectorlayerLogical");
    //   new G4PVPlacement(0,dimCalc->GetOrigin_Reflector(i),ReflectorlayerLogical,"ReflectorlayerPhysical",worldLogical,false,0,checkOverlaps);

    //   G4VSolid* mirrorSolid = new G4Box("mirrorSolid", 1.2/2. *mm, 1.2/2. *mm, reflectorT/2. *mm );
    //   ReflectorMirrorLogical_[i] = new G4LogicalVolume(mirrorSolid,FindMaterial("Aluminum"),"ReflectorMirrorLogical_");

    //   // DRsimMirrorParameterisation* mirrorParam = new DRsimMirrorParameterisation(fTowerXY.first,fTowerXY.second);
    //   DRsimMirrorParameterisation* mirrorParam = new DRsimMirrorParameterisation(fFiberX, fFiberY, fFiberWhich);
    //   G4PVParameterised* mirrorPhysical = new G4PVParameterised("mirrorPhysical",ReflectorMirrorLogical_[i],ReflectorlayerLogical,kXAxis,fTowerXY.first*fTowerXY.second/2,mirrorParam);
    //   // new G4LogicalBorderSurface("MirrorSurf",mirrorPhysical,modulePhysical,FindSurface("MirrorSurf"));
    //   new G4LogicalSkinSurface("MirrorSurf",ReflectorMirrorLogical_[i],FindSurface("MirrorSurf"));

    //   ReflectorMirrorLogical_[i]->SetVisAttributes(fVisAttrGray);
    // }
  // }
}

void DRsimDetectorConstruction::DefineCommands() {}

void DRsimDetectorConstruction::FiberImplementSingle(G4LogicalVolume* ModuleLogical__[], 
                                                    G4VSolid* fModuleSolid__[],
                                                    std::vector<float> fModuleWidth,
                                                    std::vector<float> fModuleHeight,
                                                    std::vector<float> fModuleDepth,
                                                    std::vector<G4LogicalVolume*> fiberUnitIntersection__[], 
                                                    std::vector<G4LogicalVolume*> fiberCladIntersection__[], 
                                                    std::vector<G4LogicalVolume*> fiberCoreIntersection__[]) {

  G4Tubs* tFiberClad   = new G4Tubs("tFiber",  0, clad_C_rMax, 1000. / 2. *mm, 0 *deg, 360. *deg);   // S is the same
  G4Tubs* tFiberCoreC  = new G4Tubs("tFiberC", 0, core_C_rMax, 1000. / 2. *mm, 0 *deg, 360. *deg);
  G4Tubs* tFiberCoreS  = new G4Tubs("tFiberS", 0, core_S_rMax, 1000. / 2. *mm, 0 *deg, 360. *deg);

  for (int i = 0; i < 12; i++) {
    // if (i < 9) continue;

    std::cout << i << " " << fModuleHeight.at(i) << " " << fModuleWidth.at(i) << " " << fModuleDepth.at(i) << std::endl;
    fFiberX.clear();
    fFiberY.clear();
    fFiberWhich.clear();

    int NofFiber = (int)(fModuleWidth.at(i) / 1.5);   
    int NofPlate = (int)(fModuleHeight.at(i) / 1.5);
    fLeftEdge = fmod(fModuleWidth.at(i), 1.5) / 2.;
    fBottomEdge = fmod(fModuleHeight.at(i), 1.5) / 2.;

    std::cout << "NofFiber : " << NofFiber << " | NofPlate : " << NofPlate << std::endl;

    double randDeviation = 0.; //  double randDeviation = fFiberUnitH - 1.;
    fTowerXY = std::make_pair(NofPlate,NofFiber);
    
    G4bool fWhich = false;  
    for (int k = 0; k < NofPlate; k++) {
      for (int j = 0; j < NofFiber; j++) { 
        
        //  ? fX : # of plate , fY : # of fiber in the plate

        G4float fX = fModuleHeight.at(i) *mm/2 - k*1.5 *mm - 0.75 *mm + fBottomEdge *mm;
        G4float fY = fModuleWidth.at(i) *mm/2 - j*1.5 *mm - 0.75 *mm + fLeftEdge *mm;

        G4ThreeVector tFiberCenterVector = G4ThreeVector{fX, fY, 249.9};
        fWhich = !fWhich;

        if (fModuleSolid__[i]->Inside(tFiberCenterVector) != kOutside) {
          fFiberX.push_back(fX);
          fFiberY.push_back(fY);
          fFiberWhich.push_back(fWhich);
        }
      }
      if ( NofFiber%2==0 ) { fWhich = !fWhich; }   
    }

    std::cout << i << " " << fFiberX.size() << " " << fFiberY.size() << " " << fFiberWhich.size() << std::endl;
    std::cout << i << " " << fFiberX.at(0) << " " << fFiberY.at(0) << " " << fFiberWhich.at(0) << std::endl;

    if ( doFiber ) {
      for (unsigned int j = 0; j < fFiberWhich.size(); j++) {

        if ( !fFiberWhich.at(j) ) { //c fibre

          tfiberCladIntersection = new G4IntersectionSolid("fiberClad", tFiberClad, fModuleSolid__[i], 0, G4ThreeVector(-fFiberX.at(j), -fFiberY.at(j), 0.));
          fiberCladIntersection__[i].push_back(new G4LogicalVolume(tfiberCladIntersection, FindMaterial("FluorinatedPolymer"), name));
          new G4PVPlacement(0, G4ThreeVector(fFiberX.at(j), fFiberY.at(j), 0), fiberCladIntersection__[i].at(j), name, ModuleLogical__[i], false, j, checkOverlaps);

          tfiberCoreIntersection = new G4IntersectionSolid("fiberCore", tFiberCoreC, fModuleSolid__[i], 0, G4ThreeVector(-fFiberX.at(j), -fFiberY.at(j), 0.));
          fiberCoreIntersection__[i].push_back(new G4LogicalVolume(tfiberCoreIntersection, FindMaterial("PMMA"), name));
          new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fiberCoreIntersection__[i].at(j), name, fiberCladIntersection__[i].at(j), false, j, checkOverlaps);

          fiberCladIntersection__[i].at(j)->SetVisAttributes(fVisAttrGray);
          fiberCoreIntersection__[i].at(j)->SetVisAttributes(fVisAttrBlue);
        } else { // s fibre

          tfiberCladIntersection = new G4IntersectionSolid("fiberClad", tFiberClad, fModuleSolid__[i], 0, G4ThreeVector(-fFiberX.at(j), -fFiberY.at(j), 0.));
          fiberCladIntersection__[i].push_back(new G4LogicalVolume(tfiberCladIntersection, FindMaterial("PMMA"), name));
          new G4PVPlacement(0, G4ThreeVector(fFiberX.at(j), fFiberY.at(j), 0), fiberCladIntersection__[i].at(j), name, ModuleLogical__[i], false, j, checkOverlaps);

          tfiberCoreIntersection = new G4IntersectionSolid("fiberCore", tFiberCoreS, fModuleSolid__[i], 0, G4ThreeVector(-fFiberX.at(j), -fFiberY.at(j), 0.));
          fiberCoreIntersection__[i].push_back(new G4LogicalVolume(tfiberCoreIntersection, FindMaterial("Polystyrene"), name));
          new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), fiberCoreIntersection__[i].at(j), name, fiberCladIntersection__[i].at(j), false, j, checkOverlaps);

          fiberCladIntersection__[i].at(j)->SetVisAttributes(fVisAttrGray);
          fiberCoreIntersection__[i].at(j)->SetVisAttributes(fVisAttrOrange);
        }
      }
    }
    std::cout << i << " " << ModuleLogical__[i]->GetNoDaughters() << std::endl;
  }
}

void DRsimDetectorConstruction::FiberImplement(G4int i, 
                                              G4LogicalVolume* ModuleLogical__[], 
                                              std::vector<G4LogicalVolume*> fiberUnitIntersection__[], 
                                              std::vector<G4LogicalVolume*> fiberCladIntersection__[], 
                                              std::vector<G4LogicalVolume*> fiberCoreIntersection__[]) {

  fFiberX.clear();
  fFiberY.clear();
  fFiberWhich.clear();

  int NofFiber = (int)(fModuleW / 1.5);   
  int NofPlate = (int)(fModuleH / 1.5);
  fBottomEdge = fmod(fModuleW, 1.5) / 2.;
  fLeftEdge = fmod(fModuleH, 1.5) / 2.;

  std::cout << "NofFiber : " << NofFiber << " | NofPlate : " << NofPlate << std::endl;

  double randDeviation = 0.; //  double randDeviation = fFiberUnitH - 1.;
  fTowerXY = std::make_pair(NofPlate,NofFiber);
  
  G4bool fWhich = false;  
  for (int k = 0; k < NofPlate; k++) {
    for (int j = 0; j < NofFiber; j++) { 
      /*
        ? fX : # of plate , fY : # of fiber in the plate
      */
      G4float fX = -fModuleH*mm/2 + k*1.5*mm + 0.75*mm + fBottomEdge*mm;
      G4float fY = -fModuleW*mm/2 + j*1.5*mm + 0.75*mm + fLeftEdge*mm;
      fWhich = !fWhich;
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
      fFiberWhich.push_back(fWhich);
    }
    if ( NofFiber%2==0 ) { fWhich = !fWhich; }   
  }
  
  if ( doFiber ) {
    for (unsigned int j = 0; j<fFiberX.size(); j++) {

      if ( !fFiberWhich.at(j) ) { //c fibre

        tfiberCladIntersection = new G4IntersectionSolid("fiberClad",fiberClad,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCladIntersection__[i].push_back(new G4LogicalVolume(tfiberCladIntersection,FindMaterial("FluorinatedPolymer"),name));
        new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberCladIntersection__[i].at(j),name,ModuleLogical__[i],false,j,checkOverlaps);

        tfiberCoreIntersection = new G4IntersectionSolid("fiberCore",fiberCoreC,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCoreIntersection__[i].push_back(new G4LogicalVolume(tfiberCoreIntersection,FindMaterial("PMMA"),name));
        new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberCoreIntersection__[i].at(j),name,fiberCladIntersection__[i].at(j),false,j,checkOverlaps);

        fiberCladIntersection__[i].at(j)->SetVisAttributes(fVisAttrGray);
        fiberCoreIntersection__[i].at(j)->SetVisAttributes(fVisAttrBlue);
      } else { // s fibre

        tfiberCladIntersection = new G4IntersectionSolid("fiberClad",fiberClad,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCladIntersection__[i].push_back(new G4LogicalVolume(tfiberCladIntersection,FindMaterial("PMMA"),name));
        new G4PVPlacement(0,G4ThreeVector(fFiberX.at(j),fFiberY.at(j),0),fiberCladIntersection__[i].at(j),name,ModuleLogical__[i],false,j,checkOverlaps);

        tfiberCoreIntersection = new G4IntersectionSolid("fiberCore",fiberCoreS,module,0,G4ThreeVector(-fFiberX.at(j),-fFiberY.at(j),0.));
        fiberCoreIntersection__[i].push_back(new G4LogicalVolume(tfiberCoreIntersection,FindMaterial("Polystyrene"),name));
        new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fiberCoreIntersection__[i].at(j),name,fiberCladIntersection__[i].at(j),false,j,checkOverlaps);

        fiberCladIntersection__[i].at(j)->SetVisAttributes(fVisAttrGray);
        fiberCoreIntersection__[i].at(j)->SetVisAttributes(fVisAttrOrange);
      }
    }
  }
}


