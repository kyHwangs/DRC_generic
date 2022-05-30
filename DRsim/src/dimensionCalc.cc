#include "dimensionCalc.hh"

#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"

#include <cmath>
#include <stdio.h>
#include <float.h>

using namespace std;

dimensionCalc::dimensionCalc() {

  fFrontL       = 0;
  fNofModules   = 0;
  fNofRow       = 0;
  ftower_height = 0;
  fPMTT         = 0;
  fReflectorT   = 0;
  fisModule     = false;

}

dimensionCalc::~dimensionCalc() {}

G4ThreeVector dimensionCalc::GetOrigin(G4int i) {

  double row = i/fNofRow;
  double col = i%fNofRow;

  return G4ThreeVector( -90 * (double)fNofRow/2. + row * 90 + 45, -90 * (double)fNofRow/2. + col * 90 + 45, ftower_height/2 + fFrontL);
}

G4ThreeVector dimensionCalc::GetOrigin_PMTG(G4int i) {

  double row = i/fNofRow;
  double col = i%fNofRow;

  return G4ThreeVector( -90 * (double)fNofRow/2. + row * 90 + 45, -90 * (double)fNofRow/2. + col * 90 + 45, ftower_height + fFrontL + fPMTT/2);
}

G4ThreeVector dimensionCalc::GetOrigin_Reflector(G4int i) {

  double row = i/fNofRow;
  double col = i%fNofRow;

  return G4ThreeVector( -90 * (double)fNofRow/2. + row * 90 + 45, -90 * (double)fNofRow/2. + col * 90 + 45, fFrontL - fReflectorT/2);
}

G4RotationMatrix* dimensionCalc::GetRM(G4double i) {

  G4RotationMatrix* RotMatrix = new G4RotationMatrix();
  RotMatrix->rotateX(3.14159265358979323846 * (i/180));

  return RotMatrix;
}
