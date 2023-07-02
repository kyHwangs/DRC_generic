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
  fModuleHeight = 0;
  fModuleWidth  = 0;
  fPMTT         = 0;
  fReflectorT   = 0;
  fisModule     = false;

}

dimensionCalc::~dimensionCalc() {}

G4ThreeVector dimensionCalc::GetOrigin(G4int i) {

  G4double tThetaOne = TMath::ATan(7.5 / 500.);
  G4double tThetaTwo = TMath::ATan(15. / 500.);

  if (i == 0 ) return G4ThreeVector(0, 0, 250.); // 3D Center
  
  if (i == 1 ) return G4ThreeVector(-33., 0, 250.);

  if (i == 2 ) return G4ThreeVector(7.5, -(15. + 250. * TMath::Tan(tThetaOne) + 18 / TMath::Cos(tThetaOne)), 250.);
  if (i == 3 ) return G4ThreeVector(7.5, 15. + 250. * TMath::Tan(tThetaOne) + 18 / TMath::Cos(tThetaOne), 250.);
  if (i == 4 ) return G4ThreeVector(15 + 250. * TMath::Tan(tThetaTwo) + 22.5 / TMath::Cos(tThetaTwo), 0, 250.);

  if (i == 5 ) return G4ThreeVector(-27., 34.5, 250.);
  if (i == 6 ) return G4ThreeVector(-27., -34.5, 250.);
  if (i == 7 ) return G4ThreeVector(42., -(15. + 250 * TMath::Tan(tThetaOne) + 12. / TMath::Cos(tThetaOne)), 250.);
  if (i == 8 ) return G4ThreeVector(42., 15. + 250 * TMath::Tan(tThetaOne) + 12. / TMath::Cos(tThetaOne), 250.);

  if (i == 9 ) return G4ThreeVector(0, 0, 250.); // 3D Wedge

}

G4RotationMatrix* dimensionCalc::GetRM(G4int i) {

  G4double tThetaOne = TMath::ATan(7.5 / 500.);
  G4double tThetaTwo = TMath::ATan(15. / 500.);

  if (i == 0 ) return new G4RotationMatrix(); // 3D Center
  
  if (i == 1 ) return  new G4RotationMatrix();

  if (i == 2 ) {
    G4RotationMatrix* RotMatrix = new G4RotationMatrix();

    RotMatrix->rotateX(-tThetaOne);

    return RotMatrix;
  }

  if (i == 3 ) {
    G4RotationMatrix* RotMatrix = new G4RotationMatrix();

    RotMatrix->rotateX(tThetaOne);

    return RotMatrix;
  }

  if (i == 4 ) {
    G4RotationMatrix* RotMatrix = new G4RotationMatrix();

    RotMatrix->rotateY(-tThetaTwo);

    return RotMatrix;
  }

  if (i == 5 ) return  new G4RotationMatrix();
  if (i == 6 ) return  new G4RotationMatrix();


  if (i == 7 ) {
    G4RotationMatrix* RotMatrix = new G4RotationMatrix();

    RotMatrix->rotateX(-tThetaOne);

    return RotMatrix;
  }

  if (i == 8 ) {
    G4RotationMatrix* RotMatrix = new G4RotationMatrix();

    RotMatrix->rotateX(tThetaOne);

    return RotMatrix;
  }

  if (i == 9 ) return new G4RotationMatrix(); // 3D Wedge

}

G4ThreeVector dimensionCalc::GetOrigin_PMTG(G4int i) {

  double row = i/fNofRow;
  double col = i%fNofRow;

  return G4ThreeVector( 
                        -fModuleHeight * (double)fNofRow/2. + row * fModuleHeight + fModuleHeight/2., 
                        -fModuleWidth  * (double)fNofRow/2. + col * fModuleWidth  + fModuleWidth/2., 
                        ftower_height + fFrontL + fPMTT/2.
                      );
}

G4ThreeVector dimensionCalc::GetOrigin_Reflector(G4int i) {

  double row = i/fNofRow;
  double col = i%fNofRow;

  return G4ThreeVector( 
                        -fModuleHeight * (double)fNofRow/2. + row * fModuleHeight + fModuleHeight/2., 
                        -fModuleWidth  * (double)fNofRow/2. + col * fModuleWidth  + fModuleWidth/2., 
                        fFrontL - fReflectorT/2.
                      );
}
