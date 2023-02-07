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

  if ( i == 0 ) return G4ThreeVector(22.5  , 22.5   , ftower_height/2 + fFrontL);
  if ( i == 1 ) return G4ThreeVector(22.5  , 68.25  , ftower_height/2 + fFrontL);
  if ( i == 2 ) return G4ThreeVector(-23.25, 22.5   , ftower_height/2 + fFrontL);
  if ( i == 3 ) return G4ThreeVector(-23.25, 68.25  , ftower_height/2 + fFrontL);

  if ( i == 4 ) return G4ThreeVector(30.   , -75.75 , ftower_height/2 + fFrontL);
  if ( i == 5 ) return G4ThreeVector(30.   , -45.   , ftower_height/2 + fFrontL);
  if ( i == 6 ) return G4ThreeVector(30.   , -15.   , ftower_height/2 + fFrontL);
  if ( i == 7 ) return G4ThreeVector(-0.75 , -75.75 , ftower_height/2 + fFrontL);
  if ( i == 8 ) return G4ThreeVector(-0.75 , -45.   , ftower_height/2 + fFrontL);
  if ( i == 9 ) return G4ThreeVector(-0.75 , -15.   , ftower_height/2 + fFrontL);
  if ( i == 10) return G4ThreeVector(-30.75, -75.75 , ftower_height/2 + fFrontL);
  if ( i == 11) return G4ThreeVector(-30.75, -45.   , ftower_height/2 + fFrontL);
  if ( i == 12) return G4ThreeVector(-30.75, -15.   , ftower_height/2 + fFrontL);

  // double row = i/fNofRow;
  // double col = i%fNofRow;

  // return G4ThreeVector( -90 * (double)fNofRow/2. + row * 90 + 45, -90 * (double)fNofRow/2. + col * 90 + 45, ftower_height/2 + fFrontL);
}

G4ThreeVector dimensionCalc::GetOrigin_PMTG(G4int i) {

  if ( i == 0 ) return G4ThreeVector(22.5  , 22.5   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 1 ) return G4ThreeVector(22.5  , 68.25  , ftower_height + fFrontL + fPMTT/2);
  if ( i == 2 ) return G4ThreeVector(-23.25, 22.5   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 3 ) return G4ThreeVector(-23.25, 68.25  , ftower_height + fFrontL + fPMTT/2);

  if ( i == 4 ) return G4ThreeVector(30.   , -75.75 , ftower_height + fFrontL + fPMTT/2);
  if ( i == 5 ) return G4ThreeVector(30.   , -45.   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 6 ) return G4ThreeVector(30.   , -15.   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 7 ) return G4ThreeVector(-0.75 , -75.75 , ftower_height + fFrontL + fPMTT/2);
  if ( i == 8 ) return G4ThreeVector(-0.75 , -45.   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 9 ) return G4ThreeVector(-0.75 , -15.   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 10) return G4ThreeVector(-30.75, -75.75 , ftower_height + fFrontL + fPMTT/2);
  if ( i == 11) return G4ThreeVector(-30.75, -45.   , ftower_height + fFrontL + fPMTT/2);
  if ( i == 12) return G4ThreeVector(-30.75, -15.   , ftower_height + fFrontL + fPMTT/2);

  // return G4ThreeVector( -90 * (double)fNofRow/2. + row * 90 + 45, -90 * (double)fNofRow/2. + col * 90 + 45, ftower_height + fFrontL + fPMTT/2);
}

G4ThreeVector dimensionCalc::GetOrigin_Reflector(G4int i) {

  double row = i/fNofRow;
  double col = i%fNofRow;

  return G4ThreeVector( -90 * (double)fNofRow/2. + row * 90 + 45, -90 * (double)fNofRow/2. + col * 90 + 45, fFrontL - fReflectorT/2);
}
