#include "RecoInterface.h"
#include "DRsimMirrorParameterisation.hh"
#include "DRsimCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DRsimMirrorParameterisation::DRsimMirrorParameterisation(const G4int numx, const G4int numy)
: G4VPVParameterisation()
{
  for ( G4int copyNo = 0; copyNo < numx*numy; copyNo++ ) {

    G4int column = copyNo / numy;
    G4int row = copyNo % numy;

    double colIntv = 2* 23./21;
    double rowIntv = 23./14;

    if ( RecoInterface::IsCerenkov(column,row) ) {
      fXMirror.push_back( -rowIntv*6.5 + column*rowIntv );
      fYMirror.push_back( -colIntv*5. + row*colIntv );
    }
  }
  fNumx = numx;
  fNumy = numy;
}

DRsimMirrorParameterisation::~DRsimMirrorParameterisation() {}

void DRsimMirrorParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXMirror[copyNo],fYMirror[copyNo],0.));
}
