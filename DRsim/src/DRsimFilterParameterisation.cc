#include "RecoInterface.h"
#include "DRsimFilterParameterisation.hh"
#include "DRsimCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DRsimFilterParameterisation::DRsimFilterParameterisation(const G4int numx, const G4int numy)
: G4VPVParameterisation()
{
  for ( G4int copyNo = 0; copyNo < numx*numy; copyNo++ ) {

    G4int column = copyNo / numy;
    G4int row = copyNo % numy;

    // if ( !RecoInterface::IsCerenkov(column,row) ) {
    if ( true ) {
      fXFilter.push_back( -25.*mm/2 + column*1.*mm + 0.5*mm );
      fYFilter.push_back( -25.*mm/2 + row*1.*mm + 0.5*mm );
    }
  }
  fNumx = numx;
  fNumy = numy;
}

DRsimFilterParameterisation::~DRsimFilterParameterisation() {}

void DRsimFilterParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXFilter[copyNo],fYFilter[copyNo],0.));
}
