#include "DRsimCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DRsimCellParameterisation::DRsimCellParameterisation(const G4int numx, const G4int numy)
: G4VPVParameterisation()
{
  for (G4int copyNo = 0; copyNo < numx * numy; copyNo++ ) {
    G4int column = copyNo / numy;
    G4int row = copyNo % numy;

    double colIntv = 2* 23./21;
    double rowIntv = 23./14;
    
    fXCell.push_back( -rowIntv*6.5 + column*rowIntv );
    fYCell.push_back( -colIntv*5. + row*colIntv );
  }
  fNumx = numx;
  fNumy = numy;
}

DRsimCellParameterisation::~DRsimCellParameterisation()
{}

void DRsimCellParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXCell[copyNo],fYCell[copyNo],0.));
}
