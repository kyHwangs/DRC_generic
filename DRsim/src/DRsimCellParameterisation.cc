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
    
    fXCell.push_back( -23.*mm/2 + column*1.5*mm + 1.25*mm );
    fYCell.push_back( -23.*mm/2 + row*1.5*mm + 1.*mm );
  }
  fNumx = numx;
  fNumy = numy;
}

DRsimCellParameterisation::~DRsimCellParameterisation()
{}

void DRsimCellParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXCell[copyNo],fYCell[copyNo],0.));
}
