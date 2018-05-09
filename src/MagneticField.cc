#include "MagneticField.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "math.h"

using namespace std;

MagneticField::MagneticField()
{
    //messenger = new A01MagneticFieldMessenger(this);
    G4cout << "MagneticField::MagneticField()" << G4endl;
    Bfield = 1.5*tesla;
}

MagneticField::~MagneticField()
{
    //delete messenger; 
}

void MagneticField::GetFieldValue(const double Point[4],double *Bfield) const
{
    SetField(Point, Bfield[0],Bfield[1],Bfield[2]);
}

bool MagneticField::SetField( const double *Point,  G4double & Bx, G4double & By, G4double & Bz) const {
    Bx = 0.0;
    By = 0.0;
    //Bz = Bfield;    
    Bz = 1.5*tesla;

    //G4cout<<"Bx = "<<Bx<<G4endl
    //<<"By = "<<By<<G4endl
    //<<"Bz = "<<Bz<<G4endl;

    return true;
}



