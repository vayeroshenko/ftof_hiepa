#ifndef MagneticField_H
#define MagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class MagneticField : public G4MagneticField
{
public:
  MagneticField();
  ~MagneticField();
  
  virtual void GetFieldValue( const  double Point[4],
			      double *Bfield ) const;
  
  bool SetField( const double *Point, G4double & Bx, G4double & By, G4double & Bz) const;
  void SetBz(G4double Bzval){Bfield = Bzval;};

private:
  
  G4double Bfield;
};

#endif

