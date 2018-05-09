//my
#include "fTOF_Hit.hh"
//G4
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "G4UnitsTable.hh"

G4Allocator<fTOF_Hit> fTOF_HitAllocator;

fTOF_Hit::fTOF_Hit() :
  G4VHit()
{}

fTOF_Hit::~fTOF_Hit() {}

fTOF_Hit::fTOF_Hit(const fTOF_Hit& right) : 
  G4VHit()
{
  myData = right.myData;
}

const fTOF_Hit& fTOF_Hit::operator=(const fTOF_Hit& right)
{
  myData = right.myData;
  return *this;
}

G4int fTOF_Hit::operator==(const fTOF_Hit& right) const
{
  return (this == &right) ? 1 : 0;
}

void fTOF_Hit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    G4ThreeVector pos = 
      G4ThreeVector(myData.PosX, myData.PosY, myData.PosZ);
    G4Circle circle(pos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Color color(1., 0., 0.);
    G4VisAttributes att(color);
    circle.SetVisAttributes(att);
    pVVisManager->Draw(circle);
  }
}

void fTOF_Hit::Print()
{
  /*
  G4cout << "Hit Print: "
	 << " trackID: " << trackID
	 << " energy: " << G4BestUnit(energy, "Energy")
	 << " position: " << G4BestUnit(pos, "Length")
	 << " time: " << G4BestUnit(time, "Time")
	 << G4endl;
  */  
}
