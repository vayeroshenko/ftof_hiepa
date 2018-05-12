//G4
#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Point3D.hh"
#include "G4TwoVector.hh"
#include "globals.hh"

#include "TMath.h"

//My
#include "fTOFConst.hh"
#include "crtConst.hh"


struct VolumeStruct {
  G4Material*        material;
  G4VSolid*          solid;
  G4LogicalVolume*   logical;
  G4VPhysicalVolume* physical;
  VolumeStruct() :
    material(0),
    solid(0),
    logical(0),
    physical(0)
  {;}
  ~VolumeStruct() {;}
};

struct WorldStruct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  WorldStruct() :
    sizeX(300.0*cm),
    sizeY(300.0*cm),
    sizeZ(600.0*cm)      
  {;}
};



///////////////////// fTOF bar 10.04.18 ///////////////

struct FullBarStruct : VolumeStruct {
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  FullBarStruct() :
    sizeX(fTOFConst::fullBarSizeX),
    sizeY(fTOFConst::fullBarSizeY),
    sizeZ(fTOFConst::fullBarSizeZ)
  {;}
};



struct HamPmtWindowStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  HamPmtWindowStruct() :
    sizeX(fTOFConst::hamWindowThick),
    sizeY(fTOFConst::hamSensitiveSize),
    sizeZ(fTOFConst::hamSensitiveSize)
  {;}
};

struct PlanPmtWindowStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PlanPmtWindowStruct() :
    sizeX(fTOFConst::planWindowThick),
    sizeY(fTOFConst::planSensitiveSize),
    sizeZ(fTOFConst::planSensitiveSize)
  {;}
};

struct PlanPmtChannelStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PlanPmtChannelStruct() :
    sizeX(fTOFConst::planChanThick),
    sizeY(fTOFConst::planChanSize),
    sizeZ(fTOFConst::planChanSize)
  {;}
};

struct HamPmtChannelStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  HamPmtChannelStruct() :
    sizeX(fTOFConst::hamChanThick),
    sizeY(fTOFConst::hamChanSize),
    sizeZ(fTOFConst::hamChanSize)
  {;}
};

struct HamPmtBoxStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  HamPmtBoxStruct() :
    sizeX(fTOFConst::hamBoxDepth),
    sizeY(fTOFConst::hamBoxSize),
    sizeZ(fTOFConst::hamBoxSize)
  {;}
};

struct PlanPmtBoxStruct : VolumeStruct { //////////////////// 11.04.18
  // Defined as a G4Box
  const G4double sizeX;
  const G4double sizeY;
  const G4double sizeZ;
  PlanPmtBoxStruct() :
    sizeX(fTOFConst::planBoxDepth),
    sizeY(fTOFConst::planBoxSize),
    sizeZ(fTOFConst::planBoxSize)
  {;}
};


// struct TrapezeSectorStruct: VolumeStruct {
//   const G4double shortSide;
//   const G4double longSide;
//   const G4double thickness;
//   const G4double height;
//   const G4double angle;
//   const G4double sides;
//   const G4double middleLine;
//   TrapezeSectorStruct():
//     shortSide(fTOFConst::sectorShortSide),
//     longSide(fTOFConst::sectorLongSide),
//     thickness(fTOFConst::sectorThickness),
//     height(fTOFConst::sectorHeight),
//     angle(atan((longSide-shortSide)/2./height)),
//     sides(sqrt(height*height + (longSide-shortSide)*(longSide-shortSide)/4.)),
//     middleLine((longSide+shortSide)/2.)
//   {;}
// };



struct TrapezeSectorStruct: VolumeStruct {
  const G4double shortSide;
  const G4double longSide;
  const G4double thickness;
  const G4double height;
  const G4double angle;
  const G4double sides;
  const G4double middleLine;
  TrapezeSectorStruct():
    shortSide(fTOFConst::innerSide),
    longSide(fTOFConst::outerSide),
    thickness(fTOFConst::sectorThickness),
    height(fTOFConst::outerRad - fTOFConst::innerRad),
    angle(atan((longSide-shortSide)/2./height)),
    sides(sqrt(height*height + (longSide-shortSide)*(longSide-shortSide)/4.)),
    middleLine((longSide+shortSide)/2.)
  {;}
};



struct TrapezeAbsStruct: VolumeStruct {
  const G4double shortSide;
  const G4double longSide;
  const G4double thickness;
  const G4double height;
  const G4double angle;
  const G4double sides;
  const G4double middleLine;
  TrapezeAbsStruct():
    shortSide(fTOFConst::absInnerSide),
    longSide(fTOFConst::absOuterSide),
    thickness(fTOFConst::sectorThickness),
    height(fTOFConst::outerRad - fTOFConst::innerRad),
    angle(atan((longSide-shortSide)/2./height)),
    sides(sqrt(height*height + (longSide-shortSide)*(longSide-shortSide)/4.)),
    middleLine((longSide+shortSide)/2.)
  {;}
};




// struct TrapezeMixerStruct: VolumeStruct {
//   const G4double shortSide;
//   const G4double longSide;
//   const G4double thickness;
//   const G4double height;
//   const G4double angle;
//   const G4double sides;
//   const G4double middleLine;
//   TrapezeMixerStruct():
//     shortSide(fTOFConst::sectorLongSide),
//     longSide(20.*mm),
//     thickness(fTOFConst::sectorThickness),
//     height(fTOFConst::mixerLength),
//     angle(atan(abs(longSide-shortSide)/2./height)),
//     sides(sqrt(height*height + (longSide-shortSide)*(longSide-shortSide)/4.)),
//     middleLine((longSide+shortSide)/2.)
//   {;}
// };

//////////////////////////////////////////////////////////////////
