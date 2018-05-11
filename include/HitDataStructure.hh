#ifndef HitDataStructure_h
#define HitDataStructure_h 1

#include "globals.hh"

struct HitData
{
  G4int TrackID;
  G4int ParentID;
  G4double Energy;
  G4double Wavelength;
  G4double Time;
  G4double photPathLen;
  G4int SecID;
  G4int chID;
  G4double PosX;
  G4double PosY;
  G4double PosZ;
  G4double trkMomX;
  G4double trkMomY;
  G4double trkMomZ;
  G4double trkPosX;
  G4double trkPosY;
  G4double trkPosZ;
  G4double trkT;
  G4double trkLength;

  G4double entTime;

  G4double entMomX;
  G4double entMomY;
  G4double entMomZ;
  G4double entPosX;
  G4double entPosY;
  G4double entPosZ;



  G4int trkNSideRefl;
  G4int trkSideID;

};

#endif
