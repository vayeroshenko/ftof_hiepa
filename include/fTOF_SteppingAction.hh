#ifndef fTOF_SteppingAction_h
#define fTOF_SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "fTOF_PrimaryGeneratorAction.hh"
#include "globals.hh"

class fTOF_SteppingAction : public G4UserSteppingAction
{
public:
  fTOF_SteppingAction(fTOF_PrimaryGeneratorAction* genAction);
  ~fTOF_SteppingAction();

  void UserSteppingAction(const G4Step*);

  void Reset();
  void ResetPerEvent();

private:
  void InternalReflectionProbability(G4double, G4double&);
  fTOF_PrimaryGeneratorAction* _genAction;

private:

  //Trk information  
  G4int _particleID;
  G4int _SecID;
  G4double _probOfReflection;
  G4double _trkMomX;
  G4double _trkMomY;
  G4double _trkMomZ;
  G4double _trkPosX;
  G4double _trkPosY;
  G4double _trkPosZ;
  G4double _trkT;
  G4double _trkLength;

  G4double _entMomX;
  G4double _entMomY;
  G4double _entMomZ;
  G4double _entPosX;
  G4double _entPosY;
  G4double _entPosZ;

  G4double _entTime;

  G4double _trkNSideRefl;
  G4double _trkSideID;




  //photon information
  G4int _chID;
  //G4int nKillPhot;
  //G4double  _chX;
  //G4double  _chY;
  //G4double  _chZ;

  // G4int _particleID;
  //G4double _totalPath;
  //G4double _barPath;
  //G4double _sobPath;
  //G4int _nBounceMirror1;
  //G4int _nBounceMirror2;
  //G4int _nBounceEndMirror;
  //G4int _nWedgeSide;
  //G4int _nWedgeTop;
  //G4int _nWedgeBottom;
  //G4int _nFBlockSide;
  //G4int _bar;
  //G4int _barFirst;
  //G4double _barX;
  //G4double _barY;
  //G4double _barZ;
  //G4double _barThetaX;
  //G4double _barThetaY;
  //G4double _timeLeftBar;
  //G4double _probOfReflection;
  //G4int _numberOfBounces;

};

#endif
