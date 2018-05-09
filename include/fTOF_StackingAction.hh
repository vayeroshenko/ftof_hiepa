#ifndef fTOF_StackingAction_h
#define fTOF_StackingAction_h 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class fTOF_StackingAction : public G4UserStackingAction
{
public:
  fTOF_StackingAction();
  ~fTOF_StackingAction();
  
public:
  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  void NewStage();
  void PrepareNewEvent();
  inline G4int GetTotPhotNum(){return gammaCounter;};

private:
  G4int gammaCounter;
};
#endif
