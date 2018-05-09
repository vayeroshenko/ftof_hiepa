#ifndef fTOF_TrackingAction_h
#define fTOF_TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "globals.hh"

class fTOF_TrackingAction : public G4UserTrackingAction
{
public:
  fTOF_TrackingAction();
  ~fTOF_TrackingAction();

  void PreUserTrackingAction(const G4Track*);
  void PostUserTrackingAction(const G4Track*);

private:

};
#endif
