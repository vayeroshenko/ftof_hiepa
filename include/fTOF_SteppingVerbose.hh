#ifndef fTOF_SteppingVerbose_h
#define fTOF_SteppingVerbose_h

#include "G4SteppingVerbose.hh"

class fTOF_SteppingVerbose : public G4SteppingVerbose
{
public:
  fTOF_SteppingVerbose();
  ~fTOF_SteppingVerbose();

  void StepInfo();
  void TrackingStarted();

};

#endif
