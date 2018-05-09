//My
#include "fTOF_StackingAction.hh"

//G4
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"


fTOF_StackingAction::fTOF_StackingAction() : gammaCounter(0)
{}

fTOF_StackingAction::~fTOF_StackingAction()
{}

G4ClassificationOfNewTrack
fTOF_StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{
  if (aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()
      && aTrack->GetParentID() > 0) {
    gammaCounter++;
    // Get direction
    //    const G4ThreeVector& momDir(aTrack->GetMomentumDirection());
  }

  return fUrgent;
}

void fTOF_StackingAction::NewStage()
{
  //G4cout << "Number of optical photons in this event : " << gammaCounter
  //<< G4endl;
}

void fTOF_StackingAction::PrepareNewEvent()
{
  gammaCounter = 0;
}
