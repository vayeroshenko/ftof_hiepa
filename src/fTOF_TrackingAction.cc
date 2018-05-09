#include "fTOF_TrackingAction.hh"
#include "fTOF_TrackInformation.hh"
#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4TrackingManager.hh"

fTOF_TrackingAction::fTOF_TrackingAction()
{;}

fTOF_TrackingAction::~fTOF_TrackingAction()
{;}

void fTOF_TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  if (aTrack->GetParentID() == 0 && aTrack->GetUserInformation() == 0) {
    fTOF_TrackInformation* anInfo = new fTOF_TrackInformation(aTrack);
    G4Track* theTrack = (G4Track*)aTrack;
    theTrack->SetUserInformation(anInfo);
  }
}

void fTOF_TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if (secondaries) {
    fTOF_TrackInformation* info = 
      (fTOF_TrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if (nSeco > 0) {
      for (size_t i = 0; i < nSeco; i++) {
	fTOF_TrackInformation* infoNew = new fTOF_TrackInformation(info);
	infoNew->SetMyOriginalPosition((*secondaries)[i]->GetPosition());
	infoNew->SetMyOriginalMomentum((*secondaries)[i]->GetMomentum());
	infoNew->SetMyOriginalEnergy((*secondaries)[i]->GetTotalEnergy());
	infoNew->SetMyOriginalTime((*secondaries)[i]->GetGlobalTime());
	(*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
}
