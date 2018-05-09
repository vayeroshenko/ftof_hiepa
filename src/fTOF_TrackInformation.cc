#include "fTOF_TrackInformation.hh"
#include "G4ios.hh"

G4Allocator<fTOF_TrackInformation> aTrackInformationAllocator;

fTOF_TrackInformation::fTOF_TrackInformation()
{
  originalTrackID = 0;
  particleDefinition = 0;
  originalPosition = G4ThreeVector(0.,0.,0.);
  originalMomentum = G4ThreeVector(0.,0.,0.);
  originalEnergy = 0.;
  originalTime = 0.;
  myOriginalPosition = G4ThreeVector(0.,0.,0.);
  myOriginalMomentum = G4ThreeVector(0.,0.,0.);
  myOriginalEnergy = 0.;
  myOriginalTime = 0.;
}

fTOF_TrackInformation::fTOF_TrackInformation(const G4Track* aTrack)
{
  originalTrackID = aTrack->GetTrackID();
  particleDefinition = aTrack->GetDefinition();
  originalPosition = aTrack->GetPosition();
  originalMomentum = aTrack->GetMomentum();
  originalEnergy = aTrack->GetTotalEnergy();
  originalTime = aTrack->GetGlobalTime();
  myOriginalPosition = aTrack->GetPosition();
  myOriginalMomentum = aTrack->GetMomentum();
  myOriginalEnergy = aTrack->GetTotalEnergy();
  myOriginalTime = aTrack->GetGlobalTime();
}

fTOF_TrackInformation::fTOF_TrackInformation(const fTOF_TrackInformation*
					       aTrackInfo)
{
  originalTrackID = aTrackInfo->originalTrackID;
  particleDefinition = aTrackInfo->particleDefinition;
  originalPosition = aTrackInfo->originalPosition;
  originalMomentum = aTrackInfo->originalMomentum;
  originalEnergy = aTrackInfo->originalEnergy;
  originalTime = aTrackInfo->originalTime;
  myOriginalPosition = aTrackInfo->myOriginalPosition;
  myOriginalMomentum = aTrackInfo->myOriginalMomentum;
  myOriginalEnergy = aTrackInfo->myOriginalEnergy;
  myOriginalTime = aTrackInfo->myOriginalTime;
}

fTOF_TrackInformation::~fTOF_TrackInformation() {;}

void fTOF_TrackInformation::Print() const
{
  G4cout << "Original Track ID " << originalTrackID << " at "
	 << originalPosition << G4endl;
}
