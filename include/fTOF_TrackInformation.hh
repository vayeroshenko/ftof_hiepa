#ifndef fTOF_TrackInformation_h
#define fTOF_TrackInformation_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ThreeVector.hh"

class fTOF_TrackInformation : public G4VUserTrackInformation
{
public:
  fTOF_TrackInformation();
  fTOF_TrackInformation(const G4Track* aTrack);
  fTOF_TrackInformation(const fTOF_TrackInformation* aTrackInfo);
  virtual ~fTOF_TrackInformation();

  inline void *operator new(size_t);
  inline void operator delete(void* aTrackInfo);
  inline int operator ==(const fTOF_TrackInformation& right) const
  {return (this==&right);}

  void Print() const;

private:
  G4int                 originalTrackID;
  G4ParticleDefinition* particleDefinition;
  G4ThreeVector         originalPosition;
  G4ThreeVector         originalMomentum;
  G4double              originalEnergy;
  G4double              originalTime;
  G4ThreeVector         myOriginalPosition;
  G4ThreeVector         myOriginalMomentum;
  G4double              myOriginalEnergy;
  G4double              myOriginalTime;

public:
  inline G4int                 GetOriginalTrackID() const
  {return originalTrackID;}
  inline G4ParticleDefinition* GetOriginalParticle() const
  {return particleDefinition;}
  inline G4ThreeVector         GetOriginalPosition() const 
  {return originalPosition;}
  inline G4ThreeVector         GetOriginalMomentum() const 
  {return originalMomentum;}
  inline G4double              GetOriginalEnergy() const 
  {return originalEnergy;}
  inline G4double              GetOriginalTime() const 
  {return originalTime;}
  inline G4ThreeVector         GetMyOriginalPosition() const 
  {return myOriginalPosition;}
  inline G4ThreeVector         GetMyOriginalMomentum() const 
  {return myOriginalMomentum;}
  inline G4double              GetMyOriginalEnergy() const 
  {return myOriginalEnergy;}
  inline G4double              GetMyOriginalTime() const 
  {return myOriginalTime;}

  inline void SetMyOriginalPosition(const G4ThreeVector& xyz)
  {myOriginalPosition = xyz;}
  inline void SetMyOriginalMomentum(const G4ThreeVector& xyz) 
  {myOriginalMomentum = xyz;}
  inline void SetMyOriginalEnergy(G4double e) 
  {myOriginalEnergy = e;}
  inline void SetMyOriginalTime(G4double t) 
  {myOriginalTime = t;}

};

extern G4Allocator<fTOF_TrackInformation> aTrackInformationAllocator;

inline void* fTOF_TrackInformation::operator new(size_t)
{
  void* aTrackInfo;
  aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
  return aTrackInfo;
}

inline void fTOF_TrackInformation::operator delete(void* aTrackInfo)
{ aTrackInformationAllocator.FreeSingle((fTOF_TrackInformation*)aTrackInfo);}

#endif
