#ifndef fTOF_Hit_h
#define fTOF_Hit_h 1

//G4
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

//my
#include "HitDataStructure.hh"

class fTOF_Hit : public G4VHit
{

public:
  fTOF_Hit();
  ~fTOF_Hit();
  fTOF_Hit(const fTOF_Hit&);
  const fTOF_Hit& operator=(const fTOF_Hit&);
  G4int operator==(const fTOF_Hit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  void Draw();
  void Print();

public:

  HitData myData;

private:

};

typedef G4THitsCollection<fTOF_Hit> fTOF_HitsCollection;

extern G4Allocator<fTOF_Hit> fTOF_HitAllocator;

inline void* fTOF_Hit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) fTOF_HitAllocator.MallocSingle();
  return aHit;
}

inline void fTOF_Hit::operator delete(void *aHit)
{
  fTOF_HitAllocator.FreeSingle((fTOF_Hit*) aHit);
}

#endif
