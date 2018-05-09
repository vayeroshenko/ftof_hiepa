#ifndef fTOF_SensitiveDetector_h
#define fTOF_SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "fTOF_Hit.hh"
#include "HitDataStructure.hh"

class G4HCofThisEvent;
class G4Step;
class G4TouchableHistory;

class fTOF_SensitiveDetector : public G4VSensitiveDetector
{
public:
  fTOF_SensitiveDetector(G4String);
  ~fTOF_SensitiveDetector();

  void Initialize(G4HCofThisEvent*);
  
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  G4bool ProcessHits_fTOF(const G4Step*, G4TouchableHistory*,
			   HitData hitInfo);

  void EndOfEvent(G4HCofThisEvent*);
private:

  fTOF_HitsCollection* OpticalPhotonCollection;


};

#endif
