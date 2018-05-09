#ifndef fTOF_EventAction_h
#define fTOF_EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class fTOF_RunAction;
class fTOF_SteppingAction;
class fTOF_StackingAction;
class fTOF_PrimaryGeneratorAction;

class fTOF_EventAction : public G4UserEventAction
{
public:
  fTOF_EventAction(fTOF_RunAction*, fTOF_SteppingAction*);
  ~fTOF_EventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);

  void SetStackingAction(fTOF_StackingAction* sta){_stackingAction = sta;};
  void SetPrimGenerator(fTOF_PrimaryGeneratorAction *gen){_primGenerator = gen;};
private:
  fTOF_RunAction* runAction;
  fTOF_SteppingAction* _steppingAction;
  G4int printModulo;
  G4int thePhotonCollectionID;
  
  fTOF_StackingAction* _stackingAction; 
  fTOF_PrimaryGeneratorAction* _primGenerator;
};
#endif
