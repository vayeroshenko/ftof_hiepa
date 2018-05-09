#ifndef fTOF_PhysicsList_h
#define fTOF_PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class G4Cerenkov;
class G4OpAbsorption;
class G4OpRayleigh;
class G4OpBoundaryProcess;

class fTOF_PhysicsList: public G4VUserPhysicsList
{
public:
  fTOF_PhysicsList();
  ~fTOF_PhysicsList();
  
protected:
  // Construct particle and physics process
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

protected:
  // These methods construct particles
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();

protected:
  // Construct phyisics processes and register them
  void ConstructGeneral();
  void ConstructEM();
  void ConstructOp();

protected:
  void SetNbOfPhotonsCerenkov(G4int);
  void SetVerbose(G4int verbose);
  
private:
  G4Cerenkov* theCerenkovProcess;
  G4OpAbsorption* theAbsorptionProcess;
  G4OpRayleigh* theRayleighScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;


};

#endif







