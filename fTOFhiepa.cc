#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

//#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
//#endif

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

#include "fTOF_DetectorConstruction.hh"
#include "fTOF_PhysicsList.hh"
#include "fTOF_PrimaryGeneratorAction.hh"

#include "fTOF_RunAction.hh"
#include "fTOF_EventAction.hh"
#include "fTOF_SteppingAction.hh"
#include "fTOF_StackingAction.hh"
#include "fTOF_TrackingAction.hh"
//#include "fTOF_UImessenger.hh"
#include "fTOF_SteppingVerbose.hh"

#include "Randomize.hh"

int main(int argc, char** argv)
{

  // Choose the Random Engine
  //  HepRandom::setTheEngine(new RanecuEngine);
  // Seed the random number generator manually
  //

  if(argc!=5){
    G4cout<<" ERROR of the input parameters !!! "<<G4endl
	  <<"      [0] - vis.mac or fTOF.mac or *.mac "<<G4endl
	  <<"      [1] - seed "<<G4endl
	  <<"      [2] - output root file name"<<G4endl
	  <<"      [3] - name of the particle (e+, e-, mu+, mu-, pi+, pi-, kaon+, kaon-, proton+, proton-, gamma)"<<G4endl;
    return 0;    
  }
  else{
    G4cout<<"  Input parameters         "<<G4endl
	  <<"     mac file              "<<argv[1]<<G4endl
	  <<"     seed                  "<<argv[2]<<G4endl
	  <<"     output root file name "<<argv[3]<<G4endl
	  <<"     name of the particle  "<<argv[4]<<G4endl;
  }


  G4long myseed = 345354;
  myseed = atoi(argv[2]);
  G4cout<<" myseed - "<<myseed<<G4endl;

  CLHEP::HepRandom::setTheSeed(myseed);

  // Verbose output class
  G4VSteppingVerbose::SetInstance(new fTOF_SteppingVerbose);

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager;

  // Construct mandatory initialization classes
  G4VUserDetectorConstruction* detector = new fTOF_DetectorConstruction;
  runManager->SetUserInitialization(detector);
  G4VUserPhysicsList* physics = new fTOF_PhysicsList;
  runManager->SetUserInitialization(physics);

  // Construct User Action Classes
  fTOF_RunAction* runAction = new fTOF_RunAction;
  runManager->SetUserAction(runAction);
  //
  G4String rootFileName = argv[3];
  G4cout<<rootFileName<<G4endl;
  runAction->SetOutputFileName(rootFileName);

  fTOF_PrimaryGeneratorAction* genAction = 
    new fTOF_PrimaryGeneratorAction();
  runManager->SetUserAction(genAction);
  G4String particleName = argv[4];
  G4cout<<"particle Name = "<<particleName<<G4endl;
  genAction->SetParticleName(particleName);


  fTOF_SteppingAction* steppingAction = new fTOF_SteppingAction(genAction);
  runManager->SetUserAction(steppingAction);

  fTOF_EventAction* eventAction = new fTOF_EventAction(runAction,
							 steppingAction);
  runManager->SetUserAction(eventAction);
  fTOF_StackingAction *stackingAction = new fTOF_StackingAction;
  runManager->SetUserAction(stackingAction);
  runManager->SetUserAction(new fTOF_TrackingAction);

  eventAction->SetStackingAction(stackingAction);
  eventAction->SetPrimGenerator(genAction);
  // Setup to be able to define some custom commands;
  //fTOF_UImessenger* messenger = new fTOF_UImessenger(runAction, genAction);

  // Initialize G4 kernel
  runManager->Initialize();
    // delete runManager;
    G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();
  if (argc == 1) {  // Define UI terminal for interactive mode
    UI->ApplyCommand("/run/verbose 0");
    UI->ApplyCommand("/event/verbose 0");
    UI->ApplyCommand("/tracking/verbose 0");
#if defined(G4UI_USE_TCSH)
    G4UIsession* session = new G4UIterminal(new G4UItcsh);
#else
    G4UIsession* session = new G4UIterminal;
#endif
    session->SessionStart();
    delete session;
  }
  else {

    G4String fileName = argv[1];
    if(fileName.contains("vis")){
      G4cout<<"VIS fileName "<<fileName<<G4endl;
#ifdef G4VIS_USE
      // Visualization manager
      // G4VisManager* visManager = new G4VisExecutive;
      // visManager->Initialize();
#endif    
      
      G4String command = "/control/execute ";
      UI->ApplyCommand(command+fileName);
      
#ifdef G4VIS_USE
      delete visManager;
#endif
    }
    else{
      G4cout<<"NO VIS fileName "<<fileName<<G4endl;
      G4String command = "/control/execute ";
      UI->ApplyCommand(command+fileName);    
    }

  }

  

  // Start a run
  //  G4int numberOfEvent = 3;
  //  runManager->BeamOn(numberOfEvent);

  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //delete messenger;
  //delete visManager;
  // delete runManager;
  
  return 0;
}


