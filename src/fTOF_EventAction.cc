//My
#include "fTOF_EventAction.hh"
#include "fTOF_RunAction.hh"
#include "fTOF_SteppingAction.hh"
#include "fTOF_Hit.hh"
#include "fTOF_StackingAction.hh"

//G4
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "globals.hh"

fTOF_EventAction::fTOF_EventAction(fTOF_RunAction* runact,
				     fTOF_SteppingAction* steppingAction) :
  runAction(runact), _steppingAction(steppingAction), printModulo(100)
{
  thePhotonCollectionID = -1;
}

fTOF_EventAction::~fTOF_EventAction()
{
}

void fTOF_EventAction::BeginOfEventAction(const G4Event* event)
{
  // Print number of events
  G4int eventNum = event->GetEventID();

  if (eventNum%printModulo == 0) {
    G4cout << "\n---> Begin of Event: " << eventNum << G4endl;
  }

  if (thePhotonCollectionID < 0) {
    G4String colName;
    thePhotonCollectionID = 
      G4SDManager::GetSDMpointer()->
      GetCollectionID(colName="OpticalPhotonCollection");
  }

  _steppingAction->Reset();
  _steppingAction->ResetPerEvent();
}

void fTOF_EventAction::EndOfEventAction(const G4Event* event)
{
  // Print info about end of the event
  G4int eventNum = event->GetEventID();

  if (thePhotonCollectionID < 0) return;

  // Get the Hit Collection
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  fTOF_HitsCollection * THC = 0;

  if (HCE)
    THC = (fTOF_HitsCollection*)(HCE->GetHC(thePhotonCollectionID));

  if (0 == THC) return;

  G4int nTotPhot = _stackingAction->GetTotPhotNum();
  G4int nHit = -1;
  nHit = THC->entries();

  runAction->EventInfo.EventID = eventNum;
  runAction->EventInfo.BunchXID = _primGenerator->GetBunchXID();
  runAction->EventInfo.NTotPhot = nTotPhot;
  runAction->EventInfo.Nhits = nHit;
  runAction->EventInfo.primType = event->GetPrimaryVertex()->GetPrimary()->GetPDGcode();
  //G4cout<<"primType = "<<event->GetPrimaryVertex()->GetPrimary()->GetPDGcode()<<G4endl;
  runAction->EventInfo.primMomX = event->GetPrimaryVertex()->GetPrimary()->GetPx();
  runAction->EventInfo.primMomY = event->GetPrimaryVertex()->GetPrimary()->GetPy();
  runAction->EventInfo.primMomZ = event->GetPrimaryVertex()->GetPrimary()->GetPz();
  runAction->EventInfo.primPosX = event->GetPrimaryVertex()->GetX0();
  runAction->EventInfo.primPosY = event->GetPrimaryVertex()->GetY0();
  runAction->EventInfo.primPosZ = event->GetPrimaryVertex()->GetZ0();
  runAction->EventInfo.primTime = event->GetPrimaryVertex()->GetT0();

  runAction->_nPhot = nHit;
  for (G4int i = 0; i < nHit; i++) {
    runAction->_TrackID[i] = (*THC)[i]->myData.TrackID;
    runAction->_ParentID[i] = (*THC)[i]->myData.ParentID;
    runAction->_Energy[i] = (*THC)[i]->myData.Energy;
    runAction->_Wavelength[i] = (*THC)[i]->myData.Wavelength;
    runAction->_Time[i] = (*THC)[i]->myData.Time;
    runAction->_photPathLen[i] = (*THC)[i]->myData.photPathLen;
    runAction->_SecID[i] = (*THC)[i]->myData.SecID;
    runAction->_chID[i] = (*THC)[i]->myData.chID;
    runAction->_PosX[i] = (*THC)[i]->myData.PosX;
    runAction->_PosY[i] = (*THC)[i]->myData.PosY;
    runAction->_PosZ[i] = (*THC)[i]->myData.PosZ;
    runAction->_trkMomX[i] = (*THC)[i]->myData.trkMomX;
    runAction->_trkMomY[i] = (*THC)[i]->myData.trkMomY;
    runAction->_trkMomZ[i] = (*THC)[i]->myData.trkMomZ;
    runAction->_trkPosX[i] = (*THC)[i]->myData.trkPosX;
    runAction->_trkPosY[i] = (*THC)[i]->myData.trkPosY;
    runAction->_trkPosZ[i] = (*THC)[i]->myData.trkPosZ;
    runAction->_trkT[i] = (*THC)[i]->myData.trkT;
    runAction->_trkLength[i] = (*THC)[i]->myData.trkLength;

    runAction->_entMomX[i] = (*THC)[i]->myData.entMomX;
    runAction->_entMomY[i] = (*THC)[i]->myData.entMomY;
    runAction->_entMomZ[i] = (*THC)[i]->myData.entMomZ;
    runAction->_entPosX[i] = (*THC)[i]->myData.entPosX;
    runAction->_entPosY[i] = (*THC)[i]->myData.entPosY;
    runAction->_entPosZ[i] = (*THC)[i]->myData.entPosZ;





    runAction->_entTime[i] = (*THC)[i]->myData.entTime;

    runAction->_nSideRefl[i] = (*THC)[i]->myData.trkNSideRefl;
    runAction->_sideID[i] = (*THC)[i]->myData.trkSideID;




    // Fill the structure then save
    //runAction->HitInfo = (*THC)[i]->myData;
    //runAction->tree->Fill();
  }

  runAction->tree->Fill(); 
  
}
