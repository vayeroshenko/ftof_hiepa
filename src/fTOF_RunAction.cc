//My
#include "fTOF_RunAction.hh"

//G4
#include "G4Run.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"

//root
#include "TFile.h"
#include "TTree.h"

fTOF_RunAction::fTOF_RunAction() : // @suppress("Class members should be properly initialized")
  timer(0),
  _outputFileName("fTOF.root")
  //  tree(0),
  //  hfile(0)
{
  timer = new G4Timer;
}

fTOF_RunAction::~fTOF_RunAction()
{
  delete timer;
}

void fTOF_RunAction::BeginOfRunAction(const G4Run* aRun)
{
  timer->Start();

  // Histogramming
  hfile = new TFile(_outputFileName, "RECREATE", "fTOF Simulation Data", 1);
  if (hfile->IsZombie()) exit(-1);
  tree = new TTree("T", "fTOF Data Tree");
  
  hfile->SetCompressionLevel(2);
  tree->SetAutoSave(1000000);

  // Create new event
  TTree::SetBranchStyle(0);
  //Event
  tree->Branch("EventID",  &EventInfo.EventID,  "Event/I");
  tree->Branch("BunchXID", &EventInfo.BunchXID, "BunchXID/I");
  tree->Branch("NTotPhot", &EventInfo.NTotPhot, "NTotPhot/I");
  tree->Branch("Nhits",    &EventInfo.Nhits,    "Nhits/I");
  tree->Branch("primType", &EventInfo.primType, "primType/I");
  tree->Branch("primMomX", &EventInfo.primMomX, "primMomX/D");
  tree->Branch("primMomY", &EventInfo.primMomY, "primMomY/D");
  tree->Branch("primMomZ", &EventInfo.primMomZ, "primMomZ/D");
  tree->Branch("primPosX", &EventInfo.primPosX, "primPosX/D");
  tree->Branch("primPosY", &EventInfo.primPosY, "primPosY/D");
  tree->Branch("primPosZ", &EventInfo.primPosZ, "primPosZ/D");
  tree->Branch("primTime", &EventInfo.primTime, "primTime/D");




  //Hits

  //tree->Branch("TrackID",    &HitInfo.TrackID,    "TrackID/I");
  //tree->Branch("ParentID",   &HitInfo.ParentID,   "ParentID/I");
  //tree->Branch("Energy",     &HitInfo.Energy,     "Energy/D");
  //tree->Branch("Wavelength", &HitInfo.Wavelength, "Wavelength/D");
  //tree->Branch("Time",       &HitInfo.Time,       "Time/D");
  //tree->Branch("photPathLen",&HitInfo.photPathLen,"photPathLen/D");
  //tree->Branch("SecID",      &HitInfo.SecID,      "SecID/I");
  //tree->Branch("chID",       &HitInfo.chID,       "chID/I");
  //tree->Branch("PosX",       &HitInfo.PosX,       "PosX/D");
  //tree->Branch("PosY",       &HitInfo.PosY,       "PosY/D");
  //tree->Branch("PosZ",       &HitInfo.PosZ,       "PosZ/D");
  //tree->Branch("trkMomX",    &HitInfo.trkMomX,    "trkMomX/D");
  //tree->Branch("trkMomY",    &HitInfo.trkMomY,    "trkMomY/D");
  //tree->Branch("trkMomZ",    &HitInfo.trkMomZ,    "trkMomZ/D");
  //tree->Branch("trkPosX",    &HitInfo.trkPosX,    "trkPosX/D");
  //tree->Branch("trkPosY",    &HitInfo.trkPosY,    "trkPosY/D");
  //tree->Branch("trkPosZ",    &HitInfo.trkPosZ,    "trkPosZ/D");
  //tree->Branch("trkT",       &HitInfo.trkT,       "trkT/D");
  //tree->Branch("trkLength",  &HitInfo.trkLength,  "trkLength/D");

  tree->Branch("nPhot", &_nPhot, "nPhot/I");
  tree->Branch("TrackID", _TrackID, "TrackID[nPhot]/I");
  tree->Branch("ParentID", _ParentID, "ParentID[nPhot]/I");
  tree->Branch("Energy", _Energy, "Energy[nPhot]/D");
  tree->Branch("Wavelength", _Wavelength, "Wavelength[nPhot]/D");
  tree->Branch("Time", _Time, "Time[nPhot]/D");
  tree->Branch("photPathLen", _photPathLen,"photPathLen[nPhot]/D");
  tree->Branch("SecID",_SecID, "SecID[nPhot]/I");
  tree->Branch("chID", _chID, "chID[nPhot]/I");
  tree->Branch("PosX", _PosX, "PosX[nPhot]/D");
  tree->Branch("PosY", _PosY, "PosY[nPhot]/D");
  tree->Branch("PosZ", _PosZ, "PosZ[nPhot]/D");
  tree->Branch("trkMomX", _trkMomX, "trkMomX[nPhot]/D");
  tree->Branch("trkMomY", _trkMomY, "trkMomY[nPhot]/D");
  tree->Branch("trkMomZ", _trkMomZ, "trkMomZ[nPhot]/D");
  tree->Branch("trkPosX", _trkPosX, "trkPosX[nPhot]/D");
  tree->Branch("trkPosY", _trkPosY, "trkPosY[nPhot]/D");
  tree->Branch("trkPosZ", _trkPosZ, "trkPosZ[nPhot]/D");
  tree->Branch("trkT", _trkT, "trkT[nPhot]/D");
  tree->Branch("trkLength",  _trkLength, "trkLength[nPhot]/D");


  tree->Branch("entPosX", _entPosX, "entPosX[nPhot]/D");
  tree->Branch("entPosY", _entPosY, "entPosY[nPhot]/D");
  tree->Branch("entPosZ", _entPosZ, "entPosZ[nPhot]/D");

  tree->Branch("entMomX", _entMomX, "entMomX[nPhot]/D");
  tree->Branch("entMomY", _entMomY, "entMomY[nPhot]/D");
  tree->Branch("entMomZ", _entMomZ, "entMomZ[nPhot]/D");


  tree->Branch("entTime", _entTime, "entTime[nPhot]/D");

  tree->Branch("nSideRefl", _nSideRefl, "nSideRefl[nPhot]/I");
  tree->Branch("sideId", _sideID, "SideID[nPhot]/I");

}

void fTOF_RunAction::EndOfRunAction(const G4Run* aRun)
{

  if (G4VVisManager::GetConcreteInstance()) {
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
  }

  hfile = tree->GetCurrentFile();
  hfile->Write();
  tree->Print();
  timer->Stop();

  delete tree;
  delete hfile;
  G4cout << "Time: " << *timer << G4endl;
}
