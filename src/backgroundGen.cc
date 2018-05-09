//My
#include "backgroundGen.hh"

//root
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>

//C, C++
#include <stdio.h>
#include <iostream>
#include <assert.h>

using namespace std;


backgroundGen::backgroundGen(TTree *tree)
{

  _jentry = 0;
  nbytes = 0;
  nb = 0;

  string fName = "./ana_back/Data/BackgroundWorld_SuperB_Wolf_shielded_w.root";
  //string fName = "./ana_back/Data/BackgroundWorld_SuperB_unshielded_w.root";
  cout<<fName<<endl;

  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fName.c_str());
    if (!f) {
      f = new TFile(fName.c_str());
    }
    tree = (TTree*)gDirectory->Get("tree");
  }
  Init(tree);

  nentries = fChain->GetEntriesFast();
}


backgroundGen::~backgroundGen()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}


void backgroundGen::Loop()
{

  if (fChain == 0){ 
    cout<<"fChain == 0"<<endl;
    assert(0);
  }

  //for (Long64_t jentry=0; jentry<nentries;jentry++) {
  LoadTree(_jentry);
  //if (ientry < 0) break;
  nb = fChain->GetEntry(_jentry);   nbytes += nb;
  if(_jentry%500==0){
    cout<<_jentry<<endl;
    if( nentries<_jentry)
      cout<<" nentries<_jentry "<<endl;
  }
  _jentry++;
  //}
  //GenNextPart();
  
}

Int_t backgroundGen::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t backgroundGen::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void backgroundGen::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);
   fChain->SetBranchAddress("pID", &pID, &b_pID);
   fChain->SetBranchAddress("evID", &evID, &b_evID);
   fChain->SetBranchAddress("uniqID", &uniqID, &b_uniqID);
   fChain->SetBranchAddress("baunchXID", &baunchXID, &b_baunchXID);
   fChain->SetBranchAddress("MomX", &MomX, &b_MomX);
   fChain->SetBranchAddress("MomY", &MomY, &b_MomY);
   fChain->SetBranchAddress("MomZ", &MomZ, &b_MomZ);
   fChain->SetBranchAddress("PosX", &PosX, &b_PosX);
   fChain->SetBranchAddress("PosY", &PosY, &b_PosY);
   fChain->SetBranchAddress("PosZ", &PosZ, &b_PosZ);
   fChain->SetBranchAddress("Time", &Time, &b_Time);
   Notify();
}

Bool_t backgroundGen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   return kTRUE;
}

void backgroundGen::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t backgroundGen::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
