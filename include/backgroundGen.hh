#ifndef backgroundGen_h
#define backgroundGen_h

//root
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

//C, C++
#include <stdio.h>

using namespace std;

class backgroundGen {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           pID;
   Int_t           evID;
   Int_t           uniqID;
   Int_t           baunchXID;
   Double_t        MomX;
   Double_t        MomY;
   Double_t        MomZ;
   Double_t        PosX;
   Double_t        PosY;
   Double_t        PosZ;
   Double_t        Time;

   // List of branches
   TBranch        *b_pID;   //!
   TBranch        *b_evID;   //!
   TBranch        *b_uniqID;   //!
   TBranch        *b_baunchXID;   //!
   TBranch        *b_MomX;   //!
   TBranch        *b_MomY;   //!
   TBranch        *b_MomZ;   //!
   TBranch        *b_PosX;   //!
   TBranch        *b_PosY;   //!
   TBranch        *b_PosZ;   //!
   TBranch        *b_Time;   //!

   backgroundGen(TTree *tree=0);
   virtual ~backgroundGen();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);



   Long64_t nentries;
  Long64_t nbytes;
  Long64_t nb;




  Long64_t _jentry;
  void GenNextPart(){_jentry++;}

};

#endif



