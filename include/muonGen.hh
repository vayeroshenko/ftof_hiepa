#ifndef muonGen_hh
#define muonGen_hh

#include <TROOT.h>

class crtTrk;
class plane;
class TRandom3;
class crtDetector;

class muonGen {

public :

  muonGen(Int_t mySeed);
  ~muonGen();

  //crtTrk* GenerateMuon(Double_t zz);
  void GenerateMuon(Double_t zz);
  Double_t genCos2dist();

  inline Double_t GetX(){return _x;}
  inline Double_t GetY(){return _y;}
  inline Double_t GetZ(){return _z;}
  inline Double_t GetTheta(){return _Theta;}
  inline Double_t GetPhi(){return _Phi;}

private :
  TRandom3 *_rnd;
  crtTrk *_trk; 

  Double_t _x;
  Double_t _y;
  Double_t _z;
  Double_t _Theta;//deg
  Double_t _Phi;//deg

  Bool_t Trigger();
  
  plane *_Trigg1;
  plane *_Trigg2;
  plane *_QuartzStart;
  plane *_fTOFshredTOP;
  plane *_fTOFshredBOT;
  plane *_fTOF;
};

#endif
