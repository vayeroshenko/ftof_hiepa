#ifndef crtTrk_hh
#define crtTrk_hh

//root
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TVector3.h>

using namespace std;

class crtTrk {
public :
  
  crtTrk();
  crtTrk(Double_t x1, Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2);
  crtTrk(Double_t x, Double_t y, Double_t z, Double_t theta, Double_t phi);
  ~crtTrk();
  inline Double_t getTrkTheta(){return _Theta;};
  inline Double_t getTrkPhi(){return _Phi;};
  inline TVector3 GetTrkDirectionV(){return _a;};
  inline TVector3 GetTrkPointV(){return _r;};
  
private:
  
  TVector3 _a;//direction vector 
  TVector3 _r;//point belong to the track
  Double_t _Theta;//deg
  Double_t _Phi;  //deg
  
};

#endif
