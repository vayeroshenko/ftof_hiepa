#ifndef plane_hh
#define plane_hh

//root
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TVector3.h>

using namespace std;

class crtTrk;
class TString;

class plane {
public :

  plane(TString name, Double_t dx,Double_t dy,
	Double_t x0,Double_t y0,Double_t z0, 
	Double_t angleY);
  ~plane();

  void GetIntersection(crtTrk *trk);
  
  inline TString GetName(){return _Name;};
  inline Double_t GetXint(){return _xInt;};
  inline Double_t GetYint(){return _yInt;};
  inline Double_t GetZint(){return _zInt;};
  inline Int_t GetIntersecStatus(){return _interStatus;};

private:

  //sizes
  Double_t _dx;//cm
  Double_t _dy;//cm

  //rotated angle
  Double_t _angleY;//deg

  TVector3 _n;
  TVector3 _V0;

  TString _Name;

  //parameters of trk intersection with plain 
  Double_t _xInt;//cm
  Double_t _yInt;//cm
  Double_t _zInt;//cm

  Int_t _interStatus;
  //-999 staus not defined
  //-1 no intesection
  // 0 infinit number of intersections
  // 1 one intersection with plane
  // 2 one intersection within sizes of the box

};

#endif
