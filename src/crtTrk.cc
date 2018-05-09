//my
#include "crtTrk.hh"
#include "crtConst.hh"

//root
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMath.h>
#include <TLine.h>
#include <TVector3.h>

//C, C++
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <math.h>

crtTrk::crtTrk(){
  _a.SetMagThetaPhi(0.0,0.0,0.0);
  _r.SetMagThetaPhi(0.0,0.0,0.0);
  _Theta = -999.0;
  _Phi = -999.0;
}

crtTrk::~crtTrk(){
}

crtTrk::crtTrk(Double_t x,Double_t y, Double_t z, Double_t theta, Double_t phi){
  _Theta = theta;
  _Phi = phi;
  _a.SetMagThetaPhi(1.0,_Theta*TMath::Pi()/180.0,_Phi*TMath::Pi()/180.0);
  _r.SetXYZ(x,y,z);
}

crtTrk::crtTrk(Double_t x1,Double_t y1, Double_t z1, Double_t x2, Double_t y2, Double_t z2){
  TVector3 r1(x1,y1,z1);
  TVector3 r2(x2,y2,z2);
  TVector3 dr = (r1 - r2);
  Double_t magdr = dr.Mag();
  if(magdr>0.0)
    dr.SetXYZ(dr.x()/magdr,dr.y()/magdr,dr.z()/magdr);
  else{
    cout<<" ERROR --> "<<magdr<<endl;
    assert(0);
  }
  _a = dr;
  _r = r1;
  _Theta = _a.Theta()*180.0/TMath::Pi();
  _Phi = _a.Phi()*180.0/TMath::Pi();  
}


