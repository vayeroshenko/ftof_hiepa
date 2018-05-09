//my
#include "plane.hh"
#include "crtTrk.hh"

//root
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TMath.h>
#include <TLine.h>
#include <TString.h>

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

plane::plane(TString name, Double_t dx,Double_t dy, 
	     Double_t x0,Double_t y0,Double_t z0, 
	     Double_t angleY){
  _dx = dx;
  _dy = dy;
  _angleY = angleY;
  _V0.SetXYZ(x0,y0,z0);
  _n.SetXYZ(0.0,0.0,1.0);
  if(angleY!=0.0)
    _n.RotateY(angleY/180.0*TMath::Pi());

  if(_n.Z()<0.0){
    cout<<" ERROR--> _n.Z()<0.0 "<<_n.Z()<<endl;
    assert(0);
  }
  
  _xInt = -999.0;
  _yInt = -999.0;
  _zInt = -999.0;
  _interStatus = -999;
  
  _Name = name;  
}

plane::~plane(){
}

void plane::GetIntersection(crtTrk *trk){

  Double_t t; //parameter of the trk
  _xInt = -999.0;
  _yInt = -999.0;
  _zInt = -999.0;    
  
  TVector3 trkDir = trk->GetTrkDirectionV();
  TVector3 trkPoint = trk->GetTrkPointV();

  Double_t znam = (trkDir.Dot(_n));
  Double_t chus = (_V0 - trkPoint).Dot(_n);
  
  if( znam == 0.0 && chus == 0.0 ){
    _xInt = -999.0;
    _yInt = -999.0;
    _zInt = -999.0;    
    _interStatus = 0;
  }
  else if( znam == 0.0 && chus != 0.0 ){
    _xInt = -999.0;
    _yInt = -999.0;
    _zInt = -999.0;
    _interStatus = -1;
  }
  else if(znam != 0.0){
    t = chus/znam;
    _interStatus = 1;
    TVector3 r = (trkPoint + trkDir*t);
    _xInt = r.X();
    _yInt = r.Y();
    _zInt = r.Z();
    
    if((_yInt - _V0.Y()> -1.0*_dy/2.0) && (_yInt - _V0.Y()< _dy/2.0) &&
       (_xInt - _V0.X()> -1.0*_dx*TMath::Cos(TMath::Abs(_angleY/180.0*TMath::Pi()))/2.0) && (_xInt - _V0.X()< _dx*TMath::Cos(TMath::Abs(_angleY/180.0*TMath::Pi()))/2.0)){
      _interStatus = 2;
    }

  }
  
  //cout<<endl<<"crtBox::GetIntersection   "<< _Name<<endl
  //  <<" _xInt        = "<<_xInt<<endl
  //  <<" _yInt        = "<<_yInt<<endl
  //  <<" _zInt        = "<<_zInt<<endl
  //  <<" _interStatus = "<<_interStatus<<endl;

}

