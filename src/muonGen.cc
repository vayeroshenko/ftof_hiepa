//my
#include "muonGen.hh"
#include "crtConst.hh"
#include "crtTrk.hh"
#include "plane.hh"

//root
#include <TRandom3.h>

#include <iostream>
#include <assert.h>

using namespace std;

muonGen::muonGen(Int_t mySeed){ 
  _rnd = new TRandom3(mySeed);
  _x = -999.0;
  _y = -999.0;
  _z = -999.0;
  _Theta = -999.0;
  _Phi = -999.0;

  _Trigg1 = new plane("Trigg1",(crtConst::hodoYBar_lenght),(crtConst::hodoXBar_lenght), 
		      0.0, 0.0, 0.0, 
		      0.0);

  _Trigg2 = new plane("Trigg2",(crtConst::hodoYBar_lenght),(crtConst::hodoXBar_lenght), 
		      0.0, 0.0, (crtConst::hodo_hight), 
		      0.0);


  _QuartzStart = new plane("QuartzStart",crtConst::quartzBar_lenght,crtConst::quartzBar_width,
			   4.007,0.0,(16.367+83.5),
			   -43.0);


  //_fTOFshredTOP = new plane("fTOFshred",1.0,1.0,
  //		   -6.0,0.0,(83.5 + 38.1 + 1.49 + 1.5 + 1.64 + 1.5),
  //		   0.0);

  //_fTOFshredBOT = new plane("fTOFshred",1.0,1.0,
  //		   -6.0,0.0,(83.5 + 38.1 + 1.49 + 1.5 ),
  //		   0.0);

  _fTOFshredTOP = new plane("fTOFshred",crtConst::quartzBar_lenght, crtConst::quartzBar_width,
  		   0.55,0.0,(83.5 + 38.1 + 1.49 + 1.5 + 1.64 + 1.5),
  		   0.0);

  _fTOFshredBOT = new plane("fTOFshred",crtConst::quartzBar_lenght, crtConst::quartzBar_width,
  		   0.55,0.0,(83.5 + 38.1 + 1.49 + 1.5 ),
  		   0.0);

  _fTOF = new plane("fTOF",crtConst::quartzBar_lenght,crtConst::quartzBar_width,
		    0.55,0.0,(83.5 + 38.1 + 1.49 + 1.5 + 1.64/2.0 ),
		    0.0);
  
}

muonGen::~muonGen(){
}

void muonGen::GenerateMuon(Double_t zz){
  //delete _trk;
  const Int_t nMax = 1000000;
  Bool_t trkIsOK = false;
  for(Int_t i = 0;i<nMax;i++){
    //_x = _rnd->Uniform(-1.0*crtConst::hodoYBar_lenght/2.0 - crtConst::topHodo_xShift,crtConst::hodoYBar_lenght/2.0 + crtConst::topHodo_xShift);
    //_y = _rnd->Uniform(-1.0*crtConst::hodoXBar_lenght/2.0 - crtConst::topHodo_yShift,crtConst::hodoXBar_lenght/2.0 + crtConst::topHodo_yShift);
    _x = _rnd->Uniform(-2.0*crtConst::hodoYBar_lenght, 2.0*crtConst::hodoYBar_lenght);
    _y = _rnd->Uniform(-2.0*crtConst::hodoXBar_lenght, 2.0*crtConst::hodoXBar_lenght);
    //_x = _rnd->Uniform(-6.0-20.0, -6.0+20.0);
    //_y = _rnd->Uniform(-20.0, 20.0);
    // _x = _rnd->Uniform(-35.0, 35.0);
    //_y = _rnd->Uniform(-25.0, 25.0);

    _z = zz;
    _Theta = 180.0 - genCos2dist();     //deg
    //_Theta = 177.0;     //deg
    _Phi = _rnd->Uniform(-180.0,180.0); //deg
    _trk = new crtTrk(_x,_y,_z,_Theta,_Phi);  
    if(Trigger()){
      i = nMax;
      trkIsOK = true;
      //cout<<"_x     = "<<_x<<endl
      //  <<"_y     = "<<_y<<endl
      //  <<"_z     = "<<_z<<endl
      //  <<"_Theta = "<<_Theta<<endl
      //  <<"_Phi   = "<<_Phi<<endl;
    }
    if(i != nMax)
      delete _trk;
  }

  if(!trkIsOK){
    cout<<" ERROR --> trkIsOK == false, nMax ="<<nMax<<endl
	<<"           acseptance is too small"<<endl;
    assert(0);
  }



  //return _trk;
}

Double_t muonGen::genCos2dist(){
  Double_t theta = -999.0;//deg 
  Double_t x = -999.0;
  Double_t y = -999.0;
  while(theta==-999.0){
    x = _rnd->Uniform(0.0,70.0*TMath::Pi()/180.0); //rad
    y = _rnd->Uniform(0.0,1.1); //rad
    if(TMath::Power(TMath::Cos(x),1.85)>y){
      theta = x*180.0/TMath::Pi();
    }
  }  
  return theta;
}

Bool_t muonGen::Trigger(){

  _Trigg1->GetIntersection(_trk);
  _Trigg2->GetIntersection(_trk);
  _QuartzStart->GetIntersection(_trk);
  _fTOFshredTOP->GetIntersection(_trk);
  _fTOFshredBOT->GetIntersection(_trk);
  //_fTOF->GetIntersection(_trk);

  //if(_Trigg1->GetIntersecStatus()==2 &&
  // _Trigg2->GetIntersecStatus()==2 &&
  // _QuartzStart->GetIntersecStatus()==2){
  //return true;
  //}

  if(_Trigg1->GetIntersecStatus()==2 &&
     _Trigg2->GetIntersecStatus()==2 &&
     _QuartzStart->GetIntersecStatus()==2 && 
     _fTOFshredTOP->GetIntersecStatus()==2 &&
     _fTOFshredBOT->GetIntersecStatus()==2){
    return true;
  }
  
  //if(_Trigg1->GetIntersecStatus()==2 &&
  // _Trigg2->GetIntersecStatus()==2 &&
  // _QuartzStart->GetIntersecStatus()==2 && 
  // _fTOFshredTOP->GetIntersecStatus()==2){
  //return true;
  //}
  
  //if(_QuartzStart->GetIntersecStatus()==2 && 
  //   _fTOF->GetIntersecStatus()==2){
  //  return true;
  //}
  
  //if(_fTOF->GetIntersecStatus()==2){
  //return true;
  //}
  
  return false;
}
