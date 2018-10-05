#ifndef fTOFConst_h
#define fTOFConst_h

#include "TMath.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

namespace fTOFConst{

const G4double cmMy = 10;
const G4double mmMy = 1;

//Sec
const G4double barSizeX = 20.*cm;
const G4double barSizeY = 2.0*cm;
const G4double barSizeZ = 5.0*cm;

//update from Jerry 20.11.2010
//const G4double barDeltaPosY = 5.0*mm;
const G4double barDeltaPosY = 1.8*mm;

//PmtWindow
const G4double pmtWin1SizeX = 1.93*mm;
const G4double pmtWin1SizeY = 64.0*mm;
const G4double pmtWin1SizeZ = 64.0*mm;

const G4double pmtWin2SizeX = 6.88*mm - pmtWin1SizeX;
const G4double pmtWin2SizeY = 50.0*mm;
const G4double pmtWin2SizeZ = 50.0*mm;

//channel
const G4double pmtChSizeX = 1.5*mm;
const G4double pmtChSizeY = 18.0*mm;
const G4double pmtChSizeZ = 6.0*mm;
const G4double pmtChGap = 0.1*mm;

//pmt Abs1
const G4double pmtAbs1SizeX = 1.5*mm;
const G4double pmtAbs1SizeY = 12.0*mm;
const G4double pmtAbs1SizeZ = pmtWin2SizeZ;

//pmt Box size
const G4double pmtBoxSizeX = 20.4*mm;
const G4double pmtBoxSizeY = 64.0*mm;
const G4double pmtBoxSizeZ = 64.0*mm;
const G4double pmtBoxGapp = 0.1*mm;

//abs1
const G4double abs1SizeX = 1.5*mm;
const G4double abs1SizeY = barSizeY;
const G4double abs1SizeZ = barSizeZ;

//abs2
const G4double abs2SizeX = barSizeX;
const G4double abs2SizeY = barSizeY;
const G4double abs2SizeZ = 1.5*mm;


// fTOF bar 10.04.18
const G4double barBoxSizeX = 19.8*cm;
const G4double barBoxSizeY = 1.8*cm;
const G4double barBoxSizeZ = 4.8*cm;

const G4double fullBarSizeX = 20.*cm;
const G4double fullBarSizeY = 2.*cm;
const G4double fullBarSizeZ = 5.*cm;

const G4double barChamfer = 2.*mm;

// planacon MCP PMT
const G4double planLeftDist = 3.*mm;
const G4double planRightDist = 5.*mm;
const G4double planUpDist = 20.*mm;
const G4double planDownDist = 19.*mm;

const G4double planChanThick = 1.*mm;
const G4double planChanSize = 13.25*mm;
const G4double planWindowThick = 1.5*mm;  // have no data, aproximatelly same as hamamatsu
const G4double planBoxSize = 59.*mm;
const G4double planSensitiveSize = 53.*mm;
const G4double planBoxDepth = 13.7*mm;



// hamamatsu MCP PMT
const G4double hamChanSize = 5.275*mm;
const G4double hamChanThick = 1.*mm;
const G4double hamChanGap = 0.3*mm;
const G4double hamUpDist = 5.*mm;
const G4double hamRightDist = 12.*mm;


const G4double hamSensitiveSize = 22.5*mm;
const G4double hamWindowThick = 1.5*mm;

const G4double hamBoxSize = 27.5*mm;
const G4double hamBoxDepth = 13.*mm;





// sector
// const G4double sectorLongSide = 45.*cm;
// const G4double sectorShortSide = 24.*cm;
const G4double sectorHeight = 40.*cm;
const G4double sectorThickness = 1*cm;

const G4double sectorLongSide = 3.*cm;
const G4double sectorShortSide = 2.*cm;
// const G4double sectorHeight = 0.5*cm;
// const G4double sectorThickness = 1.5*cm;


// const G4double sectorLongSide = 5.*cm;
// const G4double sectorShortSide = 2.*cm;

const G4double mixerWindowX = 20.*mm;
const G4double mixerWindowZ = 15.*mm;
const G4double mixerLength = 15.*cm;


const G4int nSec = 101;
const G4int nDrawSec = 1;

const G4int nLayers = 1;

const G4double innerRad = 47.*cm;
const G4double outerRad = 105.*cm;
const G4double centerRad = (innerRad * TMath::Cos(TMath::Pi() / nSec) +
                            outerRad * TMath::Cos(TMath::Pi() / nSec)) / 2.;

const G4double absInnerSide = 0.1*mm;
const G4double absOuterSide = absInnerSide * outerRad / innerRad;



const G4double innerSide = 2. * innerRad * TMath::Sin(TMath::Pi() / nSec) - absInnerSide;
const G4double outerSide = 2. * outerRad * TMath::Sin(TMath::Pi() / nSec) - absOuterSide;

const G4double layerDist = 1.*mm;

const G4double layerThickness = 0.1*mm;

const G4String particleName = "kaon+";
//const G4String particleName = "pi+";
const G4double particleMomentum = 2*GeV;

const G4double angle = 14 *deg;

//  const G4String particleName = "pi+";



//const G4int  N_PMT = 9;
//const G4int    N_PMT = 14;
//const G4int    N_ch  = 4; //per PMT


//time spread
//static const G4double sig_Elect = 10.0; //ps
//static const G4double sig_T0    = 15.0; //ps
//static const G4double sig_TTS   = 35.0; //ps

}

#endif
