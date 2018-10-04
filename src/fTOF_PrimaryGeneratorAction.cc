//my
#include "fTOF_PrimaryGeneratorAction.hh"
#include "fTOF_VolumeStructures.hh"
#include "backgroundGen.hh"

//G4
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

//root
#include "TMath.h"

//my
#include "crtConst.hh"
#include "muonGen.hh"

fTOF_PrimaryGeneratorAction::fTOF_PrimaryGeneratorAction() : // @suppress("Class members should be properly initialized")
  _particleGun(0),
  _particleName("pi+"),
  _particleMomentum(180.*GeV),
  _PhiAngle(0.0*deg),
  _ThetaAngle(0.0*deg),
  _singlePhoton(false)
{
  _particleGun = new G4ParticleGun(1);  
  _BunchXID = 0;

  //backGen = new backgroundGen();
  
  _muonGen = new muonGen(123112);

}

fTOF_PrimaryGeneratorAction::~fTOF_PrimaryGeneratorAction()
{
  delete _particleGun;
  //delete backGen;
}

void fTOF_PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle;
  // Correct for center of bar
  G4double xInit, yInit, zInit;
  G4double dX, dY, dZ;
  G4double Ekin, m;
  G4int pdgID;
  G4int i;

  /*
  ////////////////////////
  //_particleName = "kaon-";
  //_particleMomentum = 1500.0*MeV;
  _particleMomentum = 1.0*MeV;
  _ThetaAngle = 180*deg;
  _PhiAngle = 0.0*deg;
  xInit = 0.0*cm;
  yInit = 0.0*cm;
  //zInit = (crtConst::hodo_hight+0.5)*cm;    
  zInit = (83.5 + 38.1 + 1.49 +  1.5 + 1.64/2 + 2.5)*cm;
  //zInit = (83.5 + 38.1 + 1.49 +  1.5 + 1.64/2 + 20.0)*cm;
  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
		      + m*m) - m);
  dX = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  dY = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  dZ = std::cos(_ThetaAngle);
  //dX = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  //dZ = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  //dY = std::cos(_ThetaAngle);
  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  _particleGun->GeneratePrimaryVertex(anEvent);
  */

  ////////////////////////
  _particleName = "pi+";
  //_particleMomentum = 1500.0*MeV;
  _particleMomentum = 180.*GeV;           ////////////////////////////////// particle momentum
  //_particleMomentum = 400.0*MeV;
  _muonGen->GenerateMuon((crtConst::hodo_hight+0.5));
  //_muonGen->GenerateMuon((83.0 + 38.1 + 30.0));
  //_muonGen->GenerateMuon((84.0 + 40 + 10));
  _ThetaAngle = _muonGen->GetTheta()*deg;
  _PhiAngle = _muonGen->GetPhi()*deg;
  //xInit = _muonGen->GetX()*cm;
  //yInit = _muonGen->GetY()*cm;
  //zInit = _muonGen->GetZ()*cm;
  xInit = -999.*mm;
  yInit = -999.*mm;
  while (xInit < -6*mm || xInit > 6*mm || yInit < -6*mm || yInit > 6*mm){
	  xInit = CLHEP::RandGauss::shoot(0., 1.6*mm);
	  yInit = CLHEP::RandGauss::shoot(0., 2.1*mm);
  }


  zInit = 130*cm;             /////////////////////////////////////////////////////////////////
  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
		      + m*m) - m);
  ////dX = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  ////dY = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  ////dZ = std::cos(_ThetaAngle);
  //dX = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  //dZ = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  //dY = std::cos(_ThetaAngle);
  dX =	0;
  dZ =  -130.;
  dY =  0;
  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  _particleGun->GeneratePrimaryVertex(anEvent);


  /*
  ////////////////////////
  //_particleName = "kaon-";
  _particleMomentum = 100.0*MeV;
  _ThetaAngle = 180.0*deg;
  _PhiAngle = 180.0*deg;
  xInit = 0.0*cm;
  yInit = 0.0*cm;
  zInit = (crtConst::hodo_hight)*cm;    
  //xInit = 0.0*cm;
  //yInit = 0.0*cm;
  //zInit = 0.0*cm;    
  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
		      + m*m) - m);
  //dX = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  //dY = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  //dZ = std::cos(_ThetaAngle);
  dX = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  dZ = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  dY = std::cos(_ThetaAngle);
  //G4ThreeVector dir(dX, dY, dZ);
  G4ThreeVector dir(0.0, 0.0, -1.0);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  _particleGun->GeneratePrimaryVertex(anEvent);
  */

  /*
  ////////////////////////
  //_particleName = "kaon-";
  _particleMomentum = 3000.0*MeV;
  _ThetaAngle = 20.0*deg;
  _PhiAngle = 2.0*TMath::Pi()*(G4UniformRand())*rad;
  xInit = 0.0*cm;
  yInit = 0.0*cm;
  zInit = 0.0*cm;    
  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  pdgID = particle->GetPDGEncoding();
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
		      + m*m) - m);
  dX = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  dY = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  dZ = std::cos(_ThetaAngle);
  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  _particleGun->GeneratePrimaryVertex(anEvent);
  */

  /*
  ////////////////////////
  //_particleName = "kaon-";
  G4double  momentumMin = 0.700*GeV;
  G4double  momentumMax = 0.701*GeV;
  G4double  thetaAngleMin = 15.0*deg;
  G4double  thetaAngleMax = 25.0*deg;
  G4double  phiAngleMin = 0.0*deg;
  G4double  phiAngleMax = 360.0*deg;  
  _particleMomentum = (momentumMax - momentumMin)*G4UniformRand() + momentumMin;
  _ThetaAngle = (thetaAngleMax - thetaAngleMin)*G4UniformRand() + thetaAngleMin;
  _PhiAngle = (phiAngleMax - phiAngleMin)*G4UniformRand() + phiAngleMin;
  xInit = 0.0*cm;
  yInit = 0.0*cm;
  zInit = 0.0*cm;    
  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  pdgID = particle->GetPDGEncoding();
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
		      + m*m) - m);
  dX = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  dY = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  dZ = std::cos(_ThetaAngle);
  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  _particleGun->GeneratePrimaryVertex(anEvent);
  */

  /*
  //const G4int nMomArr   = 4;
  //const G4int nThetaArr = 3;
  //const G4int nPhiArr   = 3;
  //G4double momArr[nMomArr]     = {0.7,0.8,0.9,2.0};
  //G4double thetaArr[nThetaArr] = {17.0,20.0,24.0};
  //G4double phiArr[nPhiArr]     = {0.0,10.0,20.0};

  const G4int nMomArr   = 10;
  const G4int nThetaArr = 10;
  const G4int nPhiArr   = 10;
  G4double momArr[nMomArr]     = {0.7, 0.8, 0.9, 1.0, 1.5,
				  2.0, 2.5, 3.0, 3.5, 4.0};
  G4double thetaArr[nThetaArr] = {16.0, 17.0, 18.0, 19.0, 20.0, 
				  21.0, 22.0, 23.0, 24.0, 25.0};
  G4double phiArr[nPhiArr]     = {0.0,  3.0,  6.0,  9.0, 12.0,
				  15.0, 18.0, 21.0, 24.0, 27.0};
  
  G4int iRandMom   = 0;
  G4int iRandTheta = 0;
  G4int iRandPhi   = 0;
  ////////////////////////
  iRandMom   = GenFlatInt(0,(nMomArr-1));
  iRandTheta = GenFlatInt(0,(nThetaArr-1));
  iRandPhi   = GenFlatInt(0,(nPhiArr-1));
  //G4cout<<"iRandMom   = "<<iRandMom<<G4endl;
  //G4cout<<"iRandTheta = "<<iRandTheta<<G4endl;
  //G4cout<<"iRandPhi   = "<<iRandPhi<<G4endl;
  //_particleName = "kaon-";
  _particleMomentum = momArr[iRandMom]*GeV;
  _ThetaAngle = thetaArr[iRandTheta]*deg;
  _PhiAngle = phiArr[iRandPhi]*deg;
  xInit = 0.0*cm;
  yInit = 0.0*cm;
  zInit = 0.0*cm;    
  ///////////////////////
  _BunchXID++;
  particle = particleTable->FindParticle(_particleName);
  m = particle->GetPDGMass();
  Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
		      + m*m) - m);
  dX = std::sin(_ThetaAngle)*std::cos(_PhiAngle);
  dY = std::sin(_ThetaAngle)*std::sin(_PhiAngle);
  dZ = std::cos(_ThetaAngle);
  G4ThreeVector dir(dX, dY, dZ);
  _particleGun->SetParticleDefinition(particle);
  _particleGun->SetParticleMomentumDirection(dir);
  _particleGun->SetParticleEnergy(Ekin);  
  _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
  _particleGun->GeneratePrimaryVertex(anEvent);
  */

  /*
  for(i= 0;i<1000;i++){
    backGen->Loop();
    ////////////////////////
    if(backGen->pID==22){
      _particleName = "gamma";
    }
    if(backGen->pID==11)
      _particleName = "e-";
    if(backGen->pID==-11)
      _particleName = "e+";
    _particleMomentum = TMath::Sqrt(backGen->MomX*backGen->MomX + 
				    backGen->MomY*backGen->MomY +
				    backGen->MomZ*backGen->MomZ)*GeV;
    xInit = backGen->PosX*cm;
    yInit = backGen->PosY*cm;
    zInit = backGen->PosZ*cm;
    _BunchXID = backGen->baunchXID;
    ///////////////////////
    
    particle = particleTable->FindParticle(_particleName);
    pdgID = particle->GetPDGEncoding();
    m = particle->GetPDGMass();
    Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum 
			+ m*m) - m);
    if(_particleMomentum>0.0){
      dX = backGen->MomX/_particleMomentum;
      dY = backGen->MomY/_particleMomentum;
      dZ = backGen->MomZ/_particleMomentum;
    }
    G4ThreeVector dir(dX, dY, dZ);
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(dir);
    _particleGun->SetParticleEnergy(Ekin);  
    _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
    _particleGun->GeneratePrimaryVertex(anEvent);
  }
  */

  /*
  G4cout<<"--------------------------------------"<<G4endl
	<<" _particleMomentum "<<_particleMomentum<<G4endl
	<<" _BunchXID         "<<_BunchXID<<G4endl
	<<" _particleName     "<<_particleName<<G4endl
	<<"       dX          "<<dX<<G4endl
	<<"       dY          "<<dY<<G4endl
	<<"       dZ          "<<dZ<<G4endl 
	<<"     xInit         "<<xInit<<G4endl
	<<"     yInit         "<<yInit<<G4endl
	<<"     zInit         "<<zInit<<G4endl
	<<" _ThetaAngle       "<<_ThetaAngle<<G4endl
	<<" _PhiAngle         "<<_PhiAngle<<G4endl
	<<"       deg         "<<deg<<G4endl
	<<"++++++++++++++++++++++++++++++++++++++"<<G4endl;
  */

  /*
  G4cout<<"Momentum   = "<<_particleMomentum<<G4endl
	<<"ThetaAngle = "<<_ThetaAngle<<G4endl
	<<"PhiAngle   = "<<_PhiAngle<<G4endl
	<<"pdgID      = "<<pdgID<<G4endl;
  */

}

G4int fTOF_PrimaryGeneratorAction::GenFlatInt(G4int iMin,G4int iMax){
  G4int val;
  val = (G4int)floor((iMax - iMin + 1)*G4UniformRand() + iMin);

  return val;
}

