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
//#include "crtConst.hh"
#include "muonGen.hh"

fTOF_PrimaryGeneratorAction::fTOF_PrimaryGeneratorAction() : // @suppress("Class members should be properly initialized")
    _particleGun(0),
    _particleName("pi+"),
    _particleMomentum(0.5*GeV),
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
    G4int nParticles = 1;
    //	G4String particleName = "kaon+";
    G4String particleName = fTOFConst::particleName;
    G4double particleMomentum = fTOFConst::particleMomentum;
    //	G4double particleMomentum = 700 * MeV;

//    G4double thetaMax = TMath::Pi()/2 - TMath::ATan((fTOFConst::innerRad - 30*cm)/ (130. * cm));
//    G4double thetaMin = TMath::Pi()/2 - TMath::ATan((fTOFConst::outerRad + 30*cm)/ (130. * cm));
//    G4double phiMin = - TMath::ATan(2*fTOFConst::outerSide / (130. * cm));
//    G4double phiMax = - phiMin;




//    G4double thetaMax = 0*rad;
//    G4double thetaMin = pi * rad;
//    G4double phiMin = 0;
//    G4double phiMax = twopi*rad;

    G4double thetaMin = pi*rad / 2 - TMath::ATan((fTOFConst::centerRad +9*cm)/ (130.*cm));
    G4double thetaMax = thetaMin;
    G4double phiMin = fTOFConst::dPhiPrim;
    G4double phiMax = phiMin;

    phiMin += fTOFConst::dPhiPrim;
    phiMax += fTOFConst::dPhiPrim;


    Scan(anEvent, nParticles,
         particleName, particleMomentum,
         thetaMin, thetaMax,
         phiMin, phiMax);

}

G4int fTOF_PrimaryGeneratorAction::GenFlatInt(G4int iMin,G4int iMax){
    G4int val;
    val = (G4int)floor((iMax - iMin + 1)*G4UniformRand() + iMin);

    return val;
}

void fTOF_PrimaryGeneratorAction::Scan(	G4Event* anEvent, G4int nEv,
                                        G4String particleName, G4double particleMomentum,
                                        G4double thetaMin, G4double thetaMax,
                                        G4double phiMin, G4double phiMax)
{
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle;
    G4double xInit, yInit, zInit;
    G4double dX, dY, dZ;
    G4double Ekin, m;

    _particleName = particleName;
    _particleMomentum = particleMomentum;

    xInit = 0.0*mm;
    yInit = 0.0*cm;
    zInit = 130*cm;

    particle = particleTable->FindParticle(_particleName);
    m = particle->GetPDGMass();
    Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum
                        + m*m) - m);

    for (int i = 0; i < nEv; ++i){
        _BunchXID++;
        G4double phi = phiMin + (phiMax - phiMin) * G4UniformRand();
        G4double cosTheta = TMath::Cos(thetaMin) + (TMath::Cos(thetaMax) - TMath::Cos(thetaMin))* G4UniformRand();

        dX =   TMath::Sin(phi) * cosTheta;				// horizontal plane, normal to beam axis,
        dZ = - TMath::Sqrt(1 - cosTheta*cosTheta);		// horizontal plane, beam axis, opposite to beam
        dY = - TMath::Cos(phi) * cosTheta;				// vertical axis, upwards

        G4ThreeVector dir(dX, dY, dZ);
        _particleGun->SetParticleDefinition(particle);
        _particleGun->SetParticleMomentumDirection(dir);
        _particleGun->SetParticleEnergy(Ekin);
        _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
        _particleGun->GeneratePrimaryVertex(anEvent);
    }
}

