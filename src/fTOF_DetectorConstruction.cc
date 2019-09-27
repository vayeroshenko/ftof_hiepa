//My
#include "fTOF_DetectorConstruction.hh"
#include "fTOF_SensitiveDetector.hh"
#include "MagneticField.hh"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"

//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
// #include <G4GDMLParser.hh>

//root 
#include "TMath.h"

#include "crtConst.hh"


fTOF_DetectorConstruction::fTOF_DetectorConstruction()
{

    worldVisAtt = new G4VisAttributes();
    quartzVisAtt = new G4VisAttributes();
    sensitiveVisAtt = new G4VisAttributes();
    pmtboxVisAtt = new G4VisAttributes();

    // Define Materials to be used
    DefineMaterials();
}

fTOF_DetectorConstruction::~fTOF_DetectorConstruction()
{


    delete worldVisAtt;
    delete quartzVisAtt;
    delete sensitiveVisAtt;
    delete pmtboxVisAtt;

}

void fTOF_DetectorConstruction::DefineMaterials()
{
    G4String symbol;
    G4double a, z, density;
    G4int ncomponents, natoms;
    G4double fractionmass;


    G4Element* C =
            new G4Element("Carbon", symbol = "C", z = 6., a = 12.01*g/mole);
    G4Element* H =
            new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01*g/mole);
    G4Element* N =
            new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01*g/mole);
    G4Element* O =
            new G4Element("Oxygen", symbol = "O", z = 8., a = 16.00*g/mole);
    G4Element* Si =
            new G4Element("Silicon", symbol = "Si", z = 14., a = 28.09*g/mole);
    G4Element* Al =
            new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98*g/mole);

    // Quartz Material (SiO2)
    G4Material* SiO2 =
            new G4Material("quartz", density = 2.200*g/cm3, ncomponents = 2);
    SiO2->AddElement(Si, natoms = 1);
    SiO2->AddElement(O , natoms = 2);

    Air =
            new G4Material("Air", density = 0.000290*mg/cm3, ncomponents = 2);
    Air->AddElement(N, fractionmass = 0.7);
    Air->AddElement(O, fractionmass = 0.3);

    Bis_MSB = new G4Material("Bis_MSB",density=1.07*g/cm3,ncomponents=2);
    Bis_MSB->AddElement(H,natoms=22);
    Bis_MSB->AddElement(C,natoms=24);


    // Aluminum
    Aluminum =
            new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
    Aluminum->AddElement(Al, fractionmass = 1.0);

    // ------------- Materials -------------
    G4NistManager* nist = G4NistManager::Instance();
    Water = nist->FindOrBuildMaterial("G4_WATER");
    Glass = nist->FindOrBuildMaterial("G4_GLASS_PLATE");


    //
    // Generate and Add Material Properties Table
    //
    const G4int num = 36;
    G4double WaveLength[num];
    G4double Absorption[num];      // Default value for absorption
    G4double AirAbsorption[num];
    G4double AirRefractiveIndex[num];
    G4double WaterRefractiveIndex[num];
    G4double PhotonEnergy[num];

    // Absorption of quartz per 1m
    G4double QuartzAbsorption[num] =
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,
     0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,
     0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,
     0.990610945};

    for (int i=0; i<num; i++) {
        WaveLength[i] = (300 + i*10)*nanometer;
        Absorption[i] = 100*m;      // Fake number for no absorption
        AirAbsorption[i] = 10.*cm;   // If photon hits air, kill it
        AirRefractiveIndex[i] = 1.;
        WaterRefractiveIndex[i] = 1.3;
        PhotonEnergy[num - (i+1)] = twopi*hbarc/WaveLength[i];
        /* Absorption is given per length and G4 needs mean free path
       length, calculate it here
       mean free path length - taken as probablility equal 1/e
       that the photon will be absorbed */
        QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
        //EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*
        //epotekBarJoint.thickness;
    }

    G4double QuartzRefractiveIndex[num] =
    {1.456535,1.456812,1.4571  ,1.457399,1.457712,1.458038,
     1.458378,1.458735,1.459108,1.4595  ,1.459911,1.460344,
     1.460799,1.46128 ,1.461789,1.462326,1.462897,1.463502,
     1.464146,1.464833,1.465566,1.46635 ,1.46719 ,1.468094,
     1.469066,1.470116,1.471252,1.472485,1.473826,1.475289,
     1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

    G4MaterialPropertiesTable *WaterMPT = new G4MaterialPropertiesTable();
    WaterMPT->AddProperty("RINDEX", PhotonEnergy, WaterRefractiveIndex, num);
    Water->SetMaterialPropertiesTable(WaterMPT);


    // Wavelength shifter
    G4MaterialPropertiesTable *MPTWLS = new G4MaterialPropertiesTable();

    // BIS absorption
    G4double waveLength3[7] = {300., 340., 380., 400., 420., 460., 500};
    G4double photonEnergy5[7];
    for(int i=0; i<7; i++)
        photonEnergy5[i] = 1240./waveLength3[i]*eV;
    G4double absLen3[7] = {10*nm, 10*nm, 10*nm, 1*mm, 200*m, 200*m, 200*m};
    MPTWLS->AddProperty("WLSABSLENGTH", photonEnergy5, absLen3, 7);

    G4double ppckovEmit[8] = { 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.95 * eV, 2.6401*eV , 3.0402*eV , 3.5403*eV , 3.8404*eV};
    G4double rindexWLS[8] = { 1.5, 1.5, 1.5, 1.5, 1.504 , 1.505 , 1.515 , 1.52 };

    // BIS reemission
    G4double waveLength4[16] = {380., 390., 400., 410., 420., 430., 440.,
                                450., 460., 470., 480., 490., 500., 510., 520., 530};
    G4double photonEnergy6[16];
    for(int i=0; i<16; i++)
        photonEnergy6[i] = 1240./waveLength4[i]*eV;
    G4double reEmit4[16] = {0., 0., 0.1, 0.8, 1., 0.8, 0.5,
                            0.45, 0.3, 0.2, 0.15, 0.1, 0.05, 0.05, 0.05, 0.};
    MPTWLS->AddProperty("WLSCOMPONENT", photonEnergy6, reEmit4, 16);

    MPTWLS->AddConstProperty("WLSTIMECONSTANT", 3. * ns);
    MPTWLS-> AddConstProperty("WLSMEANNUMBERPHOTONS",1.0);
    MPTWLS->AddProperty("RINDEX", ppckovEmit, rindexWLS, 8)->SetSpline(true);
    Bis_MSB->SetMaterialPropertiesTable(MPTWLS);

    // Wavelength shifter end

    // Assign absorption and refraction to materials

    // Quartz
    G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
    QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
    QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

    SiO2->SetMaterialPropertiesTable(QuartzMPT);

    // fTOF bar 10.04.18
    bigBox.material = SiO2;


    hamWin.material = SiO2;
    planWin.material = SiO2;

    hamChan.material = Aluminum;
    planChan.material = Aluminum;

    hamBox.material = Aluminum;
    planBox.material = Aluminum;


    // Air
    G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
    AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
    AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);

    // Assign these properties to the world volume
    Air->SetMaterialPropertiesTable(AirMPT);
    world.material = Air;

    Glass->SetMaterialPropertiesTable(QuartzMPT);

}

G4VPhysicalVolume* fTOF_DetectorConstruction::Construct()
{



    //
    // Define World Volume
    //
    world.solid = new G4Box("World",
                            world.sizeX/2,
                            world.sizeY/2,
                            world.sizeZ/2);

    world.logical = new G4LogicalVolume(world.solid,
                                        world.material,
                                        "World");

    world.physical = new G4PVPlacement(0,
                                       G4ThreeVector(),
                                       world.logical,
                                       "World",
                                       0,
                                       false,
                                       0);



    G4RotationMatrix *Rot = new G4RotationMatrix();
    G4ThreeVector *Trans = new G4ThreeVector(0,0,130.*cm);


    G4VSolid* airTube = new G4Tubs("airTube",
                                   0.,
                                   fTOFConst::airTubeD/2.,
                                   fTOFConst::boxXYsize/2. - fTOFConst::airTubeD,
                                   0,
                                   twopi);

    G4VSolid* ckBoxInit = new G4Box("ckBox",
                                    fTOFConst::boxXYsize/2.,
                                    fTOFConst::boxXYsize/2.,
                                    fTOFConst::boxZsize/2.);

    G4VSolid* mirrorSolid = new G4Box("mirror",
                                      fTOFConst::boxXYsize/2. - fTOFConst::airTubeD/4.,
                                      fTOFConst::boxXYsize/2. - fTOFConst::airTubeD/4.,
                                      1.*mm);
    G4LogicalVolume* mirrorLog = new G4LogicalVolume(mirrorSolid,
                                                     Aluminum,
                                                     "mirror");

    G4VSolid* detectorSolid = new G4Tubs("detector", 0, fTOFConst::airTubeD/2., 1.*mm, 0, twopi);


    G4LogicalVolume* detectorLog = new G4LogicalVolume(detectorSolid,
                                                       Aluminum,
                                                       "detector");

    // Substraction
    // 1
    G4RotationMatrix *subsRot = new G4RotationMatrix();
    subsRot->rotateX(90.*deg);
    G4ThreeVector *subsTrans = new G4ThreeVector(fTOFConst::boxXYsize/2. - fTOFConst::airTubeD/2.,
                                                 0,
                                                 0.);

    G4VSolid *ckBoxSolid = new G4SubtractionSolid("ckBox",
                                                  ckBoxInit,
                                                  airTube,
                                                  subsRot,
                                                  *subsTrans);
    G4LogicalVolume* airTubeLogical1 = new G4LogicalVolume(airTube,
                                                           Air,
                                                           "airTube");
    G4VPhysicalVolume* airTubePhys1 = new G4PVPlacement(subsRot,
                                                        *subsTrans + G4ThreeVector(0,0,130.*cm),
                                                        airTubeLogical1,
                                                        "airTube",
                                                        world.logical,
                                                        0., false, 0);


    // 2
    subsRot = new G4RotationMatrix();
    subsRot->rotateX(90.*deg);
    subsTrans = new G4ThreeVector(- fTOFConst::boxXYsize/2. + fTOFConst::airTubeD/2.,
                                  0,
                                  0.);

    ckBoxSolid = new G4SubtractionSolid("ckBox",
                                        ckBoxSolid,
                                        airTube,
                                        subsRot,
                                        *subsTrans);
    G4LogicalVolume* airTubeLogical2 = new G4LogicalVolume(airTube,
                                                           Air,
                                                           "airTube");
    G4VPhysicalVolume* airTubePhys2 = new G4PVPlacement(subsRot,
                                                        *subsTrans + G4ThreeVector(0,0,130.*cm),
                                                        airTubeLogical2,
                                                        "airTube",
                                                        world.logical,
                                                        0., false, 0);

    // 3
    subsRot = new G4RotationMatrix();
    subsRot->rotateX(90.*deg);
    subsRot->rotateY(90.*deg);
    subsTrans = new G4ThreeVector(0,
                                  fTOFConst::boxXYsize/2. - fTOFConst::airTubeD/2.,
                                  0.);

    ckBoxSolid = new G4SubtractionSolid("ckBox",
                                        ckBoxSolid,
                                        airTube,
                                        subsRot,
                                        *subsTrans);
    G4LogicalVolume* airTubeLogical3 = new G4LogicalVolume(airTube,
                                                           Air,
                                                           "airTube");
    G4VPhysicalVolume* airTubePhys3 = new G4PVPlacement(subsRot,
                                                        *subsTrans + G4ThreeVector(0,0,130.*cm),
                                                        airTubeLogical3,
                                                        "airTube",
                                                        world.logical,
                                                        0., false, 0);

    // 4
    subsRot = new G4RotationMatrix();
    subsRot->rotateX(90.*deg);
    subsRot->rotateY(90.*deg);
    subsTrans = new G4ThreeVector(0,
                                  - fTOFConst::boxXYsize/2. + fTOFConst::airTubeD/2.,
                                  0.);

    ckBoxSolid = new G4SubtractionSolid("ckBox",
                                        ckBoxSolid,
                                        airTube,
                                        subsRot,
                                        *subsTrans);
    G4LogicalVolume* airTubeLogical4 = new G4LogicalVolume(airTube,
                                                           Air,
                                                           "airTube");
    G4VPhysicalVolume* airTubePhys4 = new G4PVPlacement(subsRot,
                                                        *subsTrans + G4ThreeVector(0,0,130.*cm),
                                                        airTubeLogical4,
                                                        "airTube",
                                                        world.logical,
                                                        0., false, 0);

    // Substraction end

    ///////////////////////////////////////////////////////////////////
    G4LogicalVolume* ckBoxLogical = new G4LogicalVolume(ckBoxSolid,
                                                        Water,
                                                        "ckBox");

    G4VPhysicalVolume* ckBoxPhys = new G4PVPlacement(Rot,
                                                     *Trans,
                                                     ckBoxLogical,
                                                     "ckBox",
                                                     world.logical,
                                                     0., false, 0);

    G4VPhysicalVolume* mirrorPhys1 = new G4PVPlacement(Rot,
                                                       *Trans - G4ThreeVector(0,0,fTOFConst::boxZsize/2. + 1.*mm),
                                                       mirrorLog,
                                                       "mirror",
                                                       world.logical,
                                                       0., false, 0);
    G4VPhysicalVolume* mirrorPhys2 = new G4PVPlacement(Rot,
                                                       *Trans + G4ThreeVector(0,0,fTOFConst::boxZsize/2. + 1.*mm),
                                                       mirrorLog,
                                                       "mirror",
                                                       world.logical,
                                                       0., false, 0);
    ////////////////////////////////////////////////////////////////////

    G4VSolid* glassTube = new G4Tubs("airTube",
                                     fTOFConst::glassTubeInD/2.,
                                     fTOFConst::glassTubeD/2.,
                                     fTOFConst::boxXYsize/2. - fTOFConst::airTubeD - 1.*mm,
                                     0,
                                     twopi);

    G4LogicalVolume* glassLogical = new G4LogicalVolume(glassTube,
                                                        Glass,
                                                        "glassTube");

    G4VSolid* wlsInTube = new G4Tubs("wlsInTube",
                                     fTOFConst::glassTubeInD/2. - fTOFConst::WLSskinThikness,
                                     fTOFConst::glassTubeD/2.,
                                     fTOFConst::boxXYsize/2. - fTOFConst::airTubeD - 1.*mm,
                                     0,
                                     twopi);
    G4LogicalVolume* wlsInLogical = new G4LogicalVolume(wlsInTube,
                                                        Bis_MSB,
                                                        "wlsIn");
    G4VSolid* wlsOutTube = new G4Tubs("wlsOutTube",
                                      fTOFConst::glassTubeD/2.,
                                      fTOFConst::glassTubeD/2.+ fTOFConst::WLSskinThikness,
                                      fTOFConst::boxXYsize/2. - fTOFConst::airTubeD - 1.*mm,
                                      0,
                                      twopi);


    G4LogicalVolume* wlsOutLogical = new G4LogicalVolume(wlsOutTube,
                                                         Bis_MSB,
                                                         "wlsOut");

    // glass tube placement
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      glassLogical,
                      "glassTube",
                      airTubeLogical1,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      glassLogical,
                      "glassTube",
                      airTubeLogical2,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      glassLogical,
                      "glassTube",
                      airTubeLogical3,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      glassLogical,
                      "glassTube",
                      airTubeLogical4,
                      0., false, 0);
    // WLS in placement
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsInLogical,
                      "wlsIn",
                      airTubeLogical1,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsInLogical,
                      "wlsIn",
                      airTubeLogical2,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsInLogical,
                      "wlsIn",
                      airTubeLogical3,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsInLogical,
                      "wlsIn",
                      airTubeLogical4,
                      0., false, 0);
    // WLS out placement
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsOutLogical,
                      "wlsOut",
                      airTubeLogical1,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsOutLogical,
                      "wlsOut",
                      airTubeLogical2,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsOutLogical,
                      "wlsOut",
                      airTubeLogical3,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(),
                      wlsOutLogical,
                      "wlsOut",
                      airTubeLogical4,
                      0., false, 0);

    // detector placements
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,fTOFConst::boxXYsize/2. - fTOFConst::airTubeD - 1.*mm),
                      detectorLog,
                      "detector 1",
                      airTubeLogical1,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,-fTOFConst::boxXYsize/2. + fTOFConst::airTubeD+1.*mm),
                      detectorLog,
                      "detector -1",
                      airTubeLogical1,
                      0., false, 0);

    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,fTOFConst::boxXYsize/2. - fTOFConst::airTubeD-1.*mm),
                      detectorLog,
                      "detector 2",
                      airTubeLogical2,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,-fTOFConst::boxXYsize/2. + fTOFConst::airTubeD+1.*mm),
                      detectorLog,
                      "detector -2",
                      airTubeLogical2,
                      0., false, 0);

    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,fTOFConst::boxXYsize/2. - fTOFConst::airTubeD-1.*mm),
                      detectorLog,
                      "detector 3",
                      airTubeLogical3,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,-fTOFConst::boxXYsize/2. + fTOFConst::airTubeD+1.*mm),
                      detectorLog,
                      "detector -3",
                      airTubeLogical3,
                      0., false, 0);

    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,fTOFConst::boxXYsize/2. - fTOFConst::airTubeD-1.*mm),
                      detectorLog,
                      "detector 4",
                      airTubeLogical4,
                      0., false, 0);
    new G4PVPlacement(new G4RotationMatrix,
                      G4ThreeVector(0,0,-fTOFConst::boxXYsize/2. + fTOFConst::airTubeD+1.*mm),
                      detectorLog,
                      "detector -4",
                      airTubeLogical4,
                      0., false, 0);






    //////////////////////////// kill volumes ////////////////////////////////








    //
    //make Imprint
    //

    //    Ra = G4RotationMatrix();
    //    Ra.rotateX(90.0*deg);
    //    Ra.rotateY(45.*deg);
    //    Ta.setX(0.);
    //    Ta.setY(0.);
    //    Ta.setZ(0.);

    //    Tr = G4Transform3D(Ra, Ta);
    //    secAssembly->MakeImprint(world.logical, Tr, 0, true);



    //-----------------------------------------------------

    //
    // Set Visualization Attributes
    //
    G4Color blue        = G4Color(0., 0., 1.);
    G4Color green       = G4Color(0., 1., 0.);
    G4Color red         = G4Color(1., 0., 0.);
    G4Color white       = G4Color(1., 1., 1.);
    G4Color cyan        = G4Color(0., 1., 1.);
    G4Color DircColor   = G4Color(0.0, 0.0, 1.0, 0.2);
    G4Color SensColor   = G4Color(0.0, 1.0, 1.0, 0.1);

    worldVisAtt->SetColor(white);
    worldVisAtt->SetVisibility(true);
    quartzVisAtt->SetColor(DircColor);
    quartzVisAtt->SetVisibility(true);
    sensitiveVisAtt->SetColor(SensColor);
    sensitiveVisAtt->SetVisibility(true);
    pmtboxVisAtt->SetColor(red);
    pmtboxVisAtt->SetVisibility(true);

    glassLogical->SetVisAttributes(quartzVisAtt);

    wlsOutLogical->SetVisAttributes(sensitiveVisAtt);
    wlsInLogical->SetVisAttributes(sensitiveVisAtt);

    //    mirrorLog->SetVisAttributes(sensitiveVisAtt);

    // fTOF vis attributes 11.04.18 //////////////////

    //    world.logical->SetVisAttributes(worldVisAtt);

    ////    fullBarLog->SetVisAttributes(quartzVisAtt);
    //    hamWin.logical->SetVisAttributes(quartzVisAtt);
    //    planWin.logical->SetVisAttributes(quartzVisAtt);

    //    hamChan.logical->SetVisAttributes(sensitiveVisAtt);
    //    planChan.logical->SetVisAttributes(sensitiveVisAtt);

    //    hamBox.logical->SetVisAttributes(pmtboxVisAtt);
    //    planBox.logical->SetVisAttributes(pmtboxVisAtt);


    //    leftLogical->SetVisAttributes(quartzVisAtt);
    //    rightLogical->SetVisAttributes(quartzVisAtt);

    //    frontLogical->SetVisAttributes(quartzVisAtt);
    //    backLogical->SetVisAttributes(quartzVisAtt);

    //////////////////////////////////////////////////



    //
    // Define Optical Borders
    //

    // Surface for killing photons at borders
    const G4int num1 = 2;
    G4double Ephoton[num1] = {1.*eV, 10.*eV};

        G4OpticalSurface* OpVolumeKillSurface =
                new G4OpticalSurface("VolumeKillSurface");
        OpVolumeKillSurface->SetType(dielectric_metal);
        OpVolumeKillSurface->SetFinish(polished);
        OpVolumeKillSurface->SetModel(glisur);

        G4double ReflectivityKill[num1] = {0., 0.};
        G4double EfficiencyKill[num1] = {1., 1.};
        G4MaterialPropertiesTable* VolumeKill = new G4MaterialPropertiesTable();
        VolumeKill->AddProperty("REFLECTIVITY", Ephoton, ReflectivityKill, num1);
        VolumeKill->AddProperty("EFFICIENCY",   Ephoton, EfficiencyKill,   num1);
        OpVolumeKillSurface->SetMaterialPropertiesTable(VolumeKill);



    G4OpticalSurface* ReflSurface = new G4OpticalSurface("ReflectiveSurface");
    G4double ReflectivityRefl[num1] = {0.95, 0.95};
    G4double EfficiencyRefl[num1] = {0., 0.};

    ReflSurface->SetType(dielectric_metal);
    ReflSurface->SetFinish(polished);
    ReflSurface->SetModel(glisur);

    G4MaterialPropertiesTable* VolumeRefl = new G4MaterialPropertiesTable();
    VolumeRefl->AddProperty("REFLECTIVITY", Ephoton, ReflectivityRefl, num1);
    VolumeRefl->AddProperty("EFFICIENCY",   Ephoton, EfficiencyRefl,   num1);
    ReflSurface->SetMaterialPropertiesTable(VolumeRefl);

    new G4LogicalSkinSurface("mirrorSurface",
                             mirrorLog, ReflSurface);

    new G4LogicalSkinSurface("detectorSurface",
                             detectorLog, OpVolumeKillSurface);




    //    new G4LogicalSkinSurface("SensitiveSurfaceLeft",
    //                             planChan.logical, OpVolumeKillSurface);
    //    new G4LogicalSkinSurface("SensitiveSurfaceRight",
    //                             hamChan.logical, OpVolumeKillSurface);

    //    //  new G4LogicalSkinSurface("SensitiveSurfaceLeft",
    //    //         leftLogical, OpVolumeKillSurface);
    //    new G4LogicalSkinSurface("SensitiveSurfaceRight",
    //                             rightLogical, OpVolumeKillSurface);
    //    new G4LogicalSkinSurface("SensitiveSurfaceLeft",
    //                             planWin.logical, OpVolumeKillSurface);
    //    new G4LogicalSkinSurface("SensitiveSurfaceRight",
    //                             hamWin.logical, OpVolumeKillSurface);



    //    new G4LogicalSkinSurface("AbsTrdSurface",
    //                             absLog, OpVolumeKillSurface);

    //    new G4LogicalSkinSurface("AbsTrdSurface",
    //                             absorber, OpVolumeKillSurface);

    //    new G4LogicalSkinSurface("AbsTrdSurface",
    //                             absLayer, OpVolumeKillSurface);


    //    new G4LogicalSkinSurface("SensitiveSurfaceLeft",
    //                             frontLogical, OpVolumeKillSurface);
    //    new G4LogicalSkinSurface("SensitiveSurfaceRight",
    //                             backLogical, OpVolumeKillSurface);



    //
    // Sensitive detector definition
    //
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    fTOF_SensitiveDetector* aSD = new fTOF_SensitiveDetector("fTOF");
    SDman->AddNewDetector(aSD);
    //    hamChan.logical->SetSensitiveDetector(aSD);
    //    planChan.logical->SetSensitiveDetector(aSD);

    detectorLog->SetSensitiveDetector(aSD);

    //  leftLogical->SetSensitiveDetector(aSD);
    //  rightLogical->SetSensitiveDetector(aSD);
    //   frontLogical->SetSensitiveDetector(aSD);
    //   backLogical->SetSensitiveDetector(aSD);





    G4cout << *(G4Material::GetMaterialTable()) << G4endl;

    printDetectorParameters();


    return world.physical;
}

void fTOF_DetectorConstruction::printDetectorParameters(){

    //fTOFConst::detTiltAngle
    //fTOFConst::det_Rmin
    //fTOFConst::det_Rmax
    //fTOFConst::det_Zmin
    //fTOFConst::det_Zmax
    //fTOFConst::N_det
    //fTOFConst::detAngleSize
    //fTOFConst::detSizeXmin
    //fTOFConst::detSizeXmax
    //fTOFConst::detlength

}
