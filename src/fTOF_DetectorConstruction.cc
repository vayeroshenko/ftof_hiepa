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

    G4Material* Air =
            new G4Material("Air", density = 0.000290*mg/cm3, ncomponents = 2);
    Air->AddElement(N, fractionmass = 0.7);
    Air->AddElement(O, fractionmass = 0.3);

    // Aluminum
    Aluminum =
            new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
    Aluminum->AddElement(Al, fractionmass = 1.0);


    //
    // Generate and Add Material Properties Table
    //
    const G4int num = 36;
    G4double WaveLength[num];
    G4double Absorption[num];      // Default value for absorption
    G4double AirAbsorption[num];
    G4double AirRefractiveIndex[num];
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

}

G4VPhysicalVolume* fTOF_DetectorConstruction::Construct()
{

    MagneticField* magField = new MagneticField();
    G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager() ->GetFieldManager();
    fieldMgr->SetDetectorField(magField,0);
    fieldMgr->CreateChordFinder(magField);

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


    //
    // MCP PMT Windows  /////////////////// 11.04.18 ////////////////////////////////////////////////
    //

    hamWin.solid = new G4Box("hamWindow",
                             hamWin.sizeX/2.0,
                             hamWin.sizeY/2.0,
                             hamWin.sizeZ/2.0);
    hamWin.logical = new G4LogicalVolume(hamWin.solid,
                                         hamWin.material,"hamWindow");

    planWin.solid = new G4Box("planWindow",
                              planWin.sizeX/2.0,
                              planWin.sizeY/2.0,
                              planWin.sizeZ/2.0);
    planWin.logical = new G4LogicalVolume(planWin.solid,
                                          planWin.material,"planWindow");

    ///////////////////////////////////////////////////////////////////////////////////////////////////


    //
    // MCP PMT Channels
    //


    hamChan.solid = new G4Box("hamChannel",
                              hamChan.sizeX/2.0,
                              hamChan.sizeY/2.0,
                              hamChan.sizeZ/2.0);
    hamChan.logical = new G4LogicalVolume(hamChan.solid,
                                          hamChan.material,"hamChannel");

    planChan.solid = new G4Box("hamChannel",
                               planChan.sizeX/2.0,
                               planChan.sizeY/2.0,
                               planChan.sizeZ/2.0);
    planChan.logical = new G4LogicalVolume(planChan.solid,
                                           planChan.material,"planChannel");


    //
    // MCP PMT Boxes
    //


    hamBox.solid = new G4Box("hamBox",
                             hamBox.sizeX/2.0,
                             hamBox.sizeY/2.0,
                             hamBox.sizeZ/2.0);
    hamBox.logical = new G4LogicalVolume(hamBox.solid,
                                         hamBox.material,"hamBox");

    planBox.solid = new G4Box("planBox",
                              planBox.sizeX/2.0,
                              planBox.sizeY/2.0,
                              planBox.sizeZ/2.0);
    planBox.logical = new G4LogicalVolume(planBox.solid,
                                          planBox.material,"planBox");




    G4RotationMatrix Ra = G4RotationMatrix();
    G4ThreeVector Ta = G4ThreeVector();
    G4Transform3D Tr;






    G4Trd *trapeze = new G4Trd(
                "trapeze",
                sector.shortSide/2.,
                sector.longSide/2.,
                sector.thickness/2.,
                sector.thickness/2.,
                sector.height/2.
                );

    G4Trd *trapeze1 = new G4Trd(
                "trapeze2",
                sector.shortSide/2.,
                sector.longSide/2.,
                fTOFConst::layerThickness/2.,
                fTOFConst::layerThickness/2.,
                sector.height/2.
                );


    G4Trd *absTrd = new G4Trd(
                "absorber",
                abs.shortSide/2.*1.1,
                abs.longSide/2.*1.1,
                abs.thickness/2. + fTOFConst::layerDist/2,
                abs.thickness/2.+ fTOFConst::layerDist/2,
                abs.height/2.
                );


    G4Trd *mixerTrd = new G4Trd(
                "mixer",
                mixer.shortSide/2.,
                mixer.longSide/2.,
                mixer.thickness/2.,
                mixer.thickness/2.,
                mixer.height/2.
                );

    G4Trd *trapAbs = new G4Trd(
                "trapeze",
                sector.shortSide/2.,
                sector.longSide/2.,
                10*um,
                10*um,
                sector.height/2.
                );


    // G4Trd *mixerTrdFront = new G4Trd(
    // "mixerFront",
    // mixer.shortSide/2.,
    // mixer.longSide/2.,
    // mixer.thickness/2.,
    // mixer.thickness/2.,
    // mixer.height/2.
    // );


    // G4Trd *mixerTrdFront = new G4Trd(
    // "mixerBack",
    // mixer.shortSide/2.,
    // mixer.longSide/2.,
    // mixer.thickness/2.,
    // mixer.thickness/2.,
    // mixer.height/2.
    // );



    G4LogicalVolume *fullBarLog = new G4LogicalVolume(trapeze,bigBox.material,"quartzBar");

    G4LogicalVolume *strAbs = new G4LogicalVolume(trapAbs,Aluminum,"strAbs");

    G4LogicalVolume *absLayer = new G4LogicalVolume(trapeze1,world.material,"absLayer");

    G4LogicalVolume *absorber = new G4LogicalVolume(absTrd,bigBox.material,"absorber");

    G4LogicalVolume *fullBarLog1 = new G4LogicalVolume(trapeze1,bigBox.material,"quartzBar1");

    G4LogicalVolume *mixerLog = new G4LogicalVolume(mixerTrd,bigBox.material,"mixer");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    G4Box *leftVol = new G4Box("leftVol",
                               sector.shortSide/2.,
                               sector.thickness/2.,
                               5.*mm
                               );

    G4LogicalVolume *leftLogical = new G4LogicalVolume(leftVol,
                                                       bigBox.material,"leftLogical");


    G4Box *rightVol = new G4Box("rightVol",
                                sector.longSide/2.,
                                sector.thickness/2.,
                                5*mm
                                );

    G4LogicalVolume *rightLogical = new G4LogicalVolume(rightVol,
                                                        bigBox.material,"rightLogical");




    //-------------------------------------------------------

    G4AssemblyVolume* secAssembly = new G4AssemblyVolume();

    //--------------------------------------------------------

    std::cout << "Short side = " << fTOFConst::innerSide/cm << " cm" << std::endl
              << "Long side = "  << fTOFConst::outerSide/cm << " cm" << std::endl;

    // one sector to be done

    G4RotationMatrix RTilt = G4RotationMatrix();
    RTilt.rotateX(fTOFConst::angle);

    for (int l = 0; l < fTOFConst::nLayers; ++l ){
        for (int j = 0; j < fTOFConst::nSec; ++j) {
            G4double i = j + l * 1./fTOFConst::nLayers;
            G4double dist = (sector.thickness + fTOFConst::layerDist + fTOFConst::layerThickness) * l;

            /////////// sector /////////////

            if (j < fTOFConst::nDrawSec){
                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();

                Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);

                Ra = Ra * RTilt;

                Ta.setX(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg));
                Ta.setY(dist);
                Ta.setZ(fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(fullBarLog,Tr);


                ////////////// Absorber for straight photons  //////////////////

                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();

                Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);
                Ta.setX(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg));
                Ta.setY((sector.thickness) * (l+0.5));
                Ta.setZ(fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta -= G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta = (Ra * (RTilt * (Ra.inverse() * Ta)));

                Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));


                Ra = Ra * RTilt;

                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(strAbs,Tr);

                ////////////// Absorber for straight photons  //////////////////

                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();

                Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);
                Ta.setX(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg));
                Ta.setY((sector.thickness) * (l-0.5));
                Ta.setZ(fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta -= G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta = (Ra * (RTilt * (Ra.inverse() * Ta)));

                Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));


                Ra = Ra * RTilt;

                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(strAbs,Tr);

                ////////////////////////////////////////////////////


            }


            ////////// absorber /////////////

            if (j < fTOFConst::nDrawSec || j == fTOFConst::nSec - 1) {
                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();

                Ra.rotateY(- 360./fTOFConst::nSec*(i+0.5) *deg + 90.*deg);

                Ra = Ra * RTilt;
                Ta.setX((fTOFConst::innerRad + fTOFConst::outerRad)/2.
                        * TMath::Cos(360./fTOFConst::nSec*(i+0.5) *deg));
                Ta.setY(dist);
                Ta.setZ((fTOFConst::innerRad + fTOFConst::outerRad)/2.
                        * TMath::Sin(360./fTOFConst::nSec*(i+0.5) *deg));

                Ta -= G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta = (Ra * (RTilt * (Ra.inverse() * Ta)));

                Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(absorber,Tr);
            }

            ///////// inner abs /////////////

            if (j < fTOFConst::nDrawSec){
                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();
                Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);



                Ta.setX((fTOFConst::innerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) - 5*mm) * TMath::Cos(360./fTOFConst::nSec*i *deg));
                Ta.setY(dist);
                Ta.setZ((fTOFConst::innerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) - 5*mm) * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta -= G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta = (Ra * (RTilt * (Ra.inverse() * Ta)));

                Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ra = Ra * RTilt;
                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(leftLogical,Tr);
            }

            /////////// outer detector ///////

            if (j < fTOFConst::nDrawSec){
                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();

                Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);


                Ta.setX((fTOFConst::outerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) + 5*mm) * TMath::Cos(360./fTOFConst::nSec*i *deg));
                Ta.setY(dist);
                Ta.setZ((fTOFConst::outerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) + 5*mm) * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta -= G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta = (Ra * (RTilt * (Ra.inverse() * Ta)));

                Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));


                Ra = Ra * RTilt;
                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(rightLogical,Tr);
            }

            ///////// layer absorber ////////

            if (l != fTOFConst::nLayers - 1){
                Ta = G4ThreeVector(0.,0.,0.);
                Ra = G4RotationMatrix();

                Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);
                Ta.setX(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg));
                Ta.setY((sector.thickness + fTOFConst::layerDist) * (l+0.5));
                Ta.setZ(fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta -= G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

                Ta = (Ra * (RTilt * (Ra.inverse() * Ta)));

                Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
                                    dist,
                                    fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));


                Ra = Ra * RTilt;

                Tr = G4Transform3D(Ra,Ta);
                secAssembly->AddPlacedVolume(absLayer,Tr);
            }

            //////////////////////////////////
        }
    }



    // for (int j = 0; j < fTOFConst::nSec; ++j) { G4double i = j + 0.5;
    //   /////////// sector /////////////


    //   Ta = G4ThreeVector(0., sector.thickness + fTOFConst::layerDist + fTOFConst::layerThickness ,0.);
    //   Ra = G4RotationMatrix();

    //   Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);
    //   Ta.setX(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg));
    //   Ta.setZ(fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));
    //   Tr = G4Transform3D(Ra,Ta);
    //   secAssembly->AddPlacedVolume(fullBarLog,Tr);

    //   ////////// absorber /////////////
    //   Ta = G4ThreeVector(0.,sector.thickness + fTOFConst::layerDist+ fTOFConst::layerThickness,0.);
    //   Ra = G4RotationMatrix();

    //   Ra.rotateY(- 360./fTOFConst::nSec*(i+0.5) *deg + 90.*deg);
    //   Ta.setX((fTOFConst::innerRad + fTOFConst::outerRad)/2.
    //           * TMath::Cos(360./fTOFConst::nSec*(i+0.5) *deg));
    //   Ta.setZ((fTOFConst::innerRad + fTOFConst::outerRad)/2.
    //           * TMath::Sin(360./fTOFConst::nSec*(i+0.5) *deg));
    //   Tr = G4Transform3D(Ra,Ta);
    //   secAssembly->AddPlacedVolume(absorber,Tr);

    //   ///////// inner abs /////////////
    //   Ta = G4ThreeVector(0.,sector.thickness + fTOFConst::layerDist+ fTOFConst::layerThickness,0.);
    //   Ra = G4RotationMatrix();

    //   Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);
    //   Ta.setX((fTOFConst::innerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) - 5*mm) * TMath::Cos(360./fTOFConst::nSec*i *deg));
    //   Ta.setZ((fTOFConst::innerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) - 5*mm) * TMath::Sin(360./fTOFConst::nSec*i *deg));
    //   Tr = G4Transform3D(Ra,Ta);
    //   secAssembly->AddPlacedVolume(leftLogical,Tr);

    //   /////////// outer detector ///////
    //   Ta = G4ThreeVector(0.,sector.thickness + fTOFConst::layerDist+ fTOFConst::layerThickness,0.);
    //   Ra = G4RotationMatrix();

    //   Ra.rotateY(- 360./fTOFConst::nSec*i *deg + 90.*deg);
    //   Ta.setX((fTOFConst::outerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) - 5*mm) * TMath::Cos(360./fTOFConst::nSec*i *deg));
    //   Ta.setZ((fTOFConst::outerRad*TMath::Cos(TMath::Pi() / fTOFConst::nSec) - 5*mm) * TMath::Sin(360./fTOFConst::nSec*i *deg));
    //   Tr = G4Transform3D(Ra,Ta);
    //   secAssembly->AddPlacedVolume(rightLogical,Tr);
    // }




    Ta = G4ThreeVector(0.,0.,0.);
    // Ta.setY(sector.thickness/2. + 0.5*mm);
    Ra = G4RotationMatrix();
    Tr = G4Transform3D(Ra,Ta);
    // secAssembly->AddPlacedVolume(fullBarLog,Tr);

    Ta = G4ThreeVector(0.,0.,0.);
    Ta.setY(-sector.thickness/2.-0.5*mm);
    Ta.setZ(0.);
    Ra = G4RotationMatrix();
    Tr = G4Transform3D(Ra,Ta);
    // secAssembly->AddPlacedVolume(fullBarLog1,Tr);

    ///////////////////////////////////////////////////////////////////////////////////////

    // mcp pmt windows 11.04.18 /////////////////////////////////////////////////

    Ta = G4ThreeVector(bigBox.sizeX/2.+hamWin.sizeX/2.,
                       0.,
                       0.);
    Ra = G4RotationMatrix();
    Tr = G4Transform3D(Ra,Ta);
    // secAssembly->AddPlacedVolume(hamWin.logical,Tr);


    Ta = G4ThreeVector(-bigBox.sizeX/2.-planWin.sizeX/2.,
                       0.,
                       0.);
    Ra = G4RotationMatrix();
    Tr = G4Transform3D(Ra,Ta);
    // secAssembly->AddPlacedVolume(planWin.logical,Tr);



    // MCP boxes

    Ta = G4ThreeVector(bigBox.sizeX/2. + hamWin.sizeX + hamBox.sizeX/2. + hamChan.sizeX,
                       0.,
                       0.);
    Ra = G4RotationMatrix();
    Tr = G4Transform3D(Ra,Ta);
    // secAssembly->AddPlacedVolume(hamBox.logical,Tr);

    Ta = G4ThreeVector(-bigBox.sizeX/2. - planWin.sizeX - planBox.sizeX/2. - planChan.sizeX,
                       0.,
                       0.);
    Ra = G4RotationMatrix();
    Tr = G4Transform3D(Ra,Ta);
    // secAssembly->AddPlacedVolume(planBox.logical,Tr);

    //

    G4int i = 0;
    G4int j = 0;
    Ra = G4RotationMatrix();




    //////////////////////////// kill volumes ////////////////////////////////






    // Ta.setZ(-sector.height/2.  - 4.999*mm);
    // Ta.setX(0.);
    // Ta.setY(sector.thickness/2. + 0.5*mm);
    // Tr = G4Transform3D(Ra, Ta);
    // secAssembly->AddPlacedVolume(leftLogical, Tr);








    G4LogicalVolume *absLog = new G4LogicalVolume(absTrd,
                                                  bigBox.material,"absTrd");
    Ra = G4RotationMatrix();
    Ta.setY(0.);
    Ta.setX(0.);
    Ta.setZ(0.);
    Tr = G4Transform3D(Ra, Ta);
    // secAssembly->AddPlacedVolume(absLog, Tr);





    Ra = G4RotationMatrix();
    Ta.setZ(sector.height/2. + 4.999*mm);
    Ta.setY(sector.thickness/2. + 0.5*mm);
    Ta.setX(0.);
    Tr = G4Transform3D(Ra, Ta);

    // secAssembly->AddPlacedVolume(rightLogical, Tr);








    // G4Box *rightVol = new G4Box("rightVol",
    //   mixer.longSide/2.,
    //   mixer.thickness/2.,
    //   5*mm
    //   );

    // Ra = G4RotationMatrix();
    // Ta.setZ(sector.height/2. + mixer.height + 4.999*mm);
    // Ta.setY(0.);
    // Ta.setX(0.);
    // Tr = G4Transform3D(Ra, Ta);
    // G4LogicalVolume *rightLogical = new G4LogicalVolume(rightVol,
    //       bigBox.material,"rightLogical");
    // secAssembly->AddPlacedVolume(rightLogical, Tr);



    Ra = G4RotationMatrix();
    Ta.setZ(sector.height/2. + mixer.height/2.);
    Ta.setY(0.);
    Ta.setX(0.);
    Tr = G4Transform3D(Ra, Ta);
    // secAssembly->AddPlacedVolume(mixerLog, Tr);






    Ra = G4RotationMatrix();

    G4Box *frontVol = new G4Box("frontVol",
                                sector.sides/2.,
                                // sector.thickness + 1*mm,
                                sector.thickness/2.,
                                5*mm
                                );
    Ra = G4RotationMatrix();
    Ta.setX(sector.middleLine/2.+4.999*mm);
    Ta.setZ(0.);
    Ta.setY(0.);
    Ra.rotateY(90*deg + sector.angle);
    Tr = G4Transform3D(Ra, Ta);

    G4LogicalVolume *frontLogical = new G4LogicalVolume(frontVol,
                                                        bigBox.material,"frontLogical");
    // secAssembly->AddPlacedVolume(frontLogical, Tr);


    G4Box *backVol = new G4Box("backVol",
                               sector.sides/2.,
                               // sector.thickness + 1*mm,
                               sector.thickness/2.,
                               5.*mm
                               );
    Ra = G4RotationMatrix();
    Ta.setZ(0.);
    Ta.setY(0.);
    Ta.setX(-sector.middleLine/2.-4.999*mm);
    Ra.rotateY(-90*deg-sector.angle);
    Tr = G4Transform3D(Ra, Ta);
    G4LogicalVolume *backLogical = new G4LogicalVolume(backVol,
                                                       bigBox.material,"backLogical");
    // secAssembly->AddPlacedVolume(backLogical, Tr);













    //
    //make Imprint
    //

    Ra = G4RotationMatrix();
    Ra.rotateY(270.0*deg);
    Ra.rotateX(90.0*deg/*+fTOFConst::angle*/);

    Ta.setX(0.);
    Ta.setY(0.);
    Ta.setZ(0.);

//    Ta += G4ThreeVector(fTOFConst::centerRad * TMath::Cos(360./fTOFConst::nSec*i *deg),
//                        0,
//                        fTOFConst::centerRad * TMath::Sin(360./fTOFConst::nSec*i *deg));

    Ta.rotateZ(270*deg);

    Tr = G4Transform3D(Ra, Ta);
    secAssembly->MakeImprint(world.logical, Tr, 0, true);



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



    // fTOF vis attributes 11.04.18 //////////////////

    world.logical->SetVisAttributes(worldVisAtt);

    fullBarLog->SetVisAttributes(quartzVisAtt);
    fullBarLog1->SetVisAttributes(quartzVisAtt);
    hamWin.logical->SetVisAttributes(quartzVisAtt);
    planWin.logical->SetVisAttributes(quartzVisAtt);

    hamChan.logical->SetVisAttributes(sensitiveVisAtt);
    planChan.logical->SetVisAttributes(sensitiveVisAtt);

    hamBox.logical->SetVisAttributes(pmtboxVisAtt);
    planBox.logical->SetVisAttributes(pmtboxVisAtt);


    leftLogical->SetVisAttributes(quartzVisAtt);
    rightLogical->SetVisAttributes(quartzVisAtt);

    frontLogical->SetVisAttributes(quartzVisAtt);
    backLogical->SetVisAttributes(quartzVisAtt);

    //////////////////////////////////////////////////



    //
    // Define Optical Borders
    //

    // Surface for killing photons at borders
    const G4int num1 = 2;
    G4double Ephoton[num1] = {1.5*eV, 5.8*eV};

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



    // G4OpticalSurface* ReflSurface = new G4OpticalSurface("ReflectiveSurface");
    // G4double ReflectivityRefl[num1] = {1., 1.};
    // G4double EfficiencyRefl[num1] = {0., 0.};

    // ReflSurface->SetType(dielectric_metal);
    // ReflSurface->SetFinish(polished);
    // ReflSurface->SetModel(glisur);

    // G4MaterialPropertiesTable* VolumeRefl = new G4MaterialPropertiesTable();
    // VolumeRefl->AddProperty("REFLECTIVITY", Ephoton, ReflectivityRefl, num1);
    // VolumeRefl->AddProperty("EFFICIENCY",   Ephoton, EfficiencyRefl,   num1);
    // ReflSurface->SetMaterialPropertiesTable(VolumeRefl);

    // new G4LogicalSkinSurface("SensitiveSurfaceLeft",
    //        mixerLog, ReflSurface);






    new G4LogicalSkinSurface("SensitiveSurfaceLeft",
                             planChan.logical, OpVolumeKillSurface);
    new G4LogicalSkinSurface("SensitiveSurfaceRight",
                             hamChan.logical, OpVolumeKillSurface);

    new G4LogicalSkinSurface("SensitiveSurfaceLeft",
                             leftLogical, OpVolumeKillSurface);
    new G4LogicalSkinSurface("SensitiveSurfaceRight",
                             rightLogical, OpVolumeKillSurface);

    new G4LogicalSkinSurface("SensitiveSurfaceRight",
                             strAbs, OpVolumeKillSurface);


    new G4LogicalSkinSurface("AbsTrdSurface",
                             absLog, OpVolumeKillSurface);

    new G4LogicalSkinSurface("AbsTrdSurface",
                             absorber, OpVolumeKillSurface);

    new G4LogicalSkinSurface("AbsTrdSurface",
                             absLayer, OpVolumeKillSurface);


    new G4LogicalSkinSurface("SensitiveSurfaceLeft",
                             frontLogical, OpVolumeKillSurface);
    new G4LogicalSkinSurface("SensitiveSurfaceRight",
                             backLogical, OpVolumeKillSurface);



    //
    // Sensitive detector definition
    //
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    fTOF_SensitiveDetector* aSD = new fTOF_SensitiveDetector("fTOF");
    SDman->AddNewDetector(aSD);
    //    hamChan.logical->SetSensitiveDetector(aSD);
    //    planChan.logical->SetSensitiveDetector(aSD);

    //    leftLogical->SetSensitiveDetector(aSD);
    //    rightLogical->SetSensitiveDetector(aSD);
    // frontLogical->SetSensitiveDetector(aSD);
    // backLogical->SetSensitiveDetector(aSD);





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
