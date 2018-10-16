#include "fTOF_SteppingAction.hh"
#include "fTOF_SensitiveDetector.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SDManager.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "fTOFConst.hh"


fTOF_SteppingAction::fTOF_SteppingAction(fTOF_PrimaryGeneratorAction* 
                                         genAction) :
    _genAction(genAction)
{
    Reset();
    ResetPerEvent();
    //nKillPhot = 0;
}

fTOF_SteppingAction::~fTOF_SteppingAction()
{ }

void fTOF_SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    G4Track* aTrack = aStep->GetTrack();
    // G4Track *aTrack = aStep->GetTrack();
    G4int trackID = aTrack->GetTrackID();
    //G4cout<<"trackID = "<<trackID<<G4endl;

    G4StepPoint* aPrePoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* aPrePV = aPrePoint->GetPhysicalVolume();
    G4StepPoint* aPostPoint = aStep->GetPostStepPoint();
    G4VPhysicalVolume* aPostPV = aPostPoint->GetPhysicalVolume();

    if (!aPostPV) return;

    //
    //----- initial Track () -----
    //
    // G4cout<<aPrePV->GetName()<<G4endl;


    _trkMomX = aStep->GetPostStepPoint()->GetMomentum().getX();
    _trkMomY = aStep->GetPostStepPoint()->GetMomentum().getY();
    _trkMomZ = aStep->GetPostStepPoint()->GetMomentum().getZ();
    _trkPosX = aStep->GetPostStepPoint()->GetPosition().getX();
    _trkPosY = aStep->GetPostStepPoint()->GetPosition().getY();
    _trkPosZ = aStep->GetPostStepPoint()->GetPosition().getZ();
    _trkT = aStep->GetPostStepPoint()->GetGlobalTime();


    // if (aPrePV->GetName().contains("mixer") &&
    if ((aPrePV->GetName().contains("quartzBar") || aPrePV->GetName().contains("quartzBar1")) &&
            aPostPV->GetName().contains("rightLogical") ) {
        G4String sdName = "fTOF";
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        fTOF_SensitiveDetector* sd =
                (fTOF_SensitiveDetector*)SDman->FindSensitiveDetector(sdName);

        HitData hitInfo;
        hitInfo.SecID = _SecID;
        hitInfo.trkMomX = _trkMomX;
        hitInfo.trkMomY = _trkMomY;
        hitInfo.trkMomZ = _trkMomZ;
        hitInfo.trkPosX = _trkPosX/mm;
        hitInfo.trkPosY = _trkPosY/mm;
        hitInfo.trkPosZ = _trkPosZ/mm;
        hitInfo.trkT = _trkT/ps;
        hitInfo.trkLength = _trkLength;
        hitInfo.chID = _chID;

        hitInfo.trkSideID = 1;

        sd->ProcessHits_fTOF(aStep, NULL, hitInfo);
        return;
        // kill
         aTrack->SetTrackStatus(fStopAndKill);
    }


    if (aPrePV->GetName().contains("World") &&
            aPostPV->GetName().contains("quartzBar"/*"strAbs"*/) && aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() == fTOFConst::particleName) {

        //    std::cout << "              ENTERING TIME               " << _trkT << std::endl;
        G4String sdName = "fTOF";
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        fTOF_SensitiveDetector* sd =
                (fTOF_SensitiveDetector*)SDman->FindSensitiveDetector(sdName);
        HitData hitInfo;
        hitInfo.entMomX = _trkMomX;
        hitInfo.entMomY = _trkMomY;
        hitInfo.entMomZ = _trkMomZ;
        hitInfo.entPosX = _trkPosX/mm;
        hitInfo.entPosY = _trkPosY/mm;
        hitInfo.entPosZ = _trkPosZ/mm;
        hitInfo.entTime = _trkT/ps;
        hitInfo.trkT = _trkT/ps;
        hitInfo.trkSideID = -500.;
        sd->ProcessHits_fTOF(aStep, NULL, hitInfo);

    }
    if (aPostPV->GetName().contains("World") &&
            aPrePV->GetName().contains("quartzBar") && aTrack->GetDynamicParticle()->GetDefinition()->GetParticleName() == fTOFConst::particleName) {  /////////////////////////////////////////////////////////////////////////////////

        // std::cout << "              ESCAPING TIME               " << _trkT << std::endl;

    }


    //
    //----- Optical photons  -----
    //
    G4ParticleDefinition* particleType = aTrack->GetDefinition();
    if (particleType != G4OpticalPhoton::OpticalPhotonDefinition())
        return;

    if (_particleID != trackID) {
        Reset();
        _particleID = trackID;
        InternalReflectionProbability(aTrack->GetTotalEnergy()/eV,
                                      _probOfReflection);
    }

    //LB
    //G4cout<<"aPrePV->GetName() = "<<aPrePV->GetName()<<G4endl
    //<<"aPostPV->GetName()= "<<aPostPV->GetName()<<G4endl;
    if (aPrePV->GetName().contains("quartzBar") &&            ////////////////////////////////////////////////////////////////////////////////
            aPostPV->GetName().contains("hamChannel")) {
        std::string nameCh = aPostPV->GetName();
        std::string numberCh = nameCh.substr(nameCh.find("_hamChannel_pv_"), 2);
        if (numberCh.find("_") != std::string::npos)
            numberCh.resize(1);
        std::istringstream stream(numberCh);
        stream >> _chID;


        //G4cout<<" _chID = "<<_chID<<G4endl;
        //G4cout<<"aPrePV->GetName()  = "<<aPrePV->GetName()<<G4endl
        //<<"aPostPV->GetName() = "<<aPostPV->GetName()<<G4endl
        //<<"_SecID             = "<<_SecID<<G4endl
        //<<"_particleID        = "<<_particleID<<G4endl
        //<<"trackID            = "<<trackID<<G4endl;
        ////_PosX = aPostPoint->GetPosition().x();
        ////_PosY = aPostPoint->GetPosition().y();
        ////_PosZ = aPostPoint->GetPosition().z();
        // Use PrePoint for direction.  PostPoint is in new media and will
        //// include refraction.
        ////const G4ThreeVector momDir = aPrePoint->GetMomentumDirection();
        ////_barThetaX = std::atan2(momDir.x(), momDir.z());
        ////_barThetaY = std::atan2(momDir.y(), momDir.z());
        ////_timeLeftBar = aPostPoint->GetGlobalTime();
    }

    G4OpBoundaryProcessStatus boundaryStatus = Undefined;
    static G4OpBoundaryProcess* boundary = NULL;

    // Find boundary process
    if (!boundary) {
        G4ProcessManager* pm =
                aStep->GetTrack()->GetDefinition()->GetProcessManager();
        G4int nprocesses = pm->GetProcessListLength();
        G4ProcessVector* pv = pm->GetProcessList();
        for (G4int i = 0; i < nprocesses; i++) {
            if ((*pv)[i]->GetProcessName() == "OpBoundary") {
                boundary = (G4OpBoundaryProcess*)(*pv)[i];
                break;
            }
        }
    }

    if (!boundary) return;

    boundaryStatus = boundary->GetStatus();

    if (aPostPoint->GetStepStatus() == fGeomBoundary) {
        G4String sdName = "fTOF";
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        fTOF_SensitiveDetector* sd =
                (fTOF_SensitiveDetector*)SDman->FindSensitiveDetector(sdName);
        G4double flat = G4UniformRand();
        switch(boundaryStatus) {
        case Absorption:
            break;
        case Detection:
            if (sd) {
                // HitData hitInfo;
                // hitInfo.SecID = _SecID;
                // hitInfo.trkMomX = _trkMomX;
                // hitInfo.trkMomY = _trkMomY;
                // hitInfo.trkMomZ = _trkMomZ;
                // hitInfo.trkPosX = _trkPosX;
                // hitInfo.trkPosY = _trkPosY;
                // hitInfo.trkPosZ = _trkPosZ;
                // hitInfo.trkT = _trkT;
                // hitInfo.trkLength = _trkLength;
                // hitInfo.chID = _chID;

                //  if (_trkPosX >  0) _trkSideID = 1;
                //  else if (_trkPosX < 0) _trkSideID = -1;


                //  hitInfo.trkNSideRefl = _trkNSideRefl;
                //  hitInfo.trkSideID = _trkSideID;


                //hitInfo.detection = true;
                //hitInfo.SecID = -3;

                //_trkMomX = -999.0;
                //_trkMomY = -999.0;
                //_trkMomZ = -999.0;
                //_trkPosX = -999.0;
                //_trkPosY = -999.0;
                //_trkPosZ = -999.0;
                //_trkT    = -999.0;
                //_trkLength = -999.0;

                //if(trackID==20)
                //G4cout<<"SecID = "<<_SecID<<G4endl;
                //if(_SecID<0 || _SecID>12)
                //hitInfo.nMirror1 = _nBounceMirror1;
                //hitInfo.nMirror2 = _nBounceMirror2;
                //hitInfo.nEndMirror = _nBounceEndMirror;
                //hitInfo.nWedgeSide = _nWedgeSide;
                //hitInfo.nWedgeTop = _nWedgeTop;
                //hitInfo.nWedgeBottom = _nWedgeBottom;
                //hitInfo.nFBlockSide = _nFBlockSide;
                //hitInfo.BarNum = _bar;
                //hitInfo.BarX = _barX;
                //hitInfo.BarY = _barY;
                //hitInfo.BarZ = _barZ;
                //hitInfo.BarThetaX = _barThetaX;
                //hitInfo.BarThetaY = _barThetaY;
                //hitInfo.TimeLeftBar = _timeLeftBar;




                // sd->ProcessHits_fTOF(aStep, NULL, hitInfo);
            }
            break;
        case FresnelReflection:
            // Reflections of surfaces of different media
            break;
        case TotalInternalReflection:
            // Add reflection probability...
            if(aTrack->GetTrackLength()>20000.0){
                G4Track* aNonConstTrack = const_cast<G4Track*>(aTrack);
                aNonConstTrack->SetTrackStatus(fStopAndKill);
                //nKillPhot++;
                //G4cout<<" fTOF_SteppingAction::UserSteppingAction;  aTrack->GetTrackLength()>20000.0 mm "<<G4endl;
            }
            if (flat > _probOfReflection) {
                G4Track* aNonConstTrack = const_cast<G4Track*>(aTrack);
                aNonConstTrack->SetTrackStatus(fStopAndKill);
            } // else if (_trkPosX > 9.99*mm || _trkPosX < -9.99*mm) _trkNSideRefl++;
            //if (aPrePoint->GetMaterial()->GetName() == "quartz"){
            //_numberOfBounces++;
            //}
            break;
        case SpikeReflection:
            break;
        default:
            break;

        }

    }
}

void fTOF_SteppingAction::ResetPerEvent(){
    _SecID   = -1;
    _trkMomX = -999.0;
    _trkMomY = -999.0;
    _trkMomZ = -999.0;
    _trkPosX = -999.0;
    _trkPosY = -999.0;
    _trkPosZ = -999.0;
    _trkT    = -999.0;
    _trkLength = -999.0;
    _trkSideID = -999.0;
    _entTime = -999.0;

    _entMomX = -999.0;
    _entMomY = -999.0;
    _entMomZ = -999.0;
    _entPosX = -999.0;
    _entPosY = -999.0;
    _entPosZ = -999.0;
}

void fTOF_SteppingAction::Reset()
{
    _particleID = -1;
    _probOfReflection = 1.0;
    _chID = -1;
    //_totalPath = 0.;
    //_barPath = 0.;
    //_sobPath = 0.;
    //_nBounceMirror1 = 0;
    //_nBounceMirror2 = 0;
    //_nBounceEndMirror = 0;
    //_nWedgeSide = 0;
    //_nWedgeTop = 0;
    //_nWedgeBottom = 0;
    //_nFBlockSide = 0;
    //_barFirst = 0;
    //_bar = -1;
    //_barX = -9999.;
    //_barY = -9999.;
    //_barZ = -9999.;
    //_barThetaX = -9999.;
    //_barThetaY = -9999.;
    //_timeLeftBar = -9999.;
    //_probOfReflection = 1.;
    //_numberOfBounces = 0;
}

void
fTOF_SteppingAction::InternalReflectionProbability(G4double energy,
                                                   G4double& probability)
{
    probability = 1.0;
    if (_genAction->SinglePhotonGenerator()) return;

    /* this function simulate the internal reflection probability per one
     bounce - each time photon bounced this function is called
     and the reflection is tested if photon reflects or disappear -
     this function estimates loss of photons during imperfections
     of bar */

    G4double opticalPhotonEnergy[36]={
        1.90744901082931,1.93725290162352,1.96800294768103,1.99974493070815,
        2.03252763449025,2.06640309506508,2.10142687633737,2.13765837420526,
        2.17516115270009,2.21400331614116,2.25425792188918,2.29600343896121,
        2.33932425856425,2.38431126353664,2.43106246478245,2.4796837140781,
        2.53028950416133,2.58300386883136,2.63796139795543,2.6953083848675,
        2.75520412675345,2.81782240236148,2.88335315590477,2.95200442152155,
        3.02400452936354,3.09960464259763,3.17908168471551,3.26274172905013,
        3.35092393794338,3.44400515844181,3.54240530582586,3.64659369717368,
        3.75709653648197,3.87450580324703,3.99948986141629,4.13280619013017};

    G4double internalReflectivity[36] =
    {0.999895281,0.999891334,0.999885743,0.999878696,0.999870426,
     0.9998612,0.999851309,0.999841055,0.999830735,0.999820635,0.999811012,
     0.999802084,0.999794018,0.999786917,0.999780807,0.999775625,
     0.999771209,0.999767282,0.999763441,0.999759146,0.999753706,
     0.999746266,0.999735798,0.999721084,0.999700708,0.99967304,
     0.999636227,0.999588178,0.999526552,0.999448747,0.999351887,
     0.999232808,0.99908805,0.998913839,0.998706078,0.998460335};

    G4int i;
    for(i = 0; i < 36;i++) {
        if(energy < opticalPhotonEnergy[i]) break;
    }


    probability = ((energy-opticalPhotonEnergy[i-1])/
            (opticalPhotonEnergy[i]-opticalPhotonEnergy[i-1]))*
            (internalReflectivity[i]-internalReflectivity[i-1]) +
            internalReflectivity[i-1];

    /* because the ratio between peak1 and peak2 did not correspond,
     the reflection probability was change to get the same
     ration 2.1:1 => the original probability is multiplied by .9992 */
    probability = probability*.9992;


    // probability = 0;
}
