// ====================================================================
//    exe_ANNRI.cc
//    ver.0 written by T. Yano, based on ex02.cc
//
//    Ref:
//    http://www-geant4.kek.jp/g4users/g4tut07/exercise-2.html
//    http://www-geant4.kek.jp/g4users/g4tut07/exercise-3.html
// ====================================================================
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#ifdef G4UI_USE_TCSH
#include "G4UItcsh.hh"
#endif
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Physics_PhysicsList.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "TrackingAction.hh"
#include "DrawMessage.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#include "G4TrajectoryDrawByParticleID.hh"

#include "fstream"
std::ofstream outf;

G4int NumGamma;
G4double GammaEnergies[15];

//===========================================================
//////////////////////// Define Structure ///////////////////
////=========================================================
////############################//############################
////++++++++++2013B0025+++++++++//+++++++++2014B0124++++++++++    
////Beam Pipe -> Al             //Beam Pipe -> Al+LiF
////best target position        //best target position
////(15.0*mm, 2.0*mm,0.0*mm)    //(1.7*mm, 3.5*mm, -3.9*mm)
////best source position        //best source position
////(15.0*mm,11.5*mm,0.0*mm)    //60Co(-1.75,3.875,-2.375)
////############################//############################
////+++++++++++++++++++++Coordinate System++++++++++++++++++++
////(POSX,POSY,POZ) = (upstream, up, left side (see on beam))
////beam from      => (30,0,0) 
////beam direction => (-1,0,0)
////upper side (0,+,0) => crystal number 1 to 7
////lower side (0,-,0) => crystal number 8 to 14
////left side  (0,0,+) => crystal number 1,2,3,11,12,13
////right side (0,0,-) => crystal number 4,5,6,8,9,10
////+++++++++++++++++++++++++++++++++++++++++

//Which experiment? 
//1=2013B0025(natural Gd), 2=2014B0124(enriched Gd) or 2015(position calculation experiment)
G4int BEAM_DUCT = 2;

//Where does the target put on?
//(POSX,POSY,POZ) = (upstream, up, left side (see on beam))
G4double SOURCE_POSX =1.7;    // -1.75;  // 1.7; 
G4double SOURCE_POSY =0.;     // 3.875;  // 3.5; 	
G4double SOURCE_POSZ =0.;   	// -2.375; //-3.9;      	  

//Where does the target holder put on?
//(POSX,POSY,POZ) = (upstream, up, left side (see on beam))
G4double HOLDER_POSX = 0.;// -1.75; 	// 1.7; 
G4double HOLDER_POSY = 0.;// 3.875; 	// 3.5; 
G4double HOLDER_POSZ = 0.;// -2.375;	// -3.9;

//What are the type of target that you use?
//1=Gd, 2=Cl, 3=Eu, 4=Co,Na,Cs, 5=None
G4int TARGET = 1;

	//When you choose "Gd target", which model do you use?
	//1:ggarnet, 2:glg4sim, 3:Geant4 default 4:ANNRI-Gd
	G4int MODEL = 4;

	//When you choose both "Gd target" and "ggarnet" or "glg4sim", you have to define the type of Gd.
	//CAPTURE = 1:natural , 2:enriched 157Gd, 3:enriched 155Gd
	//CASCADE = 1:discrete paek + continuum part, 2:discrete peaks, 3:continuum part
	G4int Gd_CAPTURE = 2;
	G4int Gd_CASCADE = 3;
	
	//When you choose "Gd target", "ggarnet" and "continuum part", you have to define the continuum parameter.
	//org:RIPL-3&&kopecky, yano_0908:newly yano tuning, yano_0909:use Iwamoto's parameter, ...
	G4String Gd157_File     = "cont_dat/Gd157_org.dat";	//"cont_dat/Gd157_yano_0909.dat";
	G4String Gd157_ROOTFile = "cont_dat/158GdContTbl__E1SLO4__HFB.root";
	G4String Gd155_File     = "cont_dat/Gd155.dat";
	G4String Gd155_ROOTFile = "cont_dat/156GdContTbl__E1SLO4__HFB.root";

//What is the type of particle that you use?
//11=ParticleGun 1 gamma, 12=ParticleGun 2 gamma, 2=GeneralParticleSource(GPS)
//If you chose Gd target, you should use GPS.
G4int PARTICLE = 2;

	//When you choose "PARTICLE=1", where from does the gamma-ray generate?
	//(POSX,POSY,POZ) = (upstream, up, left side (see on beam))
	G4double EMIT_POSX = 1.7; // -1.75;//0.0;   //	//	
	G4double EMIT_POSY = 3.5; // 3.875; // 0.0;  //	//       
	G4double EMIT_POSZ =-3.9; // -2.375;//0.0;  // 	//	

	//When you choose "PARTICLE=12", how many energies does gamma-rays have?
	//You can choose any energies.(unit is MeV)
	G4double ENERGY_1 = 5.5;   
	G4double ENERGY_2 = 2.5;  

//Do you use VETO detectors(BGO)?
//1=yes: BGO detectors work for VETO and automatically veto events
//2=No : BGO detectors do NOT work for VETO
G4int VETO = 1;

	//You should define the threshold (unit is keV)
	//BGO_thr : BGO detectors threshold (default 100keV)
	//Ene_thr : Ge detectors threshold (default 110keV)
	G4double BGO_thr = 100.;
	G4double Ene_thr = 110.;

//Visualization options
//COLOR = 1:default, 2:target and Ge, 3:Ge,Gd,BGO, 4:Al and Li, 
//COLOR = 99:View Inside <- Don't Run Using This Comand. Just Show For Check.
//INFO  = 1:full information, 2:particles, 3:hit crystals, 4:off
G4int COLOR = 1;
G4int INFO = 4;

//Definition of name of recreate file.
//Default:NewMC.root 
//You should use default, when you execute some Nrun.sh
G4String RECREATE_FILE_NAME = "NewMC.root";

DrawMessage Printing;
// ====================================================================
//     main
// ====================================================================
int main(int argc, char** argv) 
{	
	////////////////////Draw Informations///////////////////////////////
	G4String duct_info = Printing.DefDuct(BEAM_DUCT);
	G4String particle_info = Printing.DefParticle(PARTICLE);
	G4String target_info = Printing.DefTarget(TARGET,MODEL,Gd_CAPTURE,Gd_CASCADE);
	G4String veto_info = Printing.DefVETO(VETO);
	G4cout<<G4endl<<G4endl;
	G4cout<<"===============version 10 @2017.03.09================="<<G4endl;
	G4cout<<setprecision(3);
	G4cout<<"****************Detector Infomation*******************"<<G4endl;
	G4cout<<"Experiment (beam duct) : "<<duct_info<<G4endl;
	G4cout<<"      Target           : "<<target_info<<G4endl;
	G4cout<<"       VETO            : "<<veto_info<<G4endl;
	if(TARGET==1&&MODEL==1&&(Gd_CASCADE==1||Gd_CASCADE==3)){
	G4cout<<"  Continuum Component  : ["<<Gd155_File<<"] and ["<<Gd157_File<<"]"<<G4endl;
	}
	G4cout<<"Target Position(x,y,z) : "<<"("<<SOURCE_POSX<<", "<<SOURCE_POSY<<", "<<SOURCE_POSZ<<")"<<G4endl;
	G4cout<<"Holder Position(x,y,z) : "<<"("<<HOLDER_POSX<<", "<<HOLDER_POSY<<", "<<HOLDER_POSZ<<")"<<G4endl;
	G4cout<<"   Incident Particle   : "<<particle_info<<G4endl;
	if(PARTICLE==11||PARTICLE==12){
	G4cout<<"Emission Point (x,y,z) : "<<"("<<EMIT_POSX<<", "<<EMIT_POSY<<", "<<EMIT_POSZ<<")"<<G4endl;
	}
	if(PARTICLE==12){
	G4cout<<" Emission Enrgy [MeV]  : "<<ENERGY_1<<G4endl;
	G4cout<<" Emission Enrgy [MeV]  : "<<ENERGY_2<<G4endl;
	}
	G4cout<<"****************Detector Infomation*******************"<<G4endl;
	G4cout<<G4endl<<G4endl;
	sleep(3);


	////////////////////Input Log file/////////////////////////
	FILE *logf;
	logf = fopen("MC.log","a");
	time_t Starttime;
	time_t Endtime;
	srand(time(&Starttime));
	G4int seed = rand();
	if (argc >= 3) seed = atoi(argv[2]);
	CLHEP::HepRandom::setTheSeed(seed);
	fprintf(logf,"================================================================================= \n");
	fprintf(logf,"     .....:::::;;;;;##### NewMC.root has been created #####;;;;;:::::.....     \n");
	fprintf(logf,"     Random seed       : %d \n",seed);
	fprintf(logf,"      Start time       : %s \n",ctime(&Starttime));
	fprintf(logf,"        VETO           : %s \n",veto_info.c_str());
	fprintf(logf,"Experiment (beam duct) : %s \n",duct_info.c_str());
	fprintf(logf,"      Target           : %s \n",target_info.c_str());
	if(TARGET==1&&MODEL==1&&(Gd_CASCADE==1||Gd_CASCADE==3)){
		fprintf(logf,"  Continuum Component  : %s and %s \n",Gd155_File.c_str(),Gd157_File.c_str());
	}
	fprintf(logf,"Target Position(x,y,z) : (%.3f,%.3f,%.3f) \n",SOURCE_POSX,SOURCE_POSY,SOURCE_POSZ);
	fprintf(logf,"Holder Position(x,y,z) : (%.3f,%.3f,%.3f) \n \n",HOLDER_POSX,HOLDER_POSY,HOLDER_POSZ);
	fprintf(logf,"   Incident Particle   : %s \n",particle_info.c_str());
	if(PARTICLE==11||PARTICLE==12){
		fprintf(logf,"Emission Point (x,y,z) : (%.3f,%.3f,%.3f) \n",EMIT_POSX,EMIT_POSY,EMIT_POSZ);
	}
	if(PARTICLE==12){
		fprintf(logf," Emission Enrgy [MeV]  : %.3f \n",ENERGY_1);
		fprintf(logf," Emission Enrgy [MeV]  : %.3f \n",ENERGY_2);
	}
	


	// run manager
	G4RunManager* runManager= new G4RunManager;  
	G4cout << G4endl;

	// set mandatory user initialization classes...
	runManager-> SetUserInitialization(new DetectorConstruction);
	runManager-> SetUserInitialization(new Physics_PhysicsList);
	runManager-> SetUserAction(new PrimaryGeneratorAction);
	runManager-> SetUserAction(new RunAction);
	runManager-> SetUserAction(new EventAction);
	runManager-> SetUserAction(new TrackingAction);
	// user session...
	runManager-> Initialize();

#ifdef G4VIS_USE
	// initialize visualization package
	G4VisManager* visManager= new G4VisExecutive;
	visManager-> Initialize();
	G4cout << G4endl;
#endif
	G4TrajectoryDrawByParticleID *modelp = new G4TrajectoryDrawByParticleID();
	modelp->Set("neutron","green");
	modelp->Set("e-","magenta");
	modelp->Set("e+","red");
	modelp->Set("gamma","white");
	visManager->RegisterModel(modelp);
	
	if(argc==1) {
	// interactive session, if no arguments given
#ifdef G4UI_USE_TCSH
		// tcsh-like
		G4UItcsh* tcsh= new G4UItcsh("exAN(%s)[%/]:");
		G4UIterminal* session= new G4UIterminal(tcsh);
#else
		// csh-like
		G4UIterminal* session= new G4UIterminal();
		session-> SetPrompt("exAN(%s)[%/]:");
#endif
		session-> SessionStart();
		delete session;
	} else { 
		// batch mode
		G4UImanager* UImanager= G4UImanager::GetUIpointer();
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager-> ApplyCommand(command+fileName);
	}

#ifdef G4VIS_USE
	delete visManager;
#endif

	delete runManager;  
	G4cout << G4endl;
	
	// terminating...
	time(&Endtime);

	fprintf(logf,"       End time        : %s",ctime(&Endtime));
	fprintf(logf,"     Running time      : %.3fsec \n",difftime(Endtime, Starttime));
	fprintf(logf,"================================================================================= \n");

	fclose(logf);
	return 0;
}
