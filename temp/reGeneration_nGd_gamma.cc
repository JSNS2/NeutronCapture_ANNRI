
//--------------------------------------------------------------------
//Built since 2018.08.15
//
//	Author : HK_Jeon
//	Monte-Carlo Simulation Validation code.
//
//	Update Log
//	15-Aug-2018 ( V-0.1.0 ) : Construct basic code structure
//	15-Aug-2018 ( V-0.1.0 ) : Added "Output file setting"
//	16-Aug-2018 ( V-0.1.1 ) : Added processing time calculator.
//	13-Dec-2018 ( V-0.1.2 ) : Added Progress bar viualizing.
//
//	08-Mar-2019 ( V-0.2.0 ) : Initialized for ROOT ver 6
//	12-Mar-2019	( V-0.2.1 ) : Added PMT geometry Plot.
//	10-Jun-2019 ( V-0.2.2 ) : Add PMT Hitmap information for CNN study
//
//	08-Sep-2020 ( V-0.3.0 ) : New RAT version updated - Official version v2.0
//
//--------------------------------------------------------------------


#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include "TVector3.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TLorentzVector.h"


#include "RAT/DSReader.hh"
#include "RAT/DS/MC.hh"
#include "RAT/DS/MCSummary.hh"
#include "RAT/DS/Root.hh"
#include "RAT/DS/PMT.hh"
//#include "RAT/DS/EV.hh"
#include "RAT/DS/Run.hh"
#include "RAT/DS/RunStore.hh"

#include "ANNRIGd_GdNCaptureGammaGenerator.hh"
#include "ANNRIGd_GeneratorConfigurator.hh"
#include "ANNRIGd_ReactionProduct.hh"

#include "/home/mlf/zayunsna/personal_work/HKTool.h"

using namespace std;
namespace AGd = ANNRIGdGammaSpecModel;
AGd::ANNRIGd_GdNCaptureGammaGenerator * ANN_gen = 0;
const int nPMT = 120;

int NumGamma;
double GammaEnergies[15];

const double bufferRadius = 1850.;
const double bufferHeight = 1470.;
const double acrylicRadiut = 1620.;
const double acrylicHeight = 1300.;
const double TankRadius = 1570.;
const double TankHeight = 1250.;


bool CheckSinglePosForE(TVector3 start){
	bool status = false;
	double sX = start.X();
	double sY = start.Y();
	double sZ = start.Z();
	double sR = sqrt( pow( sX, 2) + pow( sY, 2) );


	if( (sR < bufferRadius && abs(sZ) < bufferHeight) ) { status = true;}
//	if( (sR < bufferRadius && abs(sZ) < bufferHeight) ) {
//		if( ( sR <= acrylicRadiut && sR >= TankRadius ) && ( abs(sZ) <= acrylicHeight && abs(sZ) >= TankHeight) ) status = false;
//		else if ( sR <= TankRadius && abs(sZ) <= TankHeight ) status = true;
//		else status = true;
//	}
	if( sX == 0 && sY == 0 && sZ == 0 ) status = false;

	return status;
}


bool CheckSinglePos(TVector3 start){
	bool status = false;
	double sX = start.X();
	double sY = start.Y();
	double sZ = start.Z();
	double sR = sqrt( pow( sX, 2) + pow( sY, 2) );


	if( (sR < TankRadius && abs(sZ) < TankHeight) ) status = true;
	//if( (abs(sX) <= TankRadius && abs(sY) <= TankRadius && abs(sZ) <= TankHeight) ) status = true;
	if( sX == 0 && sY == 0 && sZ == 0 ) status = false;

	return status;
}

TLorentzVector cal_lorentz(const TLorentzVector p1, const TLorentzVector p2){
	TVector3 p1V = p1.Vect();
	TVector3 p2V = p2.Vect();

	double amp = ( (p1V*p2V / (p2.E() + p2.Mag())) - p1.E() ) / p2.Mag();

	TVector3 ref;
	ref.SetX(p1V.X() + amp*p2V.X());
	ref.SetY(p1V.Y() + amp*p2V.Y());
	ref.SetZ(p1V.Z() + amp*p2V.Z());

	double energy = sqrt( ref.Mag2() + pow(p1.Mag() , 2) );
	TLorentzVector out(ref, energy);

	//cout << " Inner check : "  << amp << "  " << energy << endl;

	return out;
}



int main(int nargn, char ** argn){

	//---- Define the Variable

	TVector3 ini;

	//---- Select Data file address

	//string path = "/home/mlf/zayunsna/data_storage/thermal_neutron_220816/";
	string path = "/home/mlf/zayunsna/data_storage/Cf252/Cf_rat_Internal/";
	string initFile = "";
	int nFile = atoi(argn[1]);


	//---- Data file loading


	//RAT::DSReader* ds=new RAT::DSReader( path.c_str() );
	RAT::DSReader* ds=new RAT::DSReader(NULL);


	TChain runTree("runT");
	for(int i = 0; i < nFile; ++i){
		string nF = path+"Cf_"+Form("%d",i)+".root";
		cout << nF << endl;
		ds->Add(nF.c_str());
		//runTree.Add( path.c_str() );
		if( i == 0 ) initFile = nF;
	}

	long nEvent=ds->GetTotal();

	//---- Event data Initializing
	runTree.Add(initFile.c_str());
	RAT::DS::RunStore::SetReadTree(&runTree);
	cout << "\n Data File Loading ...." << endl;
	RAT::DS::Run* run = RAT::DS::RunStore::Get()->GetRun(1);


	//---- 
	if( run == 0)
	{
		cout << "\nRun / file not Found!" << endl;
		return 0;
	}
	cout << "---------------------------------- " << endl; 
	cout << "\n Data Loading Success!" <<endl;


	//---- PMT geometry information
	//

	//---- Histogram Structure

	TVector3 now;
	TVector3 prev;


	TH1D * h_multi_ori = new TH1D("h_multi_ori","; multiplicity;", 15, 0, 15);
	TH1D * h_multi_ANNRI = new TH1D("h_multi_ANNRI","; multiplicity;", 15, 0, 15);


	//---- ANNRI Definition
	const double conv_atomic_unit = 931.494; // MeV. e.g Carbon 12 = 12 * 931.494 MeV;
	const double gamma_mass = 0;
	const int Gd_CASCADE = 1; // CASCADE = 1:discrete peak + continuum part, 2:discrete peaks, 3:continuum part
	const int Gd_CAPTURE = 1; // CAPTURE = 1:natural, 2:enriched 157Gd, 3:enriched 155Gd
	string Gd155_ROOTFile = "./cont_dat/156GdContTbl__E1SLO4__HFB.root";
	string Gd157_ROOTFile = "./cont_dat/158GdContTbl__E1SLO4__HFB.root";

	//	ANN_gen = new AGd::ANNRIGd_GdNCaptureGammaGenerator();
	//	AGd::ANNRIGd_GeneratorConfigurator::Configure(*ANN_gen, Gd_CAPTURE, Gd_CASCADE, Gd155_ROOTFile, Gd157_ROOTFile); 
	//	AGd::ReactionProductVector Reaction_prod;

	ANN_gen = new AGd::ANNRIGd_GdNCaptureGammaGenerator();
	AGd::ANNRIGd_GeneratorConfigurator::Configure(*ANN_gen, Gd_CAPTURE, Gd_CASCADE, Gd155_ROOTFile, Gd157_ROOTFile); 
	AGd::ReactionProductVector Reaction_prod;


	//---- Output file set

	//TFile * out;

	TString out_file = argn[2];

	out_file = out_file;

	if( !out_file.EndsWith(".root",TString::kIgnoreCase) ){
		out_file.Append(".root");
	}
	if( !gSystem->AccessPathName(out_file) ){
		gSystem->Rename(out_file, out_file+".Previous_result");
	}

	TFile * out = new TFile(out_file,"RECREATE");

	//---- Process start

	double this_mass = 0.;

	double nE = 0.;
	TVector3 nCap_pos, nCap_mom;
	TVector3 mom;
	int ID;
	vector<TVector3> gamma_vec;
	vector<TLorentzVector> gamma_info;


	TTree * tree = new TTree("track","track");
	tree->Branch("nE", &nE, "nE/D");
	tree->Branch("nCap_pos", &nCap_pos);



	RAT::DS::MC* mc;
	cout << " Total " << nEvent << " event founded. " << endl;
	cout << "\n---------------------------------- " << endl; 
	cout << "               Start               " << endl;
	cout << "---------------------------------- " << endl; 
	cout << endl;

	//---- Processing timer set

	TStopwatch * clck = new TStopwatch();
	clck->Start();

	//---- Main Loop Start

	//nEvent =1;


	for ( long l=0 ; l < nEvent;l++ ){


		//---- Procedure status check	

		ProgressBar( l+1 , nEvent);


		//---- Event infomation loading

		RAT::DS::Root* root=ds->GetEvent(l);
		mc=root->GetMC();
		//RAT::DS::MCSummary* mcs=mc->GetMCSummary();

		gamma_vec.clear();
		gamma_info.clear();

		int nGamma = 0;
		int nGamma_ori = 0;
		int nGamma_ANNRI = 0;

		int nTrack = mc->GetMCTrackCount();
		for(int i = 0; i < nTrack; i++){
			RAT::DS::MCTrack * add_info = mc->GetMCTrack(i);
			string particle = add_info->GetParticleName();

			int nTrackStep = add_info->GetMCTrackStepCount();
			RAT::DS::MCTrackStep * trackStep_temp = add_info->GetMCTrackStep(0);
			string temp_name = trackStep_temp->GetProcess();

			cout << l << "  " << i << "  " << particle << "  " << temp_name << endl;
			
			for(int iTS = 0; iTS < nTrackStep; iTS++){
				RAT::DS::MCTrackStep * trackStep = add_info->GetMCTrackStep(iTS);
				TVector3 position = trackStep->GetEndpoint();
				string process_name = trackStep->GetProcess();
				//cout << i << "  " << particle << "  " << process_name << "  " << trackStep->GetGlobalTime() << endl;
				
				if( process_name == "nCapture"){
					if(particle == "neutron"){
						nCap_pos = trackStep->GetEndpoint();
						nCap_mom = trackStep->GetMomentum();
						nGamma_ori = 0;
						nGamma_ANNRI = 0;
					}
					if(particle == "gamma"){ 
						nGamma_ori++; 
					}
					if(particle == "Gd158"){
						Reaction_prod = ANN_gen->Generate_158Gd();
						this_mass = 158*conv_atomic_unit;
					}
					else if( particle == "Gd156"){
						Reaction_prod = ANN_gen->Generate_156Gd();
						this_mass = 156*conv_atomic_unit;
					}
					else {
						continue;
					}

					double n_energy = sqrt( nCap_mom.Mag2() + pow(this_mass,2) );
					TLorentzVector neutron(nCap_mom, n_energy);

					for(auto&& gamma_list: Reaction_prod){
						TLorentzVector gamma(gamma_list.px_, gamma_list.py_, gamma_list.pz_, gamma_list.eTot_);
						gamma_info.push_back(cal_lorentz(gamma, neutron));
						//double b = gamma.Vect().Mag();
						//double a = cal_lorentz(gamma, neutron).Vect().Mag();
						//cout << "test : " << a << "  " << b << endl;
						gamma_vec.push_back( gamma_info.back().Vect());
						nGamma_ANNRI++;
					}
					//cout << i << "  " << process_name << "  " << particle << "  " << this_mass << "  " << nCap_pos.X() << "  " << nCap_pos.Y() << "  " << nCap_pos.Z() << "  " << n_energy << "  " << gamma_vec.size() << endl;
				}
				//if(nGamma_ori != 0 || nGamma_ANNRI != 0 )cout << l << "  " << i << "  " << iTS << "  " << nGamma_ori << "  " << nGamma_ANNRI << endl;

				h_multi_ANNRI->Fill(nGamma_ANNRI);
				h_multi_ori->Fill(nGamma_ori);
			}
		}


		nGamma = (int)gamma_vec.size();
		//cout << nGamma << endl;
		for(int iGamma = 0 ; iGamma < nGamma; ++iGamma){
			ID = 22; // gamma PDG
			mom = gamma_vec.at(iGamma) * 1e-3; // Unit : GeV
/*
			cout << 1 << "  " << ID << "  " << 0 << "  " << 0;
			cout << "  " << setw(12) << mom.X() << "  " << setw(12) << mom.Y() << "  " << setw(12) << mom.Z();
			cout << setw(12) << gamma_mass;
			cout << "  " << setw(12) << 0;
			cout << "  " << setw(12) << nCap_pos.X() << "  " << setw(12) << nCap_pos.Y() << "  " << setw(12) << nCap_pos.Z();
			cout << "  " << 0 << "  " << 0 << "  " << 0 << endl;
*/
		}



		tree->Fill();
	} //---- Main loop End


	//---- timer Stop and time calculation
	clck->Stop();
	double CPU_time = clck->CpuTime();
	double Real_time = clck->RealTime();
	double Avg_time = CPU_time / nEvent;

	cout << "\n" << endl;
	cout << "----------------------------------\n " << endl; 
	cout << " Time to Validate " << nEvent << " events " << endl;
	cout << endl;
	cout << " CPU Time               = " << CPU_time << " sec" << endl;
	cout << "                          " << CPU_time/60 << " min" << endl;
	cout << " Real Time              = " << Real_time << " sec" << endl;
	cout << "                          " << Real_time/60 << " min" << endl;
	cout << " Average Time per event = " << Avg_time << " sec" << endl;
	cout << "\n---------------------------------- " << endl; 

	//---- Data Wirte

	h_multi_ori->Write();
	h_multi_ANNRI->Write();
	out->Close();

	cout << "         Processing Done           " << endl;
	cout << "---------------------------------- " << endl; 

}
