
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
#include <fstream>
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
#include "TTreeIndex.h"


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

//#include "/gpfs/home/mlf/cdshin/RAT/rat-pac-jsns2-v2.0_fixoverlapping/data/jsns2/externalVer/ext_gamma_gen/JSNS2_ANNRI_Gd/HKTool.h"

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

void GenANNRIGamma(TTree tree, vector<TVector3> mom, vector<TVector3> pos, vector<double> time){
}



int main(int nargn, char ** argn){

	//---- Define the Variable

	TVector3 ini;

	//---- Select Data file address

	string path = "/home/mlf/cdshin/Data/MC/IBDneutron_uniform/rat/";
	string datPath = "./";
	int nFile = atoi(argn[1]);

	string verboseStat = argn[3];
	bool info_print = verboseStat != "0";


	//---- Data file loading


	//RAT::DSReader* ds=new RAT::DSReader( path.c_str() );
	RAT::DSReader* ds=new RAT::DSReader(NULL);


	TChain runTree("runT");
	
	string nF = path+"IBD_"+Form("%05d",nFile)+".root";
	cout << nF << endl;
	ds->Add(nF.c_str());

/*	
	for(int i = 0; i < nFile; ++i){
		string nF = path+"thermalN_"+Form("%d",i)+".root";
		cout << nF << endl;
		ds->Add(nF.c_str());
		//runTree.Add( path.c_str() );
	}
*/	
	

	string outDATFile_name = "";
	outDATFile_name = datPath+ "./dat/ext_gamma_" + Form("%d",nFile) + ".dat";
	ofstream outDATFile;
	outDATFile.open(outDATFile_name.c_str());

	long nEvent=ds->GetTotal();

	//---- Event data Initializing
	//runTree.Add(initFile.c_str());
	//RAT::DS::RunStore::SetReadTree(&runTree);
	cout << "\n Data File Loading ...." << endl;
	//RAT::DS::Run* run = RAT::DS::RunStore::Get()->GetRun(1);


	//---- 
	//if( run == 0)
	//{
	//	cout << "\nRun / file not Found!" << endl;
	//	return 0;
	//}
	cout << "---------------------------------- " << endl; 
	cout << "\n Data Loading Success!" <<endl;


	//---- PMT geometry information
	//

	//---- Histogram Structure

	TVector3 now;
	TVector3 prev;


	TH1D * h_multi_ori = new TH1D("h_multi_ori","; multiplicity;", 15, 0, 15);
	TH1D * h_multi_ANNRI = new TH1D("h_multi_ANNRI","; multiplicity;", 15, 0, 15);

	TH1D * h_singleE_ANNRI = new TH1D("h_singleE_ANNRI","; Energy /MeV;", 50, 0, 10);
	TH1D * h_totalE_ANNRI = new TH1D("h_totalE_ANNRI","; Energy /MeV;", 50, 0, 10);

	TH1D * h_GdType_ANNRI = new TH1D("h_GdType_ANNRI",";;",4,0,4);

	TH1D * h_neutronE = new TH1D("h_neutronE","; Energy /MeV;", 50, 0, 10);


	//---- ANNRI Definition
	const double conv_atomic_unit = 931.494; // MeV. e.g Carbon 12 = 12 * 931.494 MeV;
	const int ID = 22;
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

	TVector3 nCap_pos, nCap_mom;
	TVector3 mom;
	TVector3 pos;

	TVector3 temp_mom, temp_pos;
	TString particle_save;
	TString process_save;
	double gT_save;
	int flag_save;
	int evt_save ;
	

	TTree * tree = new TTree("tree","tree");
	tree->Branch("evt_save", &evt_save, "evt_save/I");
	tree->Branch("temp_mom", &temp_mom);
	tree->Branch("temp_pos", &temp_pos);
	tree->Branch("gT_save", &gT_save, "gT_save/D");
	tree->Branch("flag_save", &flag_save, "flag_save/I");



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

	//nEvent = 10;

	cout<<"CDShin"<<endl;
	cout<<nEvent<<endl;

	for ( long l=0 ; l < nEvent;l++ ){


		//---- Procedure status check	

		//ProgressBar( l+1 , nEvent);

		//---- Event infomation loading

		RAT::DS::Root* root=ds->GetEvent(l);
		mc=root->GetMC();
		//RAT::DS::MCSummary* mcs=mc->GetMCSummary();

		flag_save = 9999;

		int nGamma_ori[3] = {0};
		double gT_ori[3] = {0.};

		tree->Reset();

		int nTrack = mc->GetMCTrackCount();

			cout<<"CDShin : "<<nTrack<<endl;
		for(int i = 0; i < nTrack; i++){
			RAT::DS::MCTrack * add_info = mc->GetMCTrack(i);
			TString particle = add_info->GetParticleName();

			RAT::DS::MCTrackStep * trackStep = add_info->GetMCTrackStep(0);
			TString processName = trackStep->GetProcess();
			double global_time = trackStep->GetGlobalTime();

			temp_mom.Clear();
			temp_pos.Clear();


			evt_save = l;
			temp_mom = trackStep->GetMomentum();
			temp_pos = trackStep->GetEndpoint();
			gT_save = global_time;

			if(particle == "deuteron" && processName == "nCapture"){
				flag_save = 0;
				gT_ori[flag_save] = global_time;
			}
			else if(particle == "Gd158" && processName == "nCapture"){
				flag_save = 1;
				gT_ori[flag_save] = global_time;
			}
			else if(particle == "Gd156" && processName == "nCapture"){
				flag_save = 2;
				gT_ori[flag_save] = global_time;
			}
			else if(particle == "gamma" && processName == "nCapture"){
				if( global_time == gT_ori[0] ) nGamma_ori[0]++;
				if( global_time == gT_ori[1] ) nGamma_ori[1]++;
				if( global_time == gT_ori[2] ) nGamma_ori[2]++;
				flag_save = 3;
			}
			else flag_save = 9999;

			if(flag_save != 9999){
				tree->Fill();
				if( flag_save != 3) h_GdType_ANNRI->Fill(flag_save);
			}
			else continue;

		}

		h_multi_ori->Fill(nGamma_ori[0]);
		h_multi_ori->Fill(nGamma_ori[1]);
		h_multi_ori->Fill(nGamma_ori[2]);

		nGamma_ori[0] = 0;
		nGamma_ori[1] = 0;
		nGamma_ori[2] = 0;

		double tot_gamma_E = 0.;
		double single_gamma_E = 0.;
		double gT = 0.;

		vector<TVector3> loop_mom, loop_pos;
		vector<double> loop_time;

		int nTreeEntry = tree->GetEntries();
		if( nTreeEntry == 0 ) continue;
		tree->BuildIndex("gT_save");
		TTreeIndex * index = (TTreeIndex*)tree->GetTreeIndex();

		//cout << " ==== " << endl;
		//cout << nTreeEntry<< endl;
		for(int iTreeEntry = 0; iTreeEntry < nTreeEntry; ++iTreeEntry){
			Long64_t local = tree->LoadTree( index->GetIndex()[iTreeEntry] );
			tree->GetEntry(local);

			//cout << iGamma << "  " << flag_save << "  " << gT_save << endl;
			if( flag_save == 1 || flag_save == 2){
				nCap_pos = temp_pos;
				nCap_mom = temp_mom;
				if( flag_save == 1){
					Reaction_prod = ANN_gen->Generate_158Gd();
					this_mass = 158*conv_atomic_unit;
					if(info_print) cout << " ++ Gd 158 " << l << " Events " << endl;
				}
				else if( flag_save == 2){
					Reaction_prod = ANN_gen->Generate_156Gd();
					this_mass = 156*conv_atomic_unit;
					if(info_print) cout << " ++ Gd 156 " << l << " Events " << endl;
				}

				double n_energy = sqrt( nCap_mom.Mag2() + pow(this_mass,2) );
				TLorentzVector neutron(nCap_mom, n_energy);

				int nGamma_ANNRI = 0;
				for(auto&& gamma_list: Reaction_prod){
					single_gamma_E = gamma_list.eTot_;
					tot_gamma_E += single_gamma_E;
					TLorentzVector gamma(gamma_list.px_, gamma_list.py_, gamma_list.pz_, gamma_list.eTot_);

					mom = cal_lorentz(gamma, neutron).Vect() * 1e-3;
					pos = nCap_pos;
					gT = gT_save;

					nGamma_ANNRI++;
					h_singleE_ANNRI->Fill(single_gamma_E);
				
					if(info_print){
						cout << flag_save << "  " << 1 << "  " << ID << "  " << 0 << "  " << 0;
						cout << "  " << setw(12) << mom.X() << "  " << setw(12) << mom.Y() << "  " << setw(12) << mom.Z();;
						cout << setw(12) << mom.Mag();;
						cout << "  " << setw(12) << gT;
						cout << "  " << setw(12) << pos.X() << "  " << setw(12) << pos.Y() << "  " << setw(12) << pos.Z();
						cout << endl;
					}

					loop_mom.push_back(mom);
					loop_pos.push_back(pos);
					loop_time.push_back(gT);
					
				/*	
					outDATFile << 1 << "  " << ID << "  " << 0 << "  " << 0;
					outDATFile << "  " << setw(12) << mom.X() << "  " << setw(12) << mom.Y() << "  " << setw(12) << mom.Z();
					outDATFile << setw(6) << gamma_mass;
					outDATFile << "  " << setw(10) << gT;
					outDATFile << "  " << setw(12) << pos.X() << "  " << setw(12) << pos.Y() << "  " << setw(12) << pos.Z();
					outDATFile << endl;
				*/	
					
				}
				h_multi_ANNRI->Fill(nGamma_ANNRI);
			}
			if( flag_save == 0){
				Long64_t local_temp = tree->LoadTree(index->GetIndex()[iTreeEntry+1]);
				tree->GetEntry(local_temp);
				mom = temp_mom * 1e-3;
				pos = temp_pos;
				gT = gT_save;

				single_gamma_E = temp_mom.Mag();
				tot_gamma_E += single_gamma_E;
				h_singleE_ANNRI->Fill(single_gamma_E);
				
				loop_mom.push_back(mom);
				loop_pos.push_back(pos);
				loop_time.push_back(gT);
				if(info_print){
					cout << " ++  Deuteron " << l << " Events " << endl;
					cout << flag_save << "  " << 1 << "  " << ID << "  " << 0 << "  " << 0;
					cout << "  " << setw(12) << mom.X() << "  " << setw(12) << mom.Y() << "  " << setw(12) << mom.Z();
					cout << setw(12) << mom.Mag();
					cout << "  " << setw(12) << gT;
					cout << "  " << setw(12) << pos.X() << "  " << setw(12) << pos.Y() << "  " << setw(12) << pos.Z();
					cout << endl;
				}
				
				/*
				outDATFile << 1 << "  " << ID << "  " << 0 << "  " << 0;
				outDATFile << "  " << setw(12) << mom.X() << "  " << setw(12) << mom.Y() << "  " << setw(12) << mom.Z();
				outDATFile << setw(6) << gamma_mass;
				outDATFile << "  " << setw(10) << gT;
				outDATFile << "  " << setw(12) << pos.X() << "  " << setw(12) << pos.Y() << "  " << setw(12) << pos.Z();
				outDATFile << endl;
				*/
				h_multi_ANNRI->Fill(1);
				
			}

			if(flag_save == 3){
				continue;
			}

		}
		h_totalE_ANNRI->Fill(tot_gamma_E);


		// dat Generation
		int nEvt = loop_mom.size();

		double init_t = 0.;
		outDATFile << nEvt << endl;
		for(int iEvt = 0; iEvt < nEvt; ++iEvt){
			//if(iEvt == 0) init_t = loop_time[iEvt];
			//double genTime = loop_time[iEvt];
			//if(genTime == loop_time[iEvt-1] - init_t) genTime = 0.;
			double genTime = 0.;
			double px = loop_mom[iEvt].X();
			double py = loop_mom[iEvt].Y();
			double pz = loop_mom[iEvt].Z();
			double x = loop_pos[iEvt].X();
			double y = loop_pos[iEvt].Y();
			double z = loop_pos[iEvt].Z();

			outDATFile << 1 << "  " << ID << "  " << 0 << "  " << 0;
			outDATFile << "  " << setw(12) << px << "  " << setw(12) << py << "  " << setw(12) << pz;
			outDATFile << setw(6) << gamma_mass;
			outDATFile << "  " << setw(10) << genTime;
			outDATFile << "  " << setw(12) << x << "  " << setw(12) << y << "  " << setw(12) << z;
			outDATFile << endl;

		}



		loop_mom.clear();
		loop_pos.clear();
		loop_time.clear();

	} //---- Main loop End
	outDATFile.close();
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
	h_singleE_ANNRI->Write();
	h_totalE_ANNRI->Write();
	h_GdType_ANNRI->Write();
	h_neutronE->Write();
	tree->Write();
	out->Close();

	cout << "         Processing Done           " << endl;
	cout << "---------------------------------- " << endl; 

}
