#include "/home/mlf/zayunsna/personal_work/HKTool.h"
using namespace std;

void draw(){

	TFile * file = new TFile("../test_hist.root");
	//TFile * file = new TFile("../test__.root");
	TH1D * h_multi_ori = (TH1D*)file->Get("h_multi_ori");
	TH1D * h_multi_ANNRI = (TH1D*)file->Get("h_multi_ANNRI");
	TH1D * h_singleE_ANNRI = (TH1D*)file->Get("h_singleE_ANNRI");
	TH1D * h_totalE_ANNRI = (TH1D*)file->Get("h_totalE_ANNRI");
	TH1D * h_GdType_ANNTI = (TH1D*)file->Get("h_GdType_ANNRI");
	TH1D * h_neutronE = (TH1D*)file->Get("h_neutronE");

	ToolHelp();

	SetPlotLabel(h_multi_ori);
	SetPlotLabel(h_multi_ANNRI);
	SetPlotLine(h_multi_ori, 1, 2);
	SetPlotLine(h_multi_ANNRI, 2, 2);

	SetPlotLine(h_singleE_ANNRI, 1, 2);
	SetPlotLine(h_totalE_ANNRI, 1, 2);

	SetPlotLabel(h_singleE_ANNRI);
	SetPlotLabel(h_totalE_ANNRI);

	string type[4] = {"Deuteron","Gd-158","Gd-156",""};

	SetPlotLine(h_GdType_ANNTI, 1, 2);
	SetPlotLabel(h_GdType_ANNTI);

	SetPlotLine(h_neutronE, 1, 2);
	SetPlotLabel(h_neutronE);

	TAxis * axis = h_GdType_ANNTI->GetXaxis();
	for(int i = 0; i < 4; ++i) axis->SetBinLabel(i+1, Form("%s", type[i].c_str()));

	h_GdType_ANNTI->GetXaxis()->SetLabelSize(0.05);

	TCanvas * can = new TCanvas("can","can",1200,500);
	can->Divide(2,1);
	can->cd(1)->SetLogy();
	InitPlot();
	h_multi_ori->Draw();
	h_multi_ANNRI->Draw("same");

	can->cd(2);
	InitPlot();
	h_GdType_ANNTI->SetStats(0);
	h_GdType_ANNTI->Draw();
//	h_neutronE->SetStats(0);
//	h_neutronE->Draw();



	TCanvas * can2 = new TCanvas("can2","can2", 1200, 500);
	can2->Divide(2,1);
	can2->cd(1)->SetLogy();
	InitPlot();
	h_singleE_ANNRI->Draw();

	can2->cd(2)->SetLogy();
	InitPlot();
	h_totalE_ANNRI->Draw();

}
