
#include <iostream>
#include <TROOT.h>
#include <TStyle.h>
#include <TPad.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1.h>
#include <TBox.h>
#include <TLatex.h>


using namespace std;


namespace Color {
	enum Code{
		FG_RED		= 31,
		FG_GREEN	= 32,
		FG_BLUE		= 34,
		FG_DEFAULT	= 39,
		BG_RED		= 41,
		BG_GREEN	= 42,
		BG_BLUE		= 44,
		BG_DEFAULT	= 49
	};
	class Modifier {
		Code code;
		public :
		Modifier(Code pCode) : code(pCode) {}
		friend std::ostream&
			operator<<(std::ostream& os, const Modifier& mod) {
				return os << "\033[" << mod.code << "m";
			}
	};
}

void ToolHelp(){

	Color::Modifier ori(Color::FG_DEFAULT);
	Color::Modifier red(Color::FG_RED);
	Color::Modifier green(Color::FG_GREEN);

	cout << "\n ------------------------------ HKTool Function List ------------------------------ " << endl;
	cout << "  -  "<<green<<"double "<<ori<<"ExpoFitWithBinReject("<<green<<"double "<<ori<<"*val, "<<green<<"double "<<ori<<"*par, "<<green<<"double"<<ori<<" start,"<<green<<" double"<<ori<<" end) " << endl;

	cout << "  -  "<<green<<"double "<<ori<<"errorMulti("<<green<<"double "<<ori<<"a, "<<green<<"double "<<ori<<"aE , "<<green<<"double "<<ori<<"b, "<<green<<"double "<<ori<<"bE) " << endl;

	cout << "  -  "<<green<<"double "<<ori<<"errorDiv("<<green<<"double "<<ori<<"a, "<<green<<"double "<<ori<<"aE , "<<green<<"double "<<ori<<"b, "<<green<<"double "<<ori<<"bE) " << endl;

	cout << "  -  "<<green<<"double "<<ori<<"CalSignificanceErr("<<green<<"double "<<ori<<"sig, "<<green<<"double "<<ori<<"sigErr, "<<green<<"double "<<ori<<"bkg, "<<green<<"double "<<ori<<"bkgErr) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"InitPlot() " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"SetPlotLabel("<<green<<"TNamed "<<ori<<"*name, "<<green<<"double "<<ori<<"yLabelSize = 0.03, "<<green<<"double "<<ori<<"yTitleOffset = 1.1, "<<green<<"double "<<ori<<"xLabelSize = 0.03, "<<green<<"double "<<ori<<"xTitleOffset = 1.1) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"SetPlotAxisTitle("<<green<<"TNamed "<<ori<<"*name, "<<green<<"TString "<<ori<<"xTitle, "<<green<<"TString "<<ori<<"yTitle) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"SetPlotAxisRange("<<green<<"TNamed "<<ori<<"*name, "<<green<<"double "<<ori<<"xMin, "<<green<<"double "<<ori<<"xMax, "<<green<<"double "<<ori<<"yMin, "<<green<<"double "<<ori<<"yMax) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"SetPlotLine("<<green<<"TNamed"<<ori<<" * name, "<<green<<"int "<<ori<<"seq = 1, "<<green<<"int "<<ori<<"width = 2) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"SetPlotNorm("<<green<<"TNamed "<<ori<<"* name, "<<green<<"double "<<ori<<"frac = 1) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"DrawWithDetectorConfigureXY("<<green<<"TH2D "<<ori<<"*inputTH2D) " << endl;

	cout << "  -  "<<green<<"void "<<ori<<"DrawWithDetectorConfigureXZ("<<green<<"TH2D "<<ori<<"*inputTH2D) " << endl;

	cout << " ---------------------------------------------------------------------------------- \n" << endl;

}

void ProgressBar(int i, int n){
	float progress = (float)i / (float)n;
	cout << " [";
	int pos = int(70*progress);
	for ( int i = 0; i < 70; ++i){
		if		( i < pos )		cout << "=";
		else if	( i == pos )	cout << ">";
		else					cout << " ";
	}
	cout << "]" << int(progress*100.0) << "%\r";
	cout.flush();
}

double ExpoFitWithBinReject(double *val, double *par, double start, double end){
	double x = val[0];
	double p1 = par[0];
	double p2 = par[1];
	//double fcn = exp(p1-x/p2);
	double fcn = p1-x/p2;

	if( x > start && x < end){
		TF1::RejectPoint();
		return 0;
	}
	return fcn;
}

double errorMulti(double a, double aE , double b, double bE){
	double tar = a*b;
	double tarSqrt = sqrt(pow(aE/a,2) + pow(bE/b, 2));
	double err = tar*tarSqrt;

	return err;
}

double errorDiv(double a, double aE , double b, double bE){
	double tar = a/b;
	double tarSqrt = sqrt(pow(aE/a,2) + pow(bE/b, 2));
	double err = tar*tarSqrt;

	return err;
}


double CalSignificanceErr(double sig, double sigErr, double bkg, double bkgErr){

	double temp = sig+bkg;
	double tempErr = sqrt(sigErr*sigErr + bkgErr*bkgErr);
	double tempErr2 = tempErr/temp * 100;

	double denomi = sqrt(temp);
	double denomiErr = sqrt(tempErr2);
	denomiErr *= denomi;

	double tar = sig/denomi;
	double tarSqrt = sqrt( pow( sigErr/sig , 2) + pow( denomiErr/denomi ,2 ) );
	double err = tar*tarSqrt;

	return err;
}


void InitPlot(){
	gStyle->SetStatColor(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetTitleH(0.05);
	gStyle->SetTitleBorderSize(1);
	gStyle->SetStatX(0.84);
	gStyle->SetStatY(0.88);
	gStyle->SetStatFontSize(0.023);
	gStyle->SetStatBorderSize(1);
	gStyle->SetStatStyle(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetFrameBorderSize(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetLegendBorderSize(1);
	gStyle->SetCanvasColor(0);
	gStyle->SetPadTopMargin(0.12);
	gStyle->SetPadBottomMargin(0.12);
	gStyle->SetPadRightMargin(0.14);
	gStyle->SetPadLeftMargin(0.12);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPadGridX(1);
	gStyle->SetPadGridY(1);
	gStyle->SetLabelSize(0.03,"X");
	gStyle->SetLabelSize(0.03,"Y");
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(0.6);
	gStyle->SetFuncWidth(1);
	gStyle->SetFuncColor(2);
	gStyle->SetEndErrorSize(0);
	gStyle->SetNdivisions(505);
	gStyle->SetPalette(55);
	gStyle->SetOptStat("neMri");
	gStyle->SetOptFit(0001);
	gStyle->SetStatW(0.3);
	gStyle->SetStatFontSize(0.03);

	gStyle->SetNumberContours(100);

	gPad->SetGrid();

	gROOT->GetColor(3)->SetRGB(0.,0.7,0.); // Green (0, 1, 0) -> (0, 0.7, 0);
	gROOT->GetColor(5)->SetRGB(1.,0.5,0.); // Yellow (1, 1, 0) -> (1, 0.5, 0)
}

double radius_sus	=	2300.;
double radius_black	=	1865.;
double radius_acryl	=	1600.;

double height_sus	=	3512.;
double height_black	=	3000.;
double height_acryl	=	2398.;


void DrawWithDetectorConfigureXY(TH2D *inputTH2D){
	TString t_name = inputTH2D->ClassName();
	if( t_name.Contains("TH2D")){
		TNamed * name = inputTH2D;
		int ID = name->Hash();
		cout << " ## \n Event vertex Drawing with Detector Configuration == X & Y axis. " << endl;
		double xMin = inputTH2D->GetXaxis()->GetXmin();
		double xMax = inputTH2D->GetXaxis()->GetXmax();
		double yMin = inputTH2D->GetYaxis()->GetXmin();
		double yMax = inputTH2D->GetYaxis()->GetXmax();
		bool conditionCheck = false;
		if( xMin <= -2400 && xMax >= 2400 && yMin <= -2400 && yMax >= 2400) conditionCheck = true;

		if( conditionCheck ){
			inputTH2D->GetXaxis()->SetRangeUser(-2400,2400);
			inputTH2D->GetYaxis()->SetRangeUser(-2400,2400);
			TEllipse * cSus = new TEllipse(0., 0., radius_sus, radius_sus);
			cSus->SetLineColor(kGray);
			cSus->SetLineWidth(4);
			cSus->SetFillStyle(0);
			TEllipse * cBlack = new TEllipse(0., 0., radius_black, radius_black);
			cBlack->SetLineColor(1);
			cBlack->SetLineWidth(4);
			cBlack->SetFillStyle(0);
			TEllipse * cAcryl = new TEllipse(0., 0., radius_acryl, radius_acryl);
			cAcryl->SetLineColor(kCyan);
			cAcryl->SetLineWidth(4);
			cAcryl->SetFillStyle(0);

			TLatex * lSus = new TLatex(); lSus->SetTextFont(132); lSus->SetTextColor(kGray);
			TLatex * lBlack = new TLatex(); lBlack->SetTextFont(132); lBlack->SetTextColor(1);
			TLatex * lAcryl = new TLatex(); lAcryl->SetTextFont(132); lAcryl->SetTextColor(kCyan);
			//TCanvas * Draw_TopView = new TCanvas("Draw_TopView","Draw_TopView",800,600);
			TCanvas * Draw_TopView = new TCanvas(Form("%d",ID),"Draw_TopView",800,600);
			Draw_TopView->SetGrid();
			inputTH2D->Draw("colz"); inputTH2D->SetStats(0);
			cSus->Draw("same");
			lSus->DrawLatex(1000, 2150, " Stainless Tank ");
			cBlack->Draw("same");
			lBlack->DrawLatex(650, 1750, " Gamma Catcher ");
			cAcryl->Draw("same");
			lAcryl->DrawLatex(900, 1350, " Acrylic Tank ");
			gStyle->SetPalette(55);
		}
		else{
			cerr << " \t Function is terminated !!! " << endl;
			cerr << " Please adjust Input histogram axis range. " << endl;
			cerr << " Input TH2D X & Y axis range should be larger than sus-tank, " << endl;
			cerr << "   xMin : -2400  || now : " << setw(6) << xMin << endl;
			cerr << "   xMax :  2400  || now : " << setw(6) << xMax << endl;
			cerr << "   yMin : -2400  || now : " << setw(6) << yMin << endl;
			cerr << "   yMax :  2400  || now : " << setw(6) << yMax << endl;
			return 0;
		}
	}
	else{
		cerr << " \t Function is terminated !!! " << endl;
		cerr << " Input Forman is not 'TH2D'. " << endl;
		return 0;
	}
}


void DrawWithDetectorConfigureXZ(TH2D *inputTH2D){
	TString t_name = inputTH2D->ClassName();
	if( t_name.Contains("TH2D")){
		TNamed * name = inputTH2D;
		int ID = name->Hash();
		cout << " ## \n Event vertex Drawing with Detector Configuration == X & Z axis. " << endl;
		double xMin = inputTH2D->GetXaxis()->GetXmin();
		double xMax = inputTH2D->GetXaxis()->GetXmax();
		double zMin = inputTH2D->GetYaxis()->GetXmin();
		double zMax = inputTH2D->GetYaxis()->GetXmax();
		bool conditionCheck = false;
		if( xMin <= -2400 && xMax >= 2400 && zMin <= -1800 && zMax >= 1800) conditionCheck = true;

		if( conditionCheck ){
			inputTH2D->GetXaxis()->SetRangeUser(-2400,2400);
			inputTH2D->GetYaxis()->SetRangeUser(-1800,1800);
			TBox * bSus = new TBox(-1*radius_sus, -1*height_sus/2., radius_sus, height_sus/2.);
			bSus->SetLineColor(kGray);
			bSus->SetLineWidth(4);
			bSus->SetFillStyle(0);
			TBox * bBlack = new TBox(-1*radius_black, -1*height_black/2., radius_black, height_black/2.);
			bBlack->SetLineColor(1);
			bBlack->SetLineWidth(4);
			bBlack->SetFillStyle(0);
			TBox * bAcryl = new TBox(-1*radius_acryl, -1*height_acryl/2., radius_acryl, height_acryl/2.);
			bAcryl->SetLineColor(kCyan);
			bAcryl->SetLineWidth(4);
			bAcryl->SetFillStyle(0);

			TLatex * lSus = new TLatex(); lSus->SetTextFont(132); lSus->SetTextColor(kGray);
			TLatex * lBlack = new TLatex(); lBlack->SetTextFont(132); lBlack->SetTextColor(1);
			TLatex * lAcryl = new TLatex(); lAcryl->SetTextFont(132); lAcryl->SetTextColor(kCyan);
			//TCanvas * Draw_SideView = new TCanvas("Draw_SideView","Draw_SideView",800,600);
			TCanvas * Draw_SideView = new TCanvas(Form("%d",ID),"Draw_SideView",800,600);
			Draw_SideView->SetGrid();
			inputTH2D->Draw("colz"); inputTH2D->SetStats(0);
			bSus->Draw("same");
			lSus->DrawLatex(-2000, 1850, " Stainless Tank ");
			bBlack->Draw("same");
			lBlack->DrawLatex(-1700, 1550, " Gamma Catcher ");
			bAcryl->Draw("same");
			lAcryl->DrawLatex(-1400, 1250, " Acrylic Tank ");
			gStyle->SetPalette(55);
		}
		else{
			cerr << " \t Function is terminated !!! " << endl;
			cerr << " Please adjust Input histogram axis range. " << endl;
			cerr << " Input TH2D X & Z axis range should be larger than sus-tank, " << endl;
			cerr << "   xMin : -2400  || now : " << setw(6) << xMin << endl;
			cerr << "   xMax :  2400  || now : " << setw(6) << xMax << endl;
			cerr << "   zMin : -1800  || now : " << setw(6) << zMin << endl;
			cerr << "   zMax :  1800  || now : " << setw(6) << zMax << endl;
			return 0;
		}
	}
	else{
		cerr << " Input Forman is not 'TH2D'. " << endl;
		cerr << " Function is terminated. " << endl;
		return 0;
	}
}


void SetPlotLabel(TNamed * name, double yLabelSize = 0.03, double yTitleOffset = 1.1, double xLabelSize = 0.03, double xTitleOffset = 1.1){
	//	InitPlot();
	TString t_name = name->ClassName();
	if( t_name.Contains("TGraph")){
		((TGraph*)name)->GetYaxis()->SetLabelSize(yLabelSize);
		((TGraph*)name)->GetYaxis()->SetTitleOffset(yTitleOffset);
		((TGraph*)name)->GetXaxis()->SetLabelSize(xLabelSize);
		((TGraph*)name)->GetXaxis()->SetTitleOffset(xTitleOffset);
	}
	if( t_name.Contains("TH1")){
		((TH1*)name)->GetYaxis()->SetLabelSize(yLabelSize);
		((TH1*)name)->GetYaxis()->SetTitleOffset(yTitleOffset);
		((TH1*)name)->GetXaxis()->SetLabelSize(xLabelSize);
		((TH1*)name)->GetXaxis()->SetTitleOffset(xTitleOffset);
	}
	if( t_name.Contains("TH2")){
		((TH2*)name)->GetYaxis()->SetLabelSize(yLabelSize);
		((TH2*)name)->GetYaxis()->SetTitleOffset(yTitleOffset);
		((TH2*)name)->GetXaxis()->SetLabelSize(xLabelSize);
		((TH2*)name)->GetXaxis()->SetTitleOffset(xTitleOffset);
	}
}

void SetPlotAxisTitle(TNamed *name, TString xTitle, TString yTitle){
	TString t_name = name->ClassName();
	if( t_name.Contains("TGraph")){
		((TGraph*)name)->GetXaxis()->SetTitle(xTitle);
		((TGraph*)name)->GetYaxis()->SetTitle(yTitle);
	}
	if( t_name.Contains("TH1")){
		((TH1*)name)->GetXaxis()->SetTitle(xTitle);
		((TH1*)name)->GetYaxis()->SetTitle(yTitle);
	}
	if( t_name.Contains("TH2")){
		((TH2*)name)->GetXaxis()->SetTitle(xTitle);
		((TH2*)name)->GetYaxis()->SetTitle(yTitle);
	}

}


void SetPlotAxisRange(TNamed *name, double xMin, double xMax, double yMin, double yMax){
	TString t_name = name->ClassName();
	if( t_name.Contains("TGraph")){
		((TGraph*)name)->GetXaxis()->SetRangeUser(xMin,xMax);
		((TGraph*)name)->GetYaxis()->SetRangeUser(yMin,yMax);
	}
	if( t_name.Contains("TH1")){
		((TH1*)name)->GetXaxis()->SetRangeUser(xMin,xMax);
		((TH1*)name)->GetYaxis()->SetRangeUser(yMin,yMax);
	}
	if( t_name.Contains("TH2")){
		((TH2*)name)->GetXaxis()->SetRangeUser(xMin,xMax);
		((TH2*)name)->GetYaxis()->SetRangeUser(yMin,yMax);
	}

}

void SetPlotLine(TNamed * name, int seq = 1, int width = 2){

	int palette[6] = {1, 1, 2, 4, 5, 7};

	TString t_name = name->ClassName();
	if( t_name.Contains("TGraph")){
		((TGraph*)name)->SetLineColor(palette[seq]);
		((TGraph*)name)->SetLineWidth(width);
		((TGraph*)name)->SetMarkerColor(palette[seq]);
	}
	if( t_name.Contains("TH1")){
		((TH1*)name)->SetLineColor(palette[seq]);
		((TH1*)name)->SetLineWidth(width);
		((TH1*)name)->SetMarkerColor(palette[seq]);
	}

}

void SetPlotNorm(TNamed * name, double frac = 1){

	TString t_name = name->ClassName();
	if( t_name.Contains("TGraph")){
		cerr << " [ERROR]:: TGraph cannot be normalized by area " << endl;
		cerr << "           Normalization has been skipped. " << endl;
	}
	if( t_name.Contains("TH1")){
		((TH1*)name)->Scale(frac/ ((TH1*)name)->Integral());
	}
}
/*
   void SetLegendInit(TLegend * name, double lineOpt = 1, double fillOpt = 0, double fillStyle = 0){
   name->SetLineColor(lineOpt);
   name->SetFillColor(fillOpt);
   name->SetFillStyle(fillStyle);
   }
   */
