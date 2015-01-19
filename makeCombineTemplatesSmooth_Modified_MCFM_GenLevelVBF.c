#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"

using namespace std;

//Initializers
int useAnomalousCouplings=kAddfLQ;
bool useDjettagging=true;
enum histtypes{kSigHist,kBkgHist,kIntHist,kBSIHist,kBSI10Hist,kBSI25Hist};
void makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics,int Djettag);
double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
void twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, TH2F* finaloutput, int outputtype, TH2F* finaloutput2=0, int output2type=-99);
void floorSignalTemplates(TH2F* hsig, TH2F* hbkg, TH2F* hinterf);
void floorBkgTemplate(TH2F* hbkg);

//Main Function, runs over all desired iterations
void makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF(){
	const int kNumSyst=5;
	int systematics[kNumSyst]={0,1,-1,2,-2};
	for(int i=0;i<kNumSyst;++i){
		for(int CoM=7;CoM<9;++CoM){
			for(int channel=0;channel<3;++channel){
				if(useDjettagging){for(int Djettag=-1;Djettag<2;++Djettag) makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],Djettag); }
				else makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],0);	
			}	
		}
	}
}

//Function to build one template
// folder = 0,1,2 (final state corresponds to 2e2mu,4mu,4e respectively)
// erg_tev = 7,8 (CoM energy)
// tFitD = [0,16] (choice of Discriminant, see FitDimensionsList.h for list; only tFitd works right now)
// Systematics = [-2,2] (Flag for systematics. 0=Nominal, +/-1=QCD, +/-2=PDF)
void makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, int Djettag){
	cout<<"SCALE "<<user_folder[folder]<<" "<<erg_tev<<" "<<Systematics<<" "<<Djettag<<endl;
	TString INPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_";
	TString INPUT_K3A_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_";
	TString INPUT_SMOOTH_NAME = "HtoZZ4l_MCFM_125p6_SmoothTemplates_";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine_";
  if (useAnomalousCouplings > 0){
    INPUT_NAME += strAnomalousType[useAnomalousCouplings];
    INPUT_K3A_NAME += strAnomalousType[useAnomalousCouplings];
    INPUT_SMOOTH_NAME += strAnomalousType[useAnomalousCouplings];
    OUTPUT_NAME += strAnomalousType[useAnomalousCouplings];
  }
  INPUT_NAME += "_Raw__GenLevelVBF_";
	INPUT_K3A_NAME += "__GenLevelVBF_";
	INPUT_SMOOTH_NAME += "__GenLevelVBF_";
	OUTPUT_NAME += "__GenLevelVBF_";		

  INPUT_NAME += TString(strFitDim[tFitD]) + "_";
	INPUT_K3A_NAME += TString(strFitDim[tFitD]) + "_";
	INPUT_SMOOTH_NAME += TString(strFitDim[tFitD]) + "_";
	OUTPUT_NAME += TString(strFitDim[tFitD]) + "_";
	TString comstring;
	comstring.Form("%i",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/",erg_tev);
	if(Djettag==0){
		if(Systematics==0) OUTPUT_NAME += "Nominal.root";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD.root";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD.root";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF.root";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF.root";

		if(Systematics==0) INPUT_SMOOTH_NAME += "Nominal.root";
		if (Systematics == 1) INPUT_SMOOTH_NAME += "SysUp_ggQCD.root";
		if (Systematics == -1) INPUT_SMOOTH_NAME += "SysDown_ggQCD.root";
		if (Systematics == 2) INPUT_SMOOTH_NAME += "SysUp_ggPDF.root";
		if (Systematics == -2) INPUT_SMOOTH_NAME += "SysDown_ggPDF.root";

		if(Systematics==0) INPUT_NAME += "Nominal.root";
		if (Systematics == 1) INPUT_NAME += "SysUp_ggQCD.root";
		if (Systematics == -1) INPUT_NAME += "SysDown_ggQCD.root";
		if (Systematics == 2) INPUT_NAME += "SysUp_ggPDF.root";
		if (Systematics == -2) INPUT_NAME += "SysDown_ggPDF.root";

		if(Systematics==0) INPUT_K3A_NAME += "Nominal.root";
		if (Systematics == 1) INPUT_K3A_NAME += "SysUp_ggQCD.root";
		if (Systematics == -1) INPUT_K3A_NAME += "SysDown_ggQCD.root";
		if (Systematics == 2) INPUT_K3A_NAME += "SysUp_ggPDF.root";
		if (Systematics == -2) INPUT_K3A_NAME += "SysDown_ggPDF.root";
	}
	else if(Djettag==-1){
		if(Systematics==0) OUTPUT_NAME += "Nominal_nonDjet.root";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD_nonDjet.root";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD_nonDjet.root";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF_nonDjet.root";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF_nonDjet.root";

		if(Systematics==0) INPUT_SMOOTH_NAME += "Nominal_nonDjet.root";
		if (Systematics == 1) INPUT_SMOOTH_NAME += "SysUp_ggQCD_nonDjet.root";
		if (Systematics == -1) INPUT_SMOOTH_NAME += "SysDown_ggQCD_nonDjet.root";
		if (Systematics == 2) INPUT_SMOOTH_NAME += "SysUp_ggPDF_nonDjet.root";
		if (Systematics == -2) INPUT_SMOOTH_NAME += "SysDown_ggPDF_nonDjet.root";

		if(Systematics==0) INPUT_NAME += "Nominal_nonDjet.root";
		if (Systematics == 1) INPUT_NAME += "SysUp_ggQCD_nonDjet.root";
		if (Systematics == -1) INPUT_NAME += "SysDown_ggQCD_nonDjet.root";
		if (Systematics == 2) INPUT_NAME += "SysUp_ggPDF_nonDjet.root";
		if (Systematics == -2) INPUT_NAME += "SysDown_ggPDF_nonDjet.root";

		if(Systematics==0) INPUT_K3A_NAME += "Nominal_nonDjet.root";
		if (Systematics == 1) INPUT_K3A_NAME += "SysUp_ggQCD_nonDjet.root";
		if (Systematics == -1) INPUT_K3A_NAME += "SysDown_ggQCD_nonDjet.root";
		if (Systematics == 2) INPUT_K3A_NAME += "SysUp_ggPDF_nonDjet.root";
		if (Systematics == -2) INPUT_K3A_NAME += "SysDown_ggPDF_nonDjet.root";
	}
	else if(Djettag==1){
		if(Systematics==0) OUTPUT_NAME += "Nominal_Djet.root";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD_Djet.root";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD_Djet.root";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF_Djet.root";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF_Djet.root";

		if(Systematics==0) INPUT_SMOOTH_NAME += "Nominal_Djet.root";
		if (Systematics == 1) INPUT_SMOOTH_NAME += "SysUp_ggQCD_Djet.root";
		if (Systematics == -1) INPUT_SMOOTH_NAME += "SysDown_ggQCD_Djet.root";
		if (Systematics == 2) INPUT_SMOOTH_NAME += "SysUp_ggPDF_Djet.root";
		if (Systematics == -2) INPUT_SMOOTH_NAME += "SysDown_ggPDF_Djet.root";

		if(Systematics==0) INPUT_NAME += "Nominal_Djet.root";
		if (Systematics == 1) INPUT_NAME += "SysUp_ggQCD_Djet.root";
		if (Systematics == -1) INPUT_NAME += "SysDown_ggQCD_Djet.root";
		if (Systematics == 2) INPUT_NAME += "SysUp_ggPDF_Djet.root";
		if (Systematics == -2) INPUT_NAME += "SysDown_ggPDF_Djet.root";

		if(Systematics==0) INPUT_K3A_NAME += "Nominal_Djet.root";
		if (Systematics == 1) INPUT_K3A_NAME += "SysUp_ggQCD_Djet.root";
		if (Systematics == -1) INPUT_K3A_NAME += "SysDown_ggQCD_Djet.root";
		if (Systematics == 2) INPUT_K3A_NAME += "SysUp_ggPDF_Djet.root";
		if (Systematics == -2) INPUT_K3A_NAME += "SysDown_ggPDF_Djet.root";
	}
	cout<<"Input Smooth: "<<INPUT_SMOOTH_NAME<<endl;
	cout<<"Input: "<<INPUT_NAME<<endl;
	cout<<"Input k3a: "<<INPUT_K3A_NAME<<endl;
	cout<<"Output: "<<OUTPUT_NAME<<endl;

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float lowside[3]={220,230,240};

	for(int lo=0;lo<1;lo++){
		TString coutput_common = user_dir + erg_dir;
		coutput_common += user_folder[folder] + "/";
		TString cinput_common = user_TemplateswithTrees_dir + erg_dir + user_folder[folder] + "/";

		TString coutput = coutput_common + OUTPUT_NAME;
		TFile* foutput = new TFile(coutput,"recreate");

		float ZZMass_cut[2]={lowside[lo],1600};
		float ZZwidth = 20.0;
		const int nbinsx=(ZZMass_cut[1]-ZZMass_cut[0])/ZZwidth;
		float kDXarray[nbinsx+1];
		for(int bin=0;bin<nbinsx+1;bin++){
			kDXarray[bin] = ZZMass_cut[0] + ZZwidth*bin;
		}
	
		int nbinsy=30;
		float kDYarray[nbinsy+1];
		float kDY_bounds[2]={0,1};
		if(tFitD==3){kDY_bounds[0]=-7.0;kDY_bounds[1]=3.0;};
		if(tFitD==4 || tFitD==5 || tFitD==9 || tFitD==10 || tFitD==11 || tFitD==15){kDY_bounds[0]=-1.0;kDY_bounds[1]=1.0;};
		for(int bin=0;bin<nbinsy+1;bin++){
			double binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		}

		const int kNumTemplates=8;
		double* overall_scale[kNumTemplates];
		double* expectedNormalizations[3][kNumTemplates];
		TH2F** D_temp_2D_expk3a[kNumTemplates];
		TH2F** D_temp_2D_exp[kNumTemplates];
		TH2F** D_temp_2D[kNumTemplates];
		TH2F* hZX_Unconditional;
		TH2F* hqqZZ_Unconditional;
		TH2F* hZX_Unconditional_exp;
		TH2F* hqqZZ_Unconditional_exp;
		TH2F* hZX_Unconditional_expk3a;
		TH2F* hqqZZ_Unconditional_expk3a;
		for (int t = 0; t < kNumTemplates; t++){
			if (t <= 2){ // gg(H)VV
				D_temp_2D[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
				D_temp_2D_exp[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
				D_temp_2D_expk3a[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];

				overall_scale[t] = new double[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
				for(int ty=0;ty<3;ty++) expectedNormalizations[ty][t] = new double[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
				for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][0]; al++){
					overall_scale[t][al] = 1;
					for(int ty=0;ty<3;ty++) expectedNormalizations[ty][t][al] = 0;
				}
			}
			else if (t >= 5 && t <= 7){ // VV(H)VV
				D_temp_2D[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
				D_temp_2D_exp[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
				D_temp_2D_expk3a[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];

				overall_scale[t] = new double[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
				for(int ty=0;ty<3;ty++) expectedNormalizations[ty][t] = new double[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
				for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++){
					overall_scale[t][al] = 1;
					for(int ty=0;ty<3;ty++) expectedNormalizations[ty][t][al] = 0;
				}
			}
			else{ // Anything else
				D_temp_2D[t] = new TH2F*[1];
				D_temp_2D_exp[t] = new TH2F*[1];
				D_temp_2D_expk3a[t] = new TH2F*[1];

				overall_scale[t] = new double[1];
				for(int ty=0;ty<3;ty++) expectedNormalizations[ty][t] = new double[1];
				for (int al = 0; al < 1; al++){
					overall_scale[t][al] = 1;
					for(int ty=0;ty<3;ty++) expectedNormalizations[ty][t][al] = 0;
				}
			}
		}

		TString cinput = cinput_common + INPUT_NAME;
		TFile* finput = new TFile(cinput,"read");
		TString cinput_k3a = cinput_common + INPUT_K3A_NAME;
		TFile* finput_k3a = new TFile(cinput_k3a,"read");

		TString cinput_smooth = cinput_common + INPUT_SMOOTH_NAME;
		TFile* finput_smooth = new TFile(cinput_smooth,"read");

		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if(t<=2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if(t>=5 && t<=7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			TString templatename_2D_exp_core;
			TString templatename_2D_core;
			TString templatename_2D_exp;
			TString templatename_2D;
			if (t < 2){
				templatename_2D_exp_core = Form("T_2D_%i", t + 1);
				templatename_2D_core = Form("T_2D_%i", t + 1);
			}
			else if (t == 2){
				templatename_2D_exp_core = Form("T_2D_%i", 4);
				templatename_2D_core = Form("T_2D_%i", 4);
			}
			else if (t == 3){
				templatename_2D_exp_core = Form("T_2D_%s", "qqZZ");
				templatename_2D_core = Form("T_2D_%s", "qqZZ");
			}
			else if (t == 4){
				templatename_2D_exp_core = Form("T_2D_%s", "ZX");
				templatename_2D_core = Form("T_2D_%s_merged", "ZX");
			}
			else if (t == 5 || t == 6){
				templatename_2D_exp_core = Form("T_2D_VBF_%i", t - 4);
				templatename_2D_core = Form("T_2D_VBF_%i", t - 4);
			}
			else if (t == 7){
				templatename_2D_exp_core = Form("T_2D_VBF_%i", 4);
				templatename_2D_core = Form("T_2D_VBF_%i", 4);
			}

			for (int al = 0; al<nAnomalousLoops; al++){
				templatename_2D = templatename_2D_core;
				templatename_2D_exp = templatename_2D_exp_core;
				if (useAnomalousCouplings == kAddfLQ && al>0){
					templatename_2D.Append(Form("_mZZ2_%i", al));
					templatename_2D_exp.Append(Form("_mZZ2_%i", al));
				}
				D_temp_2D_expk3a[t][al] = (TH2F*)finput_k3a->Get(templatename_2D_exp);
				D_temp_2D_exp[t][al] = (TH2F*)finput->Get(templatename_2D_exp);
				D_temp_2D[t][al] = (TH2F*)finput_smooth->Get(templatename_2D);
				D_temp_2D[t][al]->SetOption("colz");
				D_temp_2D_exp[t][al]->SetOption("colz");
				D_temp_2D_expk3a[t][al]->SetOption("colz");
			}

			if (t == 3){
				TString storeName = templatename_2D_core;
				storeName = storeName + "_UnConditional";
				hqqZZ_Unconditional = new TH2F(storeName, storeName, nbinsx, kDXarray, nbinsy, kDYarray);
				hqqZZ_Unconditional_exp = (TH2F*)finput->Get(storeName);
				hqqZZ_Unconditional_expk3a = (TH2F*)finput_k3a->Get(storeName);
				hqqZZ_Unconditional->SetOption("colz");
				hqqZZ_Unconditional_exp->SetOption("colz");
				hqqZZ_Unconditional_expk3a->SetOption("colz");
			}
			if (t == 4){
				TString storeName = templatename_2D_exp_core;
				storeName = storeName + "_UnConditional";
				hZX_Unconditional = new TH2F(storeName, storeName, nbinsx, kDXarray, nbinsy, kDYarray);
				hZX_Unconditional_exp = (TH2F*)finput->Get(storeName);
				hZX_Unconditional_expk3a = (TH2F*)finput_k3a->Get(storeName);
				hZX_Unconditional->SetOption("colz");
				hZX_Unconditional_exp->SetOption("colz");
				hZX_Unconditional_expk3a->SetOption("colz");

				int emptyBins[2] = { 0,0 };
				for (int binx = 1; binx <= D_temp_2D[t][0]->GetNbinsX(); binx++){
					double intBinX = D_temp_2D[t][0]->Integral(binx, binx, 0, D_temp_2D[t][0]->GetNbinsY()+1);
					if (intBinX < 5.0e-7 && emptyBins[0] != 0) emptyBins[1] = binx;
					else if (intBinX < 5.0e-7) emptyBins[0] = binx;
				}
				if (emptyBins[1] < emptyBins[0]) emptyBins[1] = emptyBins[0];
				if (emptyBins[0] != 0 && emptyBins[1] != 0){
          cout << "Empty bins: " << emptyBins[0] << '\t' << emptyBins[1] << endl;
          int minbinx = emptyBins[0];
					if (emptyBins[0]>1) minbinx = emptyBins[0]-1;
					int maxbinx = emptyBins[1];
					if (emptyBins[1]<D_temp_2D[t][0]->GetNbinsX()) maxbinx = emptyBins[1]+1;
					double nSharedBins = maxbinx - minbinx + 1;
					cout << "Min max bins: " << minbinx << '\t' << maxbinx << endl;

					for (int biny = 0; biny <= D_temp_2D[t][0]->GetNbinsY()+1; biny++){
						for (int binx = minbinx; binx <= maxbinx; binx++){
							double lowval = D_temp_2D[t][0]->GetBinContent(minbinx, biny);
							double highval = D_temp_2D[t][0]->GetBinContent(maxbinx, biny);
							cout << lowval << '\t' << highval << endl;
							double bincontent = (lowval + highval)/nSharedBins;
							D_temp_2D[t][0]->SetBinContent(binx, biny,bincontent);
						}
					}
				}
			}
			for (int al = 0; al < nAnomalousLoops; al++){
				for (int binx = 1; binx <= D_temp_2D[t][al]->GetNbinsX(); binx++){
					double binwidthx = D_temp_2D[t][al]->GetXaxis()->GetBinWidth(binx);
					for (int biny = 1; biny <= D_temp_2D[t][al]->GetNbinsY(); biny++){
						double binwidthy = D_temp_2D[t][al]->GetYaxis()->GetBinWidth(biny);
						double binwidth = binwidthx*binwidthy;
						double bincontent = D_temp_2D[t][al]->GetBinContent(binx, biny);
						if (t != 3 && t != 4) bincontent /= binwidth;
						D_temp_2D[t][al]->SetBinContent(binx, biny, bincontent);
						if (t == 3) hqqZZ_Unconditional->SetBinContent(binx, biny, bincontent / binwidth);
						if (t == 4) hZX_Unconditional->SetBinContent(binx, biny, bincontent / binwidth);
					}
				}
			}
			if (t == 4){
				hZX_Unconditional->Scale(ZX_yield[EnergyIndex][folder] / (hZX_Unconditional->Integral("width")*luminosity[EnergyIndex]));
			}
		}

		for (int t = 0; t<kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if(t<=2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if(t>=5 && t<=7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			for (int al = 0; al < nAnomalousLoops; al++){
				expectedNormalizations[0][t][al] = D_temp_2D[t][al]->Integral("width")*luminosity[EnergyIndex];
				expectedNormalizations[1][t][al] = D_temp_2D_exp[t][al]->Integral("width")*luminosity[EnergyIndex];
				expectedNormalizations[2][t][al] = D_temp_2D_expk3a[t][al]->Integral("width")*luminosity[EnergyIndex];
			}
		}

		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			for (int al = 0; al < nAnomalousLoops; al++){
				if (t == 2) expectedNormalizations[1][2][al] = expectedNormalizations[1][2][al] + expectedNormalizations[1][1][al] + expectedNormalizations[1][0][al];
				if (t == 2) expectedNormalizations[2][2][al] = expectedNormalizations[2][2][al] + expectedNormalizations[2][1][al] + expectedNormalizations[2][0][al];
				if (t < 3 || t == 6) overall_scale[t][al] = expectedNormalizations[2][t][al] / expectedNormalizations[0][t][al];
				else{
					if (t == 5){
						if (EnergyIndex == 1) overall_scale[t][al] = (expectedNormalizations[2][5][al] + expectedNormalizations[2][7][al] + expectedNormalizations[2][6][al]) / expectedNormalizations[0][t][al];
						else if (EnergyIndex == 0) overall_scale[t][al] = (25.0*expectedNormalizations[2][5][al] + 5.0*expectedNormalizations[2][7][al] + expectedNormalizations[2][6][al]) / expectedNormalizations[0][t][al];
					}
					if (t == 7) overall_scale[t][al] = (10.0*expectedNormalizations[2][5][al] + sqrt(10.0)*expectedNormalizations[2][7][al] + expectedNormalizations[2][6][al]) / expectedNormalizations[0][t][al];
				}
				cout << "Template " << t << " (anom. coupl. " << al << ") scaled rate: " << expectedNormalizations[0][t][al] * overall_scale[t][al]
					 << ", unscaled rate: " << expectedNormalizations[0][t][al]
					 << ", un-smoothened rate: " << expectedNormalizations[1][t][al] << endl;
				if (t<3 || t>4){
//					D_temp_2D[t][al]->Scale(overall_scale[t][al]);
					cout << "SCALE FOR " << D_temp_2D[t][al]->GetName() << " : " << overall_scale[t][al] << endl;
				}
				else{ // Conditionalize qqZZ and Z+X
					for (int binx = 0; binx <= D_temp_2D[t][al]->GetNbinsX()+1; binx++){
						double intBinX = D_temp_2D[t][al]->Integral(binx, binx, 0, D_temp_2D[t][al]->GetNbinsY()+1);
						for (int biny = 0; biny <= D_temp_2D[t][al]->GetNbinsY()+1; biny++){
							double bincontent = D_temp_2D[t][al]->GetBinContent(binx, biny);
							if (intBinX != 0) D_temp_2D[t][al]->SetBinContent(binx, biny, bincontent / intBinX);
						}
					}
				}
			}
		}

		for (int t = 0; t < kNumTemplates; t++){
			delete[] overall_scale[t];
			for (int ty = 0; ty < 3; ty++) delete[] expectedNormalizations[ty][t];
		}

		//Make VBF Sig/Int from linear combinations of above templates
		//5: VBF Sig
		//7: VBF Int
		//	 For 7 TeV samples, BSI25, Bkg, and BSI10 are used
		//	 For 8 TeV samples, BSI, Bkg, and BSI10 are used
		if (EnergyIndex == 0){
			if(tFitD!=0){
				for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++) twoDlinearcombination(D_temp_2D[7][al],kBSI10Hist,D_temp_2D[6][al],kBkgHist,D_temp_2D[5][al],kBSI25Hist,D_temp_2D[5][al],kSigHist,D_temp_2D[7][al],kIntHist);
			}
		}
		else if (EnergyIndex == 1){
			if(tFitD!=0){
				for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++) twoDlinearcombination(D_temp_2D[5][al],kBSIHist,D_temp_2D[6][al],kBkgHist,D_temp_2D[7][al],kBSI10Hist,D_temp_2D[5][al],kSigHist,D_temp_2D[7][al],kIntHist);
			}
		}

		if (tFitD != 0){
			for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][0]; al++){
				D_temp_2D[2][al]->Add(D_temp_2D[1][al], -1.0);
				D_temp_2D[2][al]->Add(D_temp_2D[0][al], -1.0);
			}
		}

		for (int t = 4; t < 5; t++){ // Special treatment for ZX
			if (Systematics != 0){
				for (int binx = 0; binx <= D_temp_2D[t][0]->GetNbinsX() + 1; binx++){
					double* storeOriginal = new double[D_temp_2D[t][0]->GetNbinsY() + 2];
					for (int biny = 0; biny <= D_temp_2D[t][0]->GetNbinsY() + 1; biny++){
						double bincontent = D_temp_2D[t][0]->GetBinContent(binx, biny);
						double bincontent_alt = D_temp_2D[3][0]->GetBinContent(binx, biny);
						double difference = bincontent_alt - bincontent;
						storeOriginal[biny] = bincontent;

						if (Systematics > 0) bincontent += difference;
						else bincontent -= difference;
						if (bincontent < 0) bincontent = 0;
						D_temp_2D[t][0]->SetBinContent(binx, biny, bincontent);
					}
					double intBinX = D_temp_2D[t][0]->Integral(binx, binx, 0, D_temp_2D[t][0]->GetNbinsY() + 1);
					for (int biny = 0; biny <= D_temp_2D[t][0]->GetNbinsY() + 1; biny++){
						double bincontent = D_temp_2D[t][0]->GetBinContent(binx, biny);

						if (intBinX != 0){
							D_temp_2D[t][0]->SetBinContent(binx, biny, bincontent / intBinX);
							double sysRatio = 1.0e-20;
							if (storeOriginal[biny] != 0) sysRatio = (bincontent / intBinX) / (storeOriginal[biny]);
							double unconditionalbincontent = hZX_Unconditional->GetBinContent(binx, biny);
							hZX_Unconditional->SetBinContent(binx, biny, unconditionalbincontent*sysRatio);
						}
						else{
							hZX_Unconditional->SetBinContent(binx, biny, 1.0e-20);
						}
					}
					delete[] storeOriginal;
				}
			}
			hZX_Unconditional_exp->Scale((hZX_Unconditional->Integral("width")) / (hZX_Unconditional_expk3a->Integral("width")));
			hZX_Unconditional_expk3a->Scale((hZX_Unconditional->Integral("width")) / (hZX_Unconditional_expk3a->Integral("width")));
		}

// Floor all templates since there is no class implementation for the PDF
		for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][0]; al++){
			floorSignalTemplates(D_temp_2D[0][al], D_temp_2D[1][al], D_temp_2D[2][al]);
			floorSignalTemplates(D_temp_2D_exp[0][al], D_temp_2D_exp[1][al], D_temp_2D_exp[2][al]);
			floorSignalTemplates(D_temp_2D_expk3a[0][al], D_temp_2D_expk3a[1][al], D_temp_2D_expk3a[2][al]);
		}
		for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++){
			floorSignalTemplates(D_temp_2D[5][al], D_temp_2D[6][al], D_temp_2D[7][al]);
			floorSignalTemplates(D_temp_2D_exp[5][al], D_temp_2D_exp[6][al], D_temp_2D_exp[7][al]);
			floorSignalTemplates(D_temp_2D_expk3a[5][al], D_temp_2D_expk3a[6][al], D_temp_2D_expk3a[7][al]);
		}
		for (int t=3;t<5;t++){
			floorBkgTemplate(D_temp_2D[t][0]);
			floorBkgTemplate(D_temp_2D_exp[t][0]);
			floorBkgTemplate(D_temp_2D_expk3a[t][0]);
		}
		floorBkgTemplate(hqqZZ_Unconditional);
		floorBkgTemplate(hqqZZ_Unconditional_exp);
		floorBkgTemplate(hqqZZ_Unconditional_expk3a);
		floorBkgTemplate(hZX_Unconditional);
		floorBkgTemplate(hZX_Unconditional_exp);
		floorBkgTemplate(hZX_Unconditional_expk3a);


		cout << "Integrals after everything:\nSmooth\tk3a\tRaw" << endl;
		for(int t=0;t<kNumTemplates;t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			TString oldTemplateName;
			for (int al = 0; al < nAnomalousLoops; al++){
				if (t == 4) D_temp_2D[t][al]->SetName(D_temp_2D_exp[t][al]->GetName());
				oldTemplateName = D_temp_2D_exp[t][al]->GetName();
				oldTemplateName = oldTemplateName + "_Raw";
				D_temp_2D_exp[t][al]->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t][al]->GetName();
				oldTemplateName = oldTemplateName + "_k3a";
				D_temp_2D_expk3a[t][al]->SetName(oldTemplateName);
			}

			if(t==3){
				oldTemplateName = D_temp_2D[t][0]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional";
				hqqZZ_Unconditional->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t][0]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_Raw";
				hqqZZ_Unconditional_exp->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t][0]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_k3a";
				hqqZZ_Unconditional_expk3a->SetName(oldTemplateName);
				foutput->WriteTObject(hqqZZ_Unconditional);
				foutput->WriteTObject(hqqZZ_Unconditional_exp);
				foutput->WriteTObject(hqqZZ_Unconditional_expk3a);
			}
			if(t==4){
				oldTemplateName = D_temp_2D[t][0]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional";
				hZX_Unconditional->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t][0]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_Raw";
				hZX_Unconditional_exp->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t][0]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_k3a";
				hZX_Unconditional_expk3a->SetName(oldTemplateName);
				foutput->WriteTObject(hZX_Unconditional);
				foutput->WriteTObject(hZX_Unconditional_exp);
				foutput->WriteTObject(hZX_Unconditional_expk3a);
			}

			bool* isPhysical = new bool[nAnomalousLoops];

			for (int al = 0; al < nAnomalousLoops; al++){
				isPhysical[al]=true;
				if (useAnomalousCouplings == kAddfLQ && al>0){ // Absorb coefficients due to powers of (mZZ/mH)**2 in fLQ
					double interfScale=1;

					if( (t==0||t==7) && al==1) interfScale = -2; // ggH signal or VBF interference, power 1
					else if( t==2 && al==1) interfScale = -1; // ggH interference, power 1
					else if( t==5 && (al==1 || al==3)) interfScale = -4; // VBF signal, power 1 or 3
					else if( t==5 && al==2) interfScale = 6; // VBF signal, power 2

					D_temp_2D[t][al]->Scale(interfScale);
					D_temp_2D_exp[t][al]->Scale(interfScale);
					D_temp_2D_expk3a[t][al]->Scale(interfScale);

					if(
						(t==2 && al>1) || // ggH interference does not contain terms higher than (mZZ/mH)**2
						(t==7 && al>2) || // VBF interference does not contain terms higher than pow( (mZZ/mH)**2 , 2 )
						((t==1||t==6) && al>0) // Bkg does not contain any (mZZ/mH)**2 terms
						) isPhysical[al]=false;
				}

				if(!isPhysical[al]) continue;

				foutput->WriteTObject(D_temp_2D[t][al]);
				foutput->WriteTObject(D_temp_2D_exp[t][al]);
				foutput->WriteTObject(D_temp_2D_expk3a[t][al]);

				if (t != 3 && t != 4){
					cout << "Template " << t << " (anom. coupl. " << al << "): ";
					cout << D_temp_2D[t][al]->Integral("width")*luminosity[EnergyIndex] << '\t';
					cout << D_temp_2D_expk3a[t][al]->Integral("width")*luminosity[EnergyIndex] << '\t';
					cout << D_temp_2D_exp[t][al]->Integral("width")*luminosity[EnergyIndex] << endl;
				}
				else{
					cout << "Template " << t << " (anom. coupl. " << al << "): ";
					cout << D_temp_2D[t][al]->Integral(0, D_temp_2D[t][al]->GetNbinsX()+1, 0, D_temp_2D[t][al]->GetNbinsY()+1) << '\t';
					cout << D_temp_2D_expk3a[t][al]->Integral(0, D_temp_2D_expk3a[t][al]->GetNbinsX()+1, 0, D_temp_2D_expk3a[t][al]->GetNbinsY()+1) << '\t';
					cout << D_temp_2D_exp[t][al]->Integral(0, D_temp_2D_exp[t][al]->GetNbinsX()+1, 0, D_temp_2D_exp[t][al]->GetNbinsY()+1) << endl;
				}
			}
			if(t==3){
				cout << "Template " << t << " unconditional (qqZZ): ";
				cout << hqqZZ_Unconditional->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hqqZZ_Unconditional_expk3a->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hqqZZ_Unconditional_exp->Integral("width")*luminosity[EnergyIndex] << endl;
			}
			else if(t==4){
				cout << "Template " << t << " unconditional (ZX): ";
				cout << hZX_Unconditional->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hZX_Unconditional_expk3a->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hZX_Unconditional_exp->Integral("width")*luminosity[EnergyIndex] << endl;
			}

			TString canvasname;
			TCanvas* cKD;
			TH1D* hPSmooth;
			TH1D* hPk3a;
			TH1D* hPraw;

			for (int al = 0; al < nAnomalousLoops; al++){
				if(!isPhysical[al]) continue;

				canvasname = "Compare_";
				canvasname = canvasname + D_temp_2D[t][al]->GetName();
				canvasname = canvasname + "_ZZMass";
				cKD = new TCanvas(canvasname);
				cKD->cd();
				hPSmooth = D_temp_2D[t][al]->ProjectionX();
				hPk3a = D_temp_2D_expk3a[t][al]->ProjectionX();
				hPraw = D_temp_2D_exp[t][al]->ProjectionX();
				hPSmooth->SetLineColor(kBlack);
				hPSmooth->SetLineWidth(2);
				hPSmooth->Draw("hist");
				hPraw->SetLineColor(kBlue);
				hPraw->SetLineWidth(2);
				hPraw->Draw("same");
				hPk3a->SetLineColor(kRed);
				hPk3a->SetLineWidth(2);
				hPk3a->Draw("same");
				cKD->SetLeftMargin(0.17);
				cKD->SetRightMargin(0.15);
				cKD->SetTopMargin(0.08);
				cKD->SetBottomMargin(0.17);
				cKD->Update();
				cKD->Modified();
				foutput->WriteTObject(cKD);
				delete hPSmooth;
				delete hPraw;
				delete hPk3a;
				cKD->Close();

				canvasname = "Compare_";
				canvasname = canvasname + D_temp_2D[t][al]->GetName();
				canvasname = canvasname + "_" + strFitDim[tFitD];
				cKD = new TCanvas(canvasname);
				cKD->cd();
				hPSmooth = D_temp_2D[t][al]->ProjectionY();
				hPk3a = D_temp_2D_expk3a[t][al]->ProjectionY();
				hPraw = D_temp_2D_exp[t][al]->ProjectionY();
				hPSmooth->SetLineColor(kBlack);
				hPSmooth->SetLineWidth(2);
				hPSmooth->Draw("hist");
				hPraw->SetLineColor(kBlue);
				hPraw->SetLineWidth(2);
				hPraw->Draw("same");
				hPk3a->SetLineColor(kRed);
				hPk3a->SetLineWidth(2);
				hPk3a->Draw("same");
				cKD->SetLeftMargin(0.17);
				cKD->SetRightMargin(0.15);
				cKD->SetTopMargin(0.08);
				cKD->SetBottomMargin(0.17);
				cKD->Update();
				cKD->Modified();
				foutput->WriteTObject(cKD);
				delete hPSmooth;
				delete hPraw;
				delete hPk3a;
				cKD->Close();
			}

			delete[] isPhysical;

			if(t==3 || t==4){
				canvasname = "Compare_";
				if(t==3) canvasname = canvasname + hqqZZ_Unconditional->GetName();
				if(t==4) canvasname = canvasname + hZX_Unconditional->GetName();
				canvasname = canvasname + "_ZZMass";
				cKD = new TCanvas(canvasname);
				cKD->cd();
				if(t==3){
					hPSmooth = hqqZZ_Unconditional->ProjectionX();
					hPk3a = hqqZZ_Unconditional_expk3a->ProjectionX();
					hPraw = hqqZZ_Unconditional_exp->ProjectionX();
				}
				else{
					hPSmooth = hZX_Unconditional->ProjectionX();
					hPk3a = hZX_Unconditional_expk3a->ProjectionX();
					hPraw = hZX_Unconditional_exp->ProjectionX();
				}
				hPSmooth->SetLineColor(kBlack);
				hPSmooth->SetLineWidth(2);
				hPSmooth->Draw("hist");
				hPraw->SetLineColor(kBlue);
				hPraw->SetLineWidth(2);
				hPraw->Draw("same");
				hPk3a->SetLineColor(kRed);
				hPk3a->SetLineWidth(2);
				hPk3a->Draw("same");
				cKD->SetLeftMargin(0.17);
				cKD->SetRightMargin(0.15);
				cKD->SetTopMargin(0.08);
				cKD->SetBottomMargin(0.17);
				cKD->Update();
				cKD->Modified();
				foutput->WriteTObject(cKD);
				delete hPSmooth;
				delete hPraw;
				delete hPk3a;
				cKD->Close();

				canvasname = "Compare_";
				if(t==3) canvasname = canvasname + hqqZZ_Unconditional->GetName();
				if(t==4) canvasname = canvasname + hZX_Unconditional->GetName();
				canvasname = canvasname + "_" + strFitDim[tFitD];
				cKD = new TCanvas(canvasname);
				cKD->cd();
				if(t==3){
					hPSmooth = hqqZZ_Unconditional->ProjectionY();
					hPk3a = hqqZZ_Unconditional_expk3a->ProjectionY();
					hPraw = hqqZZ_Unconditional_exp->ProjectionY();
				}
				else{
					hPSmooth = hZX_Unconditional->ProjectionY();
					hPk3a = hZX_Unconditional_expk3a->ProjectionY();
					hPraw = hZX_Unconditional_exp->ProjectionY();
				}
				hPSmooth->SetLineColor(kBlack);
				hPSmooth->SetLineWidth(2);
				hPSmooth->Draw("hist");
				hPraw->SetLineColor(kBlue);
				hPraw->SetLineWidth(2);
				hPraw->Draw("same");
				hPk3a->SetLineColor(kRed);
				hPk3a->SetLineWidth(2);
				hPk3a->Draw("same");
				cKD->SetLeftMargin(0.17);
				cKD->SetRightMargin(0.15);
				cKD->SetTopMargin(0.08);
				cKD->SetBottomMargin(0.17);
				cKD->Update();
				cKD->Modified();
				foutput->WriteTObject(cKD);
				delete hPSmooth;
				delete hPraw;
				delete hPk3a;
				cKD->Close();
			}
		}

		delete hqqZZ_Unconditional;
		delete hZX_Unconditional;
		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];
			for (int al = 0; al < nAnomalousLoops; al++){
				delete D_temp_2D_expk3a[t][al];
				delete D_temp_2D_exp[t][al];
				delete D_temp_2D[t][al];
			}
			delete[] D_temp_2D_expk3a[t];
			delete[] D_temp_2D_exp[t];
			delete[] D_temp_2D[t];
		}

		finput_smooth->Close();
		finput_k3a->Close();
		finput->Close();
		foutput->Close();
	}
}

double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype){
	double aa[3] = { first, c_first, sqrt(c_first) };
	double bb[3] = { second, c_second, sqrt(c_second) };

	double sig = (aa[0] * bb[2] - aa[2] * bb[0] + (aa[2] - bb[2])*bkg) / (aa[1] * bb[2] - aa[2] * bb[1]);
	double interf = (aa[0] * bb[1] - aa[1] * bb[0] + (aa[1] - bb[1])*bkg) / (aa[2] * bb[1] - aa[1] * bb[2]);
	if (sig <= 0){ sig = 1.0e-20; interf = 0; };

	if (outputtype == kIntHist) return interf;
	else if(outputtype==kSigHist) return sig;
	else return 0;
}

void twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, TH2F* finaloutput, int outputtype, TH2F* finaloutput2, int output2type){
	TH2F* output = (TH2F*) input->Clone();
	TH2F* output2 = (TH2F*) output->Clone();
	if(outputtype==kIntHist || outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
			if (outputtype == kSigHist){ delete output; output = (TH2F*) first->Clone(); }

			output2->Add(first, -25.0);
			output2->Add(second, -1.0);
			output2->Scale(0.2);
			if (output2type == kSigHist){ delete output2; output2 = (TH2F*) first->Clone(); }
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double bsi25 = output->GetBinContent(binx,biny);

					double weight=doLinearCombination(bsi25,25,bsi,1,bkg,outputtype);
					output->SetBinContent(binx,biny,weight);
					if (finaloutput2 != 0){
						double weight2 = doLinearCombination(bsi25, 25, bsi, 1, bkg, output2type);
						output2->SetBinContent(binx, biny, weight2);
					}
				}
			}
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			//double scaleval = 1./(-50. + 25.*sqrt(10.));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi10 = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double bsi25 = output->GetBinContent(binx,biny);

					double weight=doLinearCombination(bsi25,25,bsi10,10,bkg,outputtype);
					output->SetBinContent(binx,biny,weight);
					if (finaloutput2 != 0){
						double weight2 = doLinearCombination(bsi25, 25, bsi10, 10, bkg, output2type);
						output2->SetBinContent(binx, biny, weight2);
					}
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			//double scaleval = 1./(10 - sqrt(10));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double bsi10 = output->GetBinContent(binx,biny);

					double weight=doLinearCombination(bsi10,10,bsi,1,bkg,outputtype);
					output->SetBinContent(binx,biny,weight);
					if (finaloutput2 != 0){
						double weight2 = doLinearCombination(bsi10, 10, bsi, 1, bkg, output2type);
						output2->SetBinContent(binx, biny, weight2);
					}
				}
			}
		}
		for (int binx = 1; binx <= output->GetNbinsX(); binx++){
			for (int biny = 1; biny <= output->GetNbinsY(); biny++){
				finaloutput->SetBinContent(binx, biny, output->GetBinContent(binx, biny));
				if (finaloutput2 != 0) finaloutput2->SetBinContent(binx,biny,output2->GetBinContent(binx,biny));
			}
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	delete output;
	delete output2;
}

void floorSignalTemplates(TH2F* hsig, TH2F* hbkg, TH2F* hinterf){
	const int nbinsx = hsig->GetNbinsX();
	const int nbinsy = hsig->GetNbinsY();
	const double sign_interf = TMath::Sign(1.,hinterf->Integral());

	for (int binx = 1; binx <= nbinsx; binx++){
		for (int biny = 1; biny <= nbinsy; biny++){
			double sig = hsig->GetBinContent(binx,biny);
			double bkg = hbkg->GetBinContent(binx,biny);
			double interf = hinterf->GetBinContent(binx,biny);

			double sig_new=sig;
			double bkg_new=bkg;
			double interf_new=interf;

			if(sig<=0) sig_new=1.0e-10;
			if(bkg<=0) bkg_new=1.0e-10;
			if(sig<=0 || bkg<=0 || fabs(interf)<1.0e-10) interf_new = sign_interf*1.0e-10;
			if(pow(interf_new,2)>=(4*sig_new*bkg_new)) interf_new *= sqrt(4*sig_new*bkg_new/pow(interf_new,2))*0.9;

			hsig->SetBinContent(binx,biny,sig_new);
			hbkg->SetBinContent(binx,biny,bkg_new);
			hinterf->SetBinContent(binx,biny,interf_new);
		}
	}
}

void floorBkgTemplate(TH2F* hbkg){
	const int nbinsx = hbkg->GetNbinsX();
	const int nbinsy = hbkg->GetNbinsY();

	for (int binx = 1; binx <= nbinsx; binx++){
		for (int biny = 1; biny <= nbinsy; biny++){
			double bkg = hbkg->GetBinContent(binx,biny);
			double bkg_new=bkg;
			if(bkg<=0) bkg_new=1.0e-10;
			hbkg->SetBinContent(binx,biny,bkg_new);
		}
	}
}
