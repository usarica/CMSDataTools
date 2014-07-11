#include <iostream>
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
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"

using namespace std;

void makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool useResoVBF);

void makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF(){
	int systematics[5]={0,1,-1,2,-2};
	for(int i=0;i<5;++i){
		for(int CoM=7;CoM<9;++CoM){
			for(int channel=0;channel<3;++channel){
				makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],true);	
			}	
		}
	}
}

//Function to build one template
// folder = 0,1,2 (final state corresponds to 2e2mu,4mu,4e respectively)
// erg_tev = 7,8 (CoM energy)
// tFitD = [0,16] (choice of Discriminant, see FitDimensionsList.h for list; only tFitd works right now)
// Systematics = [-2,2] (Flag for systematics. 0=Nominal, +/-1=QCD, +/-2=PDF)
// useResoVBF = true/false (flag to use resolution smeared VBF samples, to be removed with fullsim samples)
void makeCombineTemplatesSmooth_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool useResoVBF){
	TString INPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_Raw_";
	TString INPUT_K3A_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_";
	TString INPUT_SMOOTH_NAME = "HtoZZ4l_MCFM_125p6_SmoothTemplates_";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine_";

	OUTPUT_NAME += "_GenLevelVBF_";
	if(useResoVBF) OUTPUT_NAME += "wResolution_";
	OUTPUT_NAME += TString(strFitDim[tFitD]) + "_";
	TString comstring;
	comstring.Form("%i",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/",erg_tev);
	if(Systematics==0) OUTPUT_NAME += "Nominal.root";
	if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF.root";

	INPUT_SMOOTH_NAME += "_GenLevelVBF_";
	if(useResoVBF) INPUT_SMOOTH_NAME += "wResolution_";
	INPUT_SMOOTH_NAME += TString(strFitDim[tFitD]) + "_";
	if(Systematics==0) INPUT_SMOOTH_NAME += "Nominal.root";
	if (Systematics == 1) INPUT_SMOOTH_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) INPUT_SMOOTH_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) INPUT_SMOOTH_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) INPUT_SMOOTH_NAME += "SysDown_ggPDF.root";

	INPUT_NAME += "_GenLevelVBF_";
	if(useResoVBF) INPUT_NAME += "wResolution_";
	INPUT_NAME += TString(strFitDim[tFitD]) + "_";
	if(Systematics==0) INPUT_NAME += "Nominal.root";
	if (Systematics == 1) INPUT_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) INPUT_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) INPUT_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) INPUT_NAME += "SysDown_ggPDF.root";

	INPUT_K3A_NAME += "_GenLevelVBF_";
	if(useResoVBF) INPUT_K3A_NAME += "wResolution_";
	INPUT_K3A_NAME += TString(strFitDim[tFitD]) + "_";
	if(Systematics==0) INPUT_K3A_NAME += "Nominal.root";
	if (Systematics == 1) INPUT_K3A_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) INPUT_K3A_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) INPUT_K3A_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) INPUT_K3A_NAME += "SysDown_ggPDF.root";

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float lowside[3]={220,230,240};
	double ZX_yield[2][3]={
		{ 0.1078, 0.2213, 0.3345 },
		{0.55,1.78,1.38}
	};

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
		};
	
		int nbinsy=30;
		float kDYarray[nbinsy+1];
		float kDY_bounds[2]={0,1};
		if(tFitD==3){kDY_bounds[0]=-7.0;kDY_bounds[1]=3.0;};
		if(tFitD==4 || tFitD==5 || tFitD==9 || tFitD==10 || tFitD==11 || tFitD==15){kDY_bounds[0]=-1.0;kDY_bounds[1]=1.0;};
		for(int bin=0;bin<nbinsy+1;bin++){
			double binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		};

		const int kNumTemplates=8;
		double overall_scale[kNumTemplates]={1};
		double expectedNormalizations[3][kNumTemplates]={{0}};
		TH2F* D_temp_2D_expk3a[kNumTemplates];
		TH2F* D_temp_2D_exp[kNumTemplates];
		TH2F* D_temp_2D[kNumTemplates];
		TH2F* hZX_Unconditional;
		TH2F* hqqZZ_Unconditional;
		TH2F* hZX_Unconditional_exp;
		TH2F* hqqZZ_Unconditional_exp;
		TH2F* hZX_Unconditional_expk3a;
		TH2F* hqqZZ_Unconditional_expk3a;

		TString cinput = cinput_common + INPUT_NAME;
		TFile* finput = new TFile(cinput,"read");
		TString cinput_k3a = cinput_common + INPUT_K3A_NAME;
		TFile* finput_k3a = new TFile(cinput_k3a,"read");

		TString cinput_smooth = cinput_common + INPUT_SMOOTH_NAME;
		TFile* finput_smooth = new TFile(cinput_smooth,"read");

		for (int t = 0; t < kNumTemplates; t++){
			char templatename_2D_exp[100];
			char templatename_2D[100];
			if (t < 2){
				sprintf(templatename_2D_exp, "T_2D_%i", t + 1);
				sprintf(templatename_2D, "T_2D_%i", t + 1);
			}
			else if (t == 2){
				sprintf(templatename_2D_exp, "T_2D_%i", 4);
				sprintf(templatename_2D, "T_2D_%i", 4);
			}
			else if (t == 3){
				sprintf(templatename_2D_exp, "T_2D_%s", "qqZZ");
				sprintf(templatename_2D, "T_2D_%s", "qqZZ");
			}
			else if (t == 4){
				sprintf(templatename_2D_exp, "T_2D_%s", "ZX");
				sprintf(templatename_2D, "T_2D_%s_merged", "ZX");
			}
			else if (t == 5 || t == 6){
				sprintf(templatename_2D_exp, "T_2D_VBF_%i", t - 4);
				sprintf(templatename_2D, "T_2D_VBF_%i", t - 4);
			}
			else if (t == 7){
				sprintf(templatename_2D_exp, "T_2D_VBF_%i", 4);
				sprintf(templatename_2D, "T_2D_VBF_%i", 4);
			};

			D_temp_2D_expk3a[t] = (TH2F*)finput_k3a->Get(templatename_2D_exp);
			D_temp_2D_exp[t] = (TH2F*)finput->Get(templatename_2D_exp);
			D_temp_2D[t] = (TH2F*)finput_smooth->Get(templatename_2D);
			D_temp_2D[t]->SetOption("colz");
			D_temp_2D_exp[t]->SetOption("colz");
			D_temp_2D_expk3a[t]->SetOption("colz");

			if (t == 3){
				TString storeName = templatename_2D;
				storeName = storeName + "_UnConditional";
				hqqZZ_Unconditional = new TH2F(storeName, storeName, nbinsx, kDXarray, nbinsy, kDYarray);
				hqqZZ_Unconditional_exp = (TH2F*)finput->Get(storeName);
				hqqZZ_Unconditional_expk3a = (TH2F*)finput_k3a->Get(storeName);
				hqqZZ_Unconditional->SetOption("colz");
				hqqZZ_Unconditional_exp->SetOption("colz");
				hqqZZ_Unconditional_expk3a->SetOption("colz");
			};
			if (t == 4){
				TString storeName = templatename_2D_exp;
				storeName = storeName + "_UnConditional";
				hZX_Unconditional = new TH2F(storeName, storeName, nbinsx, kDXarray, nbinsy, kDYarray);
				hZX_Unconditional_exp = (TH2F*)finput->Get(storeName);
				hZX_Unconditional_expk3a = (TH2F*)finput_k3a->Get(storeName);
				hZX_Unconditional->SetOption("colz");
				hZX_Unconditional_exp->SetOption("colz");
				hZX_Unconditional_expk3a->SetOption("colz");

				int emptyBins[2] = { 0,0 };
				for (int binx = 1; binx <= D_temp_2D[t]->GetNbinsX(); binx++){
					double intBinX = D_temp_2D[t]->Integral(binx, binx, 1, D_temp_2D[t]->GetNbinsY());
					if (intBinX < 5.0e-7 && emptyBins[0] != 0) emptyBins[1] = binx;
					else if (intBinX < 5.0e-7) emptyBins[0] = binx;
				};
				if (emptyBins[1] < emptyBins[0]) emptyBins[1] = emptyBins[0];
				cout << "Empty bins: " << emptyBins[0] << '\t' << emptyBins[1] << endl;
				if (emptyBins[0] != 0 && emptyBins[1] != 0){
					int minbinx = emptyBins[0];
					if (emptyBins[0]>1) minbinx = emptyBins[0]-1;
					int maxbinx = emptyBins[1];
					if (emptyBins[1]<D_temp_2D[t]->GetNbinsX()) maxbinx = emptyBins[1]+1;
					double nSharedBins = maxbinx - minbinx + 1;
					cout << "Min max bins: " << minbinx << '\t' << maxbinx << endl;

					for (int biny = 1; biny <= D_temp_2D[t]->GetNbinsY(); biny++){
						for (int binx = minbinx; binx <= maxbinx; binx++){
							double lowval = D_temp_2D[t]->GetBinContent(minbinx, biny);
							double highval = D_temp_2D[t]->GetBinContent(maxbinx, biny);
							cout << lowval << '\t' << highval << endl;
							double bincontent = (lowval + highval)/nSharedBins;
							D_temp_2D[t]->SetBinContent(binx, biny,bincontent);
						};
					};
				};

			};

			for (int binx = 1; binx <= D_temp_2D[t]->GetNbinsX(); binx++){
				double binwidthx = D_temp_2D[t]->GetXaxis()->GetBinWidth(binx);
				for (int biny = 1; biny <= D_temp_2D[t]->GetNbinsY(); biny++){
					double binwidthy = D_temp_2D[t]->GetYaxis()->GetBinWidth(biny);
					double binwidth = binwidthx*binwidthy;
					double bincontent = D_temp_2D[t]->GetBinContent(binx, biny);
					if (t != 3 && t != 4) bincontent /= binwidth;
					D_temp_2D[t]->SetBinContent(binx, biny, bincontent);
					if (t == 3){
						hqqZZ_Unconditional->SetBinContent(binx, biny, bincontent / binwidth);
					};
					if (t == 4){
						hZX_Unconditional->SetBinContent(binx, biny, bincontent / binwidth);
					};
				};
			};
			if (t == 4){
				hZX_Unconditional->Scale(ZX_yield[EnergyIndex][folder] / (hZX_Unconditional->Integral("width")*luminosity[EnergyIndex]));
			};
		};
		if (EnergyIndex == 0){ // USE BSI25, NO CONTAMINATION SIGNIFICANT TO CREATE STORAGE TREES
			if (folder == 2){
				D_temp_2D[2]->Add(D_temp_2D[1], -1.0);
				D_temp_2D[2]->Add(D_temp_2D[0], -25.0);
				D_temp_2D[2]->Scale(0.2);
				D_temp_2D[2]->Add(D_temp_2D[1], 1.0);
				D_temp_2D[2]->Add(D_temp_2D[0], 1.0);
			};
			if (folder == 1){
				for (int binx = 1; binx <= D_temp_2D[1]->GetNbinsX(); binx++){
					for (int biny = 1; biny <= D_temp_2D[1]->GetNbinsY(); biny++){
						double sig = D_temp_2D[0]->GetBinContent(binx, biny);
						double bsi25 = D_temp_2D[1]->GetBinContent(binx, biny);
						double bsi = D_temp_2D[2]->GetBinContent(binx, biny);

						bsi = 5.0*(bsi - sig);
						bsi25 = (bsi25 - 25.0*sig - bsi)*(-0.25);
						if (bsi25 < 0) bsi25 = 0;
						D_temp_2D[1]->SetBinContent(binx, biny, bsi25);
					};
				};
			};
		};

		for (int t = 0; t<kNumTemplates; t++){
			expectedNormalizations[0][t] = D_temp_2D[t]->Integral("width")*luminosity[EnergyIndex];
			expectedNormalizations[1][t] = D_temp_2D_exp[t]->Integral("width")*luminosity[EnergyIndex];
			expectedNormalizations[2][t] = D_temp_2D_expk3a[t]->Integral("width")*luminosity[EnergyIndex];
		};
		for(int t=0;t<kNumTemplates;t++){
			if(t==2) expectedNormalizations[1][2] = expectedNormalizations[1][2] + expectedNormalizations[1][1] + expectedNormalizations[1][0];
			if(t==2) expectedNormalizations[2][2] = expectedNormalizations[2][2] + expectedNormalizations[2][1] + expectedNormalizations[2][0];
			if(t<3 || t==6) overall_scale[t] = expectedNormalizations[2][t]/expectedNormalizations[0][t];
			else{
				if(t==5) overall_scale[t] = (25.0*expectedNormalizations[2][5] + 5.0*expectedNormalizations[2][7] + expectedNormalizations[2][6])/expectedNormalizations[0][t];
				if(t==7) overall_scale[t] = (10.0*expectedNormalizations[2][5] + sqrt(10.0)*expectedNormalizations[2][7] + expectedNormalizations[2][6])/expectedNormalizations[0][t];
			};
			cout << expectedNormalizations[0][t]*overall_scale[t] << '\t' << expectedNormalizations[0][t] << '\t' << expectedNormalizations[1][t] << endl;
			if(t<3 || t>4){
				D_temp_2D[t]->Scale(overall_scale[t]);
				cout << "SCALE FOR " << D_temp_2D[t]->GetName() << " : " << overall_scale[t] << endl;
			}
			else{
				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);
						if(intBinX!=0) D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
					};
				};
			};
		};
		for(int binx=1;binx<=nbinsx;binx++){
			for(int biny=1;biny<=nbinsy;biny++){
				double bsi25 = D_temp_2D[5]->GetBinContent(binx,biny);
				double bkg = D_temp_2D[6]->GetBinContent(binx,biny);
				double bsi10 = D_temp_2D[7]->GetBinContent(binx,biny);
				double bsi25p = bsi25 - bkg;
				double bsi10p = bsi10 - bkg;
				double signal = ( sqrt(10.0)*bsi25p - 5.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
				double interf = ( -10.0*bsi25p + 25.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
				if(signal<=0){
					signal=1.0e-20;
					interf=0;
				};
				if(bkg<0) bkg=0;

				D_temp_2D[5]->SetBinContent(binx,biny,signal);
				D_temp_2D[6]->SetBinContent(binx,biny,bkg);
				D_temp_2D[7]->SetBinContent(binx,biny,interf);
			};
		};

		if(tFitD!=0){
			D_temp_2D[2]->Add(D_temp_2D[1],-1.0);
			D_temp_2D[2]->Add(D_temp_2D[0],-1.0);
		};

		cout << "Integrals after everything:\nSmooth\tk3a\tRaw" << endl;
		for(int t=0;t<kNumTemplates;t++){
			if(Systematics!=0 && t==4){
				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double* storeOriginal = new double[D_temp_2D[t]->GetNbinsY()];
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);
						double bincontent_alt = D_temp_2D[3]->GetBinContent(binx,biny);
						double difference = bincontent_alt - bincontent;
						storeOriginal[biny-1] = bincontent;

						if(Systematics>0) bincontent += difference;
						else bincontent -= difference;
						if(bincontent <= 0) bincontent = 1.0e-20;
						D_temp_2D[t]->SetBinContent(binx,biny,bincontent);
					};
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);

						if(intBinX!=0){
							D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
							double sysRatio = 1.0e-20;
							if(storeOriginal[biny-1]!=0) sysRatio = (bincontent/intBinX)/(storeOriginal[biny-1]);
							double unconditionalbincontent = hZX_Unconditional->GetBinContent(binx,biny);
							hZX_Unconditional->SetBinContent(binx,biny,unconditionalbincontent*sysRatio);
						}
						else{
							hZX_Unconditional->SetBinContent(binx,biny,1.0e-20);
						};
					};
					delete[] storeOriginal;
				};
			};
			if(t==4){
				hZX_Unconditional_exp->Scale((hZX_Unconditional->Integral("width"))/(hZX_Unconditional_expk3a->Integral("width")));
				hZX_Unconditional_expk3a->Scale((hZX_Unconditional->Integral("width"))/(hZX_Unconditional_expk3a->Integral("width")));
			};
			if(t==4) D_temp_2D[t]->SetName(D_temp_2D_exp[t]->GetName());
			TString oldTemplateName = D_temp_2D_exp[t]->GetName();
			oldTemplateName = oldTemplateName + "_Raw";
			D_temp_2D_exp[t]->SetName(oldTemplateName);
			oldTemplateName = D_temp_2D[t]->GetName();
			oldTemplateName = oldTemplateName + "_k3a";
			D_temp_2D_expk3a[t]->SetName(oldTemplateName);

			if(t==3){
				oldTemplateName = D_temp_2D[t]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional";
				hqqZZ_Unconditional->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_Raw";
				hqqZZ_Unconditional_exp->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_k3a";
				hqqZZ_Unconditional_expk3a->SetName(oldTemplateName);
				foutput->WriteTObject(hqqZZ_Unconditional);
				foutput->WriteTObject(hqqZZ_Unconditional_exp);
				foutput->WriteTObject(hqqZZ_Unconditional_expk3a);
			};
			if(t==4){
				oldTemplateName = D_temp_2D[t]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional";
				hZX_Unconditional->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_Raw";
				hZX_Unconditional_exp->SetName(oldTemplateName);
				oldTemplateName = D_temp_2D[t]->GetName();
				oldTemplateName = oldTemplateName + "_UnConditional_k3a";
				hZX_Unconditional_expk3a->SetName(oldTemplateName);
				foutput->WriteTObject(hZX_Unconditional);
				foutput->WriteTObject(hZX_Unconditional_exp);
				foutput->WriteTObject(hZX_Unconditional_expk3a);
			};

			foutput->WriteTObject(D_temp_2D[t]);
			foutput->WriteTObject(D_temp_2D_exp[t]);
			foutput->WriteTObject(D_temp_2D_expk3a[t]);

			if(t!=3 || t!=4){
				cout << D_temp_2D[t]->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << D_temp_2D_expk3a[t]->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << D_temp_2D_exp[t]->Integral("width")*luminosity[EnergyIndex] << endl;
			}
			else{
				cout << D_temp_2D[t]->Integral(1,nbinsx,1,nbinsy) << '\t';
				cout << D_temp_2D_expk3a[t]->Integral(1,nbinsx,1,nbinsy) << '\t';
				cout << D_temp_2D_exp[t]->Integral(1,nbinsx,1,nbinsy) << endl;
			};
			if(t==3){
				cout << hqqZZ_Unconditional->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hqqZZ_Unconditional_expk3a->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hqqZZ_Unconditional_exp->Integral("width")*luminosity[EnergyIndex] << endl;
			};
			if(t==4){
				cout << hZX_Unconditional->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hZX_Unconditional_expk3a->Integral("width")*luminosity[EnergyIndex] << '\t';
				cout << hZX_Unconditional_exp->Integral("width")*luminosity[EnergyIndex] << endl;
			};


			TString canvasname = "Compare_";
			canvasname = canvasname + D_temp_2D[t]->GetName();
			canvasname = canvasname + "_ZZMass";
			TCanvas* cKD = new TCanvas(canvasname);
			cKD->cd();
			TH1D* hPSmooth = D_temp_2D[t]->ProjectionX();
			TH1D* hPk3a = D_temp_2D_expk3a[t]->ProjectionX();
			TH1D* hPraw = D_temp_2D_exp[t]->ProjectionX();
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
				};
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
			};

			canvasname = "Compare_";
			canvasname = canvasname + D_temp_2D[t]->GetName();
			canvasname = canvasname + "_" + strFitDim[tFitD];
			cKD = new TCanvas(canvasname);
			cKD->cd();
			hPSmooth = D_temp_2D[t]->ProjectionY();
			hPk3a = D_temp_2D_expk3a[t]->ProjectionY();
			hPraw = D_temp_2D_exp[t]->ProjectionY();
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

			if(t==3 || t==4){
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
				};
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
			};
		};
		delete hqqZZ_Unconditional;
		delete hZX_Unconditional;
		finput_smooth->Close();
		finput_k3a->Close();
		finput->Close();
		foutput->Close();
	};
};
