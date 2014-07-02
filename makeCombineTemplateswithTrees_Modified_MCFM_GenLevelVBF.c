#include <iostream>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"
#include "./data/HZZ4l_LeptonInterference.h"

using namespace std;

//Initializers
void makeCombineTemplates_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool useResoVBF, bool isSmooth);

//Make Lepton Interference Graph to be used later
TGraph* make_HZZ_LeptonInterferenceGraph(){
	string cinput="./data/HZZ_LeptonInterferenceGraph.root";
	float x[leptonInterf_YR3_Size];
	float y[leptonInterf_YR3_Size];
	for(int a=0;a<leptonInterf_YR3_Size;a++){
		x[a] = leptonInterf_YR3[a][0];
		y[a] = leptonInterf_YR3[a][1];
	};
	TGraph* tg = new TGraph(leptonInterf_YR3_Size,x,y);
	tg->SetName("tgHZZ_LeptonInterference");
	tg->SetTitle("H#rightarrowZZ#rightarrow4l Lepton Interference Weight on 4e, 4#mu wrt. 2e2#mu");

	return tg;
};

//Main Function, runs over all desired iterations
void makeCombineTemplates_Modified_MCFM_GenLevelVBF(){
	bool isSmooth=false;
	int systematics[5]={0,1,-1,2,-2};
	for(int i=0;i<5;++i){
		for(int usesmooth=0;usesmooth<2;++usesmooth){
			if(usesmooth==0) isSmooth=false;
			if(usesmooth==1) isSmooth=true;
			for(int CoM=7;CoM<9;++CoM){
				for(int channel=0;channel<3;++channel){
					makeCombineTemplates_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],true,isSmooth);	
				}	
			}
		}
	}
	//makeCombineTemplates_Modified_MCFM_GenLevelVBF_one(2,7,6,0,true,false);			
}

//Function to build one template
// folder = 0,1,2 (final state corresponds to 2e2mu,4mu,4e respectively)
// erg_tev = 7,8 (CoM energy)
// tFitD = [0,16] (choice of Discriminant, see FitDimensionsList.h for list; only tFitd works right now)
// Systematics = [-2,2] (Flag for systematics. 0=Nominal, +/-1=QCD, +/-2=PDF)
// useResoVBF = true/false (flag to use resolution smeared VBF samples, to be removed with fullsim samples)
// isSmooth = true/false (flag to apply smoothing at this stage, both needed for later stage)
void makeCombineTemplates_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool useResoVBF, bool isSmooth){
	char TREE_NAME[]="SelectedTree";
	TString INPUT_SM_NAME = "HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Reprocessed.root";
	TString INPUT_NAME = "HZZ4lTree_ggTo";
	TString INPUT_NAME_VBF = "HZZ4l-125_6-";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_";
	if(!isSmooth) OUTPUT_NAME += "Raw_";
	OUTPUT_NAME += "_GenLevelVBF_";
	if(useResoVBF) OUTPUT_NAME += "wResolution_";
	OUTPUT_NAME += TString(strFitDim[tFitD]) + "_";
	TString comstring;
	comstring.Form("%i",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV",erg_tev);

	if (Systematics == 0) OUTPUT_NAME += "Nominal.root";
	if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF.root";

	TString INPUT_VBFREF_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	if(!isSmooth) INPUT_VBFREF_NAME += "Raw_";
	INPUT_VBFREF_NAME += TString(strFitDim[tFitD]) + "_";
	INPUT_VBFREF_NAME += "Nominal.root";
//	if (Systematics == 0) INPUT_VBFREF_NAME += "Nominal.root";
//	if(Systematics==1) INPUT_VBFREF_NAME += "SysUp.root";
//	if(Systematics==-1) INPUT_VBFREF_NAME += "SysDown.root";

	INPUT_NAME_VBF += comstring + "TeV-";
	TString sample_VBF_suffix[4]={
		"Bkg_VBF_Phantom",
		"BSI10_VBF_Phantom",
		"BSI25_VBF_Phantom",
		"BSI_VBF_Phantom"
	};

	TString cinput_common = user_gg2VV_location + erg_dir + "/" + user_folder[folder] + "/";
	TString cinput_common_recover;
	if(folder!=2) cinput_common_recover = user_gg2VV_location + erg_dir + "/" + user_folder[2] + "/";
	TString cinput_common_qqZZ = user_gg2VV_location;
	TString cinput_common_ZX = user_gg2VV_location;
	TString cinput_common_VBF = user_gg2VV_location;

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float lowside[3]={220,230,240};
	float ZZMass_PeakCut[2]={120,130}; // Spin 0 analysis
	double ggZZ_Syst_AbsNormSyst[2][2] = { // EnergyIndex
		{ 0.0745, 0.0735 },
		{ 0.075, 0.072 }
	}; // QCD, PDF

	int isSelected;
	float templateWeight=1;
	float MC_weight;
	float MC_weight_down;
	float MC_weight_up;
	float MC_weight_Kfactor=1;
	float MC_weight_ggZZLepInt=1;
	float GenHMass=0;
	float MC_weight_Kfactor_Norm_down = 1;
	float MC_weight_Kfactor_Norm_up = 1;
	float MC_weight_PDF_Norm_down = 1;
	float MC_weight_PDF_Norm_up = 1;
	float ZZMass = 0;
	float Z1Mass=0;
	float Z2Mass=0;

	TString cinput_KDFactor = "./data/HZZ4l-KDFactorGraph";
	if (EnergyIndex == 0) cinput_KDFactor = cinput_KDFactor + "_7TeV";
	cinput_KDFactor = cinput_KDFactor + ".root";
	TFile* finput_KDFactor = new TFile(cinput_KDFactor, "read");
	TString tgkfname = "KDScale_";
	tgkfname = tgkfname + "AllFlavors_UnNormalized";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*)finput_KDFactor->Get(tgkfname);

	TString cinput_sm = cinput_common + INPUT_SM_NAME;
	TFile* fsm = new TFile(cinput_sm,"read");
	TTree* tsm = (TTree*) fsm->Get(TREE_NAME);
	tsm->SetBranchAddress("ZZMass",&ZZMass);
	tsm->SetBranchAddress("MC_weight",&MC_weight);
	double nTotalSM=0;
	for(int sev=0;sev<tsm->GetEntries();sev++){
		tsm->GetEntry(sev);
		if(ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]) nTotalSM += MC_weight;
	};
	fsm->Close();

	double luminosity[2] = { 5.051, 19.712 };
	//ggH yields
	double nSM_ScaledPeak[2][3]={
		{1.4452079,0.6087736,1.0902689},
		{6.6562629,2.6944639,5.1963998}
	}
	for (int e = 0; e < 2; e++){for (int ss = 0; ss < 3; ss++){ nSM_ScaledPeak[e][ss] /= luminosity[e];}}

	//cout << "Observed number of peak events is " << nSM_ObservedPeak*luminosity[EnergyIndex] << endl;
	//cout << "Scaled number of peak events is " << nSM_ScaledPeak*luminosity[EnergyIndex] << endl;
	TGraph* tg_interf = make_HZZ_LeptonInterferenceGraph();

	TString cinput_VBF_Sig = "./data/HZZ4l-125_6-" + comstring + "TeV-Sig_MCFM_PhantomVBF_Comparison.root";
	TFile* finput_VBF = new TFile(cinput_VBF_Sig,"read");
	TSpline3* tsp_VBF_Sig = (TSpline3*) finput_VBF->Get("Spline3");

	TH2F* hStore_ZX_Unconditional;
	TH2F* hStore_qqZZ_Unconditional;

	for(int lo=0;lo<1;lo++){
		TString user_core = user_dir;
		TString coutput_common = user_core + erg_dir;
		coutput_common = coutput_common + "/Analysis/Templates/Combine/" + user_folder[folder] + "/";

		TString coutput = coutput_common + OUTPUT_NAME;
		TFile* foutput = new TFile(coutput,"recreate");

		//Grab VBF Sig templates made in makeCombineTemplates_Modified_MCFM.c
		TString cinput_VBFRef = coutput_common + INPUT_VBFREF_NAME;
		TFile* finput_VBFRef = new TFile(cinput_VBFRef,"read");
		double nVBF_Sig_Reweighted=0;
		TH1F* hVBFRef_1D = (TH1F*) finput_VBFRef->Get("T_1D_VBF_1");
		TH2F* hVBFRef = (TH2F*) finput_VBFRef->Get("T_2D_VBF_1");
		nVBF_Sig_Reweighted = hVBFRef->Integral("width")*luminosity[EnergyIndex];
		foutput->cd();

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
		TH1F* D_temp_1D[kNumTemplates];
		TH2F* D_temp_2D[kNumTemplates];
		double overall_scale[kNumTemplates]={1};

		int nEntries;
		double nSig_Simulated=0;
		double nVBF_Sig_Simulated=0;
		float fitYval;
		isSelected=1;
		MC_weight=1;
		MC_weight_down=1;
		MC_weight_up=1;
		MC_weight_Kfactor=1;
		MC_weight_Kfactor_Norm_down = 1;
		MC_weight_Kfactor_Norm_up = 1;
		MC_weight_PDF_Norm_down = 1;
		MC_weight_PDF_Norm_up = 1;
		MC_weight_ggZZLepInt = 1;
		GenHMass=0;
		ZZMass=0;
		Z1Mass=0;
		Z2Mass=0;

		double nMZZ220[3]={0};

		//Initialize and grab each of the four Phantom trees (change to FullSim)
		TChain* tree_VBF[4];
		TH1F* h1DVBF[4];
		TH2F* h2DVBF[4];
		TH1F* h1DVBFSigRatio;
		TH2F* h2DVBFSigRatio;
		if(tFitD==0){
			h1DVBFSigRatio = new TH1F("h1DSigRatio","h1DSigRatio",nbinsx,kDXarray);
		}
		else{
			h1DVBFSigRatio = new TH1F("h1DSigRatio","h1DSigRatio",nbinsy,kDYarray);
			h2DVBFSigRatio = new TH2F("h2DSigRatio","h2DSigRatio",nbinsx,kDXarray,nbinsy,kDYarray);
		};
		for(int tr=0;tr<4;tr++){
			tree_VBF[tr] = new TChain(TREE_NAME);
		};
		TString cinput_VBF_Bkg = cinput_common_VBF + "../VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[0];
		TString cinput_VBF_BSI = cinput_common_VBF + "../VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[3];
		TString cinput_VBF_BSI10 = cinput_common_VBF + "../VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[1];
		TString cinput_VBF_BSI25 = cinput_common_VBF + "../VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[2];
		if(useResoVBF){
			cinput_VBF_Bkg += "_Reprocessed_wResolution.root";
			cinput_VBF_BSI += "_Reprocessed_wResolution.root";
			cinput_VBF_BSI10 += "_Reprocessed_wResolution.root";
			cinput_VBF_BSI25 += "_Reprocessed_wResolution.root";
		}
		else{
			cinput_VBF_Bkg += "_Reprocessed.root";
			cinput_VBF_BSI += "_Reprocessed.root";
			cinput_VBF_BSI10 += "_Reprocessed.root";
			cinput_VBF_BSI25 += "_Reprocessed.root";
		};
		tree_VBF[0]->Add(cinput_VBF_Bkg);
		tree_VBF[1]->Add(cinput_VBF_BSI10);
		tree_VBF[2]->Add(cinput_VBF_BSI25);
		tree_VBF[3]->Add(cinput_VBF_BSI);

		//Fill Template(s) from each Phantom sample with:
		// gen mZ1/mZ2>120
		// Reweighting for lepton interference (necessary)
		// Possible mZZ cuts?
		//Then write smoothed templates to file 
		for(int tr=0;tr<4;tr++){
			tree_VBF[tr]->SetBranchAddress("isSelected",&isSelected);
			tree_VBF[tr]->SetBranchAddress("GenZZMass", &ZZMass);
			tree_VBF[tr]->SetBranchAddress("GenZ1Mass", &Z1Mass);
			tree_VBF[tr]->SetBranchAddress("GenZ2Mass", &Z2Mass);
			tree_VBF[tr]->SetBranchAddress("MC_weight", &MC_weight);
			if(tFitD!=0) tree_VBF[tr]->SetBranchAddress(strFitDim[tFitD],&fitYval);

			TString templatename_1D = "htemp1DVBF_";
			templatename_1D += sample_VBF_suffix[tr];
			TString templatename_2D = "htemp2DVBF_";
			templatename_2D += sample_VBF_suffix[tr];

			if(tFitD==0){
				h1DVBF[tr] = new TH1F(templatename_1D,templatename_1D,nbinsx,kDXarray);
			}
			else{
				h1DVBF[tr] = new TH1F(templatename_1D,templatename_1D,nbinsy,kDYarray);
				h2DVBF[tr] = new TH2F(templatename_2D,templatename_2D,nbinsx,kDXarray,nbinsy,kDYarray);
			};

			int nVBFEntries = tree_VBF[tr]->GetEntries();
			cout << "Here: "  << tr << " N: " << nVBFEntries << endl;

			for(int ev=0;ev<nVBFEntries;ev++){
				tree_VBF[tr]->GetEntry(ev);
				if(isSelected!=1) continue;
				if (Z1Mass >= 120 || Z2Mass >= 120) continue;

				double weight = MC_weight;
				if(ev==0) cout << "Weight: " << weight << endl;
				if (EnergyIndex == 1){
					if (folder == 2) weight /= 2.0;
					else weight *= (tg_interf->Eval(ZZMass)) / 4.0;
				};

				//Still needed?
				if( (ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) ) continue;

				if(tFitD==0) h1DVBF[tr]->Fill(ZZMass,weight);
				else{
					h1DVBF[tr]->Fill(fitYval,weight);
					h2DVBF[tr]->Fill(ZZMass,fitYval,weight);
				};
			};
			//if(tr==3) cout << "Phantom Expected: " << nVBF_Sig_Simulated << endl;

			ZZMass = 0;
			Z1Mass = 0;
			Z2Mass = 0;
			MC_weight = 1;
			fitYval=0;

			double presmoothInt = h1DVBF[tr]->Integral();
			if(isSmooth) h1DVBF[tr]->Smooth(1,"k3a");
			double postsmoothInt = h1DVBF[tr]->Integral();
			h1DVBF[tr]->Scale(presmoothInt/postsmoothInt);

			foutput->WriteTObject(h1DVBF[tr]);
			if(tFitD!=0){
				presmoothInt = h2DVBF[tr]->Integral();
				if(isSmooth) h2DVBF[tr]->Smooth(1,"k3a");
				postsmoothInt = h2DVBF[tr]->Integral();
				h2DVBF[tr]->Scale(presmoothInt/postsmoothInt);
				foutput->WriteTObject(h2DVBF[tr]);
			};
		};

		//Make VBF Bkg/Sig/Int from linear combinations of above templates
		for(int binx=1;binx<=nbinsx;binx++){
			double bkg = h1DVBF[0]->GetBinContent(binx);
			double bsi10 = h1DVBF[1]->GetBinContent(binx);
			double bsi25 = h1DVBF[2]->GetBinContent(binx);
			double bsi25p = bsi25 - bkg;
			double bsi10p = bsi10 - bkg;
			double signal = ( sqrt(10.0)*bsi25p - 5.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
			double interf = ( -10.0*bsi25p + 25.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
			if(signal<0){
				signal=0;
				interf=0;
			};
			if(bkg<0) bkg=0;

			h1DVBF[0]->SetBinContent(binx,signal);
			h1DVBF[1]->SetBinContent(binx,bkg);
			h1DVBF[2]->SetBinContent(binx,interf);

			if(tFitD>0){
				for(int biny=1;biny<=nbinsy;biny++){
					bkg = h2DVBF[0]->GetBinContent(binx,biny);
					bsi10 = h2DVBF[1]->GetBinContent(binx,biny);
					bsi25 = h2DVBF[2]->GetBinContent(binx,biny);
					bsi25p = bsi25 - bkg;
					bsi10p = bsi10 - bkg;
					signal = ( sqrt(10.0)*bsi25p - 5.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
					interf = ( -10.0*bsi25p + 25.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
					if(signal<0){
						signal=0;
						interf=0;
					};
					if(bkg<0) bkg=0;

					h2DVBF[0]->SetBinContent(binx,biny,signal);
					h2DVBF[1]->SetBinContent(binx,biny,bkg);
					h2DVBF[2]->SetBinContent(binx,biny,interf);
				};
			};
		};

		nVBF_Sig_Simulated = h2DVBF[0]->Integral()*luminosity[EnergyIndex];
		cout << h2DVBF[0]->Integral()*luminosity[EnergyIndex] << endl;
		cout << h2DVBF[1]->Integral()*luminosity[EnergyIndex] << endl;
		cout << h2DVBF[2]->Integral()*luminosity[EnergyIndex] << endl;
		cout << h2DVBF[3]->Integral()*luminosity[EnergyIndex] << endl;

		//Template and tree filler
		for(int t=0;t<kNumTemplates;t++){
			TChain* tree = new TChain(TREE_NAME);
			TTree* templateTree;

			//Get file names, templates, and appropriate branches
			if(t==3){
				for(int b=0;b<kNumBkg;b++){
					TString cinput = cinput_common_qqZZ;
					cinput += erg_dir + "/" + user_folder[folder] + "/" + hzz4lprefix + sample_BackgroundFile[b] + "_Reprocessed.root";
					tree->Add(cinput);
					cout << cinput << endl;
				};
			}
			else if(t==4){
				TString cinput = cinput_common_ZX;
				cinput += erg_dir + "/" + user_folder[folder] + "/HZZ4lTree_DoubleOr_CRZLLTree_Reprocessed.root";
				tree->Add(cinput);
				cout << cinput << endl;
			}
			else if(t<3){
				int tp = t;
				if (EnergyIndex == 0){ // USE BSI25, NO CONTAMINATION SIGNIFICANT TO CREATE STORAGE TREES
					if (folder == 2 && tp == 2) tp = 3;
					if (folder == 1 && tp == 1) tp = 3;
				};
				TString cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				cout<<cinput_2e2mu<<endl;
				tree->Add(cinput_4e);
				cout<<cinput_4e<<endl;
				tree->Add(cinput_4mu);
				cout<<cinput_4mu<<endl;				
			};
			char templatename_1D[100];
			char templatename_2D[100];
			if(t<2){
				sprintf(templatename_1D,"T_1D_%i",t+1);
				sprintf(templatename_2D,"T_2D_%i",t+1);
			}
			else if(t==2){
				sprintf(templatename_1D,"T_1D_%i",4);
				sprintf(templatename_2D,"T_2D_%i",4);
			}
			else if(t==3){
				sprintf(templatename_1D,"T_1D_%s","qqZZ");
				sprintf(templatename_2D,"T_2D_%s","qqZZ");
			}
			else if(t==4){
				sprintf(templatename_1D,"T_1D_%s","ZX");
				sprintf(templatename_2D,"T_2D_%s","ZX");
			}
			else if(t==5 || t==6){
				sprintf(templatename_1D,"T_1D_VBF_%i",t-4);
				sprintf(templatename_2D,"T_2D_VBF_%i",t-4);
			}
			else if(t==7){
				sprintf(templatename_1D,"T_1D_VBF_%i",4);
				sprintf(templatename_2D,"T_2D_VBF_%i",4);
			};
			if(tFitD==0){
				D_temp_1D[t] = new TH1F(templatename_1D,templatename_1D,nbinsx,kDXarray);
				D_temp_1D[t]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);

				TString strTreeName = templatename_1D;
				strTreeName += "_Tree";
				templateTree = new TTree(strTreeName,strTreeName);
			}
			else{
				D_temp_1D[t] = new TH1F(templatename_1D,templatename_1D,nbinsy,kDYarray);
				D_temp_1D[t]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);

				D_temp_2D[t] = new TH2F(templatename_2D,templatename_2D,nbinsx,kDXarray,nbinsy,kDYarray);
				D_temp_2D[t]->GetXaxis()->SetTitle(strFitDim_label[0]);
				D_temp_2D[t]->GetYaxis()->SetTitle(strFitDim_label[tFitD]);

				TString strTreeName = templatename_2D;
				strTreeName += "_Tree";
				templateTree = new TTree(strTreeName,strTreeName);

				if(t==3){
					TString storeName = templatename_2D;
					storeName = storeName + "_UnConditional";
					hStore_qqZZ_Unconditional = new TH2F(storeName,storeName,nbinsx,kDXarray,nbinsy,kDYarray);
					hStore_qqZZ_Unconditional->GetXaxis()->SetTitle(strFitDim_label[0]);
					hStore_qqZZ_Unconditional->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
				};
				if(t==4){
					TString storeName = templatename_2D;
					storeName = storeName + "_UnConditional";
					hStore_ZX_Unconditional = new TH2F(storeName,storeName,nbinsx,kDXarray,nbinsy,kDYarray);
					hStore_ZX_Unconditional->GetXaxis()->SetTitle(strFitDim_label[0]);
					hStore_ZX_Unconditional->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
				};
			};
			templateTree->Branch("templateWeight",&templateWeight);
			templateTree->Branch("ZZMass",&ZZMass);
			if(tFitD!=0) templateTree->Branch(strFitDim[tFitD],&fitYval);

			//Making templates using appropriate weights
			double nTotal=0;
			if(t<5){
				tree->SetBranchAddress("GenHMass",&GenHMass);
				tree->SetBranchAddress("ZZMass",&ZZMass);
				if(t!=4) tree->SetBranchAddress("MC_weight",&MC_weight);
				else tree->SetBranchAddress("ZXfake_weightProper",&MC_weight);
				if(t<3){
					tree->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
					tree->SetBranchAddress("MC_weight_Kfactor_Norm_down", &MC_weight_Kfactor_Norm_down);
					tree->SetBranchAddress("MC_weight_Kfactor_Norm_up", &MC_weight_Kfactor_Norm_up);
					tree->SetBranchAddress("MC_weight_PDF_Norm_down", &MC_weight_PDF_Norm_down);
					tree->SetBranchAddress("MC_weight_PDF_Norm_up", &MC_weight_PDF_Norm_up);
					tree->SetBranchAddress("MC_weight_down", &MC_weight_down);
					tree->SetBranchAddress("MC_weight_up",&MC_weight_up);
				}
				else{
					MC_weight_ggZZLepInt=1;
					MC_weight_Kfactor=1;
					MC_weight_down=1;
					MC_weight_up=1;
					MC_weight_Kfactor_Norm_down = 1;
					MC_weight_Kfactor_Norm_up = 1;
					MC_weight_PDF_Norm_down = 1;
					MC_weight_PDF_Norm_up = 1;
				};
				if(t<3 || t>4) tree->SetBranchAddress("MC_weight_ggZZLepInt",&MC_weight_ggZZLepInt);
				else MC_weight_ggZZLepInt=1;

				if(tree->GetBranchStatus(strFitDim[tFitD])) tree->SetBranchAddress(strFitDim[tFitD],&fitYval);
				else if(!tree->GetBranchStatus(strFitDim[tFitD])){
					cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
					fitYval=0;
				};

				nEntries = tree->GetEntries();
				for(int ev=0;ev<nEntries;ev++){
					tree->GetEntry(ev);
					double weight = MC_weight;
					if (t<3) weight *= MC_weight_ggZZLepInt;
					if (abs(Systematics) != 1 && t<3) weight *= MC_weight_Kfactor;
					if (Systematics == -1 && t<3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_down*tgkf->Eval(125.6);
					if (Systematics == 1 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_up*tgkf->Eval(125.6);
					if (Systematics == -2 && t<3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_down;
					if (Systematics == 2 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_up;
					if (t>4) weight *= tsp_VBF_Sig->Eval(GenHMass)*MC_weight_ggZZLepInt; // REMEMBER, THIS IS FROM GGH TEMPLATES,SHOULD NOT BE THE CASE IN 7 TEV VBF PHANTOM SAMPLES

					if(t==0 && ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]) nSig_Simulated += weight;
					if(t<3 && ZZMass>=ZZMass_cut[0]) nMZZ220[t] += weight;

					if( (ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) ) continue;

					nTotal += weight;

					if(tFitD==0) D_temp_1D[t]->Fill(ZZMass,weight);
					if(fitYval!=fitYval) cout<<"Broken event! "<<ev<<"/"<<nEntries<<endl;
					else D_temp_1D[t]->Fill(fitYval,weight);
					if(tFitD>0) D_temp_2D[t]->Fill(ZZMass,fitYval,weight);
				};
			};
			cout << t << " NTotal: " << nTotal << endl;

			//Reweighting normalization
			if(t<3){
				cout << nMZZ220[t]*nSM_ScaledPeak[EnergyIndex][folder]/nSig_Simulated*luminosity[EnergyIndex] << endl;

				double myscale = nSM_ScaledPeak[EnergyIndex][folder] / nSig_Simulated;
				if (Systematics == -1 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == 1 && t < 3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == -2 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				if (Systematics == 2 && t<3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				overall_scale[t] = myscale;
				D_temp_1D[t]->Scale(myscale);
				if (tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling t: " << t << " by " << myscale << endl;
			}
			else if(t>=5){
				double myscale = (nVBF_Sig_Reweighted/nVBF_Sig_Simulated);

				cout << "Expected nVBF_peak: " << nVBF_Sig_Reweighted << endl;
				cout << "Observed nVBF_peak: " << nVBF_Sig_Simulated << endl;
				cout << "Scaling t: " << t << " by " << myscale << endl;
				overall_scale[t] = myscale;

				double integral1D=0;
				double integral2D=0;
				for(int binx=1;binx<=nbinsx;binx++){
					double bincontent = h1DVBF[t-5]->GetBinContent(binx);
					bincontent *= myscale;

					if(Systematics!=0){
						double scalingRatio = 1;
						if(t==5){
							double othercontent = hVBFRef_1D->Integral(binx,binx,"width");
							if(Systematics>0 && bincontent!=0) scalingRatio = othercontent/bincontent;
							else if(Systematics<0 && othercontent!=0) scalingRatio = bincontent/othercontent;
							else scalingRatio = 1;
							if(scalingRatio<0) scalingRatio=0;
							if(scalingRatio>4) scalingRatio=4;
							h1DVBFSigRatio->SetBinContent(binx,scalingRatio);
						}
						else scalingRatio = h1DVBFSigRatio->GetBinContent(binx);
						bincontent *= scalingRatio;
						integral1D += bincontent;
					};

					D_temp_1D[t]->SetBinContent(binx,bincontent);

					if(tFitD>0){
						for(int biny=1;biny<=nbinsy;biny++){
							bincontent = h2DVBF[t-5]->GetBinContent(binx,biny);
							bincontent *= myscale;

							if(Systematics!=0){
								double scalingRatio = 1;
								if(t==5){
									double othercontent = hVBFRef->Integral(binx,binx,biny,biny,"width");
									if(Systematics>0 && bincontent!=0) scalingRatio = othercontent/bincontent;
									else if(Systematics<0 && othercontent!=0) scalingRatio = bincontent/othercontent;
									else scalingRatio = 1;
									if(scalingRatio<0) scalingRatio=0;
									if(scalingRatio>4) scalingRatio=4;
									h2DVBFSigRatio->SetBinContent(binx,biny,scalingRatio);
								}
								else scalingRatio = h2DVBFSigRatio->GetBinContent(binx,biny);
								bincontent *= scalingRatio;
								integral2D += bincontent;
							};

							D_temp_2D[t]->SetBinContent(binx,biny,bincontent);
						};
					};
				};
				if(Systematics!=0 && t==5){
					myscale = (nVBF_Sig_Reweighted/(integral1D*luminosity[EnergyIndex]));
					cout << "VBF SYSTEMATIC 1D SCALE IS " << myscale << endl;
					D_temp_1D[t]->Scale(myscale);
					h1DVBFSigRatio->Scale(myscale);
					if(tFitD>0){
						myscale = (nVBF_Sig_Reweighted/(integral2D*luminosity[EnergyIndex]));
						cout << "VBF SYSTEMATIC 2D SCALE IS " << myscale << endl;
						D_temp_2D[t]->Scale(myscale);
						h2DVBFSigRatio->Scale(myscale);
					};
				};
			}
			else overall_scale[t]=1.0;

			//Smooths templates if desired
			if(t<5){
				double presmoothInt = D_temp_1D[t]->Integral("width");
				if(isSmooth) D_temp_1D[t]->Smooth(1,"k3a");
				double postsmoothInt = D_temp_1D[t]->Integral("width");
				D_temp_1D[t]->Scale(presmoothInt/postsmoothInt);
				if(tFitD>0){
					presmoothInt = D_temp_2D[t]->Integral("width");
					if(isSmooth) D_temp_2D[t]->Smooth(1,"k3a");
					postsmoothInt = D_temp_2D[t]->Integral("width");
					D_temp_2D[t]->Scale(presmoothInt/postsmoothInt);
					if(isSmooth) cout << "SMOOTHING SCALE: " << presmoothInt/postsmoothInt << endl;
				};
			};

			//Makes uncondtionally normalized PDF for backgrounds
			if(t>=3 && tFitD>0 && t<5){
				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);
						if(t==3) hStore_qqZZ_Unconditional->SetBinContent(binx,biny,bincontent);
						if(t==4) hStore_ZX_Unconditional->SetBinContent(binx,biny,bincontent);
						if(intBinX!=0) D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
					};
				};
			};

			//Stores total weights in tree for JB's smoother
			double nTotalRecorded=0;
			if(t<5){
				nEntries = tree->GetEntries();
				for(int ev=0;ev<nEntries;ev++){
					tree->GetEntry(ev);

					if( (ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) ) continue;

					double weight = MC_weight;
					if (t<3) weight *= MC_weight_ggZZLepInt;
					if (abs(Systematics) != 1 && t<3) weight *= MC_weight_Kfactor;
					if (Systematics == -1 && t<3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_down*tgkf->Eval(125.6);
					if (Systematics == 1 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_up*tgkf->Eval(125.6);
					if (Systematics == -2 && t<3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_down;
					if (Systematics == 2 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_up;

					templateWeight = weight * overall_scale[t];

					templateTree->Fill();
					nTotalRecorded += templateWeight;
				};
			}
			else{
				int treeCode=0;
				if(t==5) treeCode=2;
				if(t==6) treeCode=0;
				if(t==7) treeCode=1;
				nEntries = tree_VBF[treeCode]->GetEntries();

				cout<<"TEST FOR VBF "<<endl;
				for(int ev=0;ev<nEntries;ev++){
					tree_VBF[treeCode]->GetEntry(ev);
					if (isSelected != 1) continue;
					if (Z1Mass >= 120 || Z2Mass >= 120) continue;
					if ((ZZMass >= ZZMass_cut[1] || ZZMass<ZZMass_cut[0])) continue;

					double weight = MC_weight;
					if (EnergyIndex == 1){
						if (folder == 2) weight /= 2.0;
						else weight *= (tg_interf->Eval(ZZMass)) / 4.0;
					};
					weight *= overall_scale[t];

					if(Systematics!=0){
						double sysVBFScale=1;
						if(tFitD==0) sysVBFScale=h1DVBFSigRatio->GetBinContent(h1DVBFSigRatio->FindBin(ZZMass));
						else sysVBFScale=h2DVBFSigRatio->GetBinContent(h2DVBFSigRatio->FindBin(ZZMass,fitYval));
						weight *= sysVBFScale;
					};

					templateWeight = weight;

					templateTree->Fill();
					nTotalRecorded += templateWeight;
				};
			};

			if(t!=3) cout << "RECORDED YIELD IN " << templateTree->GetName() << ": " << nTotalRecorded*luminosity[EnergyIndex] << endl;
			else cout << "RECORDED YIELD IN " << templateTree->GetName() << ": " << nTotalRecorded << endl;
			foutput->WriteTObject(templateTree);
			delete templateTree;
			delete tree;
		};

		if(Systematics!=0){
			foutput->WriteTObject(h1DVBFSigRatio);
			if(tFitD!=0) foutput->WriteTObject(h2DVBFSigRatio);
		};
		delete h2DVBFSigRatio;
		delete h1DVBFSigRatio;
		for(int tr=0;tr<4;tr++){
			delete tree_VBF[tr];
			delete h1DVBF[tr];
			if(tFitD!=0) delete h2DVBF[tr];
		};

		//Makes Bkg, Sig, Int from Linear Combinations
		if (EnergyIndex == 0){ // USE BSI25, NO CONTAMINATION SIGNIFICANT TO CREATE STORAGE TREES
			if (folder == 2){
				D_temp_1D[2]->Add(D_temp_1D[1], -1.0);
				D_temp_1D[2]->Add(D_temp_1D[0], -25.0);
				D_temp_1D[2]->Scale(0.2);
				D_temp_1D[2]->Add(D_temp_1D[1], 1.0);
				D_temp_1D[2]->Add(D_temp_1D[0], 1.0);
				if (tFitD != 0){
					D_temp_2D[2]->Add(D_temp_2D[1], -1.0);
					D_temp_2D[2]->Add(D_temp_2D[0], -25.0);
					D_temp_2D[2]->Scale(0.2);
					D_temp_2D[2]->Add(D_temp_2D[1], 1.0);
					D_temp_2D[2]->Add(D_temp_2D[0], 1.0);
				};
			};
			if (folder == 1){
				for (int binx = 1; binx <= D_temp_1D[1]->GetNbinsX(); binx++){
					double sig = D_temp_1D[0]->GetBinContent(binx);
					double bsi25 = D_temp_1D[1]->GetBinContent(binx);
					double bsi = D_temp_1D[2]->GetBinContent(binx);

					bsi = 5.0*(bsi - sig);
					bsi25 = (bsi25 - 25.0*sig - bsi)*(-0.25);
					if (bsi25 < 0) bsi25 = 0;
					D_temp_1D[1]->SetBinContent(binx, bsi25);
				};

				if (tFitD != 0){
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
		};

		D_temp_1D[2]->Add(D_temp_1D[1],-1.0);
		D_temp_1D[2]->Add(D_temp_1D[0],-1.0);
		if(tFitD!=0){
			D_temp_2D[2]->Add(D_temp_2D[1],-1.0);
			D_temp_2D[2]->Add(D_temp_2D[0],-1.0);
		};

		//Divides bins by Bin Width
		cout << "Integrals after everything:\n1D\t2D" << endl;
		for(int t=0;t<kNumTemplates;t++){
			for(int binx=0;binx<nbinsx;binx++){
				double binwidthx;
				if(tFitD==0) binwidthx = kDXarray[binx+1] - kDXarray[binx];
				else binwidthx = kDYarray[binx+1] - kDYarray[binx];
				double bincontent = D_temp_1D[t]->GetBinContent(binx+1);
				if(t!=3 && t!=4) bincontent /= binwidthx;
				D_temp_1D[t]->SetBinContent(binx+1,bincontent);

				binwidthx = kDXarray[binx+1] - kDXarray[binx];

				if(tFitD!=0){
					for(int biny=0;biny<nbinsy;biny++){
						double binwidthy = kDYarray[biny+1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						bincontent = D_temp_2D[t]->GetBinContent(binx+1,biny+1);
						if(t!=3 && t!=4) bincontent /= binwidth;
						D_temp_2D[t]->SetBinContent(binx+1,biny+1,bincontent);
					};
				};
			};
			if(tFitD!=0 && t==3){
				for(int binx=0;binx<nbinsx;binx++){
					double binwidthx = kDXarray[binx+1] - kDXarray[binx];
					for(int biny=0;biny<nbinsy;biny++){
						double binwidthy = kDYarray[biny+1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						double bincontent = hStore_qqZZ_Unconditional->GetBinContent(binx+1,biny+1);
						bincontent /= binwidth;
						hStore_qqZZ_Unconditional->SetBinContent(binx+1,biny+1,bincontent);
					};
				};
			};
			if(tFitD!=0 && t==4){
				for(int binx=0;binx<nbinsx;binx++){
					double binwidthx = kDXarray[binx+1] - kDXarray[binx];
					for(int biny=0;biny<nbinsy;biny++){
						double binwidthy = kDYarray[biny+1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						double bincontent = hStore_ZX_Unconditional->GetBinContent(binx+1,biny+1);
						bincontent /= binwidth;
						hStore_ZX_Unconditional->SetBinContent(binx+1,biny+1,bincontent);
					};
				};
			};
			if(t==4 && Systematics!=0){
				double intTZX_1D = D_temp_1D[t]->Integral();
				double intTZZQQB_1D = D_temp_1D[3]->Integral();
				for(int binx=1;binx<=D_temp_1D[t]->GetNbinsX();binx++){
					double bincontent = D_temp_1D[t]->GetBinContent(binx);
					double bincontent_alt = D_temp_1D[3]->GetBinContent(binx);
					bincontent_alt *= intTZX_1D/intTZZQQB_1D;
					double difference = bincontent_alt - bincontent;

					if(Systematics>0) bincontent += difference;
					else bincontent -= difference;
					if(bincontent < 0) bincontent = 0;
					D_temp_1D[t]->SetBinContent(binx,bincontent);
				};
				D_temp_1D[t]->Scale(intTZX_1D/D_temp_1D[t]->Integral());

				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double* storeOriginal = new double[D_temp_2D[t]->GetNbinsY()];
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);
						double bincontent_alt = D_temp_2D[3]->GetBinContent(binx,biny);
						double difference = bincontent_alt - bincontent;
						storeOriginal[biny-1] = bincontent;

						if(Systematics>0) bincontent += difference;
						else bincontent -= difference;
						if(bincontent < 0) bincontent = 0;
						D_temp_2D[t]->SetBinContent(binx,biny,bincontent);
					};
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);

						if(intBinX!=0){
							D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
							double sysRatio = 0;
							if(storeOriginal[biny-1]!=0) sysRatio = (bincontent/intBinX)/(storeOriginal[biny-1]);
							double unconditionalbincontent = hStore_ZX_Unconditional->GetBinContent(binx,biny);
							hStore_ZX_Unconditional->SetBinContent(binx,biny,unconditionalbincontent*sysRatio);
						}
						else{
							hStore_ZX_Unconditional->SetBinContent(binx,biny,0);
						};
					};
					delete[] storeOriginal;
				};
				if(Systematics<0){
					double presmoothInt = D_temp_1D[t]->Integral();
					if(isSmooth) D_temp_1D[t]->Smooth(1,"k3a");
					double postsmoothInt = D_temp_1D[t]->Integral();
					D_temp_1D[t]->Scale(presmoothInt/postsmoothInt);
				};
			};
			foutput->WriteTObject(D_temp_1D[t]);
			if(tFitD!=0) foutput->WriteTObject(D_temp_2D[t]);
			if(tFitD!=0 && t==3) foutput->WriteTObject(hStore_qqZZ_Unconditional);
			if(tFitD!=0 && t==4) foutput->WriteTObject(hStore_ZX_Unconditional);

			if(t!=3 || t!=4){
				cout << D_temp_1D[t]->Integral("width")*luminosity[EnergyIndex] << '\t';
				if(tFitD!=0) cout << D_temp_2D[t]->Integral("width")*luminosity[EnergyIndex] << endl;
				else cout << endl;
			}
			else{
				cout << D_temp_1D[t]->Integral(1,nbinsx) << '\t';
				if(tFitD!=0) cout << D_temp_2D[t]->Integral(1,nbinsx,1,nbinsy) << endl;
				else cout << endl;
			};
		};
		if(tFitD!=0){
			cout << "Unconditional Integrals are: " << endl;
			cout << hStore_qqZZ_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
			cout << hStore_ZX_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
		};
		delete hStore_qqZZ_Unconditional;
		delete hStore_ZX_Unconditional;

		finput_VBFRef->Close();
		foutput->Close();
	};
	delete tgkf;
	finput_KDFactor->Close();
	finput_VBF->Close();
	delete tg_interf;
};
