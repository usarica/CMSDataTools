#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include "TChain.h"
#include "TString.h"
#include "TSpline.h"
#include "TSystem.h"
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
enum histtypes{kSigHist,kBkgHist,kIntHist,kBSIHist,kBSI10Hist,kBSI25Hist};
void makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth);
double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
double doPeakCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype);
void oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, TH1F* finaloutput, int outputtype, TH1F* finaloutput2=0, int output2type=-99);
void twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, TH2F* finaloutput, int outputtype, TH2F* finaloutput2=0, int output2type=-99);
void threeDlinearcombination(TH3F* first, int firsttype, TH3F* second, int secondtype, TH3F* input, int inputtype, TH3F* finaloutput, int outputtype, TH3F* finaloutput2=0, int output2type=-99);
void progressbar(int val, int tot);

//Make Lepton Interference Graph to be used later
TGraph* make_HZZ_LeptonInterferenceGraph(){
	float x[leptonInterf_YR3_Size];
	float y[leptonInterf_YR3_Size];
	for(int a=0;a<leptonInterf_YR3_Size;a++){
		x[a] = leptonInterf_YR3[a][0];
		y[a] = leptonInterf_YR3[a][1];
	}
	TGraph* tg = new TGraph(leptonInterf_YR3_Size,x,y);
	tg->SetName("tgHZZ_LeptonInterference");
	tg->SetTitle("H#rightarrowZZ#rightarrow4l Lepton Interference Weight on 4e, 4#mu wrt. 2e2#mu");

	return tg;
};

//Main Function, runs over all desired iterations
void makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF(){
	bool isSmooth=false;
	const int kNumSyst=5;
	int systematics[kNumSyst]={0,1,-1,2,-2};
	for(int i=0;i<kNumSyst;++i){
		for(int usesmooth=0;usesmooth<2;++usesmooth){
			if(usesmooth==0) isSmooth=false;
			if(usesmooth==1) isSmooth=true;
			for(int CoM=7;CoM<9;++CoM){
				for(int channel=0;channel<3;++channel){
					makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(channel,CoM,6,systematics[i],isSmooth);	
				}	
			}
		}
	}
}

//Function to build one template
// folder = 0,1,2 (final state corresponds to 4mu,4e,2mu2e respectively)
// erg_tev = 7,8 (CoM energy)
// tFitD = [0,16] (choice of Discriminant, see FitDimensionsList.h for list; only tFitd works right now)
// Systematics = [-2,2] (Flag for systematics. 0=Nominal, +/-1=QCD, +/-2=PDF)
// isSmooth = true/false (flag to apply smoothing at this stage, both needed for later stage)
void makeCombineTemplateswithTrees_Modified_MCFM_GenLevelVBF_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth){
	char TREE_NAME[]="SelectedTree";
	TString INPUT_NAME = "HZZ4lTree_ggTo";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine_";
	if(!isSmooth) OUTPUT_NAME += "Raw_";
	OUTPUT_NAME += "_GenLevelVBF_" + TString(strFitDim[tFitD]) + "_";
	TString comstring;
	comstring.Form("%i",erg_tev);
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/",erg_tev);
	if (Systematics == 0) OUTPUT_NAME += "Nominal.root";
	if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF.root";

	TString INPUT_VBFREF_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	if(!isSmooth) INPUT_VBFREF_NAME += "Raw_";
	INPUT_VBFREF_NAME += TString(strFitDim[tFitD]) + "_";
	INPUT_VBFREF_NAME += "Nominal.root";

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
	TString cinput_common_VBF = user_gg2VV_location; //CHANGE - made to point to correct directory

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float lowside[3]={220,230,240};
	float ZZMass_PeakCut[2]={105.6,140.6}; // Spin 0 analysis
	//use ggH yields directly
 	double nSM_ScaledPeak[2][3]={
 		{1.0902689,0.6087736,1.4452079},
 		{5.1963998,2.6944639,6.6562629}
	};
	//use VBF yields directly
	double VBF_Sig_Datacard[2][3]={
		{0.092458836,0.051755897,0.12861921},
		{0.46798807,0.24788553,0.61781689}
	};
	// disregard VH, ttH, bbH etc. at peak for the sake of offshell
	// Systematics
	double ggZZ_Syst_AbsNormSyst[2][2] = { // EnergyIndex { QCD, PDF }
		{ 0.0745, 0.0735 },
		{ 0.075, 0.072 }
	};

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

	TString cinput_KDFactor = "./data/HZZ4l-KDFactorGraph";
	if (EnergyIndex == 0) cinput_KDFactor = cinput_KDFactor + "_7TeV";
	cinput_KDFactor = cinput_KDFactor + ".root";
	TFile* finput_KDFactor = new TFile(cinput_KDFactor, "read");
	TString tgkfname = "KDScale_";
	tgkfname = tgkfname + "AllFlavors_UnNormalized";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*)finput_KDFactor->Get(tgkfname);

	double overall_VBF_scale=1;
	double nVBFPeak[4] = { 0 };
	for (int e = 0; e < 2; e++){for (int ss = 0; ss < 3; ss++){ nSM_ScaledPeak[e][ss] /= luminosity[e]; VBF_Sig_Datacard[e][ss] /= luminosity[e]; }}

	TGraph* tg_interf = make_HZZ_LeptonInterferenceGraph();

	TString cinput_VBF_Sig = "./data/HZZ4l-125_6-" + comstring + "TeV-Sig_MCFM_PhantomVBF_Comparison.root";
	TFile* finput_VBF = new TFile(cinput_VBF_Sig,"read");
	TSpline3* tsp_VBF_Sig = (TSpline3*) finput_VBF->Get("Spline3");

	TH2F* hStore_ZX_Unconditional;
	TH2F* hStore_qqZZ_Unconditional;

	for(int lo=0;lo<1;lo++){
		TString coutput_common = user_TemplateswithTrees_dir + erg_dir;
		gSystem->Exec("mkdir -p " + coutput_common);
		coutput_common += user_folder[folder] + "/";
		gSystem->Exec("mkdir -p " + coutput_common);

		cout<<"==============================="<<endl;
		cout<<"CoM Energy: "<<erg_tev<<" TeV"<<endl;
		cout<<"Decay Channel: "<<user_folder[folder]<<endl;
		cout<<endl;

		TString coutput = coutput_common + OUTPUT_NAME;
		TFile* foutput = new TFile(coutput,"recreate");

		//Grab VBF Sig templates made in makeCombineTemplates_Modified_MCFM.c
		TString cinput_VBFRef = user_dir + erg_dir + user_folder[folder] + "/" + INPUT_VBFREF_NAME;
		TFile* finput_VBFRef = new TFile(cinput_VBFRef,"read");
		cout<<cinput_VBFRef<<endl;
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
		}
	
		int nbinsy=30;
		float kDYarray[nbinsy+1];
		float kDY_bounds[2]={0,1};
		if(tFitD==3){kDY_bounds[0]=-7.0;kDY_bounds[1]=3.0;}
		if(tFitD==4 || tFitD==5 || tFitD==9 || tFitD==10 || tFitD==11 || tFitD==15){kDY_bounds[0]=-1.0;kDY_bounds[1]=1.0;}
		for(int bin=0;bin<nbinsy+1;bin++){
			double binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		}

		const int kNumTemplates=8;
		TH1F* D_temp_1D[kNumTemplates];
		TH2F* D_temp_2D[kNumTemplates];
		double overall_scale[kNumTemplates]={1};

		int nEntries;
		double nSig_Simulated=0;
		double nVBF_Sig_Simulated=0;
		float fitYval;
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
		}
		for(int tr=0;tr<4;tr++){
			tree_VBF[tr] = new TChain(TREE_NAME);
		}

		for (int decay = 0; decay < 3; decay++){
			TString channelname = user_folder[decay];
			if (decay == 2) channelname = "2e2mu";

			/*
			0: BSI
			1: BKG
			2: BSI10
			3: BSI25
			*/
			for (int type_vbf = 0; type_vbf < 4; type_vbf++){
				TString cinput_VBF_Bkg = cinput_common + "HZZ4lTree_ZZTo" + channelname + "JJ_" + sample_suffix_Phantom[type_vbf] + "_Reprocessed.root";
				tree_VBF[type_vbf]->Add(cinput_VBF_Bkg);
			}
		}

		//Fill Template(s) from each Phantom sample with:
		// gen mZ1/mZ2>120
		// Reweighting for lepton interference (necessary)
		// Possible mZZ cuts?
		//Then write smoothed templates to file 
		for(int tr=0;tr<4;tr++){
			tree_VBF[tr]->SetBranchAddress("ZZMass", &ZZMass);
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
			}

			int nVBFEntries = tree_VBF[tr]->GetEntries();
			for(int ev=0;ev<nVBFEntries;ev++){
			    progressbar(ev,tree_VBF[tr]->GetEntries());
				tree_VBF[tr]->GetEntry(ev);

				double weight = MC_weight;
				if(ev==0) cout << "Weight: " << weight << endl;
/*				if (EnergyIndex == 1){
					if (folder == 2) weight /= 2.0;
					else weight *= (tg_interf->Eval(ZZMass)) / 4.0;
				}
*/
				if (ZZMass < ZZMass_PeakCut[1] && ZZMass >= ZZMass_PeakCut[0]){
					nVBFPeak[tr] += weight;
				}
				if( ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0] ) continue;

				if(tFitD==0) h1DVBF[tr]->Fill(ZZMass,weight);
				else{
					h1DVBF[tr]->Fill(fitYval,weight);
					h2DVBF[tr]->Fill(ZZMass,fitYval,weight);
				}
			}
			cout<<endl;

			ZZMass = 0;
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
			}
		}

		//Make VBF Sig/Int from linear combinations of above templates
		//0: VBF Sig
		//2: VBF Int
		//	 For 7 TeV samples, BSI25, Bkg, and BSI10 are used
		//	 For 8 TeV samples, BSI, Bkg, and BSI10 are used
		if (EnergyIndex == 0){
			oneDlinearcombination(h1DVBF[2],kBSI10Hist,h1DVBF[1],kBkgHist,h1DVBF[3],kBSI25Hist,h1DVBF[0],kSigHist,h1DVBF[2],kIntHist);
			if(tFitD!=0){
				twoDlinearcombination(h2DVBF[2],kBSI10Hist,h2DVBF[1],kBkgHist,h2DVBF[3],kBSI25Hist,h2DVBF[0],kSigHist,h2DVBF[2],kIntHist);
			}
			double nVBFsigtemp = doPeakCombination(nVBFPeak[2],10,nVBFPeak[3],25,nVBFPeak[1],kSigHist);
			double nVBFinterftemp = doPeakCombination(nVBFPeak[2],10,nVBFPeak[3],25,nVBFPeak[1],kIntHist);
			for(int tr=0;tr<4;tr++) cout << "nVBFPeak[" << tr << "]: " << nVBFPeak[tr] << endl;
			nVBFPeak[0]=nVBFsigtemp;
			nVBFPeak[2]=nVBFinterftemp;
		}
/*		else if (folder == 2 && EnergyIndex == 1){
			oneDlinearcombination(h1DVBF[0],kBSIHist,h1DVBF[1],kBkgHist,h1DVBF[2],kBSI10Hist,h1DVBF[0],kSigHist,h1DVBF[2],kIntHist);
			if(tFitD!=0){
				twoDlinearcombination(h2DVBF[0],kBSIHist,h2DVBF[1],kBkgHist,h2DVBF[2],kBSI10Hist,h2DVBF[0],kSigHist,h2DVBF[2],kIntHist);
			}
			double nVBFsigtemp = doPeakCombination(nVBFPeak[0],1,nVBFPeak[2],10,nVBFPeak[1],kSigHist);
			double nVBFinterftemp = doPeakCombination(nVBFPeak[0],1,nVBFPeak[2],10,nVBFPeak[1],kIntHist);
			for(int tr=0;tr<4;tr++) cout << "nVBFPeak[" << tr << "]: " << nVBFPeak[tr] << endl;
			nVBFPeak[0]=nVBFsigtemp;
			nVBFPeak[2]=nVBFinterftemp;
		}
*/		else if (EnergyIndex == 1){
			oneDlinearcombination(h1DVBF[0],kBSIHist,h1DVBF[1],kBkgHist,h1DVBF[2],kBSI10Hist,h1DVBF[0],kSigHist,h1DVBF[2],kIntHist);
			if(tFitD!=0){
				twoDlinearcombination(h2DVBF[0],kBSIHist,h2DVBF[1],kBkgHist,h2DVBF[2],kBSI10Hist,h2DVBF[0],kSigHist,h2DVBF[2],kIntHist);
			}
			double nVBFsigtemp = doPeakCombination(nVBFPeak[2],10,nVBFPeak[0],1,nVBFPeak[1],kSigHist);
			double nVBFinterftemp = doPeakCombination(nVBFPeak[2],10,nVBFPeak[0],1,nVBFPeak[1],kIntHist);
			for(int tr=0;tr<4;tr++) cout << "nVBFPeak[" << tr << "]: " << nVBFPeak[tr] << endl;
			nVBFPeak[0]=nVBFsigtemp;
			nVBFPeak[2]=nVBFinterftemp;
		}

//		double vbfscale = 1;
		double vbfscale = VBF_Sig_Datacard[EnergyIndex][folder] / nVBFPeak[0];
		for (int tr = 0; tr < 4; tr++){
			h1DVBF[tr]->Scale(vbfscale);
			if (tFitD != 0) h2DVBF[tr]->Scale(vbfscale);
			if(tr==0) cout << "VBF Signal scale for syst=0: " << vbfscale << endl;
		}

		if (tFitD != 0){
			nVBF_Sig_Simulated = h2DVBF[0]->Integral()*luminosity[EnergyIndex];
			cout << h2DVBF[0]->Integral()*luminosity[EnergyIndex] << endl;
			cout << h2DVBF[1]->Integral()*luminosity[EnergyIndex] << endl;
			cout << h2DVBF[2]->Integral()*luminosity[EnergyIndex] << endl;
			cout << h2DVBF[3]->Integral()*luminosity[EnergyIndex] << endl;
		}

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
				}
			}
			else if(t==4){
				TString cinput = cinput_common_ZX;
				cinput += erg_dir + "/" + user_folder[folder] + "/HZZ4lTree_DoubleOr_CRZLLTree_Reprocessed.root";
				tree->Add(cinput);
				cout << cinput << endl;
			}
			else if(t<3){
				int tp = t;
				TString cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				cout<<cinput_2e2mu<<endl;
				tree->Add(cinput_4e);
				cout<<cinput_4e<<endl;
				tree->Add(cinput_4mu);
				cout<<cinput_4mu<<endl;
			}

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
			}
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
				}
				if(t==4){
					TString storeName = templatename_2D;
					storeName = storeName + "_UnConditional";
					hStore_ZX_Unconditional = new TH2F(storeName,storeName,nbinsx,kDXarray,nbinsy,kDYarray);
					hStore_ZX_Unconditional->GetXaxis()->SetTitle(strFitDim_label[0]);
					hStore_ZX_Unconditional->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
				}
			}
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
				}
				if(t<3 || t>4) tree->SetBranchAddress("MC_weight_ggZZLepInt",&MC_weight_ggZZLepInt);
				else MC_weight_ggZZLepInt=1;

				//ONLY CHANGE!!!
				if(tree->GetBranchStatus(strFitDim[tFitD])) tree->SetBranchAddress(strFitDim[tFitD],&fitYval);
				else if(!tree->GetBranchStatus(strFitDim[tFitD])){
					cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
					fitYval=0;
				}

				nEntries = tree->GetEntries();
				for(int ev=0;ev<nEntries;ev++){
				    progressbar(ev,tree->GetEntries());
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
					else D_temp_1D[t]->Fill(fitYval,weight);
					if(tFitD>0) D_temp_2D[t]->Fill(ZZMass,fitYval,weight);
				}
				cout<<endl;
			}
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

				cout << "Expected nVBF offshell from reweighting: " << nVBF_Sig_Reweighted << endl;
				cout << "Observed nVBF offshell from peak normalization: " << nVBF_Sig_Simulated << endl;
				cout << "VBF re-wgted/full sim. at t: " << t << " by " << myscale << endl;
				overall_scale[t] = vbfscale;

				double integral1D=0;
				double integral2D=0;
				for(int binx=1;binx<=nbinsx;binx++){
					double bincontent = h1DVBF[t-5]->GetBinContent(binx);
//					bincontent *= myscale;

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
					}

					D_temp_1D[t]->SetBinContent(binx,bincontent);

					if(tFitD>0){
						for(int biny=1;biny<=nbinsy;biny++){
							bincontent = h2DVBF[t-5]->GetBinContent(binx,biny);
//							bincontent *= myscale;

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
							}

							D_temp_2D[t]->SetBinContent(binx,biny,bincontent);
						}
					}
				}
				if(Systematics!=0 && t==5){
					if(Systematics>0) myscale = (nVBF_Sig_Reweighted/(integral1D*luminosity[EnergyIndex]));
					else  myscale = ((nVBF_Sig_Simulated/myscale)/(integral1D*luminosity[EnergyIndex]));
					cout << "VBF SYSTEMATIC 1D SCALE IS " << myscale << endl;
					D_temp_1D[t]->Scale(myscale);
					h1DVBFSigRatio->Scale(myscale);
					if(tFitD>0){
						if(Systematics>0) myscale = (nVBF_Sig_Reweighted/(integral2D*luminosity[EnergyIndex]));
						else  myscale = ((nVBF_Sig_Simulated/myscale)/(integral2D*luminosity[EnergyIndex]));
						cout << "VBF SYSTEMATIC 2D SCALE IS " << myscale << endl;
						D_temp_2D[t]->Scale(myscale);
						h2DVBFSigRatio->Scale(myscale);
					}
				}
			}
			else overall_scale[t]=1.0;

			//Smooths templates if desired
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
			}

			//Makes uncondtionally normalized PDF for backgrounds
			if(t>=3 && tFitD>0 && t<5){
				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);
						if(t==3) hStore_qqZZ_Unconditional->SetBinContent(binx,biny,bincontent);
						if(t==4) hStore_ZX_Unconditional->SetBinContent(binx,biny,bincontent);
						if(intBinX!=0) D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
					}
				}
			}

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
				}
			}
			else{
				int treeCode=0;
				if(t==5){
					//This is to account for the lack of a BSI25 sample for 8T 2e2mu
					if(EnergyIndex==0) treeCode = 3;
					else if(EnergyIndex==1) treeCode = t-5;
				}
				else if(t<8) treeCode = t-5;
				nEntries = tree_VBF[treeCode]->GetEntries();
				for(int ev=0;ev<nEntries;ev++){
					tree_VBF[treeCode]->GetEntry(ev);
					if (ZZMass >= ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) continue;

					double weight = MC_weight;
/*					if (EnergyIndex == 1){
						if (folder == 2) weight /= 2.0;
						else weight *= (tg_interf->Eval(ZZMass)) / 4.0;
					}
*/
					weight *= overall_scale[t];

					if(Systematics!=0){
						double sysVBFScale=1;
						if(tFitD==0) sysVBFScale=h1DVBFSigRatio->GetBinContent(h1DVBFSigRatio->FindBin(ZZMass));
						else sysVBFScale=h2DVBFSigRatio->GetBinContent(h2DVBFSigRatio->FindBin(ZZMass,fitYval));
						weight *= sysVBFScale;
					}

					templateWeight = weight;

					templateTree->Fill();
					nTotalRecorded += templateWeight;
				}
			}

			if(t!=3) cout << "RECORDED YIELD IN " << templateTree->GetName() << ": " << nTotalRecorded*luminosity[EnergyIndex] << endl;
			else cout << "RECORDED YIELD IN " << templateTree->GetName() << ": " << nTotalRecorded << endl;
			foutput->WriteTObject(templateTree);
			delete templateTree;
			delete tree;
		}
		if(Systematics!=0){
			foutput->WriteTObject(h1DVBFSigRatio);
			if(tFitD!=0) foutput->WriteTObject(h2DVBFSigRatio);
		}
		delete h2DVBFSigRatio;
		delete h1DVBFSigRatio;
		for(int tr=0;tr<4;tr++){
			delete tree_VBF[tr];
			delete h1DVBF[tr];
			if(tFitD!=0) delete h2DVBF[tr];
		}

		D_temp_1D[2]->Add(D_temp_1D[1],-1.0);
		D_temp_1D[2]->Add(D_temp_1D[0],-1.0);
		if(tFitD!=0){
			D_temp_2D[2]->Add(D_temp_2D[1],-1.0);
			D_temp_2D[2]->Add(D_temp_2D[0],-1.0);
		}
		cout << "Integrals after everything:\n1D\t2D" << endl;

		//Divides bins by Bin Width
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
					}
				}
			}
			if(tFitD!=0 && t==3){
				for(int binx=0;binx<nbinsx;binx++){
					double binwidthx = kDXarray[binx+1] - kDXarray[binx];
					for(int biny=0;biny<nbinsy;biny++){
						double binwidthy = kDYarray[biny+1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						double bincontent = hStore_qqZZ_Unconditional->GetBinContent(binx+1,biny+1);
						bincontent /= binwidth;
						hStore_qqZZ_Unconditional->SetBinContent(binx+1,biny+1,bincontent);
					}
				}
			}
			if(tFitD!=0 && t==4){
				for(int binx=0;binx<nbinsx;binx++){
					double binwidthx = kDXarray[binx+1] - kDXarray[binx];
					for(int biny=0;biny<nbinsy;biny++){
						double binwidthy = kDYarray[biny+1] - kDYarray[biny];
						double binwidth = binwidthx*binwidthy;
						double bincontent = hStore_ZX_Unconditional->GetBinContent(binx+1,biny+1);
						bincontent /= binwidth;
						hStore_ZX_Unconditional->SetBinContent(binx+1,biny+1,bincontent);
					}
				}
			}
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
				}
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
					}
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
						}
					}
					delete[] storeOriginal;
				}
				if(Systematics<0){
					double presmoothInt = D_temp_1D[t]->Integral();
					if(isSmooth) D_temp_1D[t]->Smooth(1,"k3a");
					double postsmoothInt = D_temp_1D[t]->Integral();
					D_temp_1D[t]->Scale(presmoothInt/postsmoothInt);
				}
			}
			foutput->WriteTObject(D_temp_1D[t]);
			if(tFitD!=0) foutput->WriteTObject(D_temp_2D[t]);
			if(tFitD!=0 && t==3) foutput->WriteTObject(hStore_qqZZ_Unconditional);
			if(tFitD!=0 && t==4) foutput->WriteTObject(hStore_ZX_Unconditional);

			if(t!=3 && t!=4){
				cout << D_temp_1D[t]->Integral("width")*luminosity[EnergyIndex] << '\t';
				if(tFitD!=0) cout << D_temp_2D[t]->Integral("width")*luminosity[EnergyIndex] << endl;
				else cout << endl;
			}
			else{
				cout << D_temp_1D[t]->Integral(1,nbinsx) << '\t';
				if(tFitD!=0) cout << D_temp_2D[t]->Integral(1,nbinsx,1,nbinsy) << endl;
				else cout << endl;
			}
		}
		if(tFitD!=0){
			cout << "Unconditional Integrals are: " << endl;
			cout << hStore_qqZZ_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
			cout << hStore_ZX_Unconditional->Integral("width")*luminosity[EnergyIndex] << endl;
		};
		delete hStore_qqZZ_Unconditional;
		delete hStore_ZX_Unconditional;

		finput_VBFRef->Close();
		foutput->Close();
	}
	delete tgkf;
	finput_KDFactor->Close();
	finput_VBF->Close();
	delete tg_interf;
}

double doLinearCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype){
	double aa[3] = { first, c_first, sqrt(c_first) };
	double bb[3] = { second, c_second, sqrt(c_second) };

	double sig = (aa[0] * bb[2] - aa[2] * bb[0] + (aa[2] - bb[2])*bkg) / (aa[1] * bb[2] - aa[2] * bb[1]);
	double interf = (aa[0] * bb[1] - aa[1] * bb[0] + (aa[1] - bb[1])*bkg) / (aa[2] * bb[1] - aa[1] * bb[2]);
	if (sig < 0){ sig = 0; interf = 0; };

	if (outputtype == kIntHist) return interf;
	else if(outputtype==kSigHist) return sig;
	else return 0;
}

double doPeakCombination(double first, double c_first, double second, double c_second, double bkg, int outputtype){
	double aa[3] = { first, 1, sqrt(c_first) };
	double bb[3] = { second, 1, sqrt(c_second) };

	double sig = (aa[0] * bb[2] - aa[2] * bb[0] + (aa[2] - bb[2])*bkg) / (aa[1] * bb[2] - aa[2] * bb[1]);
	double interf = (aa[0] * bb[1] - aa[1] * bb[0] + (aa[1] - bb[1])*bkg) / (aa[2] * bb[1] - aa[1] * bb[2]);
	if (sig < 0){ sig = 0; interf = 0; };

	cout << "VBF peak as calculated\n";
	cout << "Signal: " << sig << "\tInterf: " << interf << endl;
	
	if (outputtype == kIntHist) return interf;
	else if(outputtype==kSigHist) return sig;
	else return 0;
}


void oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, TH1F* finaloutput, int outputtype, TH1F* finaloutput2, int output2type){
	TH1F* output = (TH1F*) input->Clone();
	TH1F* output2 = (TH1F*) output->Clone();
	if(outputtype==kIntHist || outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
			if (outputtype == kSigHist){ delete output; output = (TH1F*) first->Clone(); }

			output2->Add(first, -25.0);
			output2->Add(second, -1.0);
			output2->Scale(0.2);
			if (output2type == kSigHist){ delete output2; output2 = (TH1F*) first->Clone(); }
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double bsi25 = output->GetBinContent(binx);

				double weight=doLinearCombination(bsi25,25,bsi,1,bkg,outputtype);
				output->SetBinContent(binx,weight);
				if (finaloutput2 != 0){
					double weight2 = doLinearCombination(bsi25, 25, bsi, 1, bkg, output2type);
					output2->SetBinContent(binx, weight2);
				}
			}
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi10 = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double bsi25 = output->GetBinContent(binx);

				double weight=doLinearCombination(bsi25,25,bsi10,10,bkg,outputtype);
				output->SetBinContent(binx,weight);
				if (finaloutput2 != 0){
					double weight2 = doLinearCombination(bsi25, 25, bsi10, 10, bkg, output2type);
					output2->SetBinContent(binx, weight2);
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double bsi10 = output->GetBinContent(binx);

				double weight=doLinearCombination(bsi10,10,bsi,1,bkg,outputtype);
				output->SetBinContent(binx,weight);
				if (finaloutput2 != 0){
					double weight2 = doLinearCombination(bsi10, 10, bsi, 1, bkg, output2type);
					output2->SetBinContent(binx, weight2);
				}
			}
		}
		for (int binx = 1; binx <= output->GetNbinsX(); binx++){
			finaloutput->SetBinContent(binx,output->GetBinContent(binx));
			if (finaloutput2 != 0) finaloutput2->SetBinContent(binx,output2->GetBinContent(binx));
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	delete output;
	delete output2;
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
			double scaleval = 1./(-50. + 25.*sqrt(10.));
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
			double scaleval = 1./(10 - sqrt(10));
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

void threeDlinearcombination(TH3F* first, int firsttype, TH3F* second, int secondtype, TH3F* input, int inputtype, TH3F* finaloutput, int outputtype, TH3F* finaloutput2, int output2type){
	TH3F* output = (TH3F*) input->Clone();
	TH3F* output2 = (TH3F*) output->Clone();
	if(outputtype==kIntHist || outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
			if (outputtype == kSigHist){ delete output; output = (TH3F*) first->Clone(); }

			output2->Add(first, -25.0);
			output2->Add(second, -1.0);
			output2->Scale(0.2);
			if (output2type == kSigHist){ delete output2; output2 = (TH3F*) first->Clone(); }
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
						double bsi = first->GetBinContent(binx, biny, binz);
						double bkg = second->GetBinContent(binx, biny, binz);
						double bsi25 = output->GetBinContent(binx, biny, binz);

						double weight = doLinearCombination(bsi25, 25, bsi, 1, bkg, outputtype);
						output->SetBinContent(binx, biny, binz, weight);
						if (finaloutput2 != 0){
							double weight2 = doLinearCombination(bsi25, 25, bsi, 1, bkg, output2type);
							output2->SetBinContent(binx, biny, binz, weight2);
						}
					}
				}
			}
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			double scaleval = 1./(-50. + 25.*sqrt(10.));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
						double bsi10 = first->GetBinContent(binx, biny, binz);
						double bkg = second->GetBinContent(binx, biny, binz);
						double bsi25 = output->GetBinContent(binx, biny, binz);

						double weight = doLinearCombination(bsi25, 25, bsi10, 10, bkg, outputtype);
						output->SetBinContent(binx, biny, binz, weight);
						if (finaloutput2 != 0){
							double weight2 = doLinearCombination(bsi25, 25, bsi10, 10, bkg, output2type);
							output2->SetBinContent(binx, biny, binz, weight2);
						}
					}
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			double scaleval = 1./(10 - sqrt(10));
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
						double bsi = first->GetBinContent(binx, biny, binz);
						double bkg = second->GetBinContent(binx, biny, binz);
						double bsi10 = output->GetBinContent(binx, biny, binz);

						double weight = doLinearCombination(bsi10, 10, bsi, 1, bkg, outputtype);
						output->SetBinContent(binx, biny, binz, weight);
						if (finaloutput2 != 0){
							double weight2 = doLinearCombination(bsi10, 10, bsi, 1, bkg, output2type);
							output2->SetBinContent(binx, biny, binz, weight2);
						}
					}
				}
			}
		}
		for (int binx = 1; binx <= output->GetNbinsX(); binx++){
			for (int biny = 1; biny <= output->GetNbinsY(); biny++){
				for (int binz = 1; binz <= output->GetNbinsZ(); binz++){
					finaloutput->SetBinContent(binx, biny, binz, output->GetBinContent(binx, biny, binz));
					if (finaloutput2 != 0) finaloutput2->SetBinContent(binx,biny,binz,output2->GetBinContent(binx,biny,binz));
				}
			}
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	delete output;
	delete output2;
}

void progressbar(int val, int tot){
  int percent=floor(0.01*tot);
  if(percent==0) percent=1;

  if(val%percent==0 && val!=tot){
    cout<<"[ "<<setw(3)<<val/percent<<"% |";
    for(int k=1;k<val/percent;k++) cout<<"=";
    if(val%percent!=100) cout<<">";
    for(int k=val/percent+1;k<100;k++) cout<<" ";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
  else if(val==tot){
    cout<<"[ 100% |";
    for(int k=0;k<100;k++) cout<<"=";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');        
  }
}
