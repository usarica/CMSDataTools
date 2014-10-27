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
#include "TF1.h"
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"

using namespace std;

//Initializers
int useAnomalousCouplings=kAddfLQ;
bool useDjettagging=true;
enum histtypes{kSigHist,kBkgHist,kIntHist,kBSIHist,kBSI10Hist,kBSI25Hist};
void makeCombineTemplates_Modified_MCFM_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth, int Djettag);
TH1F* oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, int outputtype);
TH2F* twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, int outputtype);
void progressbar(int val, int tot);


//Main Function, runs over all desired iterations
void makeCombineTemplates_Modified_MCFM(){
	bool isSmooth=false;
	const int kNumSyst=5;
	int systematics[kNumSyst]={0,1,-1,2,-2};
	for(int i=0;i<kNumSyst;++i){
		for(int usesmooth=0;usesmooth<2;++usesmooth){
			if(usesmooth==0) isSmooth=false;
			if(usesmooth==1) isSmooth=true;
			for(int CoM=7;CoM<9;++CoM){
				for(int channel=0;channel<3;++channel){
					if(useDjettagging){
						for(int Djettag=-1;Djettag<2;++Djettag) makeCombineTemplates_Modified_MCFM_one(channel,CoM,6,systematics[i],isSmooth,Djettag);
					}
					else makeCombineTemplates_Modified_MCFM_one(channel,CoM,6,systematics[i],isSmooth,0);	
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
// Djettag = -1,0,1 (events fail Djet > 0.5 tag, events are not tagged, events pass Djet > 0.5 tag respectively)
void makeCombineTemplates_Modified_MCFM_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth, int Djettag){
	char TREE_NAME[]="SelectedTree";
	TString INPUT_NAME = "HZZ4lTree_ggTo";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	if(useAnomalousCouplings>0) OUTPUT_NAME += strAnomalousType[useAnomalousCouplings] + "_";
	if(!isSmooth) OUTPUT_NAME += "Raw_";
	OUTPUT_NAME += TString(strFitDim[tFitD]) + "_";
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/",erg_tev);

	if (Djettag == 0){
		if (Systematics == 0) OUTPUT_NAME += "Nominal.root";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD.root";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD.root";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF.root";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF.root";
	}
	if (Djettag == -1){
		if (Systematics == 0) OUTPUT_NAME += "Nominal_nonDjet.root";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD_nonDjet.root";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD_nonDjet.root";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF_nonDjet.root";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF_nonDjet.root";
	}
	if (Djettag == 1){
		if (Systematics == 0) OUTPUT_NAME += "Nominal_Djet.root";
		if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD_Djet.root";
		if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD_Djet.root";
		if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF_Djet.root";
		if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF_Djet.root";
	}
	TString systname[5]={"SysDown_ggPDF","SysDown_ggQCD","Nominal","SysUp_ggQCD","SysUp_ggPDF"};
	TString djetname[3]={"Djet < 0.5","Nominal","Djet>=0.5"};
	TString cinput_common = user_gg2VV_location + erg_dir + "/" + user_folder[folder] + "/";
	TString cinput_common_qqZZ = user_gg2VV_location;
	TString cinput_common_ZX = user_gg2VV_location;

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float lowside[3]={220,230,240};
	float ZZMass_PeakCut[2]={105.6,140.6}; // Spin 0 analysis
	double ggZZ_Syst_AbsNormSyst[2][2] = { // EnergyIndex
		{ 0.0745, 0.0735 },
		{ 0.075, 0.072 }
	}; // QCD, PDF

	float MC_weight;
	float MC_weight_down;
	float MC_weight_up;
	float MC_weight_Kfactor=1;
	float MC_weight_ggZZLepInt=1;
	float MC_weight_Kfactor_Norm_down = 1;
	float MC_weight_Kfactor_Norm_up = 1;
	float MC_weight_PDF_Norm_down = 1;
	float MC_weight_PDF_Norm_up = 1;
	float GenHMass = 0;
	float ZZMass=0;
	float Djet=0;

	TString cinput_KDFactor = "./data/HZZ4l-KDFactorGraph";
	if (EnergyIndex == 0) cinput_KDFactor = cinput_KDFactor + "_7TeV";
	cinput_KDFactor = cinput_KDFactor + ".root";
	TFile* finput_KDFactor = new TFile(cinput_KDFactor, "read");
	TString tgkfname = "KDScale_AllFlavors_UnNormalized";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*)finput_KDFactor->Get(tgkfname);

	TString Djetcutfilename = "./data/HZZ4l-DjetCutShapes"; 
	if (EnergyIndex == 0) Djetcutfilename += "_7TeV";
	Djetcutfilename+=".root";
	TFile* Djetcutfile = new TFile(Djetcutfilename,"read");

	double mPOLE = 125.6;
 	double nSM_ScaledPeak[2][3]={
 		{1.0902689,0.6087736,1.4452079},
 		{5.1963998,2.6944639,6.6562629}
	};

	TString comstring;
	comstring.Form("%i",erg_tev);
	TString cinput_VBF_Sig = "./data/HZZ4l-125_6-" + comstring + "TeV-Sig_MCFM_PhantomVBF_Comparison.root";
	TFile* finput_VBF = new TFile(cinput_VBF_Sig,"read");
	TSpline3* tsp_VBF_Sig = (TSpline3*) finput_VBF->Get("Spline3");

	double overall_VBF_scale=1;

	//use VBF yields directly
	double VBF_Sig_Datacard[2][3]={
		{0.092458836,0.051755897,0.12861921},
		{0.46798807,0.24788553,0.61781689}
	};
	for (int e = 0; e < 2; e++){for (int ss = 0; ss < 3; ss++){ nSM_ScaledPeak[e][ss] /= luminosity[e]; VBF_Sig_Datacard[e][ss] /= luminosity[e]; }}

	for(int lo=0;lo<1;lo++){
		TString coutput_common = user_dir + erg_dir;
		gSystem->Exec("mkdir -p " + coutput_common);
		coutput_common += user_folder[folder] + "/";
		gSystem->Exec("mkdir -p " + coutput_common);


		cout<<"==============================="<<endl;
		cout<<"CoM Energy: "<<erg_tev<<" TeV"<<endl;
		cout<<"Decay Channel: "<<user_folder[folder]<<endl;
		cout<<"Systematic: "<<systname[Systematics+2]<<endl;
		cout<<"Djet cut: "<<djetname[Djettag+1]<<endl;
		cout<<endl;

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
			float binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		}

		const int kNumTemplates=8;
		TH1F** D_temp_1D[kNumTemplates];
		TH2F** D_temp_2D[kNumTemplates];
		TString templatenames[kNumTemplates]={"ggF Sig","gg Bkg","ggF Int","qqZZ","Z+X","VBF Sig","VBF Bkg","VBF Int"};
		// Build template structure, including anomalous couplings
		for (int t = 0; t < kNumTemplates; t++){
			if (t <= 2){ // gg(H)VV
				D_temp_1D[t] = new TH1F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
				D_temp_2D[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][0])];
			}
			else if (t >= 5 && t <= 7){ // VV(H)VV
				D_temp_1D[t] = new TH1F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
				D_temp_2D[t] = new TH2F*[(nAnomalousCouplingTemplates[useAnomalousCouplings][1])];
			}
			else{ // Anything else
				D_temp_1D[t] = new TH1F*[1];
				D_temp_2D[t] = new TH2F*[1];
			}
		}

		double nSig_Simulated=0;
		double nVBF_Sig_Simulated=0;
		float fitYval;
		MC_weight=1;
		MC_weight_down=1;
		MC_weight_up=1;
		MC_weight_Kfactor=1;
		MC_weight_ggZZLepInt=1;
		MC_weight_Kfactor_Norm_down = 1;
		MC_weight_Kfactor_Norm_up = 1;
		MC_weight_PDF_Norm_down = 1;
		MC_weight_PDF_Norm_up = 1;
		GenHMass = 0;
		ZZMass=0;

		double nMZZ220[3]={0};

		//Template filler
		for(int t=0;t<kNumTemplates;t++){
			TChain* tree = new TChain(TREE_NAME);
			int nAnomalousLoops = 1;
			if(t<=2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if(t>=5 && t<=7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			//Grab appropriate files for templates
			cout<<"Files: "<<endl;
			if(t==3){
				for(int b=0;b<kNumBkg;b++){
					TString cinput = cinput_common_qqZZ;
					cinput = cinput + erg_dir + "/" + user_folder[folder] + "/" + hzz4lprefix + sample_BackgroundFile[b] + "_Reprocessed.root";
					tree->Add(cinput);
					cout << cinput << endl;
				}
			}
			else if(t==4){
				TString cinput = cinput_common_ZX;
				cinput = cinput + erg_dir + "/" + user_folder[folder] + "/HZZ4lTree_DoubleOr_CRZLLTree_Reprocessed.root";
				tree->Add(cinput);
				cout << cinput << endl;
			}
			else{
				int tp;
				if(t<3) tp = t;
				if(t>4) tp = t-5;
				TString cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				tree->Add(cinput_4e);
				tree->Add(cinput_4mu);
				cout<<cinput_2e2mu<<endl;
				cout<<cinput_4e<<endl;
				cout<<cinput_4mu<<endl;
			}

			//Initialize templates
			TString templatename_1D_core;
			TString templatename_2D_core;
			TString templatename_1D;
			TString templatename_2D;
			if(t<2){
				templatename_1D_core = Form("T_1D_%i",t+1);
				templatename_2D_core = Form("T_2D_%i",t+1);
			}
			else if(t==2){
				templatename_1D_core = Form("T_1D_%i",4);
				templatename_2D_core = Form("T_2D_%i",4);
			}
			else if(t==3){
				templatename_1D_core = Form("T_1D_%s","qqZZ");
				templatename_2D_core = Form("T_2D_%s","qqZZ");
			}
			else if(t==4){
				templatename_1D_core = Form("T_1D_%s","ZX");
				templatename_2D_core = Form("T_2D_%s","ZX");
			}
			else if(t==5 || t==6){
				templatename_1D_core = Form("T_1D_VBF_%i",t-4);
				templatename_2D_core = Form("T_2D_VBF_%i",t-4);
			}
			else if(t==7){
				templatename_1D_core = Form("T_1D_VBF_%i",4);
				templatename_2D_core = Form("T_2D_VBF_%i",4);
			}
			for (int al = 0; al<nAnomalousLoops; al++){
				templatename_1D = templatename_1D_core;
				templatename_2D = templatename_2D_core;
				if (useAnomalousCouplings == kAddfLQ && al>0){
					templatename_1D.Append(Form("_mZZ2_%i", al));
					templatename_2D.Append(Form("_mZZ2_%i", al));
				}
				if (tFitD == 0){
					D_temp_1D[t][al] = new TH1F(templatename_1D, templatename_1D, nbinsx, kDXarray);
					D_temp_1D[t][al]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);
				}
				else{
					D_temp_1D[t][al] = new TH1F(templatename_1D, templatename_1D, nbinsy, kDYarray);
					D_temp_1D[t][al]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);

					D_temp_2D[t][al] = new TH2F(templatename_2D, templatename_2D, nbinsx, kDXarray, nbinsy, kDYarray);
					D_temp_2D[t][al]->GetXaxis()->SetTitle(strFitDim_label[0]);
					D_temp_2D[t][al]->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
				}
			}

			//Prepare trees
			if(tree->GetBranchStatus("GenHMass")) tree->SetBranchAddress("GenHMass",&GenHMass);
			tree->SetBranchAddress("ZZMass",&ZZMass);
			if(t!=4) tree->SetBranchAddress("MC_weight",&MC_weight);
			else tree->SetBranchAddress("ZXfake_weightProper",&MC_weight);
			if(t<3){
				tree->SetBranchAddress("MC_weight_Kfactor", &MC_weight_Kfactor);
				tree->SetBranchAddress("MC_weight_Kfactor_Norm_down", &MC_weight_Kfactor_Norm_down);
				tree->SetBranchAddress("MC_weight_Kfactor_Norm_up", &MC_weight_Kfactor_Norm_up);
				tree->SetBranchAddress("MC_weight_PDF_Norm_down", &MC_weight_PDF_Norm_down);
				tree->SetBranchAddress("MC_weight_PDF_Norm_up", &MC_weight_PDF_Norm_up);
				tree->SetBranchAddress("MC_weight_down", &MC_weight_down);
				tree->SetBranchAddress("MC_weight_up",&MC_weight_up);
			}
			else{
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
			if(t>4) tree->SetBranchAddress("Djet_VAJHU",&Djet);
			else Djet=0.;

			if(tree->GetBranchStatus(strFitDim[tFitD])) tree->SetBranchAddress(strFitDim[tFitD],&fitYval);
			else if(!tree->GetBranchStatus(strFitDim[tFitD])){
				cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
				fitYval=0;
			};
			cout << "Set variables in trees for " << templatenames[t] << endl;

			float nTotal=0;
			TF1* Djetcutshape=0;
			int nEntries = tree->GetEntries();
			if (t == 0) cout << tgkf->Eval(125.6) << endl;
			if(Djettag!=0){
				TString cutname;
				//IF USING MCFM FOR DJET CUT SHAPE
				/*if(t==0) cutname="MCFMSig_Djetcutshape";
				if(t==1) cutname="MCFMBkg_Djetcutshape";
				if(t==2) cutname="MCFMBSI_Djetcutshape";*/
				//IF USING MINLO FOR DJET CUT SHAPE
				if(t==0 || t==1 || t==2) cutname="MINLO_Djetcutshape";
				if(t==5) cutname="PhantomSig_Djetcutshape";
				if(t==6) cutname="PhantomBkg_Djetcutshape";
				if(t==7) cutname="PhantomBSI_Djetcutshape";				
				Djetcutshape = (TF1*) Djetcutfile->Get(cutname);
			}
			for(int ev=0;ev<nEntries;ev++){
				tree->GetEntry(ev);
			    progressbar(ev,tree->GetEntries());
				double weight = MC_weight;
				if (t<3) weight *= MC_weight_ggZZLepInt;
				if (abs(Systematics) != 1 && t<3) weight *= MC_weight_Kfactor;
				if (Systematics == -1 && t<3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_down*tgkf->Eval(125.6);
				if (Systematics == 1 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0])*MC_weight_Kfactor_Norm_up*tgkf->Eval(125.6);
				if (Systematics == -2 && t<3) weight *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_down;
				if (Systematics == 2 && t<3) weight *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1])*MC_weight_PDF_Norm_up;
				if (t>4) weight *= tsp_VBF_Sig->Eval(GenHMass)*MC_weight_ggZZLepInt;

				if(t==0 && ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]) nSig_Simulated += weight;
				if(t==5 && ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]) nVBF_Sig_Simulated += weight;
				if(t<3 && ZZMass>=ZZMass_cut[0]) nMZZ220[t] += weight;

				if( (ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) ) continue;

				//Code to add Djet reweighting for ggF (also reweighting for VBF using ggF samples)
				if(t<3 || t>4){
					if(Djettag==-1){
						weight *= 1-Djetcutshape->Eval(ZZMass);
					}
					if(Djettag==1){
						weight *= Djetcutshape->Eval(ZZMass);
					}
				}

				nTotal += weight; // Applies only to SM

				// Anomalous couplings loop
				for (int al = 0; al<nAnomalousLoops; al++){
					double fillWeight = weight;
					if (useAnomalousCouplings == kAddfLQ && al>0){
						fillWeight *= pow(GenHMass / mPOLE, 2 * al);
					}
					if (tFitD == 0) D_temp_1D[t][al]->Fill(ZZMass, fillWeight);
					else D_temp_1D[t][al]->Fill(fitYval, fillWeight);
					if (tFitD > 0) D_temp_2D[t][al]->Fill(ZZMass, fitYval, fillWeight);
				}
			}
			delete Djetcutshape;
			cout<<endl;
			cout << templatenames[t] << " NTotal: " << nTotal << endl;
			if(t<3) cout << nMZZ220[t]*nSM_ScaledPeak[EnergyIndex][folder]/nSig_Simulated*luminosity[EnergyIndex] << endl;
			if(t<3){
				double myscale = nSM_ScaledPeak[EnergyIndex][folder] / nSig_Simulated;
				if (Systematics == -1 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == 1 && t < 3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == -2 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				if (Systematics == 2 && t<3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				for (int al = 0; al<nAnomalousLoops; al++){
					D_temp_1D[t][al]->Scale(myscale);
					if (tFitD>0) D_temp_2D[t][al]->Scale(myscale);
					cout << "Scaling " << templatenames[t] << "(coupling:" << al << ") by " << myscale << endl;
				}
			}
			if(t==5){
				double myscale = VBF_Sig_Datacard[EnergyIndex][folder]/nVBF_Sig_Simulated;
				//myscale *= nSM_ScaledPeak[EnergyIndex][folder]/nSig_Simulated;
				overall_VBF_scale = myscale;
				for (int al = 0; al<nAnomalousLoops; al++){
					D_temp_1D[t][al]->Scale(myscale);
					if (tFitD>0) D_temp_2D[t][al]->Scale(myscale);
					cout << "Scaling " << templatenames[t] << "(coupling:" << al << ") by " << myscale << endl;
				}
			}
			if(t>5){
				double myscale = overall_VBF_scale;
				for (int al = 0; al<nAnomalousLoops; al++){
					D_temp_1D[t][al]->Scale(myscale);
					if (tFitD>0) D_temp_2D[t][al]->Scale(myscale);
					cout << "Scaling " << templatenames[t] << "(coupling:" << al << ") by " << myscale << endl;
				}
			}
			for (int al = 0; al<nAnomalousLoops; al++){
				double presmoothInt = D_temp_1D[t][al]->Integral("width");
				if (isSmooth) D_temp_1D[t][al]->Smooth(1, "k3a");
				double postsmoothInt = D_temp_1D[t][al]->Integral("width");
				D_temp_1D[t][al]->Scale(presmoothInt / postsmoothInt);
				if (tFitD>0){
					presmoothInt = D_temp_2D[t][al]->Integral("width");
					if (isSmooth) D_temp_2D[t][al]->Smooth(1, "k3a");
					postsmoothInt = D_temp_2D[t][al]->Integral("width");
					D_temp_2D[t][al]->Scale(presmoothInt / postsmoothInt);
				}

				//Conditionalize ZX and qqZZ
				if (t >= 3 && tFitD > 0 && t < 5){
					for (int binx = 0; binx <= D_temp_2D[t][al]->GetNbinsX()+1; binx++){
						double intBinX = D_temp_2D[t][al]->Integral(binx, binx, 0, D_temp_2D[t][al]->GetNbinsY()+1);
						for (int biny = 0; biny <= D_temp_2D[t][al]->GetNbinsY()+1; biny++){
							double bincontent = D_temp_2D[t][al]->GetBinContent(binx, biny);

							if (intBinX != 0) D_temp_2D[t][al]->SetBinContent(binx, biny, bincontent / intBinX);
						}
					}
				}
			}
			cout<<endl;

			delete tree;
		}

		//When pure samples aren't available, they are made via linear combinations of other samples.
		//2: ggF Int made using Sig, Bkg, and BSI samples 
		for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][0]; al++){
			D_temp_1D[2][al] = oneDlinearcombination(D_temp_1D[0][al], kSigHist, D_temp_1D[1][al], kBkgHist, D_temp_1D[2][al], kBSIHist, kIntHist);
			if (tFitD != 0) D_temp_2D[2][al] = twoDlinearcombination(D_temp_2D[0][al], kSigHist, D_temp_2D[1][al], kBkgHist, D_temp_2D[2][al], kBSIHist, kIntHist);
		}
		for (int al = 0; al < nAnomalousCouplingTemplates[useAnomalousCouplings][1]; al++){
			D_temp_1D[7][al] = oneDlinearcombination(D_temp_1D[5][al], kSigHist, D_temp_1D[6][al], kBkgHist, D_temp_1D[7][al], kBSIHist, kIntHist);
			if (tFitD != 0) D_temp_2D[7][al] = twoDlinearcombination(D_temp_2D[5][al], kSigHist, D_temp_2D[6][al], kBkgHist, D_temp_2D[7][al], kBSIHist, kIntHist);
		}
		cout << "Integrals after everything:\n1D\t2D" << endl;
		for(int t=0;t<kNumTemplates;t++){
			int nAnomalousLoops = 1;
			if(t<=2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if(t>=5 && t<=7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			for (int al = 0; al < nAnomalousLoops; al++){
				for (int binx = 0; binx < D_temp_1D[t][al]->GetNbinsX(); binx++){
					double binwidthx;
					if (tFitD == 0) binwidthx = kDXarray[binx + 1] - kDXarray[binx];
					else binwidthx = kDYarray[binx + 1] - kDYarray[binx];
					double bincontent = D_temp_1D[t][al]->GetBinContent(binx + 1);
					if (t != 3 && t != 4) bincontent /= binwidthx;
					D_temp_1D[t][al]->SetBinContent(binx + 1, bincontent);

					binwidthx = kDXarray[binx + 1] - kDXarray[binx];
				}
				if (tFitD != 0){
					for (int binx = 0; binx < D_temp_2D[t][al]->GetNbinsX(); binx++){
						for (int biny = 0; biny < D_temp_2D[t][al]->GetNbinsY(); biny++){
							double binwidthy = kDYarray[biny + 1] - kDYarray[biny];
							double binwidth = binwidthx*binwidthy;
							bincontent = D_temp_2D[t][al]->GetBinContent(binx + 1, biny + 1);
							if (t != 3 && t != 4) bincontent /= binwidth;
							D_temp_2D[t][al]->SetBinContent(binx + 1, biny + 1, bincontent);
						}
					}
				}
				if (t == 4 && Systematics != 0){
					double intTZX_1D = D_temp_1D[t][al]->Integral();
					double intTZZQQB_1D = D_temp_1D[3][al]->Integral();
					for (int binx = 1; binx <= D_temp_1D[t][al]->GetNbinsX(); binx++){
						double bincontent = D_temp_1D[t][al]->GetBinContent(binx);
						double bincontent_alt = D_temp_1D[3][al]->GetBinContent(binx);
						bincontent_alt *= intTZX_1D / intTZZQQB_1D;
						double difference = bincontent_alt - bincontent;

						if (Systematics > 0) bincontent += difference;
						else bincontent -= difference;
						if (bincontent < 0) bincontent = 0;
						D_temp_1D[t][al]->SetBinContent(binx, bincontent);
					}
					D_temp_1D[t][al]->Scale(intTZX_1D / D_temp_1D[t][al]->Integral());

					for (int binx = 1; binx <= D_temp_2D[t][al]->GetNbinsX(); binx++){
						for (int biny = 1; biny <= D_temp_2D[t][al]->GetNbinsY(); biny++){
							double bincontent = D_temp_2D[t][al]->GetBinContent(binx, biny);
							double bincontent_alt = D_temp_2D[3][al]->GetBinContent(binx, biny);
							double difference = bincontent_alt - bincontent;

							if (Systematics > 0) bincontent += difference;
							else bincontent -= difference;
							if (bincontent < 0) bincontent = 0;
							D_temp_2D[t][al]->SetBinContent(binx, biny, bincontent);
						}
						double intBinX = D_temp_2D[t][al]->Integral(binx, binx, 1, D_temp_2D[t][al]->GetNbinsY());
						for (int biny = 1; biny <= D_temp_2D[t][al]->GetNbinsY(); biny++){
							double bincontent = D_temp_2D[t][al]->GetBinContent(binx, biny);

							if (intBinX != 0) D_temp_2D[t][al]->SetBinContent(binx, biny, bincontent / intBinX);
						}
					}
					if (Systematics < 0){
						double presmoothInt = D_temp_1D[t][al]->Integral();
						if (isSmooth) D_temp_1D[t][al]->Smooth(1, "k3a");
						double postsmoothInt = D_temp_1D[t][al]->Integral();
						D_temp_1D[t][al]->Scale(presmoothInt / postsmoothInt);
					}
				}
				foutput->WriteTObject(D_temp_1D[t][al]);
				if (tFitD != 0) foutput->WriteTObject(D_temp_2D[t][al]);

				if (t != 3 && t != 4){
					cout << D_temp_1D[t][al]->Integral("width")*luminosity[EnergyIndex] << '\t';
					if (tFitD != 0) cout << D_temp_2D[t][al]->Integral("width")*luminosity[EnergyIndex] << endl;
					else cout << endl;
				}
				else{
					cout << D_temp_1D[t][al]->Integral(1, nbinsx) << '\t';
					if (tFitD != 0) cout << D_temp_2D[t][al]->Integral(1, nbinsx, 1, nbinsy) << endl;
					else cout << endl;
				}
			}
		}

		for (int t = 0; t < kNumTemplates; t++){
			int nAnomalousLoops = 1;
			if (t <= 2) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][0];
			else if (t >= 5 && t <= 7) nAnomalousLoops = nAnomalousCouplingTemplates[useAnomalousCouplings][1];

			for (int al = 0; al < nAnomalousLoops; al++){
				delete D_temp_1D[t][al];
				if (tFitD != 0) delete D_temp_2D[t][al];
			}
			delete[] D_temp_1D[t];
			if (tFitD != 0) delete[] D_temp_2D[t];
		}
		foutput->Close();
	}
	cout<<endl;
	delete tgkf;
	finput_KDFactor->Close();
	finput_VBF->Close();
};

TH1F* oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, int outputtype){
	TH1F* output = (TH1F*) input->Clone();
	if(outputtype==kIntHist){
		if(inputtype==kBSIHist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -1.0);
			output->Add(second, -1.0);
		}
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			output->Scale(-1.0);
			output->Add(first,25.0);
			output->Add(second,-24.0);
			output->Scale(0.05);
		}		
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			output->Scale(-10.0);
			output->Add(first,25.0);
			output->Add(second,-15.0);
			float scaleval = 1./(-50. + 25.*sqrt(10.));
			output->Scale(scaleval);
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			output->Scale(-1.0);
			output->Add(first,10.0);
			output->Add(second,-9.0);
			float scaleval = 1./(10 - sqrt(10));
			output->Scale(scaleval);
		}
	}
	else if(outputtype==kSigHist){
		//Need to force bins to be non-negative
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double sig = output->GetBinContent(binx);

				sig=(-5.0*bsi+4.0*bkg+sig)*(0.05);
				if (sig < 0) sig = 0;
				output->SetBinContent(binx,sig);
			}
		}
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi10 = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double sig = output->GetBinContent(binx);

				sig=(-5.*bsi10+(5.-sqrt(10))*bkg+sqrt(10)*sig)*(1./(-50+25*sqrt(10)));
				if (sig < 0) sig = 0;
				output->SetBinContent(binx,sig);
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				double bsi = first->GetBinContent(binx);
				double bkg = second->GetBinContent(binx);
				double sig = output->GetBinContent(binx);

				sig=(-sqrt(10)*bsi-(1.-sqrt(10))*bkg+sig)*(1./(10-sqrt(10)));
				if (sig < 0) sig = 0;
				output->SetBinContent(binx,sig);
			}
		}
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	return output;
}

TH2F* twoDlinearcombination(TH2F* first, int firsttype, TH2F* second, int secondtype, TH2F* input, int inputtype, int outputtype){
	TH2F* output = (TH2F*) input->Clone();
	if(outputtype==kIntHist){
		if(inputtype==kBSIHist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -1.0);
			output->Add(second, -1.0);
		}
		if(inputtype==kBSI25Hist && firsttype==kSigHist && secondtype==kBkgHist){
			output->Add(first, -25.0);
			output->Add(second, -1.0);
			output->Scale(0.2);
		}
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			output->Scale(-1.0);
			output->Add(first,25.0);
			output->Add(second,-24.0);
			output->Scale(0.05);
		}
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			output->Scale(-10.0);
			output->Add(first,25.0);
			output->Add(second,-15.0);
			float scaleval = 1./(-50.+25.*sqrt(10.));
			output->Scale(scaleval);
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			output->Scale(-1.0);
			output->Add(first,10.0);
			output->Add(second,-9.0);
			float scaleval = 1./(10-sqrt(10));
			output->Scale(scaleval);
		}
	}
	else if(outputtype==kSigHist){
		if(inputtype==kBSI25Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			//Need to force bins to be non-negative
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double sig = output->GetBinContent(binx,biny);

					sig=(-5.0*bsi+4.0*bkg+sig)*(0.05);
					if (sig < 0) sig = 0;
					output->SetBinContent(binx,biny,sig);
				}
			}
		}
		if(inputtype==kBSI25Hist && firsttype==kBSI10Hist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi10 = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double sig = output->GetBinContent(binx,biny);

					sig=(-5.*bsi10+(5.-sqrt(10))*bkg+sqrt(10)*sig)*(1./(-50+25*sqrt(10)));
					if (sig < 0) sig = 0;
					output->SetBinContent(binx,biny,sig);
				}
			}
		}
		if(inputtype==kBSI10Hist && firsttype==kBSIHist && secondtype==kBkgHist){
			for (int binx = 1; binx <= output->GetNbinsX(); binx++){
				for (int biny = 1; biny <= output->GetNbinsY(); biny++){
					double bsi = first->GetBinContent(binx,biny);
					double bkg = second->GetBinContent(binx,biny);
					double sig = output->GetBinContent(binx,biny);

					sig=(-sqrt(10)*bsi-(1.-sqrt(10))*bkg+sig)*(1./(10-sqrt(10)));
					if (sig < 0) sig = 0;
					output->SetBinContent(binx,biny,sig);
				}
			}
		}	
	}
	else{cout<<"Option not yet supported. Exiting..."<<endl; assert(0);};
	return output;
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

