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

using namespace std;

//Initializers
enum histtypes{kSigHist,kBkgHist,kIntHist,kBSIHist,kBSI10Hist,kBSI25Hist};
void makeCombineTemplates_Modified_MCFM_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth);
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
					makeCombineTemplates_Modified_MCFM_one(channel,CoM,6,systematics[i],isSmooth);	
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
void makeCombineTemplates_Modified_MCFM_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth){
	char TREE_NAME[]="SelectedTree";
	TString INPUT_NAME = "HZZ4lTree_ggTo";
	TString OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	if(!isSmooth) OUTPUT_NAME += "Raw_";
	OUTPUT_NAME += TString(strFitDim[tFitD]) + "_";
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/",erg_tev);

	if(Systematics==0) OUTPUT_NAME += "Nominal.root";
	if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD.root";
	if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD.root";
	if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF.root";
	if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF.root";
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

	TString cinput_KDFactor = "./data/HZZ4l-KDFactorGraph";
	if (EnergyIndex == 0) cinput_KDFactor = cinput_KDFactor + "_7TeV";
	cinput_KDFactor = cinput_KDFactor + ".root";
	TFile* finput_KDFactor = new TFile(cinput_KDFactor, "read");
	TString tgkfname = "KDScale_";
	tgkfname = tgkfname + "AllFlavors_UnNormalized";
	TGraphAsymmErrors* tgkf = (TGraphAsymmErrors*)finput_KDFactor->Get(tgkfname);


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
		cout<<endl;

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
			float binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
			kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
		};

		const int kNumTemplates=8;
		TH1F* D_temp_1D[kNumTemplates];
		TH2F* D_temp_2D[kNumTemplates];
		TString templatenames[kNumTemplates]={"ggF Sig","gg Bkg","ggF Int","qqZZ","Z+X","VBF Sig","VBF Bkg","VBF Int"};

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

			//Grab appropriate files for templates
			cout<<"Files: "<<endl;
			if(t==3){
				for(int b=0;b<kNumBkg;b++){
					TString cinput = cinput_common_qqZZ;
					cinput = cinput + erg_dir + "/" + user_folder[folder] + "/" + hzz4lprefix + sample_BackgroundFile[b] + "_Reprocessed.root";
					tree->Add(cinput);
					cout << cinput << endl;
				};
			}
			else if(t==4){
				TString cinput = cinput_common_ZX;
				cinput = cinput + erg_dir + "/" + user_folder[folder] + "/HZZ4lTree_DoubleOr_CRZLLTree_Reprocessed.root";
				tree->Add(cinput);
				cout << cinput << endl;
			}
			else if(t<3){
				int tp = t;
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
			else if(t>4){
				int tp = t-5;
				if(tp==2 && folder!=2 && EnergyIndex!=1) tp=3;
				if(folder==0 && tp==0) tp=2;
				cout<<user_folder[folder]<<" "<<t<<" "<<tp<<endl;
				TString cinput_2e2mu = cinput_common + "HZZ4lTree_ZZTo2e2muJJ_" + sample_suffix_Phantom[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + "HZZ4lTree_ZZTo4eJJ_" + sample_suffix_Phantom[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + "HZZ4lTree_ZZTo4muJJ_" + sample_suffix_Phantom[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				tree->Add(cinput_4e);
				tree->Add(cinput_4mu);
				cout<<cinput_2e2mu<<endl;
				cout<<cinput_4e<<endl;
				cout<<cinput_4mu<<endl;
			};

			//Initialize templates
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
			//cout<<templatename_1D<<endl;
			//cout<<templatename_2D<<endl;
			if(tFitD==0){
				D_temp_1D[t] = new TH1F(templatename_1D,templatename_1D,nbinsx,kDXarray);
				D_temp_1D[t]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);
			}
			else{
				D_temp_1D[t] = new TH1F(templatename_1D,templatename_1D,nbinsy,kDYarray);
				D_temp_1D[t]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);

				D_temp_2D[t] = new TH2F(templatename_2D,templatename_2D,nbinsx,kDXarray,nbinsy,kDYarray);
				D_temp_2D[t]->GetXaxis()->SetTitle(strFitDim_label[0]);
				D_temp_2D[t]->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
			};

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
			};
			if(t<3 || t>4) tree->SetBranchAddress("MC_weight_ggZZLepInt",&MC_weight_ggZZLepInt);
			else MC_weight_ggZZLepInt=1;

			if(tree->GetBranchStatus(strFitDim[tFitD])) tree->SetBranchAddress(strFitDim[tFitD],&fitYval);
			else if(!tree->GetBranchStatus(strFitDim[tFitD])){
				cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
				fitYval=0;
			};
			cout << "Set variables in trees for " << templatenames[t] << endl;

			float nTotal=0;
			int nEntries = tree->GetEntries();
			if (t == 0) cout << tgkf->Eval(125.6) << endl;
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

				nTotal += weight;

				if(tFitD==0) D_temp_1D[t]->Fill(ZZMass,weight);
				else D_temp_1D[t]->Fill(fitYval,weight);
				if(tFitD>0) D_temp_2D[t]->Fill(ZZMass,fitYval,weight);
			};
			cout<<endl;
			cout << templatenames[t] << " NTotal: " << nTotal << endl;
			if(t<3) cout << nMZZ220[t]*nSM_ScaledPeak[EnergyIndex][folder]/nSig_Simulated*luminosity[EnergyIndex] << endl;
			if(t<3){
				double myscale = nSM_ScaledPeak[EnergyIndex][folder] / nSig_Simulated;
				if (Systematics == -1 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == 1 && t < 3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == -2 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				if (Systematics == 2 && t<3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				D_temp_1D[t]->Scale(myscale);
				if (tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling " << templatenames[t] << " by " << myscale << endl;
			};
			if(t==5){
				double myscale = VBF_Sig_Datacard[EnergyIndex][folder]/nVBF_Sig_Simulated;
				myscale *= nSM_ScaledPeak[EnergyIndex][folder]/nSig_Simulated;
				overall_VBF_scale = myscale;
				D_temp_1D[t]->Scale(myscale);
				if(tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling " << templatenames[t] << " by " << myscale << endl;
			};
			if(t>5){
				double myscale = overall_VBF_scale;
				D_temp_1D[t]->Scale(myscale);
				if(tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling " << templatenames[t] << " by " << myscale << endl;
			};
			double presmoothInt = D_temp_1D[t]->Integral("width");
			if(isSmooth) D_temp_1D[t]->Smooth(1,"k3a");
			double postsmoothInt = D_temp_1D[t]->Integral("width");
			D_temp_1D[t]->Scale(presmoothInt/postsmoothInt);
			if(tFitD>0){
				presmoothInt = D_temp_2D[t]->Integral("width");
				if(isSmooth) D_temp_2D[t]->Smooth(1,"k3a");
				postsmoothInt = D_temp_2D[t]->Integral("width");
				D_temp_2D[t]->Scale(presmoothInt/postsmoothInt);
			};

			//Divides bins by bin width
			if(t>=3 && tFitD>0 && t<5){
				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);

						if(intBinX!=0) D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
					};
				};
			};

			cout<<endl;

			delete tree;
		};

		//When pure samples aren't available, they are made via linear combinations of other samples.
		//2: ggF Int made using Sig, Bkg, and BSI samples 
		D_temp_1D[2]=oneDlinearcombination(D_temp_1D[0],kSigHist,D_temp_1D[1],kBkgHist,D_temp_1D[2],kBSI25Hist,kIntHist);
		if(tFitD!=0) D_temp_2D[2]=twoDlinearcombination(D_temp_2D[0],kSigHist,D_temp_2D[1],kBkgHist,D_temp_2D[2],kBSI25Hist,kIntHist);
		//5: VBF Sig made using BSI, Bkg, and BSI25 samples
		//	 For 4mu samples, Bkg, BSI10, and BSI25 are used
		//	 For 8TeV 2e2mu, BSI, Bkg, and BSI10 are used
		TH1F* BSI=(TH1F*) D_temp_1D[5]->Clone();
		TH2F* BSI_2D;
		if (folder==0) D_temp_1D[5]=oneDlinearcombination(D_temp_1D[5],kBSI10Hist,D_temp_1D[6],kBkgHist,D_temp_1D[7],kBSI25Hist,kSigHist);
		else if (folder==2 && EnergyIndex==1) D_temp_1D[5]=oneDlinearcombination(D_temp_1D[5],kBSIHist,D_temp_1D[6],kBkgHist,D_temp_1D[7],kBSI10Hist,kSigHist);
		else D_temp_1D[5]=oneDlinearcombination(D_temp_1D[5],kBSIHist,D_temp_1D[6],kBkgHist,D_temp_1D[7],kBSI25Hist,kSigHist);
		D_temp_1D[5]->SetName("T_1D_VBF_1");
		if(tFitD!=0){
			BSI_2D=(TH2F*) D_temp_2D[5]->Clone();
			if (folder==0) D_temp_2D[5]=twoDlinearcombination(D_temp_2D[5],kBSI10Hist,D_temp_2D[6],kBkgHist,D_temp_2D[7],kBSI25Hist,kSigHist);
			else if (folder==2 && EnergyIndex==1) D_temp_2D[5]=twoDlinearcombination(D_temp_2D[5],kBSIHist,D_temp_2D[6],kBkgHist,D_temp_2D[7],kBSI10Hist,kSigHist);
			else D_temp_2D[5]=twoDlinearcombination(D_temp_2D[5],kBSIHist,D_temp_2D[6],kBkgHist,D_temp_2D[7],kBSI25Hist,kSigHist);
			D_temp_2D[5]->SetName("T_2D_VBF_1");
		}
		//7: VBF Int made using BSI, Bkg, and BSI25 samples
		//	 For 4mu samples, Bkg, BSI10, and BSI25 are used
		//	 For 8TeV 2e2mu, BSI, Bkg, and BSI10 are used
		if(folder==0) D_temp_1D[7]=oneDlinearcombination(BSI,kBSI10Hist,D_temp_1D[6],kBkgHist,D_temp_1D[7],kBSI25Hist,kIntHist);
		else if(folder==2 && EnergyIndex==1) D_temp_1D[7]=oneDlinearcombination(BSI,kBSIHist,D_temp_1D[6],kBkgHist,D_temp_1D[7],kBSI10Hist,kIntHist);
		else D_temp_1D[7]=oneDlinearcombination(BSI,kBSIHist,D_temp_1D[6],kBkgHist,D_temp_1D[7],kBSI25Hist,kIntHist);
		D_temp_1D[7]->SetName("T_1D_VBF_4");
		if(tFitD!=0){
			if(folder==0) D_temp_2D[7]=twoDlinearcombination(BSI_2D,kBSI10Hist,D_temp_2D[6],kBkgHist,D_temp_2D[7],kBSI25Hist,kIntHist);
			else if(folder==2 && EnergyIndex==1) D_temp_2D[7]=twoDlinearcombination(BSI_2D,kBSIHist,D_temp_2D[6],kBkgHist,D_temp_2D[7],kBSI10Hist,kIntHist);
			else D_temp_2D[7]=twoDlinearcombination(BSI_2D,kBSIHist,D_temp_2D[6],kBkgHist,D_temp_2D[7],kBSI25Hist,kIntHist);
			D_temp_2D[7]->SetName("T_2D_VBF_4");
		}

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
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);
						double bincontent_alt = D_temp_2D[3]->GetBinContent(binx,biny);
						double difference = bincontent_alt - bincontent;

						if(Systematics>0) bincontent += difference;
						else bincontent -= difference;
						if(bincontent < 0) bincontent = 0;
						D_temp_2D[t]->SetBinContent(binx,biny,bincontent);
					};
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);

						if(intBinX!=0) D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
					};
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

			if(t!=3 && t!=4){
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
		foutput->Close();
	};
	cout<<endl;
	delete tgkf;
	finput_KDFactor->Close();
	finput_VBF->Close();
};

TH1F* oneDlinearcombination(TH1F* first, int firsttype, TH1F* second, int secondtype, TH1F* input, int inputtype, int outputtype){
	TH1F* output = (TH1F*) input->Clone();
	if(outputtype==kIntHist){
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

