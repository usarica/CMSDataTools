#include <iostream>
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
void makeCombineTemplates_Modified_MCFM_one(int folder, int erg_tev, int tFitD, int Systematics, bool isSmooth);

//Main Function, runs over all desired iterations
void makeCombineTemplates_Modified_MCFM(){
	bool isSmooth=false;
	int systematics[5]={0,1,-1,2,-2};
	for(int i=0;i<5;++i){
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
	TString INPUT_SM_NAME = "HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Reprocessed.root";
	TString INPUT_SM_126_NAME = "HZZ4lTree_powheg15jhuGenV3H126_Reprocessed.root";
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
	float ZZMass_PeakCut[2]={120,130}; // Spin 0 analysis
	double ggZZ_Syst_AbsNormSyst[2][2] = { // EnergyIndex
		{ 0.0745, 0.0735 },
		{ 0.075, 0.072 }
	}; // QCD, PDF

	int genFinalState;
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
	cinput_sm = cinput_common + INPUT_SM_126_NAME;
	TFile* fsm2 = new TFile(cinput_sm,"read");
	TTree* tsm2 = (TTree*) fsm2->Get(TREE_NAME);
	tsm2->SetBranchAddress("genFinalState",&genFinalState);
	tsm2->SetBranchAddress("ZZMass",&ZZMass);
	tsm2->SetBranchAddress("MC_weight",&MC_weight);
	double nTotalSM126_notTau=0;
	double nTotalSM126_Tau=0;
	for(int sev=0;sev<tsm2->GetEntries();sev++){
		tsm2->GetEntry(sev);
		if(ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]){
			if(genFinalState!=4) nTotalSM126_notTau += MC_weight;
			else nTotalSM126_Tau += MC_weight;
		};
	};
	fsm2->Close();
	nTotalSM126_notTau *= (BR_Table[EnergyIndex][1][1]*XSEC_Table[EnergyIndex][1]) / ( BR_Table[EnergyIndex][2][1]*XSEC_Table[EnergyIndex][2] );
	nTotalSM126_Tau *= ( (BR_Table[EnergyIndex][1][0]-BR_Table[EnergyIndex][1][1])*XSEC_Table[EnergyIndex][1]) / ( (BR_Table[EnergyIndex][2][0]-BR_Table[EnergyIndex][2][1])*XSEC_Table[EnergyIndex][2] );
	double tauScale = (nTotalSM126_notTau + nTotalSM126_Tau) / nTotalSM126_notTau;
	cout << "Tau scale is " << tauScale << endl;
	double nSM_ObservedPeak = nTotalSM*tauScale;
	double nSM_ScaledPeak = nSM_ObservedPeak;
	cout << "Observed number of peak events is " << nSM_ObservedPeak*luminosity[EnergyIndex] << endl;
	cout << "Scaled number of peak events is " << nSM_ScaledPeak*luminosity[EnergyIndex] << endl;

	TString comstring;
	comstring.Form("%i",erg_tev);
	TString cinput_VBF_Sig = "./data/HZZ4l-125_6-" + comstring + "TeV-Sig_MCFM_PhantomVBF_Comparison.root";
	TFile* finput_VBF = new TFile(cinput_VBF_Sig,"read");
	TSpline3* tsp_VBF_Sig = (TSpline3*) finput_VBF->Get("Spline3");
	double VBF_Sig_Datacard_Ratio[2][3]={
		{
			0.094235231222,
			0.094839434433,
			0.097909748221
		},
		{
			0.0891967993747,
			0.0907358465455,
			0.0912144438605
		}
	};
	double overall_VBF_scale=1;

	for(int lo=0;lo<1;lo++){
		TString coutput_common = user_dir + erg_dir;
		gSystem->Exec("mkdir -p " + coutput_common);
		coutput_common += user_folder[folder] + "/";
		gSystem->Exec("mkdir -p " + coutput_common);

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

		int nEntries;
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
				if (EnergyIndex == 0){ // USE BSI25, NO CONTAMINATION SIGNIFICANT TO CREATE STORAGE TREES
					if (folder == 2 && tp == 2) tp = 3;
					if (folder == 1 && tp == 1) tp = 3;
				};
				TString cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				tree->Add(cinput_4e);
				tree->Add(cinput_4mu);
			}
			else if(t>4){
				int tp = t-5;
				if (EnergyIndex == 0){ // USE BSI25, NO CONTAMINATION SIGNIFICANT TO CREATE STORAGE TREES
					if (folder == 2 && tp == 2) tp = 3;
					if (folder == 1 && tp == 1) tp = 3;
				};
				TString cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";
				TString cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[tp] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu);
				cout << cinput_2e2mu << endl;
				tree->Add(cinput_4e);
				cout << cinput_4e << endl;
				tree->Add(cinput_4mu);
				cout << cinput_4mu << endl;
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
			cout << "Set variables in trees for " << t << endl;

			float nTotal=0;
			nEntries = tree->GetEntries();
			if (t == 0) cout << tgkf->Eval(125.6) << endl;
			for(int ev=0;ev<nEntries;ev++){
				tree->GetEntry(ev);

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
			cout << t << " NTotal: " << nTotal << endl;
			if(t<3) cout << nMZZ220[t]*nSM_ScaledPeak/nSig_Simulated*luminosity[EnergyIndex] << endl;
			if(t<3){
				double myscale = nSM_ScaledPeak / nSig_Simulated;
				if (Systematics == -1 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == 1 && t < 3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][0]);
				if (Systematics == -2 && t < 3) myscale *= (1.0 - ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				if (Systematics == 2 && t<3) myscale *= (1.0 + ggZZ_Syst_AbsNormSyst[EnergyIndex][1]);
				D_temp_1D[t]->Scale(myscale);
				if (tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling t: " << t << " by " << myscale << endl;
			};
			if(t==5){
				double myscale = VBF_Sig_Datacard_Ratio[EnergyIndex][folder];
				myscale *= nSM_ScaledPeak/nVBF_Sig_Simulated;
				overall_VBF_scale = myscale;
				D_temp_1D[t]->Scale(myscale);
				if(tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling t: " << t << " by " << myscale << endl;
			};
			if(t>5){
				double myscale = overall_VBF_scale;
				D_temp_1D[t]->Scale(myscale);
				if(tFitD>0) D_temp_2D[t]->Scale(myscale);
				cout << "Scaling t: " << t << " by " << myscale << endl;
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

			if(t>=3 && tFitD>0 && t<5){
				for(int binx=1;binx<=D_temp_2D[t]->GetNbinsX();binx++){
					double intBinX = D_temp_2D[t]->Integral(binx,binx,1,D_temp_2D[t]->GetNbinsY());
					for(int biny=1;biny<=D_temp_2D[t]->GetNbinsY();biny++){
						double bincontent = D_temp_2D[t]->GetBinContent(binx,biny);

						if(intBinX!=0) D_temp_2D[t]->SetBinContent(binx,biny,bincontent/intBinX);
					};
				};
			};

			delete tree;
		};

		//Makes Bkg, Sig, Int from linear combinations of above templates
		if (EnergyIndex == 0){ // USE BSI25, NO CONTAMINATION SIGNIFICANT TO CREATE STORAGE TREES
			if (folder == 2){
				//MAKE FUNCTION
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
				D_temp_1D[7]->Add(D_temp_1D[6], -1.0);
				D_temp_1D[7]->Add(D_temp_1D[5], -25.0);
				D_temp_1D[7]->Scale(0.2);
				D_temp_1D[7]->Add(D_temp_1D[6], 1.0);
				D_temp_1D[7]->Add(D_temp_1D[5], 1.0);
				if (tFitD != 0){
					D_temp_2D[7]->Add(D_temp_2D[6], -1.0);
					D_temp_2D[7]->Add(D_temp_2D[5], -25.0);
					D_temp_2D[7]->Scale(0.2);
					D_temp_2D[7]->Add(D_temp_2D[6], 1.0);
					D_temp_2D[7]->Add(D_temp_2D[5], 1.0);
				};
			};

			//Divides bins by bin width
			if (folder == 1){
				for (int binx = 1; binx <= D_temp_1D[1]->GetNbinsX(); binx++){
					double sig = D_temp_1D[0]->GetBinContent(binx);
					double bsi25 = D_temp_1D[1]->GetBinContent(binx);
					double bsi = D_temp_1D[2]->GetBinContent(binx);

					bsi = 5.0*(bsi-sig);
					bsi25 = (bsi25-25.0*sig-bsi)*(-0.25);
					if (bsi25 < 0) bsi25 = 0;
					D_temp_1D[1]->SetBinContent(binx,bsi25);
				};
				if (tFitD != 0){
					for (int binx = 1; binx <= D_temp_2D[1]->GetNbinsX(); binx++){
						for (int biny = 1; biny <= D_temp_2D[1]->GetNbinsY(); biny++){
							double sig = D_temp_2D[0]->GetBinContent(binx,biny);
							double bsi25 = D_temp_2D[1]->GetBinContent(binx,biny);
							double bsi = D_temp_2D[2]->GetBinContent(binx,biny);

							bsi = 5.0*(bsi - sig);
							bsi25 = (bsi25 - 25.0*sig - bsi)*(-0.25);
							if (bsi25 < 0) bsi25 = 0;
							D_temp_2D[1]->SetBinContent(binx, biny, bsi25);
						};
					};
				};

				for (int binx = 1; binx <= D_temp_1D[6]->GetNbinsX(); binx++){
					double sig = D_temp_1D[5]->GetBinContent(binx);
					double bsi25 = D_temp_1D[6]->GetBinContent(binx);
					double bsi = D_temp_1D[7]->GetBinContent(binx);

					bsi = 5.0*(bsi - sig);
					bsi25 = (bsi25 - 25.0*sig - bsi)*(-0.25);
					if (bsi25 < 0) bsi25 = 0;
					D_temp_1D[6]->SetBinContent(binx, bsi25);
				};
				if (tFitD != 0){
					for (int binx = 1; binx <= D_temp_2D[6]->GetNbinsX(); binx++){
						for (int biny = 1; biny <= D_temp_2D[6]->GetNbinsY(); biny++){
							double sig = D_temp_2D[5]->GetBinContent(binx, biny);
							double bsi25 = D_temp_2D[6]->GetBinContent(binx, biny);
							double bsi = D_temp_2D[7]->GetBinContent(binx, biny);

							bsi = 5.0*(bsi - sig);
							bsi25 = (bsi25 - 25.0*sig - bsi)*(-0.25);
							if (bsi25 < 0) bsi25 = 0;
							D_temp_2D[6]->SetBinContent(binx, biny, bsi25);
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
		D_temp_1D[7]->Add(D_temp_1D[6],-1.0);
		D_temp_1D[7]->Add(D_temp_1D[5],-1.0);
		if(tFitD!=0){
			D_temp_2D[7]->Add(D_temp_2D[6],-1.0);
			D_temp_2D[7]->Add(D_temp_2D[5],-1.0);
		};
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
	delete tgkf;
	finput_KDFactor->Close();
	finput_VBF->Close();
};
