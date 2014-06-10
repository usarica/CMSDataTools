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

void makeCombineTemplates_Modified_MCFM_GenLevelVBF(int folder, int erg_tev, int tFitD=0, int Systematics=0, bool useResoVBF=true){
	char TREE_NAME[]="SelectedTree";
	string INPUT_SM_NAME = "HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Reprocessed.root";
	string INPUT_SM_126_NAME = "HZZ4lTree_powheg15jhuGenV3H126_Reprocessed.root";
	string INPUT_NAME = "HZZ4lTree_ggTo";
	string INPUT_NAME_VBF = "HZZ4l-125_6-";
	string OUTPUT_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	OUTPUT_NAME = OUTPUT_NAME + "_GenLevelVBF_";
	if(useResoVBF) OUTPUT_NAME = OUTPUT_NAME + "wResolution_";
	OUTPUT_NAME = OUTPUT_NAME + strFitDim[tFitD] + "_";
	char cerg[1000];
	sprintf(cerg,"%iTeV",erg_tev);
	char erg_dir[1000];
	sprintf(erg_dir,"LHC_%iTeV",erg_tev);
	if(Systematics==0) OUTPUT_NAME = OUTPUT_NAME + "Nominal.root";
	if(Systematics==1) OUTPUT_NAME = OUTPUT_NAME + "SysUp.root";
	if(Systematics==-1) OUTPUT_NAME = OUTPUT_NAME + "SysDown.root";

	string INPUT_VBFREF_NAME = "HtoZZ4l_MCFM_125p6_ModifiedTemplatesForCombine_";
	INPUT_VBFREF_NAME = INPUT_VBFREF_NAME + strFitDim[tFitD] + "_";
	if(Systematics==0) INPUT_VBFREF_NAME = INPUT_VBFREF_NAME + "Nominal.root";
	if(Systematics==1) INPUT_VBFREF_NAME = INPUT_VBFREF_NAME + "SysUp.root";
	if(Systematics==-1) INPUT_VBFREF_NAME = INPUT_VBFREF_NAME + "SysDown.root";

	INPUT_NAME_VBF = INPUT_NAME_VBF + cerg + "-";
	char* sample_VBF_suffix[4]={
		"Bkg_VBF_Phantom",
		"BSI10_VBF_Phantom",
		"BSI25_VBF_Phantom",
		"BSI_VBF_Phantom"
	};

	string cinput_common = user_gg2VV_location + erg_dir + "/" + user_folder[folder] + "/";
	string cinput_common_recover;
	if(folder!=2) cinput_common_recover = user_gg2VV_location + erg_dir + "/" + user_folder[2] + "/";
	string cinput_common_qqZZ = user_gg2VV_location;
	string cinput_common_ZX = user_gg2VV_location;
	string cinput_common_VBF = user_dir;

	int EnergyIndex=1;
	if(erg_tev==7) EnergyIndex=0;
	float lowside[3]={220,230,240};
	float ZZMass_PeakCut[2]={120,130}; // Spin 0 analysis

	int genFinalState;
	int isSelected;
	float MC_weight;
	float MC_weight_down;
	float MC_weight_up;
	float MC_weight_Kfactor=1;
	float MC_weight_ggZZLepInt=1;
	float GenHMass=0;
	float ZZMass=0;
	float Z1Mass=0;
	float Z2Mass=0;

	string cinput_sm = cinput_common + INPUT_SM_NAME;
	TFile* fsm = new TFile(cinput_sm.c_str(),"read");
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
	TFile* fsm2 = new TFile(cinput_sm.c_str(),"read");
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
	double InterferenceScale=1;
	if(folder!=2) InterferenceScale = BR_Table[EnergyIndex][1][3]/(2.0*BR_Table[EnergyIndex][1][2]);
//	float nSM_ObservedPeak = nSMValues[EnergyIndex][folder];
	double nSM_ObservedPeak = nTotalSM*tauScale;
//	nSM_ObservedPeak *= (BR_Table[EnergyIndex][1][0]*XSEC_Table[EnergyIndex][1]) / ( BR_Table[EnergyIndex][2][0]*XSEC_Table[EnergyIndex][2] );
//	if(folder!=2) nSM_ObservedPeak *= (BR_Table[EnergyIndex][1][2]*XSEC_Table[EnergyIndex][1]) / ( BR_Table[EnergyIndex][2][2]*XSEC_Table[EnergyIndex][2] );
//	else nSM_ObservedPeak *= (BR_Table[EnergyIndex][1][3]*XSEC_Table[EnergyIndex][1]) / ( BR_Table[EnergyIndex][2][3]*XSEC_Table[EnergyIndex][2] );
//	float nSM_ScaledPeak = nSM_ObservedPeak*InterferenceScale;
	double nSM_ScaledPeak = nSM_ObservedPeak;
	cout << "Observed number of peak events is " << nSM_ObservedPeak*luminosity[EnergyIndex] << endl;
	cout << "Scaled number of peak events is " << nSM_ScaledPeak*luminosity[EnergyIndex] << endl;
	TGraph* tg_interf = make_HZZ_LeptonInterferenceGraph();

	char cinput_VBF_Sig[1000];
	sprintf(cinput_VBF_Sig,"%s%i%s","./data/HZZ4l-125_6-",erg_tev,"TeV-Sig_gg2VV_PhantomVBF_Comparison.root");
	TFile* finput_VBF = new TFile(cinput_VBF_Sig,"read");
	TSpline3* tsp_VBF_Sig = (TSpline3*) finput_VBF->Get("Spline3");
	double VBF_Sig_Datacard_Ratio[2][3]={
		{
			0.0891967993747, //
			0.0907358465455, //
			0.0912144438605 //
		},
		{
			0.0891967993747,
			0.0907358465455,
			0.0912144438605
		}
	};
	double VBF_Phantom_Expected[2][3]={
		{
			0.446243,
			0.217077,
			0.564907
		},
		{
			0.446243,
			0.217077,
			0.564907
		}
	};
	double overall_VBF_scale=1;

//	for(int lo=0;lo<3;lo++){
	for(int lo=0;lo<1;lo++){
		char clowside[1000];
		sprintf(clowside,"%.0f",lowside[lo]);

		string user_core = user_dir;
		string coutput_common = user_core + erg_dir;
		coutput_common = coutput_common + "/Analysis/Templates/Combine/" + user_folder[folder] + "/" + clowside + "/";

		string coutput = coutput_common + OUTPUT_NAME;
		TFile* foutput = new TFile(coutput.c_str(),"recreate");

		string cinput_VBFRef = coutput_common + INPUT_VBFREF_NAME;
		TFile* finput_VBFRef = new TFile(cinput_VBFRef.c_str(),"read");
		double nVBF_Sig_Reweighted=0;
		TH2F* hVBFRef = (TH2F*) finput_VBFRef->Get("T_2D_VBF_1");
		nVBF_Sig_Reweighted = hVBFRef->Integral("width")*luminosity[EnergyIndex];
		finput_VBFRef->Close();

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
		float nSig_Simulated=0;
		float nVBF_Sig_Simulated=0;
		float fitYval;
		isSelected=1;
		MC_weight=1;
		MC_weight_down=1;
		MC_weight_up=1;
		MC_weight_Kfactor=1;
		MC_weight_ggZZLepInt=1;
		GenHMass=0;
		ZZMass=0;
		Z1Mass=0;
		Z2Mass=0;

		double nMZZ220[3]={0};

		TChain* tree_VBF[4];
		TH1F* h1DVBF[4];
		TH2F* h2DVBF[4];
		for(int tr=0;tr<4;tr++){
			tree_VBF[tr] = new TChain(TREE_NAME);
		};
		string cinput_VBF_Bkg = cinput_common_VBF + "VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[0];
		string cinput_VBF_BSI = cinput_common_VBF + "VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[3];
		string cinput_VBF_BSI10 = cinput_common_VBF + "VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[1];
		string cinput_VBF_BSI25 = cinput_common_VBF + "VBF/" + erg_dir + "/" + user_folder[folder] + "/" + INPUT_NAME_VBF + sample_VBF_suffix[2];
		if(useResoVBF){
			cinput_VBF_Bkg = cinput_VBF_Bkg + "_Reprocessed_wResolution.root";
			cinput_VBF_BSI = cinput_VBF_BSI + "_Reprocessed_wResolution.root";
			cinput_VBF_BSI10 = cinput_VBF_BSI10 + "_Reprocessed_wResolution.root";
			cinput_VBF_BSI25 = cinput_VBF_BSI25 + "_Reprocessed_wResolution.root";
		}
		else{
			cinput_VBF_Bkg = cinput_VBF_Bkg + "_Reprocessed.root";
			cinput_VBF_BSI = cinput_VBF_BSI + "_Reprocessed.root";
			cinput_VBF_BSI10 = cinput_VBF_BSI10 + "_Reprocessed.root";
			cinput_VBF_BSI25 = cinput_VBF_BSI25 + "_Reprocessed.root";
		};
		tree_VBF[0]->Add(cinput_VBF_Bkg.c_str());
		tree_VBF[1]->Add(cinput_VBF_BSI10.c_str());
		tree_VBF[2]->Add(cinput_VBF_BSI25.c_str());
		tree_VBF[3]->Add(cinput_VBF_BSI.c_str());

		for(int tr=0;tr<4;tr++){
			tree_VBF[tr]->SetBranchAddress("isSelected",&isSelected);
			tree_VBF[tr]->SetBranchAddress("GenZZMass",&ZZMass);
			tree_VBF[tr]->SetBranchAddress("MC_weight",&MC_weight);
			if(tFitD!=0) tree_VBF[tr]->SetBranchAddress(strFitDim[tFitD],&fitYval);

			string templatename_1D = "htemp1DVBF_";
			templatename_1D = templatename_1D + sample_VBF_suffix[tr];
			string templatename_2D = "htemp2DVBF_";
			templatename_2D = templatename_2D + sample_VBF_suffix[tr];

			if(tFitD==0){
				h1DVBF[tr] = new TH1F(templatename_1D.c_str(),templatename_1D.c_str(),nbinsx,kDXarray);
			}
			else{
				h1DVBF[tr] = new TH1F(templatename_1D.c_str(),templatename_1D.c_str(),nbinsy,kDYarray);
				h2DVBF[tr] = new TH2F(templatename_2D.c_str(),templatename_2D.c_str(),nbinsx,kDXarray,nbinsy,kDYarray);
			};

			int nVBFEntries = tree_VBF[tr]->GetEntries();
			cout << "Here: "  << tr << " N: " << nVBFEntries << endl;

			for(int ev=0;ev<nVBFEntries;ev++){
				tree_VBF[tr]->GetEntry(ev);
				if(isSelected!=1) continue;

				double weight = MC_weight;
				if(ev==0) cout << "Weight: " << weight << endl;
				if(folder==2) weight /= 2.0;
				else weight *= (tg_interf->Eval(ZZMass))/4.0;

				if( (ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) ) continue;

				if(tFitD==0) h1DVBF[tr]->Fill(ZZMass,weight);
				else{
					h1DVBF[tr]->Fill(fitYval,weight);
					h2DVBF[tr]->Fill(ZZMass,fitYval,weight);
				};
			};
//			if(tr==3) cout << "Phantom Expected: " << nVBF_Sig_Simulated << endl;

			ZZMass=0;
			MC_weight=1;
			fitYval=0;
			delete tree_VBF[tr];

			double presmoothInt = h1DVBF[tr]->Integral();
			h1DVBF[tr]->Smooth(1,"k3a");
			double postsmoothInt = h1DVBF[tr]->Integral();
			h1DVBF[tr]->Scale(presmoothInt/postsmoothInt);

			foutput->WriteTObject(h1DVBF[tr]);
			if(tFitD!=0){
				presmoothInt = h2DVBF[tr]->Integral();
				h2DVBF[tr]->Smooth(1,"k3a");
				postsmoothInt = h2DVBF[tr]->Integral();
				h2DVBF[tr]->Scale(presmoothInt/postsmoothInt);
				foutput->WriteTObject(h2DVBF[tr]);
			};
		};
		for(int binx=1;binx<=nbinsx;binx++){
			double bkg = h1DVBF[0]->GetBinContent(binx);
			double bsi10 = h1DVBF[1]->GetBinContent(binx);
			double bsi25 = h1DVBF[2]->GetBinContent(binx);
			double bsi25p = bsi25 - bkg;
			double bsi10p = bsi10 - bkg;
			double signal = ( sqrt(10.0)*bsi25p - 5.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );
			double interf = ( -10.0*bsi25p + 25.0*bsi10p ) / ( 25.0*sqrt(10.0) - 50.0 );

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

		for(int t=0;t<kNumTemplates;t++){
			TChain* tree = new TChain(TREE_NAME);

			if(t==3){
				for(int b=0;b<kNumBkg;b++){
					string cinput = cinput_common_qqZZ;
					cinput = cinput + erg_dir + "/" + user_folder[folder] + "/" + hzz4lprefix + sample_BackgroundFile[b] + "_Reprocessed.root";
					tree->Add(cinput.c_str());
					cout << cinput << endl;
				};
			}
			else if(t==4){
				string cinput = cinput_common_ZX;
				cinput = cinput + erg_dir + "/" + user_folder[folder] + "/HZZ4lTree_DoubleOr_CRZLLTree_Reprocessed.root";
				tree->Add(cinput.c_str());
				cout << cinput << endl;
			}
			else if(t<3){
				string cinput_2e2mu = cinput_common + INPUT_NAME + "2e2mu_" + sample_suffix_MCFM[t] + "_Reprocessed.root";
				string cinput_4e = cinput_common + INPUT_NAME + "4e_" + sample_suffix_MCFM[t] + "_Reprocessed.root";
				string cinput_4mu = cinput_common + INPUT_NAME + "4mu_" + sample_suffix_MCFM[t] + "_Reprocessed.root";

				tree->Add(cinput_2e2mu.c_str());
				tree->Add(cinput_4e.c_str());
				tree->Add(cinput_4mu.c_str());
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
			}
			else{
				D_temp_1D[t] = new TH1F(templatename_1D,templatename_1D,nbinsy,kDYarray);
				D_temp_1D[t]->GetXaxis()->SetTitle(strFitDim_label[tFitD]);

				D_temp_2D[t] = new TH2F(templatename_2D,templatename_2D,nbinsx,kDXarray,nbinsy,kDYarray);
				D_temp_2D[t]->GetXaxis()->SetTitle(strFitDim_label[0]);
				D_temp_2D[t]->GetYaxis()->SetTitle(strFitDim_label[tFitD]);
			};

			double nTotal=0;
			if(t<5){
				tree->SetBranchAddress("GenHMass",&GenHMass);
				tree->SetBranchAddress("ZZMass",&ZZMass);
				if(t!=4) tree->SetBranchAddress("MC_weight",&MC_weight);
				else tree->SetBranchAddress("ZXfake_weightProper",&MC_weight);
				if(t<3){
					tree->SetBranchAddress("MC_weight_Kfactor",&MC_weight_Kfactor);
					tree->SetBranchAddress("MC_weight_down",&MC_weight_down);
					tree->SetBranchAddress("MC_weight_up",&MC_weight_up);
					tree->SetBranchAddress("MC_weight_ggZZLepInt",&MC_weight_ggZZLepInt);
				}
				else{
					MC_weight_ggZZLepInt=1;
					MC_weight_Kfactor=1;
					MC_weight_down=1;
					MC_weight_up=1;
				};

				int branch_set = 0;
				if(tFitD!=0) branch_set = tree->SetBranchAddress(strFitDim[tFitD],&fitYval);
				if(branch_set!=0){
					cerr << "Could NOT find branch named " << strFitDim[tFitD] << "!!! Setting strFitDim[" << tFitD << "] = 0." << endl;
					fitYval=0;
				};

				nEntries = tree->GetEntries();
				for(int ev=0;ev<nEntries;ev++){
					tree->GetEntry(ev);

					double recoveryScale=1;
	//				if(folder==0 && t==1 && (tree->GetTreeNumber() == 1) ) recoveryScale = 0.58185465;
	//				if(folder==1 && t==1 && (tree->GetTreeNumber() == 1) ) recoveryScale = 0.401331435;
					double weight = MC_weight;
					if(Systematics==-1 && t<3) weight = MC_weight_down;
					if(Systematics==1 && t<3) weight = MC_weight_up;
					if(t<3) weight *= MC_weight_Kfactor*recoveryScale*MC_weight_ggZZLepInt;;
					if(t>4) weight *= tsp_VBF_Sig->Eval(GenHMass)*MC_weight_ggZZLepInt;;

					if(t==0 && ZZMass>=ZZMass_PeakCut[0] && ZZMass<ZZMass_PeakCut[1]) nSig_Simulated += weight;
					if(t<3 && ZZMass>=ZZMass_cut[0]) nMZZ220[t] += weight;

					if( (ZZMass>=ZZMass_cut[1] || ZZMass<ZZMass_cut[0]) ) continue;

					nTotal += weight;

					if(tFitD==0) D_temp_1D[t]->Fill(ZZMass,weight);
					else D_temp_1D[t]->Fill(fitYval,weight);
					if(tFitD>0) D_temp_2D[t]->Fill(ZZMass,fitYval,weight);
				};
			};

			cout << t << " NTotal: " << nTotal << endl;
			if(t<3) cout << nMZZ220[t]*nSM_ScaledPeak/nSig_Simulated*luminosity[EnergyIndex] << endl;
			if(t<3){
				D_temp_1D[t]->Scale(nSM_ScaledPeak/nSig_Simulated);
				if(tFitD>0) D_temp_2D[t]->Scale(nSM_ScaledPeak/nSig_Simulated);
				cout << "Scaling t: " << t << " by " << nSM_ScaledPeak/nSig_Simulated << endl;
			};
			if(t>=5){
				double myscale = (nVBF_Sig_Reweighted/nVBF_Sig_Simulated);

				cout << "Expected nVBF_peak: " << nVBF_Sig_Reweighted << endl;
				cout << "Observed nVBF_peak: " << nVBF_Sig_Simulated << endl;
				cout << "Scaling t: " << t << " by " << myscale << endl;
				for(int binx=1;binx<=nbinsx;binx++){
					double bincontent = h1DVBF[t-5]->GetBinContent(binx);
					bincontent *= myscale;
					D_temp_1D[t]->SetBinContent(binx,bincontent);

					if(tFitD>0){
						for(int biny=1;biny<=nbinsy;biny++){
							bincontent = h2DVBF[t-5]->GetBinContent(binx,biny);
							bincontent *= myscale;
							D_temp_2D[t]->SetBinContent(binx,biny,bincontent);
						};
					};
				};
			};

			if(t<5){
				double presmoothInt = D_temp_1D[t]->Integral("width");
				D_temp_1D[t]->Smooth(1,"k3a");
				double postsmoothInt = D_temp_1D[t]->Integral("width");
				D_temp_1D[t]->Scale(presmoothInt/postsmoothInt);
				if(tFitD>0){
					presmoothInt = D_temp_2D[t]->Integral("width");
					D_temp_2D[t]->Smooth(1,"k3a");
					postsmoothInt = D_temp_2D[t]->Integral("width");
					D_temp_2D[t]->Scale(presmoothInt/postsmoothInt);
				};
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
		for(int tr=0;tr<4;tr++){
			delete h1DVBF[tr];
			if(tFitD!=0) delete h2DVBF[tr];
		};

		double nTotalSMSig = D_temp_1D[2]->Integral();
		D_temp_1D[2]->Add(D_temp_1D[1],-1.0);
		D_temp_1D[2]->Add(D_temp_1D[0],-1.0);
		if(tFitD!=0){
			D_temp_2D[2]->Add(D_temp_2D[1],-1.0);
			D_temp_2D[2]->Add(D_temp_2D[0],-1.0);
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
					D_temp_1D[t]->Smooth(1,"k3a");
					double postsmoothInt = D_temp_1D[t]->Integral();
					D_temp_1D[t]->Scale(presmoothInt/postsmoothInt);
				};
			};
			foutput->WriteTObject(D_temp_1D[t]);
			if(tFitD!=0) foutput->WriteTObject(D_temp_2D[t]);

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
		foutput->Close();
	};
	finput_VBF->Close();
	delete tg_interf;
};
