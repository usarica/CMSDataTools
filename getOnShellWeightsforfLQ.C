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
#include "TH1D.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "./data/ZZ4l_125_6_Samples.h"
#include "./data/FitDimensionsList.h"

using namespace std;

void getOnShellWeightsforfLQ_one(int erg_tev);

void getOnShellWeightsforfLQ(){
	getOnShellWeightsforfLQ_one(7);
	getOnShellWeightsforfLQ_one(8);
}

void getOnShellWeightsforfLQ_one(int erg_tev){
	TString inputdir = user_gg2VV_location;
	TString erg_dir;
	erg_dir.Form("LHC_%iTeV/",erg_tev);
	TString erg_name;
	erg_name.Form("%i",erg_tev);
	float Hmass=125.6;

	//GEN LEVEL
	/*TString ggH_125p6 = "HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Generated.root";
	TString VBF_125p6 = "VBFHiggs0PToZZTo4L_M-125p6_" + erg_name + "TeV-JHUGenV4_false.root";

	TChain* ggH = new TChain("GenTree");
	ggH->Add(inputdir + erg_dir + "GenSignal/" + ggH_125p6);
	TChain* VBF = new TChain("SelectedTree");
	VBF->Add(inputdir + erg_dir + "GenSignal/" + VBF_125p6);*/

	//FULL-SIM LEVEL
	TString ggH_125p6 = "HZZ4lTree_powheg15jhuGenV3-0PMH125.6_Reprocessed.root";
	TString VBF_125p6 = "HZZ4lTree_VBF0P_H125.6.root";
	TString VBFinputdir = "/scratch0/hep/ianderso/CJLST/140604/";
	TString VBFerg_dir;
	if(erg_tev==7) VBFerg_dir="PRODFSR/";
	if(erg_tev==8) VBFerg_dir="PRODFSR_8TeV/";

	TChain* ggH = new TChain("SelectedTree");
	ggH->Add(inputdir + erg_dir + "2mu2e/" + ggH_125p6);
	ggH->Add(inputdir + erg_dir + "4e/" + ggH_125p6);
	ggH->Add(inputdir + erg_dir + "4mu/" + ggH_125p6);
	TChain* VBF = new TChain("SelectedTree");
	VBF->Add(VBFinputdir + VBFerg_dir + "2mu2e/" + VBF_125p6);
	VBF->Add(VBFinputdir + VBFerg_dir + "4e/" + VBF_125p6);
	VBF->Add(VBFinputdir + VBFerg_dir + "4mu/" + VBF_125p6);

	float mZZ, mZZ_RECO;
	ggH->SetBranchAddress("GenHMass",&mZZ);
	ggH->SetBranchAddress("ZZMass",&mZZ_RECO);
	VBF->SetBranchAddress("GenHMass",&mZZ);
	VBF->SetBranchAddress("ZZMass",&mZZ_RECO);

	double ggH_0=1.;
	double ggH_1=0.;
	double ggH_2=0.;
	double VBF_0=1.;
	double VBF_1=0.;
	double VBF_2=0.;
	double VBF_3=0.;
	double VBF_4=0.;

	for(int iEvt=0;iEvt<ggH->GetEntries();++iEvt){
		ggH->GetEntry(iEvt);
		if(mZZ_RECO<105.6 || mZZ_RECO>141.6) continue;
		ggH_1+=(double)(pow(mZZ,2.)/pow(Hmass,2.));
		ggH_2+=(double)pow(pow(mZZ,2.)/pow(Hmass,2.),2.);
	}
	ggH_1/=ggH->GetEntries();
	ggH_2/=ggH->GetEntries();

	for(int iEvt=0;iEvt<VBF->GetEntries();++iEvt){
		VBF->GetEntry(iEvt);
		if(mZZ_RECO<105.6 || mZZ_RECO>141.6) continue;
		VBF_1+=(double)(pow(mZZ,2.)/pow(Hmass,2.));
		VBF_2+=(double)pow(pow(mZZ,2.)/pow(Hmass,2.),2.);
		VBF_3+=(double)pow(pow(mZZ,2.)/pow(Hmass,2.),3.);
		VBF_4+=(double)pow(pow(mZZ,2.)/pow(Hmass,2.),4.);
	}
	VBF_1/=VBF->GetEntries();
	VBF_2/=VBF->GetEntries();
	VBF_3/=VBF->GetEntries();
	VBF_4/=VBF->GetEntries();

	cout<<erg_tev<<endl;
	cout<<"ggH"<<endl;
	cout<<setprecision(14)<<ggH_0<<endl;
	cout<<setprecision(14)<<ggH_1<<endl;
	cout<<setprecision(14)<<ggH_2<<endl;
	cout<<"VBF"<<endl;
	cout<<setprecision(14)<<VBF_0<<endl;
	cout<<setprecision(14)<<VBF_1<<endl;
	cout<<setprecision(14)<<VBF_2<<endl;
	cout<<setprecision(14)<<VBF_3<<endl;
	cout<<setprecision(14)<<VBF_4<<endl;

	TH1D* ggH_ratios = new TH1D("ggH_ratios","ggH_ratios",3,0.,3.);
	TH1D* VBF_ratios = new TH1D("VBF_ratios","VBF_ratios",5,0.,5.);
	ggH_ratios->Fill(0.,ggH_0);
	ggH_ratios->Fill(1.,ggH_1);
	ggH_ratios->Fill(2.,ggH_2);
	VBF_ratios->Fill(0.,VBF_0);
	VBF_ratios->Fill(1.,VBF_1);
	VBF_ratios->Fill(2.,VBF_2);
	VBF_ratios->Fill(3.,VBF_3);
	VBF_ratios->Fill(4.,VBF_4);

	TString outputname = "./data/m4l_ratios_" + erg_name + "TeV_RECO.root";  
	TFile* output = new TFile(outputname,"recreate");
	output->WriteTObject(ggH_ratios);
	output->WriteTObject(VBF_ratios);
	output->Close();
}