#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TROOT.h"
#include "TF1.h"
#include <string>
#include <sstream>
#include <iostream>
#include <vector>

using namespace std;

TString writeloc = "Plots/Ratios/";
TString templatenames[10]={"T_2D_1","T_2D_2","T_2D_4","T_2D_qqZZ_UnConditional","T_2D_qqZZ","T_2D_ZX_UnConditional","T_2D_ZX","T_2D_VBF_1","T_2D_VBF_2","T_2D_VBF_4"};
TString filenames[5]={
	"HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysUp_ggQCD.root",
	"HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysUp_ggPDF.root",
	"HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysDown_ggQCD.root",
	"HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysDown_ggPDF.root",
	"HtoZZ4l_MCFM_125p6_ModifiedSmoothTemplatesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_Nominal.root"
};
TString treefilenames[5]={
	"HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_Nominal.root",
	"HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysDown_ggPDF.root",
	"HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysDown_ggQCD.root",
	"HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysUp_ggPDF.root",
	"HtoZZ4l_MCFM_125p6_ModifiedTemplateswithTreesForCombine__GenLevelVBF_wResolution_D_Gamma_gg_r10_SysUp_ggQCD.root"
};

void compareonesmoothedtemplate(int CoM, int channel, int numfile, int numtemplate);
void compareonenormalization(int CoM, int channel, int numfile, int numtemplate);


void comparesmoothedtemplates(){
	for(int CoM=7;CoM<9;++CoM){
		for(int channel=0;channel<3;++channel){
			for(int file=0;file<5;++file){
				for(int temp=0;temp<10;++temp){
					compareonesmoothedtemplate(CoM,channel,file,temp);
				}
			}
		}
	}
}

void compareonesmoothedtemplate(int CoM, int channel, int numfile, int numtemplate){
	TString newfileloc="/home/ianderso/Work/HighMass/AnomalousCouplings/HiggsWidth_OLD/UsingLastProduction/";
	TString oldfileloc="/home/ianderso/Work/HighMass/AnomalousCouplings/HiggsWidth_PostICHEP/FromUlascan_020414/";
	if(CoM==7){
		newfileloc+="LHC_7TeV/";
		oldfileloc+="7TeV/";
	} else if(CoM==8){
		newfileloc+="LHC_8TeV/";
		oldfileloc+="8TeV/";
	}
	if(channel==0){ newfileloc+="2mu2e/250314/"; oldfileloc+="2mu2e/";}
	if(channel==1){ newfileloc+="4e/250314/"; oldfileloc+="4e/";}
	if(channel==2){ newfileloc+="4mu/250314/"; oldfileloc+="4mu/";}

	TFile* newfile = new TFile(newfileloc + filenames[numfile],"read");
	TFile* oldfile = new TFile(oldfileloc + filenames[numfile],"read");

	TH2F* newtemp = (TH2F*) newfile->Get(templatenames[numtemplate]);
	TH2F* oldtemp = (TH2F*) oldfile->Get(templatenames[numtemplate]);

	TH2F* ratio = (TH2F*) newtemp->Clone("ratio");
	ratio->Divide(oldtemp);

	TCanvas* c = new TCanvas("c","c",800,800);
	c->cd();
	ratio->Draw("COLZ");

	TString filename="Ratio_";
	if(CoM==7) filename+="7TeV_";
	if(CoM==8) filename+="8TeV_";
	if(channel==0) filename+="2e2mu_";
	if(channel==1) filename+="4e_";
	if(channel==2) filename+="4mu_";
	if(numfile==0) filename+="SysUp_QCD_";
	if(numfile==1) filename+="SysUp_PDF_";
	if(numfile==2) filename+="SysDn_QCD_";
	if(numfile==3) filename+="SysDn_PDF_";
	if(numfile==4) filename+="Nominal_";
	TString tempstring;
	tempstring.Form("%i",numtemplate);
	filename+=tempstring;
	
	cout<<filename<<" normalization: "<<newtemp->Integral("width")<<"/"<<oldtemp->Integral("width")<<"="<<newtemp->Integral("width")/oldtemp->Integral("width")<<endl;

	c->SaveAs(writeloc + filename + ".eps");
	c->SaveAs(writeloc + filename + ".pdf");
	c->SaveAs(writeloc + filename + ".png");
	c->SaveAs(writeloc + filename + ".root");
	c->SaveAs(writeloc + filename + ".C");
}

void comparenormalizations(){
	for(int CoM=7;CoM<9;++CoM){
		for(int channel=0;channel<2;++channel){
			for(int file=0;file<5;++file){
				for(int temp=0;temp<10;++temp){
					if(temp==3 || temp==5) continue;
					compareonenormalization(CoM,channel,file,temp);
				}
			}
		}
	}
}

void compareonenormalization(int CoM, int channel, int numfile, int numtemplate){
	TString newfileloc="/home/ianderso/Work/HighMass/AnomalousCouplings/HiggsWidth_PostICHEP/UsingLastProduction/";
	TString oldfileloc="/home/ianderso/Work/HighMass/AnomalousCouplings/HiggsWidth_PostICHEP/FromUlascan_020414/";
	if(CoM==7){
		newfileloc+="LHC_7TeV/";
		oldfileloc+="7TeV/";
	} else if(CoM==8){
		newfileloc+="LHC_8TeV/";
		oldfileloc+="8TeV/";
	}
	if(channel==0){ newfileloc+="2mu2e/"; oldfileloc+="2mu2e/";}
	if(channel==1){ newfileloc+="4e/"; oldfileloc+="4e/";}
	if(channel==2){ newfileloc+="4mu/"; oldfileloc+="4mu/";}

	TString onetreename = templatenames[numtemplate] + "_Tree";
	TChain* newfile = new TChain(onetreename);
	TChain* oldfile = new TChain(onetreename);
	newfile->Add(newfileloc+treefilenames[numfile]);
	oldfile->Add(oldfileloc+treefilenames[numfile]);

	float wt;
	float totold=0.;
	float totnew=0.;
	newfile->SetBranchAddress("templateWeight",&wt);
	oldfile->SetBranchAddress("templateWeight",&wt);

	for(int i=0;i<oldfile->GetEntries();++i) {oldfile->GetEntry(i); totold+=wt;}
	for(int i=0;i<newfile->GetEntries();++i) {newfile->GetEntry(i); totnew+=wt;}

	cout<<totnew/totold<<endl;
}
