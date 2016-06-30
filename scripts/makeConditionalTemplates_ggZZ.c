#include <iostream>
#include <ctime>
#include <fstream>
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
#include "../data/ZZ4l_Samples.h"
#include <ZZMatrixElement/MELA/interface/Mela.h>

using namespace std;

// Constants to affect the template code
const TString fixedDate="290616";

// Initializers
void makeConditionalTemplates_ggZZ_one(int channel, int erg_tev, int Systematics);
void progressbar(int val, int tot);
TString todaysdate();
SimpleParticleCollection_t constructFourVectors(short GenLepId[4], float GenLepPtEtaPhi[3][4]);
void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray);

// Main Function, runs over all desired iterations
void makeConditionalTemplates_ggZZ(){
  for (int CoM=13; CoM<=13; CoM++){ for (int channel=0; channel<3; channel++){ for (int syst=0; syst<=0; syst++) makeConditionalTemplates_ggZZ_one(channel, CoM, syst); } }
}

// Function to build one template
// folder = 0,1,2 (final state corresponds to 4mu, 4e, 2mu2e respectively)
// erg_tev = 13 (CoM energy)
void makeConditionalTemplates_ggZZ_one(int folder, int erg_tev, int Systematics){
  float mPOLE=125.;
  float wPOLE=4.07e-3;
  int EnergyIndex=(erg_tev==13 ? 0 : 0);

  const float ZZMass_cut[2]={ 70., 3010. };
  const int nbinsy=30;
  float kDYarray[nbinsy+1];
  const float kDY_bounds[2]={ 0, 1 };
  for (int bin=0; bin<nbinsy+1; bin++){
    float binwidth = (kDY_bounds[1] - kDY_bounds[0])/nbinsy;
    kDYarray[bin] = kDY_bounds[0] + binwidth*bin;
  }


  TString erg_dir = Form("LHC_%iTeV/", erg_tev);
  TString strdate = todaysdate();
  if (fixedDate!="") strdate=fixedDate;
  cout << "Today's date: " << strdate << endl;

  TString INPUT_NAME = user_treefile;
  TString OUTPUT_NAME = Form("HtoZZ4l_ggTo%s_ConditionalTemplatesForCombine_", user_folder[folder].Data());
  if (Systematics == 0) OUTPUT_NAME += "Nominal";
  else if (Systematics == 1) OUTPUT_NAME += "SysUp_ggQCD";
  else if (Systematics == -1) OUTPUT_NAME += "SysDown_ggQCD";
  else if (Systematics == 2) OUTPUT_NAME += "SysUp_ggPDF";
  else if (Systematics == -2) OUTPUT_NAME += "SysDown_ggPDF";
  else{ cerr << "Invalid systematics " << Systematics << endl; return; }

  TString OUTPUT_LOG_NAME = OUTPUT_NAME;
  OUTPUT_NAME += ".root";
  OUTPUT_LOG_NAME += ".log";

  TString cinput_common = user_gg2VV_location + "/";
  TString coutput_common = user_dir + erg_dir + "Templates/" + strdate + "/gg/";
  gSystem->Exec("mkdir -p " + coutput_common);

  float varKD=0;
  float weight=1;
  float weightHypo=1;
  float MC_weight=1;
  float MC_weight_norm=1;
  float MC_weight_xsec=1;
  float MC_weight_Kfactor=1;
  short Z1Flav, Z2Flav;
  int ZZFlav;
  float ZZMass, p0plus_VAJHU, bkg_VAMCFM;
  float GenHMass;
  short GenLepId[4];
  float GenLepPtEtaPhi[3][4];

  TVar::VerbosityLevel verbosity = TVar::SILENT;
  Mela* mela = new Mela(erg_tev, mPOLE, verbosity);

  TString coutput = coutput_common + OUTPUT_NAME;
  TString coutput_log = coutput_common + OUTPUT_LOG_NAME;
  ofstream tout(coutput_log.Data(), ios::out);
  cout << "Opened file " << coutput_log << endl;
  TFile* foutput = new TFile(coutput, "recreate");

  cout << endl;
  cout << "===============================" << endl;
  cout << "CoM Energy: " << erg_tev << " TeV" << endl;
  cout << "Decay Channel: " << user_folder[folder] << endl;
  cout << "===============================" << endl;
  cout << endl;
  tout << "===============================" << endl;
  tout << "CoM Energy: " << erg_tev << " TeV" << endl;
  tout << "Decay Channel: " << user_folder[folder] << endl;
  tout << "===============================" << endl;
  tout << endl;

  const unsigned int kNumTemplates=4;
  std::vector<TH2F*> D_temp_2D[kNumTemplates];
  TTree* newtree[kNumTemplates]={ 0 };
  TString templatenames[kNumTemplates]={
    "ggF Sig",
    "gg Bkg",
    "ggF Int",
    "ggF Int Perp"
  };

  //Template filler
  for (unsigned int t=0; t<kNumTemplates; t++){
    //Grab appropriate files for templates
    cout << "Files: " << endl;
    const unsigned int nTrees=6;
    std::vector<std::pair<TTree*, TH1F*>> tree;

    //Initialize templates
    TString templatename_2D;
    // Build template structure, including anomalous couplings
    if (t<2) templatename_2D = Form("T_2D_%i", t+1);
    else if (t==2) templatename_2D = Form("T_2D_12%i", 4);
    else if (t==3) templatename_2D = Form("T_2D_12%i_perp", 4);
    newtree[t] = new TTree(Form("%s_Tree", templatename_2D.Data()), "");
    newtree[t]->Branch("ZZMass", &ZZMass);
    newtree[t]->Branch("GenHMass", &GenHMass);
    newtree[t]->Branch("KD", &varKD);
    newtree[t]->Branch("weight", &weight);
    newtree[t]->Branch("weightHypo", &weightHypo);

    const unsigned int ngenchannels=6;
    TString genchannel[ngenchannels]={
      "4mu",
      "4e",
      "2e2mu",
      "2e2tau",
      "2mu2tau",
      "4tau"
    };
    TFile* fin[ngenchannels]={ 0 };
    std::vector<std::pair<unsigned int, float>> massArray;
    for (unsigned int igch=0; igch<ngenchannels; igch++){
      TString cinput = cinput_common;
      if (t<3) cinput = cinput + sample_suffix_MCFM[t] + genchannel[igch] + "/" + INPUT_NAME;
      else cinput = cinput + sample_suffix_MCFM[2] + genchannel[igch] + "/" + INPUT_NAME;
      cout << "Opening file " << cinput << endl;
      tout << "Opening file " << cinput << endl;
      fin[igch] = TFile::Open(cinput, "read");
      if (fin[igch]!=0 && !fin[igch]->IsZombie() && fin[igch]->IsOpen()){
        // Prepare trees
        cout << "Set variables in trees for " << templatenames[t] << " / " << genchannel[igch] << endl;
        tout << "Set variables in trees for " << templatenames[t] << " / " << genchannel[igch] << endl;
        foutput->cd();
        TTree* theTree = (TTree*)fin[igch]->Get(user_treename);
        TH1F* theCounters = (TH1F*)fin[igch]->Get(user_countersname);
        if (theTree->GetBranchStatus("GenHMass")){
          theTree->SetBranchAddress("GenHMass", &GenHMass);
          for (unsigned int ilep=0; ilep<4; ilep++){
            theTree->SetBranchAddress(Form("GenLep%iId", ilep+1), &(GenLepId[ilep]));
            theTree->SetBranchAddress(Form("GenLep%iPt", ilep+1), &(GenLepPtEtaPhi[0][ilep]));
            theTree->SetBranchAddress(Form("GenLep%iEta", ilep+1), &(GenLepPtEtaPhi[1][ilep]));
            theTree->SetBranchAddress(Form("GenLep%iPhi", ilep+1), &(GenLepPtEtaPhi[2][ilep]));
          }
        }
        theTree->SetBranchAddress("ZZMass", &ZZMass);
        theTree->SetBranchAddress("Z1Flav", &Z1Flav);
        theTree->SetBranchAddress("Z2Flav", &Z2Flav);
        theTree->SetBranchAddress("overallEventWeight", &MC_weight_norm);
        theTree->SetBranchAddress("xsec", &MC_weight_xsec);
        theTree->SetBranchAddress("KFactorggZZ", &MC_weight_Kfactor);
        theTree->SetBranchAddress("p0plus_VAJHU", &p0plus_VAJHU);
        theTree->SetBranchAddress("bkg_VAMCFM", &bkg_VAMCFM);
        tree.push_back(std::pair<TTree*, TH1F*>(theTree, theCounters));
      }
      else cerr << "Could not open file!" << endl;
    }

    TVar::Process sampleProcess=TVar::bkgZZ;
    TVar::Process targetProcess=TVar::bkgZZ;
    TVar::Production sampleProduction = TVar::ZZGG;
    TVar::MatrixElement sampleME=TVar::MCFM;
    if ((t<3 && sample_suffix_MCFM[t]==sample_suffix_MCFM[kBkg]) || (t==3 && sample_suffix_MCFM[t-1]==sample_suffix_MCFM[kBkg])) sampleProcess=TVar::bkgZZ;
    else if ((t<3 && sample_suffix_MCFM[t]==sample_suffix_MCFM[kSig]) || (t==3 && sample_suffix_MCFM[t-1]==sample_suffix_MCFM[kSig])) sampleProcess=TVar::HSMHiggs;
    else if ((t<3 && sample_suffix_MCFM[t]==sample_suffix_MCFM[kBSI]) || (t==3 && sample_suffix_MCFM[t-1]==sample_suffix_MCFM[kBSI])) sampleProcess=TVar::bkgZZ_SMHiggs;
    if (t==0) targetProcess=TVar::HSMHiggs;
    else if (t==1) targetProcess=TVar::bkgZZ;
    else targetProcess=TVar::bkgZZ_SMHiggs;

    unsigned int ctr=0;
    float nTotal=0;
    for (unsigned int itr=0; itr<tree.size(); itr++){
      cout << "First pass over tree " << itr << endl;
      tout << "First pass over tree " << itr << endl;
      float wSumGenWeights=1./tree.at(itr).second->GetBinContent(40); // HARDCODED NUMBER!!!! 
      const unsigned int nEntries = tree.at(itr).first->GetEntries();
      cout << "Number of entries before filtering: " << nEntries << endl;
      tout << "Number of entries before filtering: " << nEntries << endl;
      for (unsigned int ev=0; ev<nEntries; ev++){
        tree.at(itr).first->GetEntry(ev);
        progressbar(ev, tree.at(itr).first->GetEntries());

        ZZFlav=Z1Flav*Z2Flav;
        if (folder==0 && abs(ZZFlav)!=pow(13, 4)) continue;
        else if (folder==1 && abs(ZZFlav)!=pow(11, 4)) continue;
        else if (folder==2 && abs(ZZFlav)!=pow(13*11, 2)) continue;

        varKD = p0plus_VAJHU/(p0plus_VAJHU + bkg_VAMCFM/(1./0.45-1.));
        if (varKD!=varKD) continue;

        MC_weight = MC_weight_norm*MC_weight_xsec*MC_weight_Kfactor*wSumGenWeights;
        // Protect against any KD exceeding boundaries
        if (varKD>=kDYarray[nbinsy]){
          cout << "Found varKD == " << varKD;
          varKD = kDYarray[nbinsy] - (kDYarray[nbinsy]-kDYarray[nbinsy-1])*0.1*(ev+1.)/(nEntries+1.);
          cout << ". Corrected to have " << varKD << endl;
        }
        if (varKD<kDYarray[0]) varKD = kDYarray[0] + (kDYarray[1]-kDYarray[0])*0.1*ev/nEntries;
        if (varKD>=kDYarray[nbinsy] || varKD<kDYarray[0]) cout << "Fix has been numerically unsuccessful for " << tree.at(itr).first->GetName() << endl;

        weight = MC_weight;

        SimpleParticleCollection_t daughters = constructFourVectors(GenLepId, GenLepPtEtaPhi);
        mela->setInputEvent(&daughters, (SimpleParticleCollection_t*)0, (SimpleParticleCollection_t*)0, false);
        GenHMass = mela->getCurrentCandidate()->m(); // Since we add artificial masses, better to recompute m.

        float sampleWeight=1;
        mela->setProcess(sampleProcess, sampleME, sampleProduction);
        mela->setMelaHiggsMassWidth(mPOLE, wPOLE, 0);
        if (sampleProcess!=TVar::bkgZZ) mela->selfDHzzcoupl[0][0][0] = sqrt(MCFM_widthrescale);
        mela->computeP(sampleWeight, false);

        float targetWeight=1;
        mela->setProcess(targetProcess, sampleME, sampleProduction);
        mela->setMelaHiggsMassWidth(mPOLE, wPOLE, 0);
        if (t!=1){
          double g1Re = pow(GenHMass, 2)-pow(mPOLE, 2);
          double g1Im = mPOLE*wPOLE;
          if (t!=3){
            mela->selfDHzzcoupl[0][0][0] = g1Re/pow(GenHMass, 2);
            mela->selfDHzzcoupl[0][0][1] = g1Im/pow(GenHMass, 2);
          }
          else{
            mela->selfDHzzcoupl[0][0][0] = -g1Im/pow(GenHMass, 2);
            mela->selfDHzzcoupl[0][0][1] = g1Re/pow(GenHMass, 2);
          }
        }
        mela->computeP(targetWeight, false);

        mela->resetInputEvent();

        weightHypo = targetWeight/sampleWeight;
        if (std::isnan(weightHypo) || std::isinf(weightHypo)) continue;
        weight *= weightHypo;
        nTotal += weight;

        addByLowest(ctr, GenHMass, massArray);
        newtree[t]->Fill();
        ctr++;
      }
      cout << endl;
      tout << endl;
    }

    // Close the input file
    for (unsigned int igch=0; igch<ngenchannels; igch++){ if (fin[igch]!=0 && fin[igch]->IsOpen()) fin[igch]->Close(); }

    cout << templatenames[t] << " Total Simulated: " << nTotal*luminosity[EnergyIndex] << '\n' << endl;
    tout << templatenames[t] << " Total Simulated: " << nTotal*luminosity[EnergyIndex] << '\n' << endl;

    const unsigned int nFilteredEntries = newtree[t]->GetEntries();
    newtree[t]->GetEntry(massArray.at(0).first);
    float firstXVal=GenHMass;
    newtree[t]->GetEntry(massArray.at(nFilteredEntries-1).first);
    float lastXVal=GenHMass;
    float infimum = (float)((int)firstXVal); infimum -= (float)(((int)infimum)%10);
    float supremum = (float)((int)(lastXVal+0.5)); supremum += (float)(10-((int)supremum)%10);
    cout << "Nentries = " << nFilteredEntries << " | mzz = " << firstXVal << " - " << lastXVal << "(" << infimum << ", " << supremum << ")" << endl;

    int nbinsx=0;
    int divisor=6000;
    const int nbins_th=50;
    while (nbinsx<nbins_th){
      if (divisor>1000) divisor -= 1000;
      else if (divisor>100) divisor -= 100;
      else break;
      nbinsx=nFilteredEntries/divisor+1;
    }
    cout << "nbinsx=" << nbinsx << endl;
    if (nbinsx<3) cerr << "Not enough bins!" << endl;
    float* kDXarray = new float[nbinsx+1];
    kDXarray[0]=infimum;
    kDXarray[nbinsx]=supremum;
    int ev_stepsize = nFilteredEntries/nbinsx;
    cout << "Event step size: " << ev_stepsize << endl;
    cout << "Boundary (" << 0 << ") = " << kDXarray[0] << endl;
    for (int ix=1; ix<nbinsx; ix++){
      int ev = massArray.at(ix*ev_stepsize).first;
      newtree[t]->GetEntry(ev);
      float bhigh = GenHMass;
      ev = massArray.at(ix*ev_stepsize-1).first;
      float blow = GenHMass;
      kDXarray[ix]=(bhigh+blow)*0.5;
      cout << "Boundary (" << ix << ")= " << kDXarray[ix] << " [event " << ev << ", step " << ix*ev_stepsize << "]" << endl;
    }
    cout << "Boundary (" << nbinsx << ") = " << kDXarray[nbinsx] << endl;

    for (unsigned int al = 0; al<1; al++){
      D_temp_2D[t].push_back(new TH2F(templatename_2D, templatename_2D, nbinsx, kDXarray, nbinsy, kDYarray));
      D_temp_2D[t].at(al)->GetXaxis()->SetTitle("D^{kin}_{bkg}");
      D_temp_2D[t].at(al)->GetYaxis()->SetTitle("Events");
    }

    cout << "Number of entries after filtering: " << nFilteredEntries << endl;
    tout << "Number of entries after filtering: " << nFilteredEntries << endl;
    for (unsigned int ev=0; ev<nFilteredEntries; ev++){
      newtree[t]->GetEntry(ev);
      progressbar(ev, newtree[t]->GetEntries());
      // Anomalous couplings loop
      for (unsigned int al=0; al<D_temp_2D[t].size(); al++){
        double fillWeight = weight;
        D_temp_2D[t].at(al)->Fill(GenHMass, varKD, fillWeight);
      }
    }

    delete[] kDXarray;
  }


  for (unsigned int t=0; t<kNumTemplates; t++){
    for (unsigned int al=0; al<D_temp_2D[t].size(); al++){
      foutput->WriteTObject(D_temp_2D[t].at(al));
      delete D_temp_2D[t].at(al);
    }
    foutput->WriteTObject(newtree[t]);
    delete newtree[t];
  }
  foutput->Close();
  tout.close();
  delete mela;
}

void progressbar(int val, int tot){
  int percent=floor(0.01*tot);
  if (percent==0) percent=1;

  if (val%percent==0 && val!=tot){
    cout<<"[ "<<setw(3)<<val/percent<<"% |";
    for (int k=1; k<val/percent; k++) cout<<"=";
    if (val%percent!=100) cout<<">";
    for (int k=val/percent+1; k<100; k++) cout<<" ";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
  else if (val==tot){
    cout<<"[ 100% |";
    for (int k=0; k<100; k++) cout<<"=";
    cout<<"| ]";
    fflush(stdout);
    putchar('\r');
  }
}

TString todaysdate(){
  TString result;
  time_t now = time(0);
  tm* ltm = localtime(&now);
  int year = 1900 + ltm->tm_year-2000;
  int month = 1 + ltm->tm_mon;
  int day = ltm->tm_mday;
  if(month>=10) result = Form("%i%i%i", day, month, year);
  else result = Form("%i%i%i%i", day, 0, month, year);
  return result;
}


void addByLowest(unsigned int ev, float val, std::vector<std::pair<unsigned int, float>>& valArray){
  bool inserted = false;
  for (std::vector<std::pair<unsigned int, float>>::iterator it = valArray.begin(); it<valArray.end(); it++){
    if ((*it).second>=val){
      inserted=true;
      valArray.insert(it, std::pair<unsigned int, float>(ev, val));
      break;
    }
  }
  if (!inserted) valArray.push_back(std::pair<unsigned int, float>(ev, val));
}


SimpleParticleCollection_t constructFourVectors(short GenLepId[4], float GenLepPtEtaPhi[3][4]){
  SimpleParticleCollection_t result;
  for (int idau=0; idau<4; idau++){
    float mass=0;
    if (abs(GenLepId[idau])==11) mass=0.000511;
    else if (abs(GenLepId[idau])==13) mass=0.10566;
    else if (abs(GenLepId[idau])==15) mass=1.7768;
    else if (abs(GenLepId[idau])==5) mass=4.75;
    else if (abs(GenLepId[idau])==4) mass=1.275;
    TLorentzVector mom; mom.SetPtEtaPhiM((double)GenLepPtEtaPhi[0][idau], (double)GenLepPtEtaPhi[1][idau], (double)GenLepPtEtaPhi[2][idau], (double)mass);
    result.push_back(SimpleParticle_t((int)GenLepId[idau], mom));
  }
  return result;
}



