#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <unordered_map>
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TF1.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TRandom.h"
#include "interface/CalcHelpers.h"
#include "interface/Samples.h"
#include "interface/SampleHelpers.h"


#ifndef doDebugKD
#define doDebugKD false
#endif
#ifndef doDebugKDExt
#define doDebugKDExt false
#endif

using namespace std;
using namespace CalcHelpers;
using namespace SampleHelpers;


//////////////////
// KD functions //
//////////////////
float constructSimpleKD(vector<float>& vars, float constant=1){
  assert(!checkNonZero(vars) || vars.size()!=2);
  if (!checkNanInf(vars)) return -999;
  float res = vars[0]/(vars[0]+constant*vars[1]);
  return res;
}
float constructPA2PB1KD(vector<float>& vars, float constant=1){
  assert(!checkNonZero(vars) || vars.size()!=3);
  if (!checkNanInf(vars)) return -999;
  return vars[0]*vars[1]/(vars[0]*vars[1]+constant*vars[2]);
}
float constructPA1PB1PBp1KD(vector<float>& vars, float constant=1){
  assert(!checkNonZero(vars) || vars.size()!=3);
  if (!checkNanInf(vars)) return -999;
  return vars[0]/(vars[0]+constant*(vars[1]+vars[2]));
}
float constructPA2PB2KD(vector<float>& vars, float constant=1){
  assert(!checkNonZero(vars) || vars.size()!=4);
  if (!checkNanInf(vars)) return -999;
  return vars[0]*vars[1]/(vars[0]*vars[1]+constant*vars[2]*vars[3]);
}

float constructDjjVH(vector<float>& vars, float constant=1){
  return constructPA2PB1KD(vars, constant);
}
float constructDjVBF(vector<float>& vars, float constant=1){
  return constructPA2PB1KD(vars, constant);
}
float constructDjjVBF(vector<float>& vars, float constant=1){
  return constructSimpleKD(vars, constant);
}
float constructDbkgkin(vector<float>& vars, float constant=1){
  return constructSimpleKD(vars, constant);
}
float constructDbkgdec(vector<float>& vars, float constant=1){
  return constructPA1PB1PBp1KD(vars, constant);
}
float constructDjjbkgEWQCD(vector<float>& vars, float constant=1){
  // vars[0]: JJVBF sig
  // vars[1]: JJZH sig
  // vars[2]: JJWH sig

  // vars[3]: JJVBF bkg
  // vars[4]: JJZH bkg
  // vars[5]: JJWH bkg
  // vars[6]: JJQCD bkg

  // vars[7]: JJZH P(mjj)
  // vars[8]: JJWH P(mjj)

  // vars[9]: JJVBF sig const
  // vars[10]: JJZH sig const
  // vars[11]: JJWH sig const

  // vars[12]: JJVBF bkg const
  // vars[13]: JJZH bkg const
  // vars[14]: JJWH bkg const
  // vars[15]: JJQCD bkg const
  assert(!checkNonZero(vars) || vars.size()!=17);
  if (!checkNanInf(vars)) return -999;

  float const_ew_sig = 1./(1./vars[9]+1./vars[10]+1./vars[11]);
  float const_ew_bkg = 1./(1./vars[12]+1./vars[13]+1./vars[14]);
  float const_ewqcd_bkg = 1./(1./const_ew_bkg+1./vars[15]);

  //float ewsig = vars[0]/const_ew_sig;
  float vbf = vars[0]/vars[9];
  float zh = vars[1]/vars[10];
  float wh = vars[2]/vars[11];
  //float vbfzh_vbfwh = (ewsig - vbf - zh - wh);

  //float ewbkg = vars[4]/const_ew_bkg;
  float vbs = vars[3]/vars[12];
  float zzz = vars[4]/vars[13];
  float wzz = vars[5]/vars[14];
  float qcdzz = vars[6]/vars[15];
  //float vbszzz_vbswzz = (ewbkg - vbs - zzz - wzz);

  float sjj = pow(vars[16], 2);
  float scale_Pmjj_vb = 1;
  float scale_Pmjj_z = vars[7]*(pow(sjj-pow(91.1876, 2), 2)+pow(91.1876*2.4952, 2));
  float scale_Pmjj_w = vars[8]*(pow(sjj-pow(80.399, 2), 2)+pow(80.399*2.085, 2));
  //float scale_Pmjj_int = sqrt(scale_Pmjj_z + scale_Pmjj_w)*sqrt(scale_Pmjj_vb);

  vbf *= scale_Pmjj_vb;
  vbs *= scale_Pmjj_vb;

  zh *= scale_Pmjj_z;
  zzz *= scale_Pmjj_z;

  wh *= scale_Pmjj_w;
  wzz *= scale_Pmjj_w;

  //vbfzh_vbfwh *= scale_Pmjj_int;
  //vbszzz_vbswzz *= scale_Pmjj_int;

  // Ignore VBF-VH or VBS-VZZ interference terms
  float PA = (vbf + zh + wh)*const_ew_sig;
  float PB = (vbs + zzz + wzz + qcdzz)*const_ewqcd_bkg;
  return PA/(PA+constant*PB);
}


///////////////////
// Event helpers //
///////////////////
void getSamplesList(float sqrts, vector<TString> s, vector<TString>& vs){
  for (auto& ss : s){
    vector<TString> dumappend = constructSamplesList(ss, sqrts);
    appendVector<TString>(vs, dumappend);
  }
}
void getSamplePairs(float sqrts, vector<TString> s1, vector<TString> s2, vector<TString>& vs1, vector<TString>& vs2){
  getSamplesList(sqrts, s1, vs1);
  getSamplesList(sqrts, s2, vs2);
}
void bookSampleTrees(
  float sqrts,
  const vector<TString>& strSamples,

  short* Z1Id, 
  short* Z2Id,
  float* DiJetMass,

  TString strTrueBranch,
  float* var_true,
  vector<pair<TString, float*>> str_val_reco,
  vector<pair<TString, float*>> str_weight,

  int& nEntries,
  vector<TFile*>& finputList,
  vector<pair<TTree*, TH1F*>>& treeList
  ){
  TString cinput_main = CJLSTsamplesdir;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";

  for (int is=0; is<(int)strSamples.size(); is++){
    TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[is]).Data());
    TFile* finput = TFile::Open(cinput, "read");
    cout << "Opening file " << cinput << "..." << endl;
    TTree* tree=0;
    if (finput!=0){
      if (finput->IsOpen() && !finput->IsZombie()){
        cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
        tree = (TTree*)finput->Get(TREE_NAME);
        if (tree!=0){
          cout << TREE_NAME << " is found." << endl;
          vector<unsigned int> weightExists;
          for (auto& p : str_weight){ if (tree->GetBranchStatus(p.first)!=0) weightExists.push_back(1); else weightExists.push_back(0); }
          tree->SetBranchStatus("*", 0);

          bookBranch(tree, "Z1Flav", Z1Id);
          bookBranch(tree, "Z2Flav", Z2Id);
          bookBranch(tree, "DiJetMass", DiJetMass);
          bookBranch(tree, strTrueBranch, var_true);

          for (auto& p : str_val_reco){ if (p.first!="DiJetMass") bookBranch(tree, p.first, p.second); }
          for (unsigned int ip=0; ip<str_weight.size(); ip++){
            pair<TString, float*>& p = str_weight.at(ip);
            if (weightExists[ip]==1) bookBranch(tree, p.first, p.second);
            else{ cout << " - Weight " << p.first << " is not applicable" << endl; }
          }
          nEntries += tree->GetEntries();
          treeList.push_back(pair<TTree*, TH1F*>(tree, (TH1F*)finput->Get(COUNTERS_NAME)));
          finputList.push_back(finput);
        }
        else finput->Close();
      }
      else if (finput->IsOpen()) finput->Close();
    }
  }
  cout << "NEntries = " << nEntries << " over " << treeList.size() << " trees." << endl;
}
void getEvents(
  vector<pair<TTree*, TH1F*>>& treeList,
  int& nEntries,
  short& Z1Id, short& Z2Id, float& DiJetMass,
  float& varTrack,
  vector<float>& valReco,
  vector<TString>& nameReco,
  vector<float>& weights,
  vector<SimpleEntry>& index,
  float(*KDfcn)(vector<float>&, float),
  TString strcustomselection
  ){
  vector<short> matchdecid;
  if (strcustomselection.Contains("2l2l") || strcustomselection.Contains("2e2mu")) matchdecid.push_back(121*169);
  else if (strcustomselection.Contains("4l")){ matchdecid.push_back(121*121); matchdecid.push_back(169*169); }
  else if (strcustomselection.Contains("4mu")){ matchdecid.push_back(169*169); }
  else if (strcustomselection.Contains("4e")){ matchdecid.push_back(121*121); }

  pair<float, float> mjjcut(-1, -1);
  if (strcustomselection.Contains("VBFLike")){ mjjcut.first=150; }
  else if (strcustomselection.Contains("VHLike")){ mjjcut.first=40; mjjcut.second=150; }

  pair<float, float> m4lcut(-1, -1);
  if (strcustomselection.Contains("LowMass")){ m4lcut.first=100; m4lcut.second=160; }
  else if (strcustomselection.Contains("HighMass")){ m4lcut.first=160; }

  if (matchdecid.size()>0){
    cout << "Matching events to either of ";
    for (auto& mid : matchdecid) cout << mid << " ";
    cout << endl;
  }
  if (mjjcut.first>=0. || mjjcut.second>=0.) cout << "Matching events to an mjj cut of [ " << mjjcut.first << " , " << mjjcut.second << " ]" << endl;

  unsigned ev_acc=0;
  for (int ev=0; ev<nEntries; ev++){
    for (auto& w : weights) w=1;
    float scale = getEntry(treeList, ev);
    bool recoPos=false;
    unsigned int itv=0;
    for (auto& v:valReco){
      if (nameReco[itv]=="DiJetMass") v = DiJetMass;
      if (std::isnan(v) || std::isinf(v) || v<0.){
        if (std::isnan(v) || std::isinf(v)) cerr << "Invalid valReco[" << itv << "] = " << v << " is being discarded at mass " << varTrack << endl;
        recoPos=false;
        break;
      }
      else if (v>0.) recoPos=true;
      itv++;
    }
    bool doProcess = recoPos;
    
    bool testMatchDec=(matchdecid.size()==0);
    for (short& testid : matchdecid){ if (testid == Z1Id*Z2Id){ testMatchDec=true; break; } }
    doProcess &= testMatchDec;

    bool testMatchRecoM4L = (m4lcut.first<0. || m4lcut.first<=varTrack) && (m4lcut.second<0. || varTrack<=m4lcut.second);
    doProcess &= testMatchRecoM4L;

    bool testMatchRecoMJJ = (mjjcut.first<0. || (mjjcut.first<=DiJetMass && DiJetMass>=0.)) && (mjjcut.second<0. || (DiJetMass<=mjjcut.second && DiJetMass>=0.));
    doProcess &= testMatchRecoMJJ;

    if (!doProcess) continue;

    // Weight has to be positive-definite for TProfile later on
    float wgt = scale;
    for (auto const& w : weights) wgt *= w;
    if (std::isnan(wgt) || std::isinf(wgt) || wgt<=0.){
      // If weight is NaN, it is a big problem.
      if (std::isnan(wgt) || std::isinf(wgt)) cerr << "Invalid weight " << wgt << " is being discarded at mass " << varTrack << endl;
      continue;
    }

    // Need to also test KD explicitly
    float KD = KDfcn(valReco, 1.);
    if (std::isnan(KD) || std::isinf(KD) || KD<0.) continue;

    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    if (doDebugKD && ev_acc==100000) break;

    SimpleEntry theEntry(ev, varTrack, valReco, wgt);
    if (doDebugKDExt & doDebugKD) theEntry.print();
    addByLowest(index, theEntry, false);
    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;
}
void LoopForConstant(
  vector<SimpleEntry>(&index)[2],
  vector<unsigned int>(&indexboundaries)[2],
  float(*KDfcn)(vector<float>&, float),
  TProfile* px,
  TH1F* hrec,
  unsigned int nstepsiter=100
  ){
  int nbins = indexboundaries[0].size()-1;

  for (int bin=0; bin<nbins; bin++){
    cout << "Bin " << bin << " / " << nbins << " is now being scrutinized..." << endl;

    // First find the average KD in each sample
    float sumKD[2]={ 0 }; float sumWgt[2]={ 0 }; float avgKD[2]={ 0 };
    for (unsigned int ih=0; ih<2; ih++){
      unsigned int& evlow = indexboundaries[ih].at(bin);
      unsigned int& evhigh = indexboundaries[ih].at(bin+1);
      cout << " - Scanning events [ " << evlow << " , " << evhigh << " ) for sample set " << ih << endl;
      for (unsigned int ev=evlow; ev<evhigh; ev++){
        float KD = KDfcn(index[ih].at(ev).recoval, 1.);
        if (std::isnan(KD) || std::isinf(KD)){
          cerr << "Something went terribly wrong! KD is " << KD << endl;
          for (auto& v : index[ih].at(ev).recoval) cerr << v << " ";
          cerr << endl;
          continue;
        }
        else if (KD<0.){
          cerr << "KD(ev=" << ev << ") is invalid (" << KD << ")" << endl;
          continue;
        }
        float& varTrack = index[ih].at(ev).trackingval;
        float& weight = index[ih].at(ev).weight;
        if (std::isnan(weight) || std::isinf(weight) || weight<=0.){ cerr << "Invalid weight " << weight << " is being discarded." << endl; continue; }
        sumKD[ih] += KD*weight;
        sumWgt[ih] += weight;

        if (px->GetXaxis()->GetBinLowEdge(bin+1)>varTrack || px->GetXaxis()->GetBinLowEdge(bin+2)<=varTrack) cerr
          << "Something terrible happened! " << varTrack << " is outside ["
          << px->GetXaxis()->GetBinLowEdge(bin+1) << " , " << px->GetXaxis()->GetBinLowEdge(bin+2)
          << "]" << endl;

        px->Fill(varTrack, varTrack, weight);
      }
      cout << " - Sum of weights for sample set " << ih << " = " << sumWgt[ih] << endl;
      avgKD[ih]=sumKD[ih]/sumWgt[ih];
      cout << " - Average KD for sample set " << ih << " = " << avgKD[ih] << endl;
    }

    float marginlow=1;
    float marginhigh=20;
    float Cfound=0;
    float centralConstant = (avgKD[0]+avgKD[1])*0.5; centralConstant = 1./(1./centralConstant-1.);
    unsigned int it=0;
    while (true){
      cout << " - Iteration " << it << " with margins = " << marginlow << ", " << marginhigh << endl;
      cout << "   - Checking c = [ " << centralConstant*(1.-marginlow) << " , " << centralConstant*(1.+marginhigh) << " ]" << endl;
      float mindiff=1;
      float finalFraction[2]={ 2, 2 };
      unsigned int nsteps = nstepsiter;
      if (it==0) nsteps*=1000;
      else if (it==1) nsteps*=100;
      else if (it==2) nsteps*=10;
      for (unsigned int step=0; step<=nsteps; step++){
        float testC = centralConstant*((1.-marginlow) + (marginhigh+marginlow)*(float(step))/((float)nsteps));

        float sumWgtAll[2]={ 0 };
        float sumWgtHalf[2]={ 0 };
        for (unsigned int ih=0; ih<2; ih++){
          unsigned int& evlow = indexboundaries[ih].at(bin);
          unsigned int& evhigh = indexboundaries[ih].at(bin+1);
          for (unsigned int ev=evlow; ev<evhigh; ev++){
            float KD = KDfcn(index[ih].at(ev).recoval, testC);
            if (std::isnan(KD) || std::isinf(KD)){
              cerr << "Something went terribly wrong! KD is " << KD << endl;
              for (auto& v : index[ih].at(ev).recoval) cerr << v << " ";
              cerr << endl;
              continue;
            }
            else if (KD<0.){
              cerr << "KD(ev=" << ev << ") is invalid (" << KD << ")" << endl;
              continue;
            }
            float& weight = index[ih].at(ev).weight;
            if (std::isnan(weight) || std::isinf(weight) || weight<=0.){ cerr << "Invalid weight " << weight << " is being discarded." << endl; continue; }
            if (KD==-999.) continue;
            sumWgtAll[ih] += weight;
            if (
              (KD>=0.5 && ih==0)
              ||
              (KD<0.5 && ih==1)
              ) sumWgtHalf[ih] += weight;
          }
          sumWgtHalf[ih]=sumWgtHalf[ih]/sumWgtAll[ih];
        }

        if (mindiff>fabs(sumWgtHalf[0]-sumWgtHalf[1])){
          finalFraction[0] = sumWgtHalf[0];
          finalFraction[1] = sumWgtHalf[1];
          mindiff=fabs(sumWgtHalf[0]-sumWgtHalf[1]);
          Cfound=testC;
        }
      }
      cout << "  - New c found = " << Cfound << " (old was " << centralConstant << ")" << endl;
      cout << "  - Final fractions were = " << finalFraction[0] << " , " << finalFraction[1] << endl;
      if (fabs(Cfound/centralConstant-1)<1e-4 && it>0) break;
      centralConstant=Cfound;
      if (it>2){
        marginhigh /= float(nsteps)/10.;
        marginlow /= float(nsteps)/10.;
      }
      else{
        marginhigh /= 5.;
        marginlow /= 2.;
      }
      Cfound=0;
      it++;
    }

    hrec->SetBinContent(bin+1, centralConstant);
  }
}
void getKDConstantByMass(
  float sqrts, TString strname,
  vector<TString>& strRecoBranch,
  vector<TString>(&strSamples)[2],
  vector<TString>(&strExtraWeights)[2],
  float(*KDfcn)(vector<float>&, float),
  unsigned int divisor,
  const bool writeFinalTree,
  TString strcustomselection="",
  vector<pair<vector<float>, pair<float, float>>>* manualboundary_validity_pairs=0
  ){
  short Z1Id, Z2Id;
  float DiJetMass;

  const TString strTrueBranch = "ZZMass";

  vector<TString> strWeights; strWeights.push_back(TString("overallEventWeight")); strWeights.push_back(TString("xsec"));
  strWeights.push_back(TString("KFactor_QCD_ggZZ_Nominal")); strWeights.push_back(TString("KFactor_EW_qqZZ")); strWeights.push_back(TString("KFactor_QCD_qqZZ_M"));

  float varTrack;
  vector<float> valReco; for (unsigned int i=0; i<strRecoBranch.size(); i++) valReco.push_back(float(0));
  vector<TString> strAllWeights[2];
  for (unsigned int ih=0; ih<2; ih++){
    for (auto& s : strWeights) strAllWeights[ih].push_back(s);
    for (auto& s : strExtraWeights[ih]) strAllWeights[ih].push_back(s);
  }

  vector<SimpleEntry> index[2];

  gSystem->Exec("mkdir -p ./output/KDConstants");
  TString coutput = Form("KDConstant_m4l_%s", strname.Data());
  if (strcustomselection!="") coutput += Form("_%s", strcustomselection.Data());
  if (sqrts>0.) coutput += Form("%.0fTeV", sqrts);
  TFile* foutput = TFile::Open(Form("./output/KDConstants/%s%s", coutput.Data(), ".root"), "recreate");

  int nEntries[2]={ 0 };
  float infimum=0;
  float supremum=sqrts*1000.;
  for (unsigned int ih=0; ih<2; ih++){
    // Keep track of TFile and TTree objects
    vector<TFile*> finputList;
    vector<pair<TTree*, TH1F*>> treeList;

    vector<pair<TString, float*>> str_weight;
    vector<float> weights; for (unsigned int i=0; i<strAllWeights[ih].size(); i++) weights.push_back(float(1)); // overallEventWeight, xsec, extras
    for (unsigned int i=0; i<strAllWeights[ih].size(); i++) str_weight.push_back(pair<TString, float*>(strAllWeights[ih].at(i), &(weights.at(i))));

    vector<pair<TString, float*>> str_val_reco;
    for (unsigned int i=0; i<strRecoBranch.size(); i++) str_val_reco.push_back(pair<TString, float*>(strRecoBranch.at(i), &(valReco.at(i))));

    bookSampleTrees(
      sqrts, strSamples[ih],
      &Z1Id, &Z2Id, &DiJetMass,
      strTrueBranch, &varTrack, str_val_reco, str_weight,
      nEntries[ih], finputList, treeList
      );
    getEvents(
      treeList,
      nEntries[ih],
      Z1Id, Z2Id, DiJetMass,
      varTrack,
      valReco, strRecoBranch,
      weights,
      index[ih],
      KDfcn,
      strcustomselection
      );
    for (unsigned int f=0; f<finputList.size(); f++) finputList.at(f)->Close();

    float firstVal=index[ih].at(0).trackingval; firstVal = (float)((int)firstVal); firstVal -= (float)(((int)firstVal)%10);
    float lastVal=index[ih].at(index[ih].size()-1).trackingval; lastVal = (float)((int)(lastVal+0.5)); lastVal += (float)(10-((int)lastVal)%10);
    infimum = max(firstVal, infimum);
    supremum = min(lastVal, supremum);
  }
  for (unsigned int ih=0; ih<2; ih++){
    SimpleEntry::cropByTrueVal(index[ih], infimum, supremum);
    nEntries[ih]=index[ih].size();
    cout << "Nentries remaining = " << nEntries[ih] << " | truth = [ " << infimum << " , " << supremum << " ]" << endl;
  }
  
  vector<unsigned int> indexboundaries[2];
  for (unsigned int ih=0; ih<2; ih++) indexboundaries[ih].push_back(0);
  {
    unsigned int iit[2]={ divisor, divisor };
    while (iit[0]<index[0].size() && iit[1]<index[1].size()){
      if (
        (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
        <
        (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
        ){
        while (
          (iit[0]+1)<index[0].size() &&
          (
          (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          <
          (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
          )
          &&
          (
          (index[0].at(iit[0]+1).trackingval+index[0].at(iit[0]).trackingval)*0.5
          <
          (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
          )
          ) iit[0]++;
      }
      else if (
        (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
        >
        (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
        ){
        while (
          (iit[1]+1)<index[1].size() &&
          (
          (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          >
          (index[1].at(iit[1]).trackingval+index[1].at(iit[1]-1).trackingval)*0.5
          )
          &&
          (
          (index[0].at(iit[0]).trackingval+index[0].at(iit[0]-1).trackingval)*0.5
          >
          (index[1].at(iit[1]+1).trackingval+index[1].at(iit[1]).trackingval)*0.5
          )
          ) iit[1]++;
      }

      if (
        (index[0].size()-iit[0])<divisor
        ||
        (index[1].size()-iit[1])<divisor
        ) break;

      for (unsigned int ih=0; ih<2; ih++){ indexboundaries[ih].push_back(iit[ih]); iit[ih] += divisor; }
    }
  }
  for (unsigned int ih=0; ih<2; ih++) indexboundaries[ih].push_back(index[ih].size());
  for (unsigned int ih=0; ih<2; ih++) cout << "Final size of indexboundaries[" << ih << "] = " << indexboundaries[ih].size() << endl;

  unsigned int nbins=indexboundaries[0].size()-1;
  vector<float> binboundarylist;
  binboundarylist.push_back(infimum);
  for (unsigned int ix=1; ix<nbins; ix++){
    float binboundary = max(
      (index[0].at(indexboundaries[0].at(ix)-1).trackingval+index[0].at(indexboundaries[0].at(ix)).trackingval)*0.5
      ,
      (index[1].at(indexboundaries[1].at(ix)-1).trackingval+index[1].at(indexboundaries[1].at(ix)).trackingval)*0.5
      );

    cout << "Initial bin boundary for bin " << ix << ": " << binboundary << endl;

    bool skip=false;
    if (manualboundary_validity_pairs!=0){
      for (auto const& mbvpair : (*manualboundary_validity_pairs)){
        const pair<float, float>& valrange = mbvpair.second;
        skip = (
          (valrange.first<0 || valrange.first<=binboundary)
          &&
          (valrange.second<0 || valrange.second>binboundary)
          );
        if (skip) break;
      }
    }
    if (!skip) binboundarylist.push_back(binboundary);
  }
  binboundarylist.push_back(supremum);

  if (manualboundary_validity_pairs!=0){
    for (auto& mbvpair : (*manualboundary_validity_pairs)){
      vector<float>& bvals = mbvpair.first;
      for (unsigned int ib=0; ib<bvals.size(); ib++){
        float& bval = bvals.at(ib);
        cout << "Adding manual boundary " << bval << endl;
        addByLowest<float>(binboundarylist, bval, true);
      }
    }
  }

  nbins = binboundarylist.size()-1;
  float* binning = new float[nbins+1];
  for (unsigned int ix=0; ix<=nbins; ix++){
    binning[ix] = binboundarylist[ix];
    cout << "Boundary (" << ix << ")= " << binning[ix] << endl;
  }

  // Recalibrate index boundaries
  if (manualboundary_validity_pairs!=0){
    for (unsigned int ih=0; ih<2; ih++){
      indexboundaries[ih].clear();
      indexboundaries[ih].push_back(0);
      unsigned int ix=1;
      for (unsigned int ev=1; ev<index[ih].size(); ev++){
        if (ix==nbins) break;
        if (index[ih].at(ev).trackingval>=binning[ix]){
          indexboundaries[ih].push_back(ev);
          ix++;
        }
      }
      indexboundaries[ih].push_back(index[ih].size());
    }
  }

  foutput->cd();

  TH1F* h_varTrack_Constant = new TH1F("varTrue_Constant", "", nbins, binning); h_varTrack_Constant->Sumw2();
  TProfile* p_varTrack = new TProfile("avg_varTrue", "", nbins, binning); p_varTrack->Sumw2();
  delete[] binning;

  LoopForConstant(
    index, indexboundaries,
    KDfcn,
    p_varTrack,
    h_varTrack_Constant,
    100
    );

  TGraphErrors* gr = makeGraphFromTH1(p_varTrack, h_varTrack_Constant, "gr_varTrue_Constant");
  foutput->WriteTObject(p_varTrack);
  foutput->WriteTObject(h_varTrack_Constant);
  foutput->WriteTObject(gr);
  delete gr;
  delete h_varTrack_Constant;
  delete p_varTrack;
  foutput->Close();
}

/*
SPECIFIC COMMENT:
- Multiplies by Pmjj
*/
void getKDConstant_DjjVH(TString strprod, float sqrts=13){
  const bool writeFinalTree=false;
  float divisor=20000;

  vector<TString> extraweights[2];

  vector<TString> strRecoBranch;
  if (strprod=="ZH_Sig_POWHEG"){
    strRecoBranch.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    strRecoBranch.push_back("p_HadZH_mavjj_JECNominal");
  }
  else if (strprod=="WH_Sig_POWHEG"){
    strRecoBranch.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    strRecoBranch.push_back("p_HadWH_mavjj_JECNominal");
  }
  strRecoBranch.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
  vector<TString> strSamples[2];
  if (strprod=="ZH_Sig_POWHEG"){
    vector<TString> s1; s1.push_back("ZH_Sig_POWHEG");
    vector<TString> s2; s2.push_back("gg_Sig_POWHEG");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }
  else if (strprod=="WH_Sig_POWHEG"){
    vector<TString> s1; s1.push_back("WH_Sig_POWHEG");
    vector<TString> s2; s2.push_back("gg_Sig_POWHEG");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }
  else{
    cerr << "Production " << strprod << " is unknown." << endl;
    assert(0);
  }

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  if (strprod=="ZH_Sig_POWHEG"){
    pair<float, float> valrange(230, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(245);
    manualboundaries.push_back(300);
    manualboundaries.push_back(500);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  else if (strprod=="WH_Sig_POWHEG"){
    pair<float, float> valrange(195, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(196);
    manualboundaries.push_back(209);
    manualboundaries.push_back(300);
    manualboundaries.push_back(500);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }

  getKDConstantByMass(
    sqrts, Form("Djj%s", strprod.Data()),
    strRecoBranch, strSamples, extraweights,
    constructDjjVH, divisor, writeFinalTree, "",
    &manualboundary_validity_pairs
    );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_DjjVBF(float sqrts=13){
  const bool writeFinalTree=false;
  TString strprod="VBF";
  float divisor=40000;

  vector<TString> extraweights[2];

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal");
  strRecoBranch.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
  vector<TString> strSamples[2];
  {
    vector<TString> s1; s1.push_back("VBF_Sig_POWHEG");
    vector<TString> s2; s2.push_back("gg_Sig_POWHEG");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  {
    pair<float, float> valrange(750, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(770);
    manualboundaries.push_back(850);
    manualboundaries.push_back(1100);
    manualboundaries.push_back(1400);
    //manualboundaries.push_back(1700);
    manualboundaries.push_back(2100);
    manualboundaries.push_back(2650);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }

  getKDConstantByMass(
    sqrts, Form("Djj%s", strprod.Data()),
    strRecoBranch, strSamples, extraweights,
    constructDjjVBF, divisor, writeFinalTree, "",
    &manualboundary_validity_pairs
    );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_DjVBF(float sqrts=13){
  const bool writeFinalTree=false;
  TString strprod="VBF";
  float divisor=40000;

  vector<TString> extraweights[2];

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_JVBF_SIG_ghv1_1_JHUGen_JECNominal");
  strRecoBranch.push_back("pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal");
  strRecoBranch.push_back("p_JQCD_SIG_ghg2_1_JHUGen_JECNominal");
  vector<TString> strSamples[2];
  {
    vector<TString> s1; s1.push_back("VBF_Sig_POWHEG");
    vector<TString> s2; s2.push_back("gg_Sig_POWHEG");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  {
    pair<float, float> valrange(1000, 2000);
    vector<float> manualboundaries;
    manualboundaries.push_back(1400);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }

  getKDConstantByMass(
    sqrts, Form("Dj%s", strprod.Data()),
    strRecoBranch, strSamples, extraweights,
    constructDjVBF, divisor, writeFinalTree, "",
    &manualboundary_validity_pairs
    );
}

/*
SPECIFIC COMMENT:
Add bin boundaries 75, 105 and 120 manually
*/
void getKDConstant_Dbkgkin(TString strchannel, float sqrts=13){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;

  const bool writeFinalTree=false;
  //float divisor=25000;
  //if (strchannel=="2l2l" || strchannel=="2e2mu") divisor = 50000;
  float divisor=21000;
  if (strchannel=="2l2l" || strchannel=="2e2mu") divisor = 50000;

  vector<TString> extraweights[2];
  extraweights[1].push_back(TString("p_Gen_QQB_BKG_MCFM"));

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen");
  strRecoBranch.push_back("p_QQB_BKG_MCFM");
  vector<TString> strSamples[2];
  {
    vector<TString> s1; s1.push_back("gg_Sig_POWHEG"); s1.push_back("VBF_Sig_POWHEG"); s1.push_back("gg_Sig_SM_MCFM");
    //vector<TString> s2; s2.push_back("qq_Bkg_Combined");
    vector<TString> s2; s2.push_back("qq_Bkg_Combined"); s2.push_back("gg_Bkg_MCFM");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }

  // Define manual bin boundaries
  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  {
    pair<float, float> valrange(70, 142);
    vector<float> manualboundaries;
    manualboundaries.push_back(75); manualboundaries.push_back(85);
    manualboundaries.push_back(89); manualboundaries.push_back(93); manualboundaries.push_back(96);
    manualboundaries.push_back(100); manualboundaries.push_back(105); manualboundaries.push_back(110); manualboundaries.push_back(115);
    manualboundaries.push_back(120);  manualboundaries.push_back(123); manualboundaries.push_back(135);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  {
    pair<float, float> valrange(600, 2500);
    vector<float> manualboundaries;
    manualboundaries.push_back(700); manualboundaries.push_back(900); manualboundaries.push_back(1100); manualboundaries.push_back(1400); manualboundaries.push_back(1900);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }

  getKDConstantByMass(
    sqrts, "Dbkgkin",
    strRecoBranch, strSamples, extraweights,
    constructDbkgkin, divisor, writeFinalTree,
    strchannel,
    &manualboundary_validity_pairs
    );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_Dbkgdec(TString strchannel, float sqrts=13){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;

  const bool writeFinalTree=false;
  float divisor=25000;
  if (strchannel=="2l2l" || strchannel=="2e2mu") divisor = 40000;

  vector<TString> extraweights[2];

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_GG_SIG_kappaTopBot_1_ghz1_1_MCFM");
  strRecoBranch.push_back("p_QQB_BKG_MCFM");
  strRecoBranch.push_back("p_GG_BKG_MCFM");
  vector<TString> strSamples[2];
  {
    vector<TString> s1; s1.push_back("gg_Sig_POWHEG"); s1.push_back("VBF_Sig_POWHEG"); s1.push_back("gg_Sig_SM_MCFM");
    vector<TString> s2; s2.push_back("qq_Bkg_Combined"); s2.push_back("gg_Bkg_MCFM");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }
  getKDConstantByMass(
    sqrts, "Dbkgdec",
    strRecoBranch, strSamples, extraweights,
    constructDbkgdec, divisor, writeFinalTree,
    strchannel
    );
}

/* SPECIFIC COMMENT: NONE */
void getKDConstant_DbkgjjEWQCD(TString strchannel, float sqrts=13){
  if (!(strchannel.Contains("2l2l") || strchannel.Contains("4l"))) return;
  if (!(strchannel.Contains("VBFLike") || strchannel.Contains("VHLike_LowMass"))) return;

  const bool writeFinalTree=false;
  float divisor=40000;
  if (strchannel.Contains("VHLike_LowMass")) divisor = 1000;

  vector<TString> extraweights[2];

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_JJVBF_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("p_HadZH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("p_Had_WH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("p_JJVBF_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("p_HadZH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("p_Had_WH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("p_JJQCD_BKG_MCFM_JECNominal");

  strRecoBranch.push_back("p_HadZH_mavjj_JECNominal");
  strRecoBranch.push_back("p_HadWH_mavjj_JECNominal");

  strRecoBranch.push_back("pConst_JJVBF_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_HadZH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_Had_WH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_JJVBF_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_HadZH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_Had_WH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_JJQCD_BKG_MCFM_JECNominal");

  strRecoBranch.push_back("DiJetMass");

  vector<TString> strSamples[2];
  {
    vector<TString> s1; s1.push_back("VBF_Sig_POWHEG"); s1.push_back("ZH_Sig_POWHEG"); s1.push_back("WH_Sig_POWHEG");
    vector<TString> s2; s2.push_back("qq_Bkg_Combined"); s2.push_back("VV_Bkg_Phantom");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }

  vector<pair<vector<float>, pair<float, float>>> manualboundary_validity_pairs;
  /*
  {
    pair<float, float> valrange(70, 120);
    vector<float> manualboundaries;
    manualboundaries.push_back(105);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  */
  if (strchannel.Contains("VBFLike")){
    pair<float, float> valrange(750, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(770);
    manualboundaries.push_back(850);
    manualboundaries.push_back(1100);
    manualboundaries.push_back(1400);
    //manualboundaries.push_back(1700);
    manualboundaries.push_back(2100);
    manualboundaries.push_back(2650);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }

  getKDConstantByMass(
    sqrts, "DbkgjjEWQCD",
    strRecoBranch, strSamples, extraweights,
    constructDbkgkin, divisor, writeFinalTree, strchannel,
    &manualboundary_validity_pairs
    );
}


void generic_SmoothKDConstantProducer(
  float sqrts, TString strname, TString strcustomselection,
  TF1* (*lowf)(TSpline3*, double, double, bool),
  TF1* (*highf)(TSpline3*, double, double, bool),
  bool useFaithfulSlopeFirst, bool useFaithfulSlopeSecond,
  vector<pair<pair<double, double>, unsigned int>>* addpoints=0
  ){
  const double xmin=0;
  const double xmax=(sqrts>0 ? (double)sqrts*1000. : 15000.);

  gSystem->Exec("mkdir -p ./output/KDConstants");
  TString cinput = Form("KDConstant_m4l_%s", strname.Data());
  if (strcustomselection!="") cinput += Form("_%s", strcustomselection.Data());
  if (sqrts>0.) cinput += Form("%.0fTeV", sqrts);

  TFile* finput = TFile::Open(Form("./output/KDConstants/%s%s", cinput.Data(), ".root"), "read");
  TFile* foutput = TFile::Open(Form("./output/KDConstants/Smooth%s%s", cinput.Data(), ".root"), "recreate");
  foutput->cd();

  //TH1F* h_varTrack_Constant = (TH1F*)finput->Get("varTrack_Constant");
  //TProfile* p_varTrack = (TProfile*)finput->Get("avg_varTrack");
  TGraphErrors* tg = (TGraphErrors*)finput->Get("gr_varTrack_Constant");
  foutput->WriteTObject(tg);

  if (addpoints!=0){ for (auto& prange : *addpoints) addPointsBetween(tg, prange.first.first, prange.first.second, prange.second); }

  int n = tg->GetN();
  double* xx = tg->GetX();
  double* ex = tg->GetEX();
  double* yy = tg->GetY();
  double* ey = tg->GetEY();

  TSpline3* sp;
  sp = convertGraphToSpline3_MightFaithfulSlopes(tg, useFaithfulSlopeFirst, useFaithfulSlopeSecond);
  double tglow = xx[0];
  double tghigh = xx[tg->GetN()-1];
  TF1* lowFcn = lowf(sp, xmin, tglow, true);
  TF1* highFcn = highf(sp, tghigh, xmax, false);
  lowFcn->SetNpx((int)(tglow-xmin)*5);
  highFcn->SetNpx((int)(xmax-tghigh)*5);

  vector<pair<double, double>> points;
  for (double xval=xmin; xval<tglow; xval+=1){
    double yval = lowFcn->Eval(xval);
    addByLowest<double, double>(points, xval, yval);
  }
  for (int ix=0; ix<n; ix++){
    addByLowest<double, double>(points, xx[ix], yy[ix]);
  }
  int tghigh_int = ((int)((tghigh+1.)/100.+0.5))*100;
  if (tghigh>=(double)tghigh_int) tghigh_int+=100;
  for (double xval=tghigh_int; xval<=xmax; xval+=100){
    double yval = highFcn->Eval(xval);
    addByLowest<double, double>(points, xval, yval);
  }

  int nn_new = points.size();
  cout << "Number of new points: " << nn_new-n << endl;
  double* xy_new[2];
  for (unsigned int i=0; i<2; i++) xy_new[i] = new double[nn_new];
  for (int ix=0; ix<nn_new; ix++){
    xy_new[0][ix] = points.at(ix).first;
    xy_new[1][ix] = points.at(ix).second;
  }

  delete highFcn;
  delete lowFcn;
  delete sp;

  TGraph* tg_updated = new TGraph(nn_new, xy_new[0], xy_new[1]);
  tg_updated->SetName(Form("%s_Smooth", tg->GetName()));
  foutput->WriteTObject(tg_updated);

  sp = convertGraphToSpline3(tg_updated);
  foutput->WriteTObject(sp);
  delete sp;
  delete tg_updated;
  for (unsigned int i=0; i<2; i++) delete[] xy_new[i];

  foutput->Close();
  finput->Close();
}

void SmoothKDConstantProducer_DjjVH(TString strprod){
  generic_SmoothKDConstantProducer(
    13, Form("Djj%s", strprod.Data()), "",
    &getFcn_a0plusa1timesX,
    &getFcn_a0timesexpa1X,
    false, false
    );
}

void SmoothKDConstantProducer_DjjVBF(){
  TString strprod="VBF";

  vector<pair<pair<double, double>, unsigned int>> addpoints;
  {
    pair<double, double> xminmax(1000, 3000); unsigned int nadd=5;
    pair<pair<double, double>, unsigned int> addsingle(xminmax, nadd); addpoints.push_back(addsingle);
  }

  generic_SmoothKDConstantProducer(
    13, Form("Djj%s", strprod.Data()), "",
    &getFcn_a0plusa1timesX,
    &getFcn_a0timesexpa1X,
    //&getFcn_a0plusa1overX,
    true, true,
    &addpoints
    );
}

void SmoothKDConstantProducer_DjVBF(){
  TString strprod="VBF";
  generic_SmoothKDConstantProducer(
    13, Form("Dj%s", strprod.Data()), "",
    &getFcn_a0plusa1timesX,
    &getFcn_a0timesexpa1X,
    true, true
    );
}

void SmoothKDConstantProducer_Dbkgkin(TString strchannel){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;
  generic_SmoothKDConstantProducer(
    13, Form("Dbkgkin_%s", strchannel.Data()), "",
    &getFcn_a0plusa1timesX,
    &getFcn_a0timesexpa1X,
    //&getFcn_a0plusa1overXN<6>,
    true, false
    );
}


/*
g-constants for AC
*/

TGraph* getSingleTGraph(TString fname){
  TDirectory* curdir = gDirectory;

  cout << "Opening file " << fname << endl;

  TFile* finput = TFile::Open(Form("JHUGenXsec/%s%s", fname.Data(), ".root"), "read");
  TGraph* tgold = (TGraph*)finput->Get("Graph");
  double* xx = tgold->GetX();
  double* yy = tgold->GetY();

  vector<pair<double, double>> xyinterm;
  for (int ip=0; ip<tgold->GetN(); ip++){
    if (std::isnan(xx[ip]) || std::isinf(xx[ip])) continue;
    if (std::isnan(yy[ip]) || std::isinf(yy[ip])) continue;
    xyinterm.push_back(pair<double, double>(xx[ip], yy[ip]));
  }
  TGraph* tginterm = makeGraphFromPair(xyinterm, "tginterm");
  TSpline3* spinterm = convertGraphToSpline3(tginterm);
  delete tginterm;

  vector<pair<double, double>> xyfinal;
  for (int ip=0; ip<tgold->GetN(); ip++){
    if (std::isnan(xx[ip]) || std::isinf(xx[ip])) continue;
    if (std::isnan(yy[ip]) || std::isinf(yy[ip])) yy[ip] = spinterm->Eval(xx[ip]);
    xyfinal.push_back(pair<double, double>(xx[ip], yy[ip]));
  }

  curdir->cd();
  TGraph* result = makeGraphFromPair(xyfinal, Form("tgfinal_%s", fname.Data()));

  finput->Close();
  return result;
}

void generic_gConstantProducer(TString strprod, TString strhypo, bool useproddec=false){
  TString strSM = "SM";
  if (strhypo=="L1Zgs") strSM += "_photoncut";

  TFile* foutput = TFile::Open(Form("gConstant_%s%s_%s%s", strprod.Data(), (useproddec ? "_ProdDec" : ""), strhypo.Data(), ".root"), "recreate");
  TGraph* tgSMhypo = getSingleTGraph(Form("%s_%s", strprod.Data(), strSM.Data()));
  TGraph* tgBSMhypo = getSingleTGraph(Form("%s_%s", strprod.Data(), strhypo.Data()));

  if (useproddec){
    TGraph* tgSMhypodec = getSingleTGraph(Form("HZZ2e2mu_%s", strSM.Data()));
    TGraph* tgBSMhypodec = getSingleTGraph(Form("HZZ2e2mu_%s", strhypo.Data()));

    TGraph* tgSMhyponew = multiplyTGraphs(tgSMhypo, tgSMhypodec);
    TGraph* tgBSMhyponew = multiplyTGraphs(tgBSMhypo, tgBSMhypodec);

    delete tgSMhypo; tgSMhypo=tgSMhyponew;
    delete tgBSMhypo; tgBSMhypo=tgBSMhyponew;
  }

  TGraph* result = divideTGraphs(tgSMhypo, tgBSMhypo, 0.5, 0.5);
  foutput->WriteTObject(result);
  TSpline3* spresult = convertGraphToSpline3(result);
  foutput->WriteTObject(spresult);
  cout << "Result of spline at 125 GeV = " << spresult->Eval(125) << endl;
  delete spresult; delete result;
  foutput->Close();
}

void gConstantProducer(){
  generic_gConstantProducer("WH_Sig_POWHEG", "g2");
  generic_gConstantProducer("WH_Sig_POWHEG", "g4");
  generic_gConstantProducer("WH_Sig_POWHEG", "L1");
  generic_gConstantProducer("ZH_Sig_POWHEG", "g2");
  generic_gConstantProducer("ZH_Sig_POWHEG", "g4");
  generic_gConstantProducer("ZH_Sig_POWHEG", "L1");
  generic_gConstantProducer("ZH_Sig_POWHEG", "L1Zgs");
  generic_gConstantProducer("VBF", "g2");
  generic_gConstantProducer("VBF", "g4");
  generic_gConstantProducer("VBF", "L1");
  generic_gConstantProducer("VBF", "L1Zgs");
  generic_gConstantProducer("HZZ2e2mu", "g2");
  generic_gConstantProducer("HZZ2e2mu", "g4");
  generic_gConstantProducer("HZZ2e2mu", "L1");
  generic_gConstantProducer("HZZ2e2mu", "L1Zgs");
}


void testDbkgkinGGZZvsQQZZ(){
  const unsigned int nsamples=3;
  vector<TString> strSamples[nsamples];
  vector<TString> s1; s1.push_back("qq_Bkg_Combined"); getSamplesList(13, s1, strSamples[0]);
  vector<TString> s2; s2.push_back("gg_Bkg_MCFM"); getSamplesList(13, s2, strSamples[1]);
  vector<TString> s3; s3.push_back("gg_Sig_POWHEG"); s3.push_back("VBF_Sig_POWHEG"); s3.push_back("gg_Sig_SM_MCFM"); getSamplesList(13, s3, strSamples[2]);

  TString strchannel[3]={ "4e", "4mu", "2e2mu" };

  float ZZMass;
  vector<float> vars;
  float extrawgt=1;

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_GG_SIG_ghg2_1_ghz1_1_JHUGen"); vars.push_back(float(0));
  strRecoBranch.push_back("p_QQB_BKG_MCFM"); vars.push_back(float(0));

  TString cinput_main = CJLSTsamplesdir;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";

  TH2F* hist[nsamples];
  TH1F* hist1d[nsamples];
  TFile* foutput = TFile::Open("ftmp.root", "recreate");

  TFile* fcKD[3];
  TSpline3* spcKD[3];
  for (unsigned int ic=0; ic<3; ic++){
    fcKD[ic] = TFile::Open(Form("SmoothKDConstant_m4l_Dbkgkin_%s13TeV.root", strchannel[ic].Data()), "read");
    spcKD[ic] = (TSpline3*)fcKD[ic]->Get("sp_gr_varTrue_Constant_Smooth");
  }

  for (unsigned int ih=0; ih<nsamples; ih++){
    foutput->cd();
    hist[ih]=new TH2F(Form("hist%i", ih), "", 68, 600, 4000, 20, 0, 1);
    hist1d[ih]=new TH1F(Form("hist1d%i", ih), "", 68, 600, 4000);
    //hist[ih]=new TH2F(Form("hist%i", ih), "", 50, 70, 170, 30, 0, 1);
    //hist1d[ih]=new TH1F(Form("hist1d%i", ih), "", 50, 70, 170);
    //hist[ih]=new TH2F(Form("hist%i", ih), "", 145, 100, 3000, 20, 0, 1);
    //hist1d[ih]=new TH1F(Form("hist1d%i", ih), "", 145, 100, 3000);
    for (int is=0; is<(int)strSamples[ih].size(); is++){
      TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[ih][is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;

            tree->SetBranchStatus("*", 0);
            float* ZZMassptr = &ZZMass;
            short Z1Flav, Z2Flav;
            short* Z1Flavptr = &Z1Flav;
            short* Z2Flavptr = &Z2Flav;
            bookBranch(tree, "ZZMass", ZZMassptr);
            bookBranch(tree, "Z1Flav", Z1Flavptr);
            bookBranch(tree, "Z2Flav", Z2Flavptr);
            float* extrawgtptr = &extrawgt;
            if (ih==1) bookBranch(tree, "p_Gen_QQB_BKG_MCFM", extrawgtptr);
            for (unsigned int iv=0; iv<vars.size(); iv++){ float* varptr = &(vars[iv]); bookBranch(tree, strRecoBranch[iv], varptr); }

            for (int ev=0; ev<tree->GetEntries(); ev++){
              extrawgt=1;
              tree->GetEntry(ev);
              
              unsigned int ic =
                0*(Z1Flav*Z2Flav==pow(11*11, 2))
                +1*(Z1Flav*Z2Flav==pow(13*13, 2))
                +2*(Z1Flav*Z2Flav==pow(11*13, 2));
              float cKDval = spcKD[ic]->Eval(ZZMass);

              float KD = constructDbkgkin(vars, cKDval);
              hist[ih]->Fill(ZZMass, KD, extrawgt);
              hist1d[ih]->Fill(ZZMass, extrawgt);
            }

            finput->Close();
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    for (int binx=1; binx<=hist[ih]->GetNbinsX(); binx++){
      double integral = hist[ih]->Integral(binx, binx, 1, hist[ih]->GetNbinsY());
      if (integral!=0.){
        for (int biny=1; biny<=hist[ih]->GetNbinsY(); biny++){
          hist[ih]->SetBinContent(binx, biny, hist[ih]->GetBinContent(binx, biny)/integral);
        }
      }
    }
    foutput->WriteTObject(hist[ih]);
    foutput->WriteTObject(hist1d[ih]);
    delete hist[ih];
    delete hist1d[ih];
  }
  for (unsigned int ic=0; ic<3; ic++) fcKD[ic]->Close();
  foutput->Close();
}

void testDjjEWQCDEWvsQCD(TString strcustomselection=""){
  const unsigned int nsamples=3;
  vector<TString> strSamples[nsamples];
  vector<TString> s1; s1.push_back("JJVBF"); s1.push_back("ZH"); s1.push_back("WH");
  vector<TString> s2; s2.push_back("qq_Bkg_Combined"); s2.push_back("VV_Bkg_Phantom");
  vector<TString> s3; s3.push_back("JJQCD"); s3.push_back("gg_Sig_MCFM");
  getSamplesList(13, s1, strSamples[0]);
  getSamplesList(13, s2, strSamples[1]);
  getSamplesList(13, s3, strSamples[2]);

  TString strchannel[2]={ "4l", "2l2l" };

  pair<float, float> mjjcut(-1, -1);
  if (strcustomselection.Contains("VBFLike")){ mjjcut.first=150; }
  else if (strcustomselection.Contains("VHLike")){ mjjcut.first=40; mjjcut.second=150; }

  pair<float, float> m4lcut(-1, -1);
  if (strcustomselection.Contains("LowMass")){ m4lcut.first=100; m4lcut.second=160; }
  else if (strcustomselection.Contains("HighMass")){ m4lcut.first=160; }

  float ZZMass;
  vector<float> vars;
  float extrawgt=1;

  vector<TString> strRecoBranch;
  strRecoBranch.push_back("p_JJVBF_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("p_HadZH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("p_Had_WH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("p_JJVBF_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("p_HadZH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("p_Had_WH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("p_JJQCD_BKG_MCFM_JECNominal");

  strRecoBranch.push_back("p_HadZH_mavjj_JECNominal");
  strRecoBranch.push_back("p_HadWH_mavjj_JECNominal");

  strRecoBranch.push_back("pConst_JJVBF_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_HadZH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_Had_WH_SIG_ghv1_1_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_JJVBF_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_HadZH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_Had_WH_BKG_MCFM_JECNominal");
  strRecoBranch.push_back("pConst_JJQCD_BKG_MCFM_JECNominal");

  strRecoBranch.push_back("DiJetMass");

  for (auto const& b:strRecoBranch) vars.push_back(float(0));

  TString cinput_main = CJLSTsamplesdir;
  TString TREE_NAME = "ZZTree/candTree";
  TString COUNTERS_NAME = "ZZTree/Counters";

  TH2F* hist[nsamples];
  TH1F* hist1d[nsamples];
  TFile* foutput = TFile::Open("ftmp.root", "recreate");

  TFile* fcKD[3];
  TSpline3* spcKD[3];
  //for (unsigned int ic=0; ic<3; ic++){
  //  fcKD[ic] = TFile::Open(Form("SmoothKDConstant_m4l_Dbkgkin_%s13TeV.root", strchannel[ic].Data()), "read");
  //  spcKD[ic] = (TSpline3*)fcKD[ic]->Get("sp_gr_varTrue_Constant_Smooth");
  //}

  for (unsigned int ih=0; ih<nsamples; ih++){
    foutput->cd();
    float xlow = 70;
    float xhigh = 4000;
    if (m4lcut.first>=0.) xlow = m4lcut.first;
    if (m4lcut.second>=0.) xhigh = m4lcut.second;
    float xwidth=2;
    if (xhigh>1000.) xwidth=50;
    int nxbins = (xhigh-xlow)/xwidth;
    hist[ih]=new TH2F(Form("hist%i", ih), "", nxbins, xlow, xhigh, 20, 0, 1);
    hist1d[ih]=new TH1F(Form("hist1d%i", ih), "", nxbins, xlow, xhigh);
    for (int is=0; is<(int)strSamples[ih].size(); is++){
      TString cinput = Form("%s/%s/ZZ4lAnalysis.root", cinput_main.Data(), (strSamples[ih][is]).Data());
      TFile* finput = TFile::Open(cinput, "read");
      cout << "Opening file " << cinput << "..." << endl;
      TTree* tree=0;
      if (finput!=0){
        if (finput->IsOpen() && !finput->IsZombie()){
          cout << cinput << " opened. Extracting tree " << TREE_NAME << "..." << endl;
          tree = (TTree*)finput->Get(TREE_NAME);
          if (tree!=0){
            cout << TREE_NAME << " is found." << endl;

            tree->SetBranchStatus("*", 0);
            float* ZZMassptr = &ZZMass;
            short Z1Flav, Z2Flav;
            short* Z1Flavptr = &Z1Flav;
            short* Z2Flavptr = &Z2Flav;
            bookBranch(tree, "ZZMass", ZZMassptr);
            bookBranch(tree, "Z1Flav", Z1Flavptr);
            bookBranch(tree, "Z2Flav", Z2Flavptr);
            float* extrawgtptr = &extrawgt;
            //if (ih==1) bookBranch(tree, "p_Gen_QQB_BKG_MCFM", extrawgtptr);
            for (unsigned int iv=0; iv<vars.size(); iv++){ float* varptr = &(vars[iv]); bookBranch(tree, strRecoBranch[iv], varptr); }

            for (int ev=0; ev<tree->GetEntries(); ev++){
              extrawgt=1;
              tree->GetEntry(ev);

              float& DiJetMass = vars[16];
              bool testMatchRecoMJJ = (mjjcut.first<0. || (mjjcut.first<=DiJetMass && DiJetMass>=0.)) && (mjjcut.second<0. || (DiJetMass<=mjjcut.second && DiJetMass>=0.));
              if (!testMatchRecoMJJ || DiJetMass<0.) continue;

              bool testMatchRecoM4L = (m4lcut.first<0. || m4lcut.first<=ZZMass) && (m4lcut.second<0. || ZZMass<=m4lcut.second);
              if (!testMatchRecoM4L) continue;

              unsigned int ic =
                0*(Z1Flav*Z2Flav==pow(11*11, 2) || Z1Flav*Z2Flav==pow(13*13, 2))
                +1*(Z1Flav*Z2Flav==pow(11*13, 2));
              float cKDval = /*spcKD[ic]->Eval(ZZMass)*/1;

              float KD = constructDjjbkgEWQCD(vars, cKDval);
              hist[ih]->Fill(ZZMass, KD, extrawgt);
              hist1d[ih]->Fill(ZZMass, extrawgt);
            }

            finput->Close();
          }
          else finput->Close();
        }
        else if (finput->IsOpen()) finput->Close();
      }
    }
    for (int binx=1; binx<=hist[ih]->GetNbinsX(); binx++){
      double integral = hist[ih]->Integral(binx, binx, 1, hist[ih]->GetNbinsY());
      if (integral!=0.){
        for (int biny=1; biny<=hist[ih]->GetNbinsY(); biny++){
          hist[ih]->SetBinContent(binx, biny, hist[ih]->GetBinContent(binx, biny)/integral);
        }
      }
    }
    foutput->WriteTObject(hist[ih]);
    foutput->WriteTObject(hist1d[ih]);
    delete hist[ih];
    delete hist1d[ih];
  }
  //for (unsigned int ic=0; ic<3; ic++) fcKD[ic]->Close();
  foutput->Close();
}
