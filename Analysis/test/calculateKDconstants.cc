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
#include "HelperFunctions.h"
#include "Samples.h"
#include "SampleHelpers.h"
#include "CJLSTSet.h"
#include "DiscriminantClasses.h"


#ifndef doDebugKD
#define doDebugKD true
#endif
#ifndef doDebugKDExt
#define doDebugKDExt false
#endif

using namespace std;
using namespace HelperFunctions;
using namespace SampleHelpers;
using namespace DiscriminantClasses;


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
CJLSTSet* bookSampleTrees(
  const vector<TString>& strSamples,

  const TString& strTrackVar,
  const vector<TString>& strExtraWeightsList,
  const vector<TString>& strKDVarsList,

  const vector<TString>& strExtraFloatVarsList,
  const vector<TString>& strExtraShortVarsList
  ){
  cout << "Begin bookSampleTrees" << endl;

  CJLSTSet* theSampleSet = new CJLSTSet(strSamples);

  cout << "bookSampleTrees: Calling theSampleSet->bookOverallEventWgt" << endl;
  theSampleSet->bookOverallEventWgt();
  cout << "bookSampleTrees: Calling theSampleSet->setPermanentWeights" << endl;
  theSampleSet->setPermanentWeights(true, true, true, false);

  int nEntries = 0;
  cout << "bookSampleTrees: Looping over the trees in theSampleSet" << endl;
  for (auto& tree:theSampleSet->getCJLSTTreeList()){
    tree->bookBranch<float>(strTrackVar, 0);
    for (auto& s:strExtraWeightsList) tree->bookBranch<float>(s, 1);
    for (auto& s:strKDVarsList) tree->bookBranch<float>(s, 0);
    for (auto& s:strExtraFloatVarsList) tree->bookBranch<float>(s, 0);
    for (auto& s:strExtraShortVarsList) tree->bookBranch<short>(s, 0);
    tree->silenceUnused();
    nEntries += tree->getSelectedNEvents();
    cout << "Events accumulated: " << nEntries << endl;
  }

  cout << "NEntries = " << nEntries << " over " << theSampleSet->getCJLSTTreeList().size() << " trees." << endl;

  cout << "End bookSampleTrees" << endl;
  return theSampleSet;
}
void getEvents(
  CJLSTSet* theSet,

  const TString& strTrackVar,
  const vector<TString>& strExtraWeightsList,
  const vector<TString>& strKDVarsList,

  vector<SimpleEntry>& index,
  Discriminant* const& KDbuilder,
  TString strcustomselection
  ){
  cout << "Begin getKDConstantByMass" << endl;

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

  int ev=0, ev_acc=0;
  CJLSTTree* tree=nullptr;
  while ((tree = theSet->getSelectedEvent(ev))){
    ev++;

    float scale = theSet->getPermanentWeight(tree)*theSet->getOverallEventWgt(tree);

    bool doProcess = true;
    { // Flavor check
      short Z1Id, Z2Id;
      tree->getVal("Z1Flav", Z1Id);
      tree->getVal("Z2Flav", Z2Id);
      bool testMatchDec=(matchdecid.size()==0);
      for (short& testid : matchdecid){ if (testid == Z1Id*Z2Id){ testMatchDec=true; break; } }
      doProcess &= testMatchDec;
    }

    vector<float> valReco;
    for (auto const& s : strKDVarsList){
      float tmp;
      tree->getVal(s, tmp);
      valReco.push_back(tmp);
    }

    float varTrack;
    tree->getVal(strTrackVar, varTrack);
    { // Check m4l cut
      bool testMatchRecoM4L = (m4lcut.first<0. || m4lcut.first<=varTrack) && (m4lcut.second<0. || varTrack<=m4lcut.second);
      doProcess &= testMatchRecoM4L;
    }

    { // Check mJJ cut
      float DiJetMass=-1;
      tree->getVal("DiJetMass", DiJetMass);
      bool testMatchRecoMJJ = (mjjcut.first<0. || (mjjcut.first<=DiJetMass && DiJetMass>=0.)) && (mjjcut.second<0. || (DiJetMass<=mjjcut.second && DiJetMass>=0.));
      doProcess &= testMatchRecoMJJ;
    }

    if (!doProcess) continue;

    // Weight has to be positive-definite for TProfile later on
    float wgt = scale;
    for (auto const& w : strExtraWeightsList){
      float wval;
      tree->getVal(w, wval);
      wgt *= wval;
    }
    if (std::isnan(wgt) || std::isinf(wgt) || wgt<=0.){
      // If weight is NaN, it is a big problem.
      if (std::isnan(wgt) || std::isinf(wgt)) cerr << "Invalid weight " << wgt << " is being discarded at mass " << varTrack << endl;
      continue;
    }

    // Need to also test KD explicitly
    float KD = KDbuilder->update(valReco);
    if (std::isnan(KD) || std::isinf(KD) || KD<0.) continue;

    if (ev_acc%10000==0) cout << "Pre-processing event " << ev << endl;
    if (doDebugKD && ev_acc==20000) break;

    SimpleEntry theEntry(ev, varTrack, valReco, wgt);
    if (doDebugKDExt & doDebugKD) theEntry.print();
    addByLowest(index, theEntry, false);

    ev_acc++;
  }
  cout << "Number of valid entries: " << ev_acc << endl;
  cout << "End getKDConstantByMass" << endl;

}
void LoopForConstant(
  vector<SimpleEntry>(&index)[2],
  vector<unsigned int>(&indexboundaries)[2],
  Discriminant* const& KDbuilder,
  TProfile* px,
  TH1F* hrec,
  unsigned int nstepsiter=100
  ){
  cout << "Begin getKDConstantByMass" << endl;

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
        float KD = KDbuilder->update(index[ih].at(ev).recoval);
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
            KDbuilder->update(index[ih].at(ev).recoval);
            float KD = KDbuilder->applyAdditionalC(testC);
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

  cout << "End getKDConstantByMass" << endl;
}
void getKDConstantByMass(
  float sqrts, TString strname,
  vector<TString>& strRecoBranch,
  vector<TString>& strExtraRecoBranches,
  vector<TString>(&strSamples)[2],
  vector<TString>(&strExtraWeights)[2],
  unsigned int divisor,
  const bool writeFinalTree,
  TString strcustomselection="",
  vector<pair<vector<float>, pair<float, float>>>* manualboundary_validity_pairs=0
  ){
  cout << "Begin getKDConstantByMass" << endl;

  Discriminant* KDbuilder = constructKDFromType(strname);
  if (!KDbuilder) return;

  const TString strTrackVar = "ZZMass";
  vector<TString> strWeights; strWeights.push_back(TString("KFactor_QCD_ggZZ_Nominal")); strWeights.push_back(TString("KFactor_EW_qqZZ")); strWeights.push_back(TString("KFactor_QCD_qqZZ_M"));

  vector<TString> strAllWeights[2];
  for (unsigned int ih=0; ih<2; ih++){
    for (auto& s : strWeights) strAllWeights[ih].push_back(s);
    for (auto& s : strExtraWeights[ih]) strAllWeights[ih].push_back(s);
  }

  vector<TString> strExtraShortVars; strExtraShortVars.push_back("Z1Flav"); strExtraShortVars.push_back("Z2Flav");

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
    CJLSTSet* theSet = bookSampleTrees(
      strSamples[ih],

      strTrackVar,
      strAllWeights[ih],
      strRecoBranch,

      strExtraRecoBranches,
      strExtraShortVars
      );
    getEvents(
      theSet,

      strTrackVar,
      strAllWeights[ih],
      strRecoBranch,

      index[ih],
      KDbuilder,
      strcustomselection
      );

    delete theSet;

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

  TH1F* h_varTrack_Constant = new TH1F("varReco_Constant", "", nbins, binning); h_varTrack_Constant->Sumw2();
  TProfile* p_varTrack = new TProfile("avg_varReco", "", nbins, binning); p_varTrack->Sumw2();
  delete[] binning;

  LoopForConstant(
    index, indexboundaries,
    KDbuilder,
    p_varTrack,
    h_varTrack_Constant,
    100
    );

  TGraphErrors* gr = makeGraphFromTH1(p_varTrack, h_varTrack_Constant, "gr_varReco_Constant");
  foutput->WriteTObject(p_varTrack);
  foutput->WriteTObject(h_varTrack_Constant);
  foutput->WriteTObject(gr);
  delete gr;
  delete h_varTrack_Constant;
  delete p_varTrack;
  foutput->Close();

  delete KDbuilder;

  cout << "End getKDConstantByMass" << endl;
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
  if (strprod=="ZHG"){
    strRecoBranch.push_back("p_HadZH_SIG_ghz1_1_JHUGen_JECNominal");
    strRecoBranch.push_back("p_HadZH_mavjj_JECNominal");
  }
  else if (strprod=="WH"){
    strRecoBranch.push_back("p_HadWH_SIG_ghw1_1_JHUGen_JECNominal");
    strRecoBranch.push_back("p_HadWH_mavjj_JECNominal");
  }
  strRecoBranch.push_back("p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal");
  vector<TString> strExtraRecoBranches;

  vector<TString> strSamples[2];
  if (strprod=="ZH"){
    vector<TString> s1; s1.push_back("ZH_Sig_POWHEG");
    vector<TString> s2; s2.push_back("gg_Sig_POWHEG");
    getSamplePairs(sqrts, s1, s2, strSamples[0], strSamples[1]);
  }
  else if (strprod=="WH"){
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
  if (strprod=="ZH"){
    pair<float, float> valrange(230, 3500);
    vector<float> manualboundaries;
    manualboundaries.push_back(245);
    manualboundaries.push_back(300);
    manualboundaries.push_back(500);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(
      manualboundaries, valrange
      ));
  }
  else if (strprod=="WH"){
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
    strRecoBranch, strExtraRecoBranches, strSamples, extraweights,
    divisor, writeFinalTree, "",
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
  vector<TString> strExtraRecoBranches;

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
    strRecoBranch, strExtraRecoBranches, strSamples, extraweights,
    divisor, writeFinalTree, "",
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
  vector<TString> strExtraRecoBranches;

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
    strRecoBranch, strExtraRecoBranches, strSamples, extraweights,
    divisor, writeFinalTree, "",
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
  vector<TString> strExtraRecoBranches;

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
    manualboundaries.push_back(120); manualboundaries.push_back(123); manualboundaries.push_back(135);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(manualboundaries, valrange));
  }
  {
    pair<float, float> valrange(600, 2500);
    vector<float> manualboundaries;
    manualboundaries.push_back(700); manualboundaries.push_back(900); manualboundaries.push_back(1100); manualboundaries.push_back(1400); manualboundaries.push_back(1900);
    manualboundary_validity_pairs.push_back(pair<vector<float>, pair<float, float>>(manualboundaries, valrange));
  }

  getKDConstantByMass(
    sqrts, "Dbkgkin",
    strRecoBranch, strExtraRecoBranches, strSamples, extraweights,
    divisor, writeFinalTree,
    strchannel,
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

  TGraphErrors* tg = (TGraphErrors*)finput->Get("gr_varReco_Constant");
  foutput->WriteTObject(tg);

  if (addpoints!=0){ for (auto& prange : *addpoints) addPointsBetween(tg, prange.first.first, prange.first.second, prange.second); }

  int n = tg->GetN();
  double* xx = tg->GetX();
  double* ex = tg->GetEX();
  double* yy = tg->GetY();
  double* ey = tg->GetEY();

  TSpline3* sp;
  sp = convertGraphToSpline3(tg, useFaithfulSlopeFirst, useFaithfulSlopeSecond);
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
    &getFcn_a0plusa1timesXN<1>,
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
    &getFcn_a0plusa1timesXN<1>,
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
    &getFcn_a0plusa1timesXN<1>,
    &getFcn_a0timesexpa1X,
    true, true
    );
}

void SmoothKDConstantProducer_Dbkgkin(TString strchannel){
  if (strchannel!="2e2mu" && strchannel!="4e" && strchannel!="4mu") return;
  generic_SmoothKDConstantProducer(
    13, Form("Dbkgkin_%s", strchannel.Data()), "",
    &getFcn_a0plusa1timesXN<1>,
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

