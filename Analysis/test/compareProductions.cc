#include "common_includes.h"


using namespace std;

void compareProductions(TString strOldProd, TString strNewProd, TString strSample, int maxmatch=-1){
  TString strinput[2];
  strinput[0] = "root://lxcms03.cern.ch//data3/Higgs/" + strOldProd + "/" + strSample + "/ZZ4lAnalysis.root";
  strinput[1] = "root://lxcms03.cern.ch//data3/Higgs/" + strNewProd + "/" + strSample + "/ZZ4lAnalysis.root";
  TFile* finput[2]={ nullptr };
  TTree* tree[2]={ nullptr };
  unordered_map<TString, pair<short, short>> sbranches;
  unordered_map<TString, pair<float, float>> fbranches;
  unordered_map<TString, pair<vector<short>*, vector<short>*>> vsrbranches;
  unordered_map<TString, pair<vector<short>, vector<short>>> vsbranches;
  unordered_map<TString, pair<vector<float>*, vector<float>*>> vfrbranches;
  unordered_map<TString, pair<vector<float>, vector<float>>> vfbranches;
  sbranches["genFinalState"]=pair<short, short>(0, 0);
  sbranches["nCleanedJetsPt30"]=pair<short, short>(0, 0);
  sbranches["Z1Flav"]=pair<short, short>(0, 0);
  sbranches["Z2Flav"]=pair<short, short>(0, 0);
  fbranches["GenHMass"]=pair<float, float>(0, 0);
  fbranches["GenZ1Mass"]=pair<float, float>(0, 0);
  fbranches["GenZ2Mass"]=pair<float, float>(0, 0);
  fbranches["ZZMass"]=pair<float, float>(0, 0);
  fbranches["Z1Mass"]=pair<float, float>(0, 0);
  fbranches["Z2Mass"]=pair<float, float>(0, 0);
  fbranches["DiJetMass"]=pair<float, float>(0, 0);
  fbranches["p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal"]=pair<float, float>(0, 0);
  fbranches["p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal"]=pair<float, float>(0, 0);
  fbranches["pConst_JJVBF_SIG_ghv1_1_JHUGen_JECNominal"]=pair<float, float>(0, 0);
  fbranches["pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal"]=pair<float, float>(0, 0);
  fbranches["p_JJVBF_SIG_ghv1_1_JHUGen_JECUp"]=pair<float, float>(0, 0);
  fbranches["p_JJQCD_SIG_ghg2_1_JHUGen_JECUp"]=pair<float, float>(0, 0);
  fbranches["pConst_JJVBF_SIG_ghv1_1_JHUGen_JECUp"]=pair<float, float>(0, 0);
  fbranches["pConst_JJQCD_SIG_ghg2_1_JHUGen_JECUp"]=pair<float, float>(0, 0);
  fbranches["p_JJVBF_SIG_ghv1_1_JHUGen_JECDn"]=pair<float, float>(0, 0);
  fbranches["p_JJQCD_SIG_ghg2_1_JHUGen_JECDn"]=pair<float, float>(0, 0);
  fbranches["pConst_JJVBF_SIG_ghv1_1_JHUGen_JECDn"]=pair<float, float>(0, 0);
  fbranches["pConst_JJQCD_SIG_ghg2_1_JHUGen_JECDn"]=pair<float, float>(0, 0);
  vfrbranches["LepPt"]=pair<vector<float>*, vector<float>*>(0, 0);
  vfrbranches["LepEta"]=pair<vector<float>*, vector<float>*>(0, 0);
  vfrbranches["LepPhi"]=pair<vector<float>*, vector<float>*>(0, 0);
  vsrbranches["LepLepId"]=pair<vector<short>*, vector<short>*>(0, 0);
  vfrbranches["JetPt"]=pair<vector<float>*, vector<float>*>(0, 0);
  vfrbranches["JetEta"]=pair<vector<float>*, vector<float>*>(0, 0);
  vfrbranches["JetPhi"]=pair<vector<float>*, vector<float>*>(0, 0);
  vfrbranches["JetMass"]=pair<vector<float>*, vector<float>*>(0, 0);
  for (auto it=vsrbranches.begin(); it!=vsrbranches.end(); it++) vsbranches[it->first]=pair<vector<short>, vector<short>>(vector<short>(), vector<short>());
  for (auto it=vfrbranches.begin(); it!=vfrbranches.end(); it++) vfbranches[it->first]=pair<vector<float>, vector<float>>(vector<float>(), vector<float>());
  for (unsigned int f=0; f<2; f++){
    finput[f] = TFile::Open(strinput[f], "read");
    cout << "File[" << f << "] = " << finput[f]->GetName() << " open!" << endl;
    tree[f] = (TTree*) finput[f]->Get("ZZTree/candTree");
    cout << "Tree[" << f << "] = " << tree[f]->GetName() << " taken! Nevents: " << tree[f]->GetEntries() << endl;
    tree[f]->SetBranchStatus("*", 0);
    for (auto it=sbranches.begin(); it!=sbranches.end(); it++) bookBranch(tree[f], it->first, (f==0 ? &(it->second.first) : &(it->second.second)));
    for (auto it=fbranches.begin(); it!=fbranches.end(); it++) bookBranch(tree[f], it->first, (f==0 ? &(it->second.first) : &(it->second.second)));
    for (auto it=vsrbranches.begin(); it!=vsrbranches.end(); it++) bookBranch(tree[f], it->first, (f==0 ? &(it->second.first) : &(it->second.second)));
    for (auto it=vfrbranches.begin(); it!=vfrbranches.end(); it++) bookBranch(tree[f], it->first, (f==0 ? &(it->second.first) : &(it->second.second)));
  }
  
  TFile* foutput = TFile::Open("comparison.root", "recreate");
  TTree* compTree = new TTree("compTree", "");
  for (unsigned int f=0; f<2; f++){
    TString strOldNew = (f==0 ? "Old" : "New");
    for (auto it=sbranches.begin(); it!=sbranches.end(); it++) compTree->Branch(Form("%s%s", it->first.Data(), strOldNew.Data()), (f==0 ? &(it->second.first) : &(it->second.second)));
    for (auto it=fbranches.begin(); it!=fbranches.end(); it++) compTree->Branch(Form("%s%s", it->first.Data(), strOldNew.Data()), (f==0 ? &(it->second.first) : &(it->second.second)));
    for (auto it=vsbranches.begin(); it!=vsbranches.end(); it++) compTree->Branch(Form("%s%s", it->first.Data(), strOldNew.Data()), (f==0 ? &(it->second.first) : &(it->second.second)));
    for (auto it=vfbranches.begin(); it!=vfbranches.end(); it++) compTree->Branch(Form("%s%s", it->first.Data(), strOldNew.Data()), (f==0 ? &(it->second.first) : &(it->second.second)));
  }
  int nmatched=0;
  int croffset=0;
  for (int ev=0; ev<tree[1]->GetEntries(); ev++){
    if (ev%1000==0){
      cout << "Trying to match event " << ev << " in new production..." << endl;
      cout << "So far matched: " << nmatched << endl;
    }
    tree[1]->GetEntry(ev);
    bool hasMatch=false;
    if (ev+croffset<tree[0]->GetEntries()){
      tree[0]->GetEntry(ev+croffset);
      hasMatch = (
        fabs(fbranches["GenHMass"].second/fbranches["GenHMass"].first-1.)<1e-5
        &&
        fabs(fbranches["GenZ1Mass"].second/fbranches["GenZ1Mass"].first-1.)<1e-5
        &&
        fabs(fbranches["GenZ2Mass"].second/fbranches["GenZ2Mass"].first-1.)<1e-5
        );
    }
    if (!hasMatch){
      for (int cr=0; cr<tree[0]->GetEntries(); cr++){
        tree[0]->GetEntry(cr);
        if (
          fabs(fbranches["GenHMass"].second/fbranches["GenHMass"].first-1.)<1e-5
          &&
          fabs(fbranches["GenZ1Mass"].second/fbranches["GenZ1Mass"].first-1.)<1e-5
          &&
          fabs(fbranches["GenZ2Mass"].second/fbranches["GenZ2Mass"].first-1.)<1e-5
          ){ // Found match
          hasMatch=true;
          croffset = cr-ev;
          break;
        }
      }
    }
    if (hasMatch){
      for (auto it=vsbranches.begin(); it!=vsbranches.end(); it++){
        it->second.first = *(vsrbranches[it->first].first);
        it->second.second = *(vsrbranches[it->first].second);
      }
      for (auto it=vfbranches.begin(); it!=vfbranches.end(); it++){
        it->second.first = *(vfrbranches[it->first].first);
        it->second.second = *(vfrbranches[it->first].second);
      }
      compTree->Fill();
      nmatched++;
    }
    else cout << "Event " << ev << " does not have a match!" << endl;
    if (maxmatch>0 && nmatched==maxmatch){
      cout << "Breaking at event " << ev << " since maxmatch=" << maxmatch << " is reached." << endl;
      break;
    }
  }
  foutput->WriteTObject(compTree);
  delete compTree;
  foutput->Close();
  for (unsigned int f=0; f<2; f++) finput[f]->Close();
}

void computeMELA(){
  TFile* finput = TFile::Open("comparison.root", "read");
  TTree* tree = (TTree*) finput->Get("compTree");

  float ZZMass, p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal, pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
  vector<float>* LepPt=nullptr;
  vector<float>* LepEta=nullptr;
  vector<float>* LepPhi=nullptr;
  vector<short>* LepLepId=nullptr;
  vector<float>* JetPt=nullptr;
  vector<float>* JetEta=nullptr;
  vector<float>* JetPhi=nullptr;
  vector<float>* JetMass=nullptr;

  bookBranch(tree, "ZZMassNew", &ZZMass);
  bookBranch(tree, "p_JJQCD_SIG_ghg2_1_JHUGen_JECNominalNew", &p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);
  bookBranch(tree, "pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominalNew", &pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal);

  bookBranch(tree, "LepPtNew", &LepPt);
  bookBranch(tree, "LepEtaNew", &LepEta);
  bookBranch(tree, "LepPhiNew", &LepPhi);
  bookBranch(tree, "LepLepIdNew", &LepLepId);

  bookBranch(tree, "JetPtNew", &JetPt);
  bookBranch(tree, "JetEtaNew", &JetEta);
  bookBranch(tree, "JetPhiNew", &JetPhi);
  bookBranch(tree, "JetMassNew", &JetMass);

  Mela mela(13, 125, TVar::DEBUG);

  vector<MELAParticle*> particleList;
  vector<MELACandidate*> candList;
  MELACandidate* cands[3]={ nullptr };
  float mevals[3][2]={ { 0 } };
  float constvals[3][2]={ { 0 } };

  for (int ev=0; ev<tree->GetEntries(); ev++){
    tree->GetEntry(ev);

    if (cands[0] && cands[1] && cands[2]) break;
    if (ZZMass<120 || ZZMass>130) continue;
    int icand = -1
      + (pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal>0. && pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal<2.)*1
      + (pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal>2. && pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal<10.)*2
      + (pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal>10. && pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal<18.)*3;
    if (icand<0 || cands[icand]) continue;
    TLorentzVector jet1(0, 0, 1e-3, 1e-3), jet2(0, 0, 1e-3, 1e-3);
    jet1.SetPtEtaPhiM(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetMass->at(0));
    jet2.SetPtEtaPhiM(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetMass->at(1));

    SimpleParticleCollection_t associated;
    associated.push_back(SimpleParticle_t(0, jet1));
    associated.push_back(SimpleParticle_t(0, jet2));

    SimpleParticleCollection_t daughters;
    for (unsigned int idau=0; idau<4; idau++){
      TLorentzVector ptmp; ptmp.SetPtEtaPhiM(LepPt->at(idau), LepEta->at(idau), LepPhi->at(idau), 0.);
      daughters.push_back(SimpleParticle_t(int(LepLepId->at(idau)), ptmp));
    }

    MELACandidate* cand = TUtil::ConvertVectorFormat(
      &daughters,
      &associated,
      nullptr,
      false,
      &particleList, &candList
    );
    mela.setCurrentCandidate(cand);
    cands[icand]=cand;
    mela.setProcess(TVar::HSMHiggs, TVar::JHUGen, TVar::JJQCD);
    float p_hjj_VAJHU, pconst_hjj_VAJHU;
    mela.computeProdP(p_hjj_VAJHU, true);
    mela.getConstant(pconst_hjj_VAJHU);
    mevals[icand][0]=p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
    mevals[icand][1]=p_hjj_VAJHU;
    constvals[icand][0]=pConst_JJQCD_SIG_ghg2_1_JHUGen_JECNominal;
    constvals[icand][1]=pconst_hjj_VAJHU;

    MELAout << "Candidate " << icand << ": " << endl;
    MELAout << cands[icand] << endl;
    MELAout << "Tree me / const = " << mevals[icand][0] << " / " << constvals[icand][0] << " = " << mevals[icand][0]/constvals[icand][0] << endl;
    MELAout << "Computed me / const = " << mevals[icand][1] << " / " << constvals[icand][1] << " = " << mevals[icand][1]/constvals[icand][1] << endl;
  }

  for (auto& cand:candList) delete cand;
  for (auto& part:particleList) delete part;
  finput->Close();
}