#include "CJLSTSet.h"
#include <algorithm>


using namespace std;


CJLSTSet::CJLSTSet(const TString& strname){
  addCJLSTTree(strname);
}
CJLSTSet::CJLSTSet(const std::vector<TString>& strlist){
  addCJLSTTreeList(strlist);
}
CJLSTSet::~CJLSTSet(){
  for (auto& tree:treeList) delete tree;
  treeList.clear();
}
bool CJLSTSet::addCJLSTTree(const TString& strname){
  CJLSTTree* tree = new CJLSTTree(strname);
  if (tree->isValid()) treeList.push_back(tree);
  else{ delete tree; tree=nullptr; }
  if (!tree) cerr << "CJLSTSet::addCJLSTTree(" << strname << ") is invalid!" << endl;
  else cout << "CJLSTSet::addCJLSTTree(" << strname << ") is successful!" << endl;
  if (tree) tree->setAssociatedSet(this);
  return (tree!=nullptr);
}
bool CJLSTSet::addCJLSTTreeList(const std::vector<TString>& strlist){
  bool res=true;
  for (auto const& s:strlist) res &= addCJLSTTree(s);
  return res;
}
bool CJLSTSet::dissociateCJLSTTree(CJLSTTree*& tree){
  if (!tree) return false;
  auto it = std::find(treeList.begin(), treeList.end(), tree);
  if (it!=treeList.end()){ treeList.erase(it); tree->setAssociatedSet(nullptr); return true; }
  return false;
}
bool CJLSTSet::associateCJLSTTree(CJLSTTree*& tree){
  if (!tree) return false;
  CJLSTSet* theSet = tree->getAssociatedSet();
  if (theSet==this) return true;
  else if (theSet){
    theSet->dissociateCJLSTTree(tree);
    treeList.push_back(tree);
    tree->setAssociatedSet(this);
    return true;
  }
  return false;
}


void CJLSTSet::bookXS(){
  std::vector<TString> xsecvars=SampleHelpers::getXsecBranchNames();
  for (auto& tree:treeList){ for (auto& v:xsecvars) tree->bookBranch<float>(v, 1); }
}
void CJLSTSet::bookOverallEventWgt(){
  for (auto& tree:treeList){
    // Get these individually instead of using the overall evt. wgt. directly
    // so that failed tree also gets a meaningful weight
    // overallEventWeight = PUWeight * genHEPMCweight * dataMCWeight * trigEffWeight
    tree->bookBranch<float>("PUWeight", 1);
    tree->bookBranch<float>("genHEPMCweight", 1);
    tree->bookBranch<float>("dataMCWeight", 1);
    tree->bookBranch<float>("trigEffWeight", 1);
  }
}

CJLSTTree* CJLSTSet::getCJLSTTree(TString sampleid) const{
  for (auto tree:treeList){ if (tree->sampleIdentifier==sampleid) return tree; }
  CJLSTTree* res=nullptr;
  return res;
}
const std::vector<CJLSTTree*>& CJLSTSet::getCJLSTTreeList() const{ return treeList; }
std::vector<CJLSTTree*>& CJLSTSet::getCJLSTTreeList(){ return treeList; }

float CJLSTSet::getOverallEventWgt(CJLSTTree* const sample) const{
  if (!sample) return 0;
  float PUWeight=1;
  float genHEPMCweight=1;
  float dataMCWeight=1;
  float trigEffWeight=1;
  sample->getVal("PUWeight", PUWeight);
  sample->getVal("genHEPMCweight", genHEPMCweight);
  sample->getVal("dataMCWeight", dataMCWeight);
  sample->getVal("trigEffWeight", trigEffWeight);
  return (PUWeight*genHEPMCweight*dataMCWeight*trigEffWeight);
}
float CJLSTSet::getOverallEventWgt(TString sampleid) const{ return getOverallEventWgt(getCJLSTTree(sampleid)); }

float CJLSTSet::getPermanentWeight(CJLSTTree* const sample) const{
  if (!sample) return 0;
  auto it = permanentWeights.find(sample);
  if (it==permanentWeights.cend()) return 1;
  return it->second;
}
float CJLSTSet::getPermanentWeight(TString sampleid) const{ return getPermanentWeight(getCJLSTTree(sampleid)); }
void CJLSTSet::setPermanentWeights(const CJLSTSet::NormScheme scheme, const bool useNormPerMass, const bool useNgenWPU){
  cout << "CJLSTSet::setPermanentWeights(" << scheme << "," << useNormPerMass << "," << useNgenWPU << ") is called." << endl;

  if (scheme==NormScheme_NgenOverNgenWPU && (!useNgenWPU || useNormPerMass)) return; // permanentWeight=1 would be the case

  const bool renormalizeWeights = !( // !!!!!
    scheme==NormScheme_None
    ||
    scheme==NormScheme_OneOverNgen
    ||
    scheme==NormScheme_NgenOverNgenWPU
    ||
    scheme==NormScheme_XsecOnly
    ||
    scheme==NormScheme_XsecOverNgen
    );
  const bool useXS = (
    scheme==NormScheme_XsecOnly
    ||
    scheme==NormScheme_XsecOverNgen
    ||
    scheme==NormScheme_XsecOverNgen_RenormBySumXsecOverNgen
    ||
    scheme==NormScheme_XsecOverNgen_RelRenormToSumNgen
    );
  const bool useNgen = (
    scheme==NormScheme_OneOverNgen
    ||
    scheme==NormScheme_OneOverNgen_RenormBySumOneOverNgen
    ||
    scheme==NormScheme_OneOverNgen_RelRenormToSumNgen
    ||
    scheme==NormScheme_NgenOverNgenWPU
    ||
    scheme==NormScheme_XsecOverNgen
    ||
    scheme==NormScheme_XsecOverNgen_RenormBySumXsecOverNgen
    ||
    scheme==NormScheme_XsecOverNgen_RelRenormToSumNgen
    );
  cout << "- CJLSTSet::setPermanentWeights: (renormalizeWeights,useXS,useNgen) = (" << renormalizeWeights << "," << useXS << "," << useNgen << ")" << endl;

  std::vector<TString> xsecvars=SampleHelpers::getXsecBranchNames();
  vector<CJLSTTree*> extraBookXS;
  if (useXS){
    for (auto& tree:treeList){
      for (auto& v:xsecvars){ if (!tree->branchExists(v)){ extraBookXS.push_back(tree); break; } }
    }
    bookXS();
  }

  vector<pair<float, vector<CJLSTTree*>>> massgrouping;
  float mh=-1;
  if (useNormPerMass){ // Group everything by mass
    for (auto& tree:treeList){
      mh = tree->MHVal;
      cout << "- MH(" << tree->sampleIdentifier << ") = " << mh << endl;
      int whichgroup=-1;
      for (unsigned int ig=0; ig<massgrouping.size(); ig++){ // Need to keep track of group number
        if (massgrouping.at(ig).first==mh){ whichgroup=int(ig); break; }
      }
      if (whichgroup>=0) massgrouping.at(whichgroup).second.push_back(tree);
      else{ massgrouping.push_back(pair<float, vector<CJLSTTree*>>(mh, vector<CJLSTTree*>())); massgrouping.back().second.push_back(tree); }
    }
  }
  else{ // Everything goes into the mH=-1 (front) group
    massgrouping.push_back(pair<float, vector<CJLSTTree*>>(mh, vector<CJLSTTree*>()));
    for (auto& tree:treeList) massgrouping.front().second.push_back(tree);
  }
  for (auto& mg:massgrouping){
    float sumwgt=0;
    float sumxsec=0;
    float sumngen=0;
    float sumNonZero=0;
    for (auto& tree:mg.second){
      float xsec=1;
      if (useXS){
        if (!(tree->getSelectedEvent(0) || tree->getFailedEvent(0))) xsec=0; // Kill contribution to sum of weights from this tree
        else{ for (auto& v:xsecvars){ float xv; tree->getVal(v, xv); xsec *= xv; } }
      }
      if (xsec<=0.){ cerr << "- XSec=" << xsec << " is invalid for sample " << tree->sampleIdentifier << ". Setting to 1" << endl; xsec=1; }

      float ngen_nopuhep=tree->getNGenNoPU();
      float ngen_wpuhep=tree->getNGenWithPU();
      float ngen=1;
      if (useNgen){
        if (!useNgenWPU) ngen = ngen_nopuhep;
        else ngen = ngen_wpuhep;
        if (scheme==NormScheme_NgenOverNgenWPU) ngen /= ngen_nopuhep;
      }
      if (ngen==0.) ngen=1;

      permanentWeights[tree] = xsec/ngen;

      cout << "- " << tree->sampleIdentifier << " specifics: "
        << "XSec: " << xsec << ", "
        << "Ngen: " << ngen << ", "
        << (renormalizeWeights ? "permanent weight before renormalization: " : "permanent weight: ") << permanentWeights[tree]
        << endl;

      sumwgt += permanentWeights[tree];
      if (permanentWeights[tree]!=0.){
        sumxsec += xsec;
        sumngen += ngen_nopuhep;
        sumNonZero += 1;
      }
    }
    if (renormalizeWeights){
      float renormMult = 0;
      if (scheme==NormScheme_OneOverNgen_RelRenormToSumNgen && sumNonZero!=0.) renormMult = sumngen/sumNonZero;
      else if (scheme==NormScheme_XsecOverNgen_RelRenormToSumNgen && sumxsec!=0.) renormMult = sumngen/sumxsec;
      else if (sumwgt!=0.) renormMult = 1./sumwgt;

      cout << "- Renormalizing the mass(" << mg.first << ") grouping by x" << renormMult << endl;
      for (auto& tree:mg.second){
        permanentWeights[tree] = permanentWeights[tree]*renormMult;
      }
    }
  }

  for (auto& tree:extraBookXS){ for (auto& v:xsecvars) tree->releaseBranch(v); }
}

CJLSTTree* CJLSTSet::getSelectedEvent(const int evid){
  int ev = evid;
  CJLSTTree* isfound=nullptr;
  for (auto& tree:treeList){
    int nevts = tree->getSelectedNEvents();
    if (ev<nevts){ tree->getSelectedEvent(ev); isfound = tree; break; }
    else ev -= nevts;
  }
  return isfound;
}
CJLSTTree* CJLSTSet::getFailedEvent(const int evid){
  int ev = evid;
  CJLSTTree* isfound=nullptr;
  for (auto& tree:treeList){
    int nevts = tree->getFailedNEvents();
    if (ev<nevts){ tree->getFailedEvent(ev); isfound = tree; break; }
    else ev -= nevts;
  }
  return isfound;
}


