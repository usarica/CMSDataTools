#include "CJLSTSet.h"


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
  if (tree==nullptr) std::cerr << "CJLSTSet::addCJLSTTree(" << strname << ") is invalid!" << std::endl;
  else std::cout << "CJLSTSet::addCJLSTTree(" << strname << ") is successful!" << std::endl;
  if (tree) tree->setAssociatedSet(this);
  return (tree!=nullptr);
}
bool CJLSTSet::addCJLSTTreeList(const std::vector<TString>& strlist){
  bool res=true;
  for (auto const& s:strlist) res &= addCJLSTTree(s);
  return res;
}


void CJLSTSet::bookXS(){
  for (auto& tree:treeList) tree->bookBranch<float>("xsec", 1);
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

CJLSTTree* CJLSTSet::getCJLSTTree(TString sampleid){
  for (auto tree:treeList){ if (tree->sampleIdentifier==sampleid) return tree; }
  CJLSTTree* res=0;
  return res;
}
const std::vector<CJLSTTree*>& CJLSTSet::getCJLSTTreeList() const{ return treeList; }
std::vector<CJLSTTree*>& CJLSTSet::getCJLSTTreeList(){ return treeList; }

float CJLSTSet::getOverallEventWgt(CJLSTTree* sample){
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
float CJLSTSet::getOverallEventWgt(TString sampleid){ return getOverallEventWgt(getCJLSTTree(sampleid)); }

float CJLSTSet::getPermanentWeight(CJLSTTree* sample){
  if (permanentWeights.empty() || permanentWeights.find(sample)==permanentWeights.end()) return 1;
  return permanentWeights[sample];
}
float CJLSTSet::getPermanentWeight(TString sampleid){ return getPermanentWeight(getCJLSTTree(sampleid)); }
void CJLSTSet::setPermanentWeights(bool useXS, bool useNormPerMass, bool useNgen, bool renormalizeWeights){
  cout << "CJLSTSet::setPermanentWeights is called." << endl;

  vector<CJLSTTree*> extraBookXS;
  if (useXS){
    for (auto& tree:treeList){ if (!tree->branchExists("xsec")) extraBookXS.push_back(tree); }
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
    for (auto& tree:mg.second){
      float xsec=1;
      if (useXS){
        if (!tree->getSelectedEvent(0)) xsec=0; // Kill contribution to sum of weights from this tree
        tree->getVal("xsec", xsec);
      }
      if (xsec<=0.){ cerr << "- XSec=" << xsec << " is invalid for sample " << tree->sampleIdentifier << ". Setting to 1" << endl; xsec=1; }

      float ngen=1;
      if (useNgen) ngen = tree->getNGenWithPU();
      if (ngen==0.) ngen=1;

      permanentWeights[tree] = xsec/ngen;

      cout << "- " << tree->sampleIdentifier << " specifics: "
        << "XSec: " << xsec << ", "
        << "Ngen: " << ngen << ", "
        << (renormalizeWeights ? "permanent weight before renormalization: " : "permanent weight: ") << permanentWeights[tree]
        << endl;

      sumwgt += permanentWeights[tree];
    }
    if (renormalizeWeights && sumwgt!=0.){
      cout << "- Renormalizing the mass(" << mg.first << ") grouping by 1./" << sumwgt << endl;
      for (auto& tree:mg.second) permanentWeights[tree] = permanentWeights[tree]/sumwgt;
    }
  }

  for (auto& tree:extraBookXS) tree->releaseBranch("xsec");
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


