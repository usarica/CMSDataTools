#include <sys/types.h>
#include <dirent.h>
#include "SampleHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


namespace SampleHelpers{
  shared_ptr<Mela> GlobalMELA;
}

void SampleHelpers::makeGlobalMELA(int CoM, TVar::VerbosityLevel verbosity){ if (!GlobalMELA) GlobalMELA.reset(new Mela(CoM, 125, verbosity)); }

std::vector<TString> SampleHelpers::lsdir(TString indir){
  std::vector<TString> res;

  struct dirent* ep;
  DIR* dp = opendir(indir.Data());
  if (dp != NULL){
    while ((ep = readdir(dp))) res.push_back(ep->d_name);
    closedir(dp);
  }
  else MELAerr << "Couldn't open the directory" << endl;

  return res;
}


float SampleHelpers::findPoleMass(const TString samplename){
  float mass = -1;
  if (samplename=="") return mass;
  else if (samplename.Contains("_M125")) return 125; // JHUGen samples
  else if (samplename.Contains("_M125p6")) return 125.6; // JHUGen samples
  else if (samplename.Contains("_M126")) return 126; // JHUGen samples
  std::string strtmp = samplename.Data();
  std::size_t extpos = strtmp.find(".root");
  if (extpos!=std::string::npos) strtmp.erase(extpos, 5);
  std::vector<std::string> strsplit;
  HelperFunctions::splitOptionRecursive(strtmp, strsplit, 'H');
  if (strsplit.size()>1){
    std::string strmass = strsplit.at(1);
    if (strmass=="f05ph0" || strmass=="f05ph90") strmass = strsplit.at(2);
    strsplit.clear();
    HelperFunctions::splitOptionRecursive(strmass, strsplit, '_');
    strmass = strsplit.at(0);
    mass = std::stod(strmass);
  }
  return mass;
}
TTree* SampleHelpers::findTree(std::vector<TTree*> const& treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t);
    int nevts = tree->GetEntries();
    if (ev<nevts) return tree;
    else ev -= nevts;
    if (ev<0) cerr << "findTree::ERROR: Could not find the event " << evid << endl;
  }
  return 0;
}
bool SampleHelpers::branchExists(TTree* tree, TString strname){
  if (!tree) return false;
  bool found=false;
  const TList* blist = (const TList*) tree->GetListOfBranches();
  for (int ib=0; ib<blist->GetSize(); ib++){
    TObject* branch = blist->At(ib);
    if (!branch) continue;
    TString bname = branch->GetName();
    if (strname==bname){ found=true; break; }
  }
  return found;
}
bool SampleHelpers::aliasExists(TTree* tree, TString strname){
  if (!tree) return false;
  bool found=false;
  const TList* blist = (const TList*) tree->GetListOfAliases();
  for (int ib=0; ib<blist->GetSize(); ib++){
    TObject* branch = blist->At(ib);
    if (!branch) continue;
    TString bname = branch->GetName();
    if (strname==bname){ found=true; break; }
  }
  return found;
}
void SampleHelpers::getEntry(std::vector<TTree*>& treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t);
    int nevts = tree->GetEntries();
    if (ev<nevts){ tree->GetEntry(ev); break; }
    else ev -= nevts;
    if (ev<0) cerr << "getEntry::ERROR: Could not find the event " << evid << endl;
  }
}
float SampleHelpers::getEntry(std::vector<std::pair<TTree*, TH1F*>>& treeList, int evid){
  int ev = evid;
  for (unsigned int t=0; t<treeList.size(); t++){
    TTree* tree = treeList.at(t).first;
    int nevts = tree->GetEntries();
    if (ev<nevts){ tree->GetEntry(ev); return float(1./treeList.at(t).second->GetBinContent(40)); }
    else ev -= nevts;
    if (ev<0) cerr << "getEntry::ERROR: Could not find the event " << evid << endl;
  }
  return 0;
}
