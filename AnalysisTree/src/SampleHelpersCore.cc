#include <sys/types.h>
#include <dirent.h>
#include <regex>
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
    while ((ep = readdir(dp))){
      TString strdir = ep->d_name;
      if (strdir == "." || strdir == "..") continue;
      res.push_back(strdir);
    }
    closedir(dp);
  }
  else MELAerr << "SampleHelpers::lsdir: Could not open the directory " << indir << "." << endl;

  return res;
}


float SampleHelpers::findPoleMass(const TString samplename){
  float mass = -1;
  if (samplename=="") return mass;

  std::string strSampleName = samplename.Data();
  std::regex strMatch(".*_M([0-9]*p*[0-9]*)_.*");
  std::smatch matches_MXYZ;
  std::regex_match(strSampleName, matches_MXYZ, strMatch);
  if (matches_MXYZ.size()==2){
    std::string strmass = matches_MXYZ[1].str();
    HelperFunctions::replaceString(strmass, "p", ".");
    try{
      mass = std::stof(strmass);
    }
    catch (std::invalid_argument& e){
      MELAerr << "SampleHelpers::findPoleMass: Sample '" << samplename << "' contains the mass string '" << strmass << "', but parsing of this mass string failed!" << endl;
      assert(0);
    }
    return mass;
  }
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
  const TList* blist = nullptr;
  // First search all branches
  blist = (const TList*) tree->GetListOfBranches();
  for (int ib=0; ib<blist->GetSize(); ib++){
    TObject* branch = blist->At(ib);
    if (!branch) continue;
    TString bname = branch->GetName();
    if (strname==bname){ found=true; break; }
  }
  // It is possible that the branch is more composite, so search every leaf in that case.
  if (!found){
    blist = (const TList*) tree->GetListOfLeaves();
    for (int ib=0; ib<blist->GetSize(); ib++){
      TObject* branch = blist->At(ib);
      if (!branch) continue;
      TString bname = branch->GetName();
      if (strname==bname){ found=true; break; }
    }
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
