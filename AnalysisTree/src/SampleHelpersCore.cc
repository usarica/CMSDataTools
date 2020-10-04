#include <sys/types.h>
#include <dirent.h>
#include <regex>
#include "SampleHelpersCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


namespace SampleHelpers{
  shared_ptr<Mela> GlobalMELA;
  TDirectory* const rootTDirectory = gDirectory;
  volatile sig_atomic_t doSignalInterrupt = 0;
}

void SampleHelpers::makeGlobalMELA(int CoM, TVar::VerbosityLevel verbosity){ if (!GlobalMELA) GlobalMELA.reset(new Mela(CoM, 125, verbosity)); }

void SampleHelpers::setSignalInterrupt(int snum){
  if (snum==SIGINT) doSignalInterrupt = 1;
}

std::vector<TString> SampleHelpers::lsdir(TString indir, HostHelpers::Hosts const* target_host){
  std::vector<TString> res;

  if (!target_host || *target_host == HostHelpers::GetHostLocation()){
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
  }
  else{
    TString localredirector = HostHelpers::GetHostLocalRedirector(*target_host, true);
    TString pathToStore = HostHelpers::GetHostPathToStore(*target_host);
    if (localredirector=="") MELAerr << "SampleHelpers::lsdir: Could not locate the site local director to run ls." << endl;
    else{
      TString tmpfname = Form("tmpfile_SampleHelpers_lsdir_%s_%i.tmp", indir.Data(), *target_host);
      while (tmpfname.Contains("/")) HelperFunctions::replaceString(tmpfname, "/", "_");

      TString indir_eff = indir;
      if (pathToStore!="" && indir_eff.Contains(pathToStore)) HelperFunctions::replaceString<TString, const TString>(indir_eff, pathToStore, "");
      while (indir_eff.Contains("//")) HelperFunctions::replaceString(indir_eff, "//", "/");

      TString x509_proxy = "/tmp/x509up_u$(id -u)";
      if (getenv("X509_USER_PROXY")!=nullptr) x509_proxy = getenv("X509_USER_PROXY");

      TString strcmd;
      if (*target_host == HostHelpers::kUCSDT2) strcmd = Form("env -i X509_USER_PROXY=%s gfal-ls root://%s/%s", x509_proxy.Data(), localredirector.Data(), indir_eff.Data());
      else strcmd = Form("env -i X509_USER_PROXY=%s xrdfs %s ls %s", x509_proxy.Data(), localredirector.Data(), indir_eff.Data());
      strcmd = strcmd + " > " + tmpfname;
      int status_cmd = HostHelpers::ExecuteCommand(strcmd);
      if (status_cmd==0){
        TString strRemove = indir_eff;
        if (!strRemove.EndsWith("/")) strRemove = strRemove + '/';

        ifstream fin;
        fin.open(tmpfname.Data(), std::ios_base::in);
        if (fin.good()){
          std::string str_in;
          while (!fin.eof()){
            getline(fin, str_in);
            if (str_in!=""){
              HelperFunctions::replaceString<std::string, const char*>(str_in, strRemove.Data(), "");
              while (str_in.find("/")==0) str_in = str_in.substr(1);
              res.push_back(str_in.data());
            }
          }
        }
        fin.close();
        std::remove(tmpfname.Data());
      }
      else{
        MELAerr << "SampleHelpers::lsdir: Command '" << strcmd << "' returned exit status " << status_cmd << "." << endl;
        assert(0);
      }
    }
  }

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
