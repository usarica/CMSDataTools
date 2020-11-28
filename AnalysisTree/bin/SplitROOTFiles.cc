#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "HostHelpersCore.h"
#include "HelperFunctions.h"


using namespace std;
using namespace HostHelpers;
using namespace HelperFunctions;


int main(int argc, char** argv){
  constexpr int iarg_offset=1; // argv[0]==[Executable name]

  constexpr int istatus_unreadable = 1;
  constexpr int istatus_invalid = 2;
  constexpr int istatus_invalid_recovered = 3;
  constexpr int istatus_writefail = 4;

  bool print_help=false, has_help=false;
  string indir="";
  vector<string> filenames;
  unsigned int nchunks=0;

  for (int iarg=iarg_offset; iarg<argc; iarg++){
    string strarg = argv[iarg];
    string wish, value;
    splitOption(strarg, wish, value, '=');

    if (wish.empty()){
      if (value.find(".root")!=string::npos) filenames.push_back(value);
      else if (value=="help"){ print_help=has_help=true; }
    }
    else if (wish=="indir") indir = value;
    else if (wish=="nchunks"){ nchunks = stoi(value); }
    else{
      cerr << "SplitROOTFiles: Unknown argument " << wish << " = " << value << endl;
      print_help=true;
    }
  }
  if (indir!="" && *(indir.rbegin())!='/') indir += '/';
  if (filenames.empty()){
    if (!has_help) cerr << "SplitROOTFiles: No input file is specified." << endl;
    print_help |= true;
  }
  else if (nchunks==0){
    if (!has_help) cerr << "SplitROOTFiles: Number of chunks needs to be greater than 0." << endl;
    print_help |= true;
  }
  else{
    for (auto& s:filenames){
      if (s.find('/')==0){
        if (indir!="/" && indir!=""){
          s = s.substr(1);
          s = indir + s;
        }
      }
      else s = indir + s;
    }
  }
  if (print_help){
    cout << "SplitROOTFiles options:\n\n";
    cout << "- No option specifier: Input files with extension .root. Multiple input files can be passed as different arguments.\n\n";
    cout << "- indir: Location of input files. Default=\"\"\n\n";
    cout << "- nchunks: Number of final chunks for each file. Default=\"\"\n\n";

    cout << endl;
    return 0;
  }

  for (auto const& fname:filenames){
    if (!HostHelpers::FileReadable(fname.data())){
      cerr << "SplitROOTFiles: Input file " << fname << " is not readable!" << endl;
      return istatus_unreadable;
    }

    bool valid = false;
    bool isRecovered = false;

    TDirectory* curdir = gDirectory;

    TFile* finput = TFile::Open(fname.data(), "read");
    if (finput){
      if (finput->IsOpen() && !finput->IsZombie()){
        if (!finput->TestBit(TFile::kRecovered)) valid = true; // Note that the file is not closed.
        else{
          cout << "SplitROOTFiles: Input file " << fname << " went through a recovery, but we would not know if everything is fine. Will assume as a failure." << endl;
          isRecovered = true;
          finput->Close();
        }
      }
      else{
        if (finput->IsOpen()) finput->Close();
        else delete finput;
        finput = nullptr;
      }
    }

    if (!valid){
      cerr << "SplitROOTFiles: Input file " << fname << " is invalid!" << endl;
      if (isRecovered) return istatus_invalid_recovered;
      else return istatus_invalid;
    }

    cout << "SplitROOTFiles: Splitting input file " << fname << "..." << endl;
    std::vector<TFile*> foutputlist; foutputlist.reserve(nchunks);
    for (unsigned int ichunk=0; ichunk<nchunks; ichunk++){
      std::string stroutput = fname;
      std::string strchunk = Form("_chunk_%u_of_%u%s", ichunk, nchunks, ".root");
      HelperFunctions::replaceString<std::string, std::string const>(stroutput, ".root", strchunk);
      TFile* foutput = TFile::Open(stroutput.data(), "recreate");
      if (foutput) foutputlist.push_back(foutput);
      else{
        cerr << "SplitROOTFiles: Failed to create new file " << stroutput << endl;
        for (auto& ff:foutputlist) ff->Close();
        finput->Close();
        return istatus_writefail;
      }
    }

    std::vector<TDirectory*> outputdirs; outputdirs.reserve(foutputlist.size());
    for (auto& ff:foutputlist) outputdirs.push_back(ff);
    HelperFunctions::distributeObjects(finput, outputdirs);

    for (auto& ff:foutputlist) ff->Close();
    finput->Close();

    curdir->cd();
  }

  return 0;
}
