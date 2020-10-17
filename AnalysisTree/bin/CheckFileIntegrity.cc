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

  bool print_help=false, has_help=false;
  string indir="";
  vector<string> filenames;

  for (int iarg=iarg_offset; iarg<argc; iarg++){
    string strarg = argv[iarg];
    string wish, value;
    splitOption(strarg, wish, value, '=');

    if (wish.empty()){
      if (value.find(".root")!=string::npos) filenames.push_back(value);
      else if (value=="help"){ print_help=has_help=true; }
    }
    else if (wish=="indir") indir = value;
    else{
      cerr << "CheckFileIntegrity: Unknown argument " << wish << " = " << value << endl;
      print_help=true;
    }
  }
  if (indir!="" && *(indir.rbegin())!='/') indir += '/';
  if (filenames.empty()){
    if (!has_help) cerr << "CheckFileIntegrity: No input file is specified." << endl;
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
    cout << "CheckFileIntegrity options:\n\n";
    cout << "- No option specifier: Input files with extension .root. Multiple input files can be passed as different arguments.\n\n";
    cout << "- indir: Location of input files. Default=\"\"\n\n";

    cout << endl;
    return 0;
  }

  for (auto const& fname:filenames){
    if (!HostHelpers::FileReadable(fname.data())){
      cerr << "CheckFileIntegrity: Input file " << fname << " is not readable!" << endl;
      return istatus_unreadable;
    }

    bool valid = false;
    bool isRecovered = false;

    TFile* finput = TFile::Open(fname.data(), "read");
    if (finput){
      if (finput->IsOpen() && !finput->IsZombie()){
        if (!finput->TestBit(TFile::kRecovered)) valid = true;
        else{
          cout << "CheckFileIntegrity: Input file " << fname << " went through a recovery, but we would not know if everything is fine. Will assume as a failure." << endl;
          isRecovered = true;
        }
        finput->Close();
      }
      else{
        if (finput->IsOpen()) finput->Close();
        else delete finput;
        finput=nullptr;
      }
    }

    if (!valid){
      cerr << "CheckFileIntegrity: Input file " << fname << " is invalid!" << endl;
      if (isRecovered) return istatus_invalid_recovered;
      else return istatus_invalid;
    }
  }

  return 0;
}
