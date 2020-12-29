#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
#include "HostHelpersCore.h"
#include "SampleHelpersCore.h"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


int main(int argc, char** argv){
  constexpr int iarg_offset=1; // argv[0]==[Executable name]

  int exit_status = 1;
  bool print_help = false;
  TString strCmd;
  std::vector<TString> strArgs;
  for (int iarg=iarg_offset; iarg<argc; iarg++){
    TString strarg = argv[iarg];
    if (strarg == "help"){ print_help = true; exit_status = 0; break; }
    else if (strCmd == "") strCmd = strarg;
    else strArgs.push_back(strarg);
  }

  if (strCmd == "") print_help = true;
  if (!print_help){
    TString strCmdLower = strCmd; strCmdLower.ToLower();
    if (strCmdLower == "lstrip" || strCmdLower == "rstrip" || strCmdLower == "lrstrip"){
      if (strArgs.empty() || strArgs.size()>2) print_help = true;
      else{
        auto res = strArgs.front();
        const char* strip_chars = nullptr;
        if (strArgs.size()==2) strip_chars = strArgs.back().Data();
        if (strCmdLower == "lstrip") HelperFunctions::lstrip(res, strip_chars);
        else if (strCmdLower == "rstrip") HelperFunctions::rstrip(res, strip_chars);
        else HelperFunctions::lrstrip(res, strip_chars);
        cout << res << endl;
        exit_status = 0;
      }
    }
    else if (strCmdLower == "lsdir"){
      if (strArgs.size()!=1) print_help = true;
      else{
        auto vres = SampleHelpers::lsdir(strArgs.front());
        for (auto const& res:vres) cout << res << endl;
        exit_status = 0;
      }
    }
    else if (strCmdLower == "gethostlocalredirector" || strCmdLower == "gethostpathtostore"){
      if (strArgs.size()!=1) print_help = true;
      else{
        HostHelpers::Hosts theHost = HostHelpers::GetHostLocation(strArgs.front());
        if (strCmdLower == "gethostlocalredirector") cout << HostHelpers::GetHostLocalRedirector(theHost, true) << endl;
        else cout << HostHelpers::GetHostPathToStore(theHost) << endl;
        exit_status = 0;
      }
    }
    else if (strCmdLower == "getstandardhostpathtostore"){
      if (strArgs.size()!=2) print_help = true;
      else{
        HostHelpers::Hosts theHost = HostHelpers::GetHostLocation(strArgs.back());
        TString res = HostHelpers::GetStandardHostPathToStore(strArgs.front(), theHost);
        cout << res << endl;
        exit_status = 0;
      }
    }
    else{
      MELAerr << "Command " << strCmd << " with arguments " << strArgs << " is not supported." << endl;
      print_help = true;
    }
  }

  if (print_help){
    cout << "ExecuteCompiledCommand usage:" << endl;
    cout << "  help: Print this help message." << endl;
    cout << "  lstrip [line] [characters (optional)]: Strip line off a set of characters from the beginning. If character string is empty, strip off whitespace." << endl;
    cout << "  rstrip [line] [characters (optional)]: Strip line off a set of characters from the end. If character string is empty, strip off whitespace." << endl;
    cout << "  lrstrip [line] [characters (optional)]: lstrip and rstrip combined." << endl;
    cout << "  lsdir [directory]: The directory should already contain local redirectors if needed by an ordinary gfal-ls. Local redirectors should not be put for local files." << endl;
    cout << "  GetHostLocalRedirector [host identifier]: The host identifier should be a common string such as 't2.ucsd.edu', 'iihe.ac.be', or 'eos.cms.cern.ch'." << endl;
    cout << "  GetHostPathToStore [host identifier]: Same as above." << endl;
    cout << "  GetStandardHostPathToStore [path] [host identifier]: The path can start with '/store', or it can be an absolute path. Same as above for the host identifier." << endl;
    cout << "  GetX509Proxy: This is a call to get the best x509 user proxy, collecting a bunch of checks together." << endl;
  }

  return exit_status;
}
