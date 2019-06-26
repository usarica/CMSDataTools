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

  bool print_help=false, has_help=false;
  string outfile, indir="./", outdir="./";
  vector<string> filenames;

  for (int iarg=iarg_offset; iarg<argc; iarg++){
    string strarg = argv[iarg];
    string wish, value;
    splitOption(strarg, wish, value, '=');

    if (wish.empty()){
      if (value.find(".lhe")!=string::npos) filenames.push_back(value);
      else if (value=="help"){ print_help=has_help=true; }
    }
    else if (wish=="indir") indir = value;
    else if (wish=="outdir") outdir = value;
    else if (wish=="outfile") outfile = value;
    else{
      cerr << "MergeLHEFiles: Unknown argument " << wish << " = " << value << endl;
      print_help=true;
    }
  }
  if (*(indir.rbegin())!='/') indir += '/';
  if (*(outdir.rbegin())!='/') outdir += '/';
  if (outfile.empty()){
    if (!has_help) cerr << "MergeLHEFiles: No output file is specified. This argument is mandatory." << endl;
    print_help |= true;
  }
  else outfile = outdir + outfile;
  if (filenames.empty()){
    if (!has_help) cerr << "MergeLHEFiles: No input file is specified." << endl;
    print_help |= true;
  }
  else{ for (auto& s:filenames) s = indir + s; }
  if (print_help){
    cout << "MergeLHEFiles options:\n\n";
    cout << "- No option specifier: Input files with extension .lhe. Multiple input files can be passed as different arguments.\n\n";
    cout << "- outfile: Output file name (mandatory).\n\n";
    cout << "- indir: Location of input files. Default=\"./\"\n\n";
    cout << "- outdir: Location of the output file. Default=\"./\"\n\n";

    cout << endl;
    return 0;
  }

  string const strMatch_lhebegin = "<LesHouchesEvents";
  string const strMatch_lheend = "</LesHouchesEvents>";
  string const strMatch_lhecommentbegin = "<!--";
  string const strMatch_lhecommentend = "-->";
  string const strMatch_initbegin = "<init>";
  string const strMatch_initend = "</init>";
  string const strMatch_eventbegin = "<event>";
  string const strMatch_eventend = "</event>";

  string const strMatch_lhecomment_leptonfilter = "Lepton filter information";

  size_t ninputfiles=filenames.size();
  vector<ifstream> inputfiles(ninputfiles);
  for (size_t ifile=0; ifile<ninputfiles; ifile++){
    string const& fname = filenames.at(ifile);
    ifstream& fin = inputfiles.at(ifile);
    fin.open(fname.c_str());
    if (fin.good()) cout << "MergeLHEFiles: Input file " << fname << " is opened." << endl;
    else{ cerr << "MergeLHEFiles: Input file " << fname << " cannot be opened." << endl; fin.close(); }
  }

  ofstream fout(outfile.c_str());
  if (fout.good()) cout << "MergeLHEFiles: Output file " << outfile << " is opened." << endl;
  else{ cerr << "MergeLHEFiles: Output file " << outfile << " cannot be opened." << endl; fout.close(); }

  if (fout.is_open()){
    bool hasLeptonFilterComment=false;
    unsigned int nLeptonFilterEventsProcessed=0;
    unsigned int nLeptonFilterEventsAccepted=0;

    bool write_lhebegin=false;
    bool write_lhebegincomments=false;
    bool write_lhebegincomments_progress=false;
    bool write_init=false;
    bool write_init_progress=false;
    //bool write_lheendcomments=false;
    //bool write_lheend=false;

    // Initial LHE begin, comments and the init block
    // Should be written only once per LHE output
    for (ifstream& fin:inputfiles){
      if (!fin.is_open()) continue;
      while (!fin.eof()){
        string strin="";
        getline(fin, strin);

        if (!write_lhebegin && strin.find(strMatch_lhebegin.c_str())!=string::npos){
          fout << strin << endl;
          write_lhebegin=true;
        }
        if (!write_lhebegincomments && (write_lhebegincomments_progress || strin.find(strMatch_lhecommentbegin.c_str())!=string::npos)){
          fout << strin << endl;
          write_lhebegincomments_progress=true;
          if (strin.find(strMatch_lhecommentend.c_str())!=string::npos){ // Stop writing LHE begin comments
            write_lhebegincomments_progress=false;
            write_lhebegincomments=true;
          }
        }
        if (!write_init && (write_init_progress || strin.find(strMatch_initbegin.c_str())!=string::npos)){
          fout << strin << endl;
          write_init_progress=true;
        }

        if (strin.find(strMatch_initend.c_str())!=string::npos){ // Stop writing the init block or looping over the file
          write_init_progress=false;
          write_init=true;
          break; // Break from the input file for now
        }
      }
    }

    // Event writing step
    for (ifstream& fin:inputfiles){
      if (!fin.is_open()) continue;
      if (fin.eof()) cerr << "File is at the end! Cannot proceed to event writing." << endl;
      else cout << "Processing event writing step..." << endl;

      bool write_event=false;
      bool write_event_progress=false;
      bool accumulate_leptonefficiency=false;

      while (!fin.eof()){
        string strin="";
        getline(fin, strin);

        if (strin.find(strMatch_lheend.c_str())!=string::npos) break;

        if (!write_event_progress){
          if (strin.find(strMatch_lhecommentbegin.c_str())!=string::npos) write_event=true;
          if (strin.find(strMatch_lhecomment_leptonfilter.c_str())!=string::npos){
            accumulate_leptonefficiency=true;
          }
          if (strin.find(strMatch_lhecommentend.c_str())!=string::npos){
            accumulate_leptonefficiency=false;
          }
        }

        if (accumulate_leptonefficiency){
          hasLeptonFilterComment=true;
          if (strin.find("events processed")!=string::npos){
            while (strin.find(' ')!=string::npos) replaceString<std::string, const std::string>(strin, " ", "");
            string wish, value;
            splitOption(strin, wish, value, ':');
            unsigned int ntmp = atoi(value.c_str());
            nLeptonFilterEventsProcessed += ntmp;
          }
          else if (strin.find("events accepted")!=string::npos){
            while (strin.find(' ')!=string::npos) replaceString<std::string, const std::string>(strin, " ", "");
            string wish, value;
            splitOption(strin, wish, value, ':');
            unsigned int ntmp = atoi(value.c_str());
            nLeptonFilterEventsAccepted += ntmp;
          }
        }

        if (!write_event && (write_event_progress || strin.find(strMatch_eventbegin.c_str())!=string::npos)){
          fout << strin << endl;
          write_event_progress=true;
          if (strin.find(strMatch_eventend.c_str())!=string::npos) write_event_progress=false;
        }
      }
    }

    // Special lines for the end of the output LHE file
    if (hasLeptonFilterComment){
      fout << "<!-- Filter information:" << endl;
      fout << "     Events processed: " << nLeptonFilterEventsProcessed << endl;
      fout << "     Events accepted: " << nLeptonFilterEventsAccepted << endl;
      fout << "     Filter efficiency: " << double(nLeptonFilterEventsAccepted)/double(nLeptonFilterEventsProcessed) << ' ' << strMatch_lhecommentend << endl; // Writes in fraction rather than percentage as done in original files
    }
    fout << strMatch_lheend << endl; // Finalize the output LHE file
    fout.close();
  }

  // Close the input LHE files
  for (ifstream& fin:inputfiles){
    if (!fin.is_open()) continue;
    fin.close();
  }
}
