#include "FileTransferHelpers.h"
#include "TSystem.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


void FileTransferHelpers::InitiateCondorFileTransfer(TString indir, TString fname, TString outsite, TString outdir, TString renamefile, int ntries){
  if (renamefile=="") renamefile = fname;
  if (indir.EndsWith("/")) indir.Remove(indir.Length()-1);
  if (outdir.EndsWith("/")) outdir.Remove(outdir.Length()-1);
  if (indir == "") indir = ".";

  TString strcmd = Form("copyFromCondorToSite.sh %s %s %s %s %s", indir.Data(), fname.Data(), outsite.Data(), outdir.Data(), renamefile.Data());
  MELAout << "FileTransferHelpers::InitiateCondorFileTransfer: Running command \"" << strcmd << "\"" << endl;
  int copy_status = HostHelpers::ExecuteCommand(strcmd);
  if (copy_status!=0){
    MELAerr << "FileTransferHelpers::InitiateCondorFileTransfer: Copy status is " << copy_status << " != 0, so something went wrong!";
    if (copy_status==-1){
      MELAerr << " The problem most likely occured due to the call to the system function, so no re-trials will be performed.";
      ntries = 0;
    }
    else if (ntries>0) MELAerr << " Retrying at most " << ntries << " times.";
    MELAerr << endl;
    if (ntries>0){
      int itry=0;
      while (copy_status!=0){
        MELAout << "FileTransferHelpers::InitiateCondorFileTransfer: Trial " << itry+1 << " after status " << copy_status << "..." << endl;
        copy_status = HostHelpers::ExecuteCommand(strcmd);
        itry++;
        if (itry>=ntries) break;
      }
    }
  }
  MELAout << "FileTransferHelpers::InitiateCondorFileTransfer: Final copy status is " << copy_status << endl;
  if (copy_status == -1){ // Workaround in order to ensure file transfer outside the script
    ofstream olf("EXTERNAL_TRANSFER_CMD_LIST.LST", ios_base::app);
    olf << strcmd << endl;
    olf.close();
  }
}
