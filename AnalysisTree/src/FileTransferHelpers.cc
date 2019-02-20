#include "FileTransferHelpers.h"
#include "TSystem.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


void FileTransferHelpers::InitiateCondorFileTransfer(TString indir, TString fname, TString outsite, TString outdir, TString renamefile, int ntries){
  if (renamefile=="") renamefile = fname;
  if (indir.EndsWith("/")) indir.Remove(indir.Length()-1);
  if (outdir.EndsWith("/")) outdir.Remove(outdir.Length()-1);

  int copy_status = gSystem->Exec(Form("copyFromCondorToSite.sh \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"", indir.Data(), fname.Data(), outsite.Data(), outdir.Data(), renamefile.Data()));
  if (copy_status!=0 && ntries>0){
    MELAerr << "FileTransferHelpers::InitiateCondorFileTransfer: Copy status is " << copy_status << " != 0, so something went wrong! Retrying at most " << ntries << " times." << endl;
    int itry=0;
    while (copy_status!=0){
      copy_status = gSystem->Exec(Form("copyFromCondorToSite.sh \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"", indir.Data(), fname.Data(), outsite.Data(), outdir.Data(), renamefile.Data()));
      itry++;
      if (itry==ntries) break;
    }
  }
  MELAout << "FileTransferHelpers::InitiateCondorFileTransfer: Final copy status is " << copy_status << endl;
}
