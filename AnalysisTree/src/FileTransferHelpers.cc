#include "FileTransferHelpers.h"
#include "TSystem.h"


void FileTransferHelpers::InitiateFileTransfer(TString indir, TString fname, TString outsite, TString outdir, TString renamefile){
  if (renamefile=="") renamefile = fname;
  if (indir.EndsWith("/")) indir.Remove(indir.Length()-1);
  if (outdir.EndsWith("/")) outdir.Remove(outdir.Length()-1);

  gSystem->Exec(Form("copyFromCondorToSite.sh \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"", indir.Data(), fname.Data(), outsite.Data(), outdir.Data(), renamefile.Data()));
}
