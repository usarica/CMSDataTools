#ifndef FILETRANSFERHELPERS_H
#define FILETRANSFERHELPERS_H

#include "HostHelpersCore.h"

namespace FileTransferHelpers{
  void InitiateCondorFileTransfer(TString indir, TString fname, TString outsite, TString outdir, TString renamefile="", int ntries=-1);
}


#endif
