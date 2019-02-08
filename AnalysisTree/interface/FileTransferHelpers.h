#ifndef FILETRANSFERHELPERS_H
#define FILETRANSFERHELPERS_H

#include "HostHelpersCore.h"

namespace FileTransferHelpers{
  void InitiateFileTransfer(TString indir, TString fname, TString outsite, TString outdir, TString renamefile="");
}


#endif
