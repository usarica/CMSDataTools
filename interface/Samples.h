#ifndef SAMPLES_H
#define SAMPLES_H

#include "HostHelpers.h"

// CJLST samples directory
const TString CJLSTdate = "170623";
const TString CJLSTrootdir = HostHelpers::GetCJLSTSamplesDirectory(CJLSTdate);
const TString CJLSTsamplesdir = CJLSTrootdir + "/" + CJLSTdate;


#endif
