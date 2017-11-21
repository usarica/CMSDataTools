#ifndef SAMPLES_H
#define SAMPLES_H

#include "HostHelpers.h"

// CJLST samples directory
const TString CJLSTdate = "171104";
const TString CJLSTrootdir = HostHelpers::GetCJLSTSamplesDirectory(CJLSTdate);
const TString CJLSTsamplesdir = CJLSTrootdir + "/" + CJLSTdate;

const TString TREE_NAME = "ZZTree/candTree";
const TString TREE_FAILED_NAME = "ZZTree/candTree_failed";
const TString COUNTERS_NAME = "ZZTree/Counters";

const unsigned int theSqrts = 13;

#endif
