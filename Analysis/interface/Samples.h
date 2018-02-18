#ifndef SAMPLES_H
#define SAMPLES_H

#include "HostHelpers.h"
#include <string>

// CJLST samples directory
const unsigned int CJLSTversion = 180212;
//const unsigned int CJLSTversion = 180121;
const TString CJLSTdate = std::to_string(CJLSTversion);
const TString CJLSTrootdir = HostHelpers::GetCJLSTSamplesDirectory(CJLSTdate);
const TString CJLSTsamplesdir = CJLSTrootdir + "/" + CJLSTdate;

const TString TREE_NAME = "ZZTree/candTree";
const TString TREE_FAILED_NAME = "ZZTree/candTree_failed";
const TString TREE_CRZLL_NAME = "CRZLLTree/candTree";
const TString TREE_CRZL_NAME = "CRZLTree/candTree";
const TString COUNTERS_NAME = "ZZTree/Counters";
const TString COUNTERS_CRZLL_NAME = "CRZLLTree/Counters";
const TString COUNTERS_CRZL_NAME = "CRZLTree/Counters";

const unsigned int theSqrts = 13;
const float xsecScale = 1e3;

#endif
