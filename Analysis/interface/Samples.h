#ifndef SAMPLES_H
#define SAMPLES_H

#include "HostHelpers.h"
#include <string>

// Package directory
#ifndef xstr_lit
#define xstr_lit(s) str_lit(s)
#define str_lit(s) #s
#endif
#ifndef _higgswidthpkgpathstr_
#ifndef _higgswidthpkgpath_
#define _higgswidthpkgpath_ ./
#endif
#define _higgswidthpkgpathstr_ xstr_lit(_higgswidthpkgpath_)
#endif
const TString HIGGSWIDTHPKGPATH = _higgswidthpkgpathstr_;

// LHC sqrts and data period
constexpr unsigned int theSqrts = 13;
const TString theDataPeriod = "2016";

// MH and GH reference values
constexpr float MHRefVal=125;
constexpr float GHRefVal=4.07e-3;

// CJLST samples directory
constexpr unsigned int CJLSTversion = 180721;
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

// CJLST xsec*BR units
constexpr float xsecScale = 1e3;

// ZZMass infimum and supremum
constexpr float ZZMass_Infimum = 70.; // Analyzer cut
constexpr float ZZMass_Supremum = theSqrts*1000.;

#endif
