#ifndef BASEEDMINPUTTREE_HPP
#define BASEEDMINPUTTREE_HPP

#include "SampleHelpersCore.h"
#include "BaseEDMInputTree.h"


template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_TBits_t>(){ BaseTree::resetBranch<BaseTree::BranchType_TBits_t>(); for (auto& it:bridgeTBitss){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_bool_t>(){ BaseTree::resetBranch<BaseTree::BranchType_bool_t>(); for (auto& it:bridgebools){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_short_t>(){ BaseTree::resetBranch<BaseTree::BranchType_short_t>(); for (auto& it:bridgeshorts){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_uint_t>(){ BaseTree::resetBranch<BaseTree::BranchType_uint_t>(); for (auto& it:bridgeuints){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_int_t>(){ BaseTree::resetBranch<BaseTree::BranchType_int_t>(); for (auto& it:bridgeints){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_ulong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_ulong_t>(); for (auto& it:bridgeulongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_long_t>(){ BaseTree::resetBranch<BaseTree::BranchType_long_t>(); for (auto& it:bridgelongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_ulonglong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_ulonglong_t>(); for (auto& it:bridgeulonglongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_longlong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_longlong_t>(); for (auto& it:bridgelonglongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_float_t>(){ BaseTree::resetBranch<BaseTree::BranchType_float_t>(); for (auto& it:bridgefloats){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_double_t>(){ BaseTree::resetBranch<BaseTree::BranchType_double_t>(); for (auto& it:bridgedoubles){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_string_t>(){ BaseTree::resetBranch<BaseTree::BranchType_string_t>(); for (auto& it:bridgestrings){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_TString_t>(){ BaseTree::resetBranch<BaseTree::BranchType_TString_t>(); for (auto& it:bridgeTStrings){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_CMSLorentzVector_t>(){ BaseTree::resetBranch<BaseTree::BranchType_CMSLorentzVector_t>(); for (auto& it:bridgeCMSLorentzVectors){ if (it.second) it.second->reset(); } }

template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vbool_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vbool_t>(); for (auto& it:bridgeVbools){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vshort_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vshort_t>(); for (auto& it:bridgeVshorts){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vuint_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vuint_t>(); for (auto& it:bridgeVuints){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vint_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vint_t>(); for (auto& it:bridgeVints){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vulong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vulong_t>(); for (auto& it:bridgeVulongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vlong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vlong_t>(); for (auto& it:bridgeVlongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vulonglong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vulonglong_t>(); for (auto& it:bridgeVulonglongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vlonglong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vlonglong_t>(); for (auto& it:bridgeVlonglongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vfloat_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vfloat_t>(); for (auto& it:bridgeVfloats){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vdouble_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vdouble_t>(); for (auto& it:bridgeVdoubles){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vstring_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vstring_t>(); for (auto& it:bridgeVstrings){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vTString_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vTString_t>(); for (auto& it:bridgeVTStrings){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vCMSLorentzVector_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vCMSLorentzVector_t>(); for (auto& it:bridgeVCMSLorentzVectors){ if (it.second) it.second->reset(); } }

template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvbool_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvbool_t>(); for (auto& it:bridgeVVbools){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvshort_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvshort_t>(); for (auto& it:bridgeVVshorts){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvuint_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvuint_t>(); for (auto& it:bridgeVVuints){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvint_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvint_t>(); for (auto& it:bridgeVVints){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvulong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvulong_t>(); for (auto& it:bridgeVVulongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvlong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvlong_t>(); for (auto& it:bridgeVVlongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvulonglong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvulonglong_t>(); for (auto& it:bridgeVVulonglongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvlonglong_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvlonglong_t>(); for (auto& it:bridgeVVlonglongs){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvfloat_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvfloat_t>(); for (auto& it:bridgeVVfloats){ if (it.second) it.second->reset(); } }
template<> void BaseEDMInputTree::resetEDMBranch<BaseTree::BranchType_vvdouble_t>(){ BaseTree::resetBranch<BaseTree::BranchType_vvdouble_t>(); for (auto& it:bridgeVVdoubles){ if (it.second) it.second->reset(); } }


template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_TBits_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_TBits_t>(branchname); for (auto& it:bridgeTBitss){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeTBitss.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_bool_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_bool_t>(branchname); for (auto& it:bridgebools){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgebools.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_short_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_short_t>(branchname); for (auto& it:bridgeshorts){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeshorts.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_uint_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_uint_t>(branchname); for (auto& it:bridgeuints){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeuints.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_int_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_int_t>(branchname); for (auto& it:bridgeints){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeints.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_ulong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_ulong_t>(branchname); for (auto& it:bridgeulongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeulongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_long_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_long_t>(branchname); for (auto& it:bridgelongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgelongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_ulonglong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_ulonglong_t>(branchname); for (auto& it:bridgeulonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeulonglongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_longlong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_longlong_t>(branchname); for (auto& it:bridgelonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgelonglongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_float_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_float_t>(branchname); for (auto& it:bridgefloats){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgefloats.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_double_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_double_t>(branchname); for (auto& it:bridgedoubles){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgedoubles.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_string_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_string_t>(branchname); for (auto& it:bridgestrings){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgestrings.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_TString_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_TString_t>(branchname); for (auto& it:bridgeTStrings){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeTStrings.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_CMSLorentzVector_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_CMSLorentzVector_t>(branchname); for (auto& it:bridgeCMSLorentzVectors){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeCMSLorentzVectors.erase(branchname); }

template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vbool_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vbool_t>(branchname); for (auto& it:bridgeVbools){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVbools.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vshort_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vshort_t>(branchname); for (auto& it:bridgeVshorts){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVshorts.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vuint_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vuint_t>(branchname); for (auto& it:bridgeVuints){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVuints.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vint_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vint_t>(branchname); for (auto& it:bridgeVints){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVints.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vulong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vulong_t>(branchname); for (auto& it:bridgeVulongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVulongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vlong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vlong_t>(branchname); for (auto& it:bridgeVlongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVlongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vulonglong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vulonglong_t>(branchname); for (auto& it:bridgeVulonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVulonglongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vlonglong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vlonglong_t>(branchname); for (auto& it:bridgeVlonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVlonglongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vfloat_t>(branchname); for (auto& it:bridgeVfloats){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVfloats.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vdouble_t>(branchname); for (auto& it:bridgeVdoubles){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVdoubles.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vstring_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vstring_t>(branchname); for (auto& it:bridgeVstrings){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVstrings.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vTString_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vTString_t>(branchname); for (auto& it:bridgeVTStrings){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVTStrings.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vCMSLorentzVector_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vCMSLorentzVector_t>(branchname); for (auto& it:bridgeVCMSLorentzVectors){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVCMSLorentzVectors.erase(branchname); }

template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvbool_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvbool_t>(branchname); for (auto& it:bridgeVVbools){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVbools.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvshort_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvshort_t>(branchname); for (auto& it:bridgeVVshorts){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVshorts.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvuint_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvuint_t>(branchname); for (auto& it:bridgeVVuints){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVuints.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvint_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvint_t>(branchname); for (auto& it:bridgeVVints){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVints.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvulong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvulong_t>(branchname); for (auto& it:bridgeVVulongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVulongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvlong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvlong_t>(branchname); for (auto& it:bridgeVVlongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVlongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvulonglong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvulonglong_t>(branchname); for (auto& it:bridgeVVulonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVulonglongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvlonglong_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvlonglong_t>(branchname); for (auto& it:bridgeVVlonglongs){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVlonglongs.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvfloat_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvfloat_t>(branchname); for (auto& it:bridgeVVfloats){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVfloats.erase(branchname); }
template<> void BaseEDMInputTree::removeEDMBranch<BaseTree::BranchType_vvdouble_t>(TString branchname){ BaseTree::removeBranch<BaseTree::BranchType_vvdouble_t>(branchname); for (auto& it:bridgeVVdoubles){ if (it.first==branchname){ delete it.second; it.second=0; } } bridgeVVdoubles.erase(branchname); }


template<> bool BaseEDMInputTree::bookEDMBranch<TBits>(TString branchname, TBits valdef){
  if (valTBitss.find(branchname)==valTBitss.end()) valTBitss[branchname] = new std::pair<TBits, TBits>(valdef, valdef);
  else{ valTBitss[branchname]->first=valdef; valTBitss[branchname]->second=valdef; }
  if (bridgeTBitss.find(branchname)==bridgeTBitss.end()) bridgeTBitss[branchname] = new CMSEDMWrapperLinker<TBits>(&(valTBitss[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeTBitss[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeTBitss[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<bool>(TString branchname, bool valdef){
  if (valbools.find(branchname)==valbools.end()) valbools[branchname] = new std::pair<bool, bool>(valdef, valdef);
  else{ valbools[branchname]->first=valdef; valbools[branchname]->second=valdef; }
  if (bridgebools.find(branchname)==bridgebools.end()) bridgebools[branchname] = new CMSEDMWrapperLinker<bool>(&(valbools[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgebools[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgebools[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<short>(TString branchname, short valdef){
  if (valshorts.find(branchname)==valshorts.end()) valshorts[branchname] = new std::pair<short, short>(valdef, valdef);
  else{ valshorts[branchname]->first=valdef; valshorts[branchname]->second=valdef; }
  if (bridgeshorts.find(branchname)==bridgeshorts.end()) bridgeshorts[branchname] = new CMSEDMWrapperLinker<short>(&(valshorts[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeshorts[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeshorts[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<unsigned int>(TString branchname, unsigned int valdef){
  if (valuints.find(branchname)==valuints.end()) valuints[branchname] = new std::pair<unsigned int, unsigned int>(valdef, valdef);
  else{ valuints[branchname]->first=valdef; valuints[branchname]->second=valdef; }
  if (bridgeuints.find(branchname)==bridgeuints.end()) bridgeuints[branchname] = new CMSEDMWrapperLinker<unsigned int>(&(valuints[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeuints[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeuints[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<int>(TString branchname, int valdef){
  if (valints.find(branchname)==valints.end()) valints[branchname] = new std::pair<int, int>(valdef, valdef);
  else{ valints[branchname]->first=valdef; valints[branchname]->second=valdef; }
  if (bridgeints.find(branchname)==bridgeints.end()) bridgeints[branchname] = new CMSEDMWrapperLinker<int>(&(valints[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeints[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeints[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<unsigned long>(TString branchname, unsigned long valdef){
  if (valulongs.find(branchname)==valulongs.end()) valulongs[branchname] = new std::pair<unsigned long, unsigned long>(valdef, valdef);
  else{ valulongs[branchname]->first=valdef; valulongs[branchname]->second=valdef; }
  if (bridgeulongs.find(branchname)==bridgeulongs.end()) bridgeulongs[branchname] = new CMSEDMWrapperLinker<unsigned long>(&(valulongs[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeulongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeulongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<long>(TString branchname, long valdef){
  if (vallongs.find(branchname)==vallongs.end()) vallongs[branchname] = new std::pair<long, long>(valdef, valdef);
  else{ vallongs[branchname]->first=valdef; vallongs[branchname]->second=valdef; }
  if (bridgelongs.find(branchname)==bridgelongs.end()) bridgelongs[branchname] = new CMSEDMWrapperLinker<long>(&(vallongs[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgelongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgelongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<unsigned long long>(TString branchname, unsigned long long valdef){
  if (valulonglongs.find(branchname)==valulonglongs.end()) valulonglongs[branchname] = new std::pair<unsigned long long, unsigned long long>(valdef, valdef);
  else{ valulonglongs[branchname]->first=valdef; valulonglongs[branchname]->second=valdef; }
  if (bridgeulonglongs.find(branchname)==bridgeulonglongs.end()) bridgeulonglongs[branchname] = new CMSEDMWrapperLinker<unsigned long long>(&(valulonglongs[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeulonglongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeulonglongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<long long>(TString branchname, long long valdef){
  if (vallonglongs.find(branchname)==vallonglongs.end()) vallonglongs[branchname] = new std::pair<long long, long long>(valdef, valdef);
  else{ vallonglongs[branchname]->first=valdef; vallonglongs[branchname]->second=valdef; }
  if (bridgelonglongs.find(branchname)==bridgelonglongs.end()) bridgelonglongs[branchname] = new CMSEDMWrapperLinker<long long>(&(vallonglongs[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgelonglongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgelonglongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<float>(TString branchname, float valdef){
  if (valfloats.find(branchname)==valfloats.end()) valfloats[branchname] = new std::pair<float, float>(valdef, valdef);
  else{ valfloats[branchname]->first=valdef; valfloats[branchname]->second=valdef; }
  if (bridgefloats.find(branchname)==bridgefloats.end()) bridgefloats[branchname] = new CMSEDMWrapperLinker<float>(&(valfloats[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgefloats[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgefloats[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<double>(TString branchname, double valdef){
  if (valdoubles.find(branchname)==valdoubles.end()) valdoubles[branchname] = new std::pair<double, double>(valdef, valdef);
  else{ valdoubles[branchname]->first=valdef; valdoubles[branchname]->second=valdef; }
  if (bridgedoubles.find(branchname)==bridgedoubles.end()) bridgedoubles[branchname] = new CMSEDMWrapperLinker<double>(&(valdoubles[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgedoubles[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgedoubles[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::string>(TString branchname, std::string valdef){
  if (valstrings.find(branchname)==valstrings.end()) valstrings[branchname] = new std::pair<std::string, std::string>(valdef, valdef);
  else{ valstrings[branchname]->first=valdef; valstrings[branchname]->second=valdef; }
  if (bridgestrings.find(branchname)==bridgestrings.end()) bridgestrings[branchname] = new CMSEDMWrapperLinker<std::string>(&(valstrings[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgestrings[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgestrings[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<TString>(TString branchname, TString valdef){
  if (valTStrings.find(branchname)==valTStrings.end()) valTStrings[branchname] = new std::pair<TString, TString>(valdef, valdef);
  else{ valTStrings[branchname]->first=valdef; valTStrings[branchname]->second=valdef; }
  if (bridgeTStrings.find(branchname)==bridgeTStrings.end()) bridgeTStrings[branchname] = new CMSEDMWrapperLinker<TString>(&(valTStrings[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeTStrings[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeTStrings[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<CMSLorentzVector>(TString branchname, CMSLorentzVector valdef){
  if (valCMSLorentzVectors.find(branchname)==valCMSLorentzVectors.end()) valCMSLorentzVectors[branchname] = new std::pair<CMSLorentzVector, CMSLorentzVector>(valdef, valdef);
  else{ valCMSLorentzVectors[branchname]->first=valdef; valCMSLorentzVectors[branchname]->second=valdef; }
  if (bridgeCMSLorentzVectors.find(branchname)==bridgeCMSLorentzVectors.end()) bridgeCMSLorentzVectors[branchname] = new CMSEDMWrapperLinker<CMSLorentzVector>(&(valCMSLorentzVectors[branchname]->first));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeCMSLorentzVectors[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeCMSLorentzVectors[branchname]->getWrapperRef()));
  return true;
}

template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<bool>*>(TString branchname, std::vector<bool>*/* valdef*/){
  valVbools[branchname] = nullptr;
  if (bridgeVbools.find(branchname)==bridgeVbools.end()) bridgeVbools[branchname] = new CMSEDMWrapperLinker<std::vector<bool>, std::vector<bool>*>(&(valVbools[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVbools[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVbools[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<short>*>(TString branchname, std::vector<short>*/* valdef*/){
  valVshorts[branchname] = nullptr;
  if (bridgeVshorts.find(branchname)==bridgeVshorts.end()) bridgeVshorts[branchname] = new CMSEDMWrapperLinker<std::vector<short>, std::vector<short>*>(&(valVshorts[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVshorts[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVshorts[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<unsigned int>*>(TString branchname, std::vector<unsigned int>*/* valdef*/){
  valVuints[branchname] = nullptr;
  if (bridgeVuints.find(branchname)==bridgeVuints.end()) bridgeVuints[branchname] = new CMSEDMWrapperLinker<std::vector<unsigned int>, std::vector<unsigned int>*>(&(valVuints[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVuints[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVuints[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<int>*>(TString branchname, std::vector<int>*/* valdef*/){
  valVints[branchname] = nullptr;
  if (bridgeVints.find(branchname)==bridgeVints.end()) bridgeVints[branchname] = new CMSEDMWrapperLinker<std::vector<int>, std::vector<int>*>(&(valVints[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVints[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVints[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<unsigned long>*>(TString branchname, std::vector<unsigned long>*/* valdef*/){
  valVulongs[branchname] = nullptr;
  if (bridgeVulongs.find(branchname)==bridgeVulongs.end()) bridgeVulongs[branchname] = new CMSEDMWrapperLinker<std::vector<unsigned long>, std::vector<unsigned long>*>(&(valVulongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVulongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVulongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<long>*>(TString branchname, std::vector<long>*/* valdef*/){
  valVlongs[branchname] = nullptr;
  if (bridgeVlongs.find(branchname)==bridgeVlongs.end()) bridgeVlongs[branchname] = new CMSEDMWrapperLinker<std::vector<long>, std::vector<long>*>(&(valVlongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVlongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVlongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<unsigned long long>*>(TString branchname, std::vector<unsigned long long>*/* valdef*/){
  valVulonglongs[branchname] = nullptr;
  if (bridgeVulonglongs.find(branchname)==bridgeVulonglongs.end()) bridgeVulonglongs[branchname] = new CMSEDMWrapperLinker<std::vector<unsigned long long>, std::vector<unsigned long long>*>(&(valVulonglongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVulonglongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVulonglongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<long long>*>(TString branchname, std::vector<long long>*/* valdef*/){
  valVlonglongs[branchname] = nullptr;
  if (bridgeVlonglongs.find(branchname)==bridgeVlonglongs.end()) bridgeVlonglongs[branchname] = new CMSEDMWrapperLinker<std::vector<long long>, std::vector<long long>*>(&(valVlonglongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVlonglongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVlonglongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<float>*>(TString branchname, std::vector<float>*/* valdef*/){
  valVfloats[branchname] = nullptr;
  if (bridgeVfloats.find(branchname)==bridgeVfloats.end()) bridgeVfloats[branchname] = new CMSEDMWrapperLinker<std::vector<float>, std::vector<float>*>(&(valVfloats[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVfloats[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVfloats[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<double>*>(TString branchname, std::vector<double>*/* valdef*/){
  valVdoubles[branchname] = nullptr;
  if (bridgeVdoubles.find(branchname)==bridgeVdoubles.end()) bridgeVdoubles[branchname] = new CMSEDMWrapperLinker<std::vector<double>, std::vector<double>*>(&(valVdoubles[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVdoubles[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVdoubles[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::string>*>(TString branchname, std::vector<std::string>*/* valdef*/){
  valVstrings[branchname] = nullptr;
  if (bridgeVstrings.find(branchname)==bridgeVstrings.end()) bridgeVstrings[branchname] = new CMSEDMWrapperLinker<std::vector<std::string>, std::vector<std::string>*>(&(valVstrings[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVstrings[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVstrings[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<TString>*>(TString branchname, std::vector<TString>*/* valdef*/){
  valVTStrings[branchname] = nullptr;
  if (bridgeVTStrings.find(branchname)==bridgeVTStrings.end()) bridgeVTStrings[branchname] = new CMSEDMWrapperLinker<std::vector<TString>, std::vector<TString>*>(&(valVTStrings[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVTStrings[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVTStrings[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<CMSLorentzVector>*>(TString branchname, std::vector<CMSLorentzVector>*/* valdef*/){
  valVCMSLorentzVectors[branchname] = nullptr;
  if (bridgeVCMSLorentzVectors.find(branchname)==bridgeVCMSLorentzVectors.end()) bridgeVCMSLorentzVectors[branchname] = new CMSEDMWrapperLinker<std::vector<CMSLorentzVector>, std::vector<CMSLorentzVector>*>(&(valVCMSLorentzVectors[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVCMSLorentzVectors[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVCMSLorentzVectors[branchname]->getWrapperRef()));
  return true;
}

template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<bool>>*>(TString branchname, std::vector<std::vector<bool>>*/* valdef*/){
  valVVbools[branchname] = nullptr;
  if (bridgeVVbools.find(branchname)==bridgeVVbools.end()) bridgeVVbools[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<bool>>, std::vector<std::vector<bool>>*>(&(valVVbools[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVbools[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVbools[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<short>>*>(TString branchname, std::vector<std::vector<short>>*/* valdef*/){
  valVVshorts[branchname] = nullptr;
  if (bridgeVVshorts.find(branchname)==bridgeVVshorts.end()) bridgeVVshorts[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<short>>, std::vector<std::vector<short>>*>(&(valVVshorts[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVshorts[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVshorts[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<unsigned int>>*>(TString branchname, std::vector<std::vector<unsigned int>>*/* valdef*/){
  valVVuints[branchname] = nullptr;
  if (bridgeVVuints.find(branchname)==bridgeVVuints.end()) bridgeVVuints[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>*>(&(valVVuints[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVuints[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVuints[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<int>>*>(TString branchname, std::vector<std::vector<int>>*/* valdef*/){
  valVVints[branchname] = nullptr;
  if (bridgeVVints.find(branchname)==bridgeVVints.end()) bridgeVVints[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<int>>, std::vector<std::vector<int>>*>(&(valVVints[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVints[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVints[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<unsigned long>>*>(TString branchname, std::vector<std::vector<unsigned long>>*/* valdef*/){
  valVVulongs[branchname] = nullptr;
  if (bridgeVVulongs.find(branchname)==bridgeVVulongs.end()) bridgeVVulongs[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<unsigned long>>, std::vector<std::vector<unsigned long>>*>(&(valVVulongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVulongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVulongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<long>>*>(TString branchname, std::vector<std::vector<long>>*/* valdef*/){
  valVVlongs[branchname] = nullptr;
  if (bridgeVVlongs.find(branchname)==bridgeVVlongs.end()) bridgeVVlongs[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<long>>, std::vector<std::vector<long>>*>(&(valVVlongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVlongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVlongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<unsigned long long>>*>(TString branchname, std::vector<std::vector<unsigned long long>>*/* valdef*/){
  valVVulonglongs[branchname] = nullptr;
  if (bridgeVVulonglongs.find(branchname)==bridgeVVulonglongs.end()) bridgeVVulonglongs[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<unsigned long long>>, std::vector<std::vector<unsigned long long>>*>(&(valVVulonglongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVulonglongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVulonglongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<long long>>*>(TString branchname, std::vector<std::vector<long long>>*/* valdef*/){
  valVVlonglongs[branchname] = nullptr;
  if (bridgeVVlonglongs.find(branchname)==bridgeVVlonglongs.end()) bridgeVVlonglongs[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<long long>>, std::vector<std::vector<long long>>*>(&(valVVlonglongs[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVlonglongs[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVlonglongs[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<float>>*>(TString branchname, std::vector<std::vector<float>>*/* valdef*/){
  valVVfloats[branchname] = nullptr;
  if (bridgeVVfloats.find(branchname)==bridgeVVfloats.end()) bridgeVVfloats[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<float>>, std::vector<std::vector<float>>*>(&(valVVfloats[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVfloats[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVfloats[branchname]->getWrapperRef()));
  return true;
}
template<> bool BaseEDMInputTree::bookEDMBranch<std::vector<std::vector<double>>*>(TString branchname, std::vector<std::vector<double>>*/* valdef*/){
  valVVdoubles[branchname] = nullptr;
  if (bridgeVVdoubles.find(branchname)==bridgeVVdoubles.end()) bridgeVVdoubles[branchname] = new CMSEDMWrapperLinker<std::vector<std::vector<double>>, std::vector<std::vector<double>>*>(&(valVVdoubles[branchname]));
  SampleHelpers::bookEDMBranch(tree, branchname, &(bridgeVVdoubles[branchname]->getWrapperRef()));
  SampleHelpers::bookEDMBranch(failedtree, branchname, &(bridgeVVdoubles[branchname]->getWrapperRef()));
  return true;
}


template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_TBits_t>(TString branchname){ return this->bookEDMBranch<TBits>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_bool_t>(TString branchname){ return this->bookEDMBranch<bool>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_short_t>(TString branchname){ return this->bookEDMBranch<short>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_uint_t>(TString branchname){ return this->bookEDMBranch<unsigned int>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_int_t>(TString branchname){ return this->bookEDMBranch<int>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_ulong_t>(TString branchname){ return this->bookEDMBranch<unsigned long>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_long_t>(TString branchname){ return this->bookEDMBranch<long>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_ulonglong_t>(TString branchname){ return this->bookEDMBranch<unsigned long long>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_longlong_t>(TString branchname){ return this->bookEDMBranch<long long>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_float_t>(TString branchname){ return this->bookEDMBranch<float>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_double_t>(TString branchname){ return this->bookEDMBranch<double>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_string_t>(TString branchname){ return this->bookEDMBranch<std::string>(branchname, ""); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_TString_t>(TString branchname){ return this->bookEDMBranch<TString>(branchname, ""); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_CMSLorentzVector_t>(TString branchname){ return this->bookEDMBranch<CMSLorentzVector>(branchname, CMSLorentzVector(0, 0, 0, 0)); }

template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vbool_t>(TString branchname){ return this->bookEDMBranch<std::vector<bool>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vshort_t>(TString branchname){ return this->bookEDMBranch<std::vector<short>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vuint_t>(TString branchname){ return this->bookEDMBranch<std::vector<unsigned int>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vint_t>(TString branchname){ return this->bookEDMBranch<std::vector<int>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vulong_t>(TString branchname){ return this->bookEDMBranch<std::vector<unsigned long>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vlong_t>(TString branchname){ return this->bookEDMBranch<std::vector<long>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vulonglong_t>(TString branchname){ return this->bookEDMBranch<std::vector<unsigned long long>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vlonglong_t>(TString branchname){ return this->bookEDMBranch<std::vector<long long>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vfloat_t>(TString branchname){ return this->bookEDMBranch<std::vector<float>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vdouble_t>(TString branchname){ return this->bookEDMBranch<std::vector<double>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vstring_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::string>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vTString_t>(TString branchname){ return this->bookEDMBranch<std::vector<TString>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vCMSLorentzVector_t>(TString branchname){ return this->bookEDMBranch<std::vector<CMSLorentzVector>*>(branchname, 0); }

template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvbool_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<bool>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvshort_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<short>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvuint_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<unsigned int>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvint_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<int>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvulong_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<unsigned long>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvlong_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<long>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvulonglong_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<unsigned long long>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvlonglong_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<long long>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvfloat_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<float>>*>(branchname, 0); }
template<> bool BaseEDMInputTree::bookEDMBranch<BaseTree::BranchType_vvdouble_t>(TString branchname){ return this->bookEDMBranch<std::vector<std::vector<double>>*>(branchname, 0); }


template<> void CMSEDMWrapperLinker<TBits>::assignProductToTarget(CMSEDMWrapperLinker<TBits>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<bool>::assignProductToTarget(CMSEDMWrapperLinker<bool>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<short>::assignProductToTarget(CMSEDMWrapperLinker<short>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<unsigned int>::assignProductToTarget(CMSEDMWrapperLinker<unsigned int>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<int>::assignProductToTarget(CMSEDMWrapperLinker<int>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<unsigned long>::assignProductToTarget(CMSEDMWrapperLinker<unsigned long>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<long>::assignProductToTarget(CMSEDMWrapperLinker<long>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<unsigned long long>::assignProductToTarget(CMSEDMWrapperLinker<unsigned long long>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<long long>::assignProductToTarget(CMSEDMWrapperLinker<long long>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<float>::assignProductToTarget(CMSEDMWrapperLinker<float>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<double>::assignProductToTarget(CMSEDMWrapperLinker<double>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<std::string>::assignProductToTarget(CMSEDMWrapperLinker<std::string>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<TString>::assignProductToTarget(CMSEDMWrapperLinker<TString>::Wrapped_t& product){ *targetVal = product; }
template<> void CMSEDMWrapperLinker<CMSLorentzVector>::assignProductToTarget(CMSEDMWrapperLinker<CMSLorentzVector>::Wrapped_t& product){ *targetVal = product; }

template<> void CMSEDMWrapperLinker<std::vector<bool>, std::vector<bool>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<bool>, std::vector<bool>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<short>, std::vector<short>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<short>, std::vector<short>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<unsigned int>, std::vector<unsigned int>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<unsigned int>, std::vector<unsigned int>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<int>, std::vector<int>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<int>, std::vector<int>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<unsigned long>, std::vector<unsigned long>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<unsigned long>, std::vector<unsigned long>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<long>, std::vector<long>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<long>, std::vector<long>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<unsigned long long>, std::vector<unsigned long long>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<unsigned long long>, std::vector<unsigned long long>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<long long>, std::vector<long long>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<long long>, std::vector<long long>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<float>, std::vector<float>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<float>, std::vector<float>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<double>, std::vector<double>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<double>, std::vector<double>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::string>, std::vector<std::string>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::string>, std::vector<std::string>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<TString>, std::vector<TString>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<TString>, std::vector<TString>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<CMSLorentzVector>, std::vector<CMSLorentzVector>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<CMSLorentzVector>, std::vector<CMSLorentzVector>*>::Wrapped_t& product){ *targetVal = &product; }

template<> void CMSEDMWrapperLinker<std::vector<std::vector<bool>>, std::vector<std::vector<bool>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<bool>>, std::vector<std::vector<bool>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<short>>, std::vector<std::vector<short>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<short>>, std::vector<std::vector<short>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<int>>, std::vector<std::vector<int>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<int>>, std::vector<std::vector<int>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<unsigned long>>, std::vector<std::vector<unsigned long>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<unsigned long>>, std::vector<std::vector<unsigned long>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<long>>, std::vector<std::vector<long>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<long>>, std::vector<std::vector<long>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<unsigned long long>>, std::vector<std::vector<unsigned long long>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<unsigned long long>>, std::vector<std::vector<unsigned long long>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<long long>>, std::vector<std::vector<long long>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<long long>>, std::vector<std::vector<long long>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<float>>, std::vector<std::vector<float>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<float>>, std::vector<std::vector<float>>*>::Wrapped_t& product){ *targetVal = &product; }
template<> void CMSEDMWrapperLinker<std::vector<std::vector<double>>, std::vector<std::vector<double>>*>::assignProductToTarget(CMSEDMWrapperLinker<std::vector<std::vector<double>>, std::vector<std::vector<double>>*>::Wrapped_t& product){ *targetVal = &product; }


template<> void CMSEDMWrapperLinker<TBits, TBits>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: ";
  if (targetVal){ for (unsigned int ibit=0; ibit<targetVal->GetNbits(); ibit++) MELAout << targetVal->TestBitNumber(ibit); }
  else MELAout << "null";
  MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<bool, bool>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<short, short>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<unsigned int, unsigned int>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<int, int>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<unsigned long, unsigned long>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<long, long>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<unsigned long long, unsigned long long>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<long long, long long>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<float, float>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<double, double>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::string, std::string>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<TString, TString>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<CMSLorentzVector, CMSLorentzVector>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << " (address: " << targetVal << ")" << std::endl;
}

template<> void CMSEDMWrapperLinker<std::vector<bool>, std::vector<bool>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<short>, std::vector<short>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<unsigned int>, std::vector<unsigned int>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<int>, std::vector<int>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<unsigned long>, std::vector<unsigned long>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<long>, std::vector<long>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<unsigned long long>, std::vector<unsigned long long>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<long long>, std::vector<long long>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<float>, std::vector<float>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<double>, std::vector<double>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::string>, std::vector<std::string>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<TString>, std::vector<TString>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<CMSLorentzVector>, std::vector<CMSLorentzVector>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}

template<> void CMSEDMWrapperLinker<std::vector<std::vector<bool>>, std::vector<std::vector<bool>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<short>>, std::vector<std::vector<short>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<unsigned int>>, std::vector<std::vector<unsigned int>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<int>>, std::vector<std::vector<int>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<unsigned long>>, std::vector<std::vector<unsigned long>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<long>>, std::vector<std::vector<long>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<unsigned long long>>, std::vector<std::vector<unsigned long long>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<long long>>, std::vector<std::vector<long long>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<float>>, std::vector<std::vector<float>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}
template<> void CMSEDMWrapperLinker<std::vector<std::vector<double>>, std::vector<std::vector<double>>*>::print() const{
  using MELAStreamHelpers::MELAout;
  Wrapped_t const* product = nullptr; if (var) product = var->product();
  MELAout << "\t\t- edm product: "; if (product) MELAout << *product; else MELAout << "null"; MELAout << " (address: " << product << ")" << std::endl;
  MELAout << "\t\t- Target: "; if (targetVal && *targetVal) MELAout << **targetVal; else MELAout << "null"; MELAout << " (address: "; if (targetVal) MELAout << *targetVal; else MELAout << "null"; MELAout << ")" << std::endl;
}

#endif
