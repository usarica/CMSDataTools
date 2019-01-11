#include "TSystem.h"
#include "TDirectory.h"
#include "BaseEDMInputTree.h"
#include "BaseEDMInputTree.hpp"
#include "HelperFunctions.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;
using namespace HelperFunctions;


BaseEDMInputTree::BaseEDMInputTree() :
  BaseTree()
{}
BaseEDMInputTree::BaseEDMInputTree(const TString cinput, const TString treename, const TString failedtreename, const TString countersname) :
  BaseTree(cinput, treename, failedtreename, countersname)
{}
BaseEDMInputTree::~BaseEDMInputTree(){
  HelperFunctions::cleanUnorderedMap(bridgeTBitss);
  HelperFunctions::cleanUnorderedMap(bridgebools);
  HelperFunctions::cleanUnorderedMap(bridgeshorts);
  HelperFunctions::cleanUnorderedMap(bridgeuints);
  HelperFunctions::cleanUnorderedMap(bridgeints);
  HelperFunctions::cleanUnorderedMap(bridgeulongs);
  HelperFunctions::cleanUnorderedMap(bridgelongs);
  HelperFunctions::cleanUnorderedMap(bridgeulonglongs);
  HelperFunctions::cleanUnorderedMap(bridgelonglongs);
  HelperFunctions::cleanUnorderedMap(bridgefloats);
  HelperFunctions::cleanUnorderedMap(bridgedoubles);
  HelperFunctions::cleanUnorderedMap(bridgestrings);
  HelperFunctions::cleanUnorderedMap(bridgeTStrings);
  HelperFunctions::cleanUnorderedMap(bridgeCMSLorentzVectors);

  HelperFunctions::cleanUnorderedMap(bridgeVbools);
  HelperFunctions::cleanUnorderedMap(bridgeVshorts);
  HelperFunctions::cleanUnorderedMap(bridgeVuints);
  HelperFunctions::cleanUnorderedMap(bridgeVints);
  HelperFunctions::cleanUnorderedMap(bridgeVulongs);
  HelperFunctions::cleanUnorderedMap(bridgeVlongs);
  HelperFunctions::cleanUnorderedMap(bridgeVulonglongs);
  HelperFunctions::cleanUnorderedMap(bridgeVlonglongs);
  HelperFunctions::cleanUnorderedMap(bridgeVfloats);
  HelperFunctions::cleanUnorderedMap(bridgeVdoubles);
  HelperFunctions::cleanUnorderedMap(bridgeVstrings);
  HelperFunctions::cleanUnorderedMap(bridgeVTStrings);
  HelperFunctions::cleanUnorderedMap(bridgeVCMSLorentzVectors);

  HelperFunctions::cleanUnorderedMap(bridgeVVbools);
  HelperFunctions::cleanUnorderedMap(bridgeVVshorts);
  HelperFunctions::cleanUnorderedMap(bridgeVVuints);
  HelperFunctions::cleanUnorderedMap(bridgeVVints);
  HelperFunctions::cleanUnorderedMap(bridgeVVulongs);
  HelperFunctions::cleanUnorderedMap(bridgeVVlongs);
  HelperFunctions::cleanUnorderedMap(bridgeVVulonglongs);
  HelperFunctions::cleanUnorderedMap(bridgeVVlonglongs);
  HelperFunctions::cleanUnorderedMap(bridgeVVfloats);
  HelperFunctions::cleanUnorderedMap(bridgeVVdoubles);
}

void BaseEDMInputTree::synchronizeEDMBranches(){
  for (auto& it:bridgeTBitss){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgebools){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeshorts){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeuints){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeints){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeulongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgelongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeulonglongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgelonglongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgefloats){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgedoubles){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgestrings){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeTStrings){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeCMSLorentzVectors){ if (it.second) it.second->synchronize(); }

  for (auto& it:bridgeVbools){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVshorts){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVuints){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVints){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVulongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVlongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVulonglongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVlonglongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVfloats){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVdoubles){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVstrings){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVTStrings){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVCMSLorentzVectors){ if (it.second) it.second->synchronize(); }

  for (auto& it:bridgeVVbools){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVshorts){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVuints){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVints){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVulongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVlongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVulonglongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVlonglongs){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVfloats){ if (it.second) it.second->synchronize(); }
  for (auto& it:bridgeVVdoubles){ if (it.second) it.second->synchronize(); }
}

void BaseEDMInputTree::print() const{
  BaseTree::print();

  for (auto const& it:bridgeTBitss){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgebools){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeshorts){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeuints){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeints){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeulongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgelongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeulonglongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgelonglongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgefloats){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgedoubles){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgestrings){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeTStrings){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeCMSLorentzVectors){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }

  for (auto const& it:bridgeVbools){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVshorts){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVuints){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVints){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVulongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVlongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVulonglongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVlonglongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVfloats){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVdoubles){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVstrings){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVTStrings){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVCMSLorentzVectors){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }

  for (auto const& it:bridgeVVbools){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVshorts){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVuints){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVints){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVulongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVlongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVulonglongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVlonglongs){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVfloats){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
  for (auto const& it:bridgeVVdoubles){ if (it.second){ MELAout << "\t- " << it.first << " details:" << endl; it.second->print(); } }
}


bool BaseEDMInputTree::getSelectedEvent(int ev){
  bool result = BaseTree::getSelectedEvent(ev);
  synchronizeEDMBranches();
  return result;
}
bool BaseEDMInputTree::getFailedEvent(int ev){
  bool result = BaseTree::getFailedEvent(ev);
  synchronizeEDMBranches();
  return result;
}

void BaseEDMInputTree::resetBranches(){
  BaseTree::resetBranches();

  this->resetEDMBranch<BaseTree::BranchType_TBits_t>();
  this->resetEDMBranch<BaseTree::BranchType_bool_t>();
  this->resetEDMBranch<BaseTree::BranchType_short_t>();
  this->resetEDMBranch<BaseTree::BranchType_uint_t>();
  this->resetEDMBranch<BaseTree::BranchType_int_t>();
  this->resetEDMBranch<BaseTree::BranchType_ulong_t>();
  this->resetEDMBranch<BaseTree::BranchType_long_t>();
  this->resetEDMBranch<BaseTree::BranchType_ulonglong_t>();
  this->resetEDMBranch<BaseTree::BranchType_longlong_t>();
  this->resetEDMBranch<BaseTree::BranchType_float_t>();
  this->resetEDMBranch<BaseTree::BranchType_double_t>();
  this->resetEDMBranch<BaseTree::BranchType_string_t>();
  this->resetEDMBranch<BaseTree::BranchType_TString_t>();
  this->resetEDMBranch<BaseTree::BranchType_CMSLorentzVector_t>();

  this->resetEDMBranch<BaseTree::BranchType_vbool_t>();
  this->resetEDMBranch<BaseTree::BranchType_vshort_t>();
  this->resetEDMBranch<BaseTree::BranchType_vuint_t>();
  this->resetEDMBranch<BaseTree::BranchType_vint_t>();
  this->resetEDMBranch<BaseTree::BranchType_vulong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vlong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vulonglong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vlonglong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vfloat_t>();
  this->resetEDMBranch<BaseTree::BranchType_vdouble_t>();
  this->resetEDMBranch<BaseTree::BranchType_vstring_t>();
  this->resetEDMBranch<BaseTree::BranchType_vTString_t>();
  this->resetEDMBranch<BaseTree::BranchType_vCMSLorentzVector_t>();

  this->resetEDMBranch<BaseTree::BranchType_vvbool_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvshort_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvuint_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvint_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvulong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvlong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvulonglong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvlonglong_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvfloat_t>();
  this->resetEDMBranch<BaseTree::BranchType_vvdouble_t>();
}
