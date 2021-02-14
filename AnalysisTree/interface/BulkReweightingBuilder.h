#ifndef REWEIGHTINGBUILDER_H
#define REWEIGHTINGBUILDER_H

#include "BaseTree.h"
#include "ExtendedBinning.h"
#include "ReweightingFunctions.h"


class BulkReweightingBuilder{
protected:
  std::vector<BaseTree*> registeredTrees;

  ExtendedBinning binning;
  std::vector<TString> const strBinningVars;
  ReweightingFunctions::ReweightingVariableBinFunction_t rule_binningVar;
  std::unordered_map< BaseTree*, std::vector<float*> > binningVarRefs;

  std::vector<TString> const strNominalWeights;
  ReweightingFunctions::ReweightingFunction_t rule_nominalweights;
  std::unordered_map< BaseTree*, std::vector<float*> > componentRefs_nominalweights;

  std::vector<TString> const strCrossSectionWeights;
  ReweightingFunctions::ReweightingFunction_t rule_xsecweights;
  std::unordered_map< BaseTree*, std::vector<float*> > componentRefs_xsecweights;

  std::vector<std::vector<TString>> strReweightingWeightsList;
  std::vector<ReweightingFunctions::ReweightingFunction_t> rule_reweightingweights_list;
  std::vector< std::pair<double, double> > reweightingweights_frac_tolerance_pair_list;
  std::unordered_map< BaseTree*, std::vector<std::vector<float*> > > componentRefsList_reweightingweights;

  std::unordered_map<BaseTree*, double> normWeights;

  // Derived variables
  std::unordered_map< BaseTree*, std::vector< std::vector<float> > > absWeightThresholdsPerBinList;
  std::unordered_map< BaseTree*, std::vector<double> > sum_normwgts_all;
  std::unordered_map< BaseTree*, std::vector<double> > sum_normwgts_nonzerorewgt;
  std::unordered_map< BaseTree*, std::vector< std::vector< std::pair<double, double> > > > sum_wgts_withrewgt;
  std::unordered_map< BaseTree*, std::vector<double> > NeffsPerBin;
  std::unordered_map< BaseTree*, std::vector<double> > sampleNormalization;
  std::unordered_map< BaseTree*, std::vector<double> > sampleZeroMECompensation;
  std::unordered_map< BaseTree*, double > samplePairwiseNormalization;

public:
  BulkReweightingBuilder(
    ExtendedBinning const& binning_,
    std::vector<TString> const& strBinningVars_,
    std::vector<TString> const& strNominalWeights_,
    std::vector<TString> const& strCrossSectionWeights_,
    ReweightingFunctions::ReweightingVariableBinFunction_t rule_binningVar_,
    ReweightingFunctions::ReweightingFunction_t rule_nominalweights_,
    ReweightingFunctions::ReweightingFunction_t rule_xsecweights_
  );
  virtual ~BulkReweightingBuilder(){}

  void addReweightingWeights(
    std::vector<TString> const& strReweightingWeights_,
    ReweightingFunctions::ReweightingFunction_t rule_reweightingweights_,
    double thr_wgt=0.9995, double tolerance=5.
  );
  void registerTree(BaseTree* tree, double extNorm=1);

  // thr_wgt and tol_wgt adjust the large weight removal threshold and tolerance.
  // thr_wgt=-1 defaults to the default of ReweightingFunctions::getAbsWeightThresholdsPerBinByFixedFractionalThreshold.
  // thr_frac_Neff with value >0 sets the contribution of samples with Neff fraction<thr_frac_Neff to 0.
  // This is primarily to ensure that narrow-width samples cannot enter into bins far away from their pole.
  void setup(
    int ihypo_Neff, std::vector<std::pair<BaseTree*, BaseTree*>> const* tree_normTree_pairs=nullptr,
    float thr_frac_Neff=-1.
  );
  // This function sets the vectors by reading from a file
  void setupFromFile(TString cinput);

  double getOverallReweightingNormalization(BaseTree* tree) const;
  double getSamplePairwiseNormalization(BaseTree* tree) const;
  bool checkWeightsBelowThreshold(BaseTree* tree) const;

  // A few get functions to acquire the actual set of variables used
  std::vector<TString> const& getBinningVars() const{ return strBinningVars; }
  std::vector<TString> const& getNominalWeightVars() const{ return strNominalWeights; }
  std::vector<TString> const& getCrossSectionWeightVars() const{ return strCrossSectionWeights; }
  std::vector<std::vector<TString>> const& getReweightingWeightVarList() const{ return strReweightingWeightsList; }

  // This is to write the reweigthing specifications to a file
  void writeToFile(TFile* foutput) const;

  // Print the existing reweighting information
  void print() const;

};


#endif
