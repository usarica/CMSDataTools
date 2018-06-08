#ifndef TEMPLATESEVENTANALYZER_H
#define TEMPLATESEVENTANALYZER_H

#include "common_includes.h"

// Analyzer for template production
class TemplatesEventAnalyzer : public BaseTreeLooper{
protected:
  Channel channel;
  Category category;

  bool allowSSChannel;
  bool recordCategorizationKDs;
  bool recordKDVariables;

  std::vector<std::pair<float, float>> mass_boundaries;

  bool runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  TemplatesEventAnalyzer(Channel channel_, Category category_);
  TemplatesEventAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_);
  TemplatesEventAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_);
  TemplatesEventAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_);

  // Record categorization discriminants
  void setAllowSSChannel(bool flag){ allowSSChannel=flag; }
  // Record categorization discriminants
  void setRecordCategorizationKDs(bool flag){ recordCategorizationKDs=flag; }
  // Record discriminant input MEs
  void setRecordKDVariables(bool flag){ recordKDVariables=flag; }
  // Reco. mass windows
  void addMassWindow(std::pair<float, float> const boundaries){ mass_boundaries.push_back(boundaries); }
};

TemplatesEventAnalyzer::TemplatesEventAnalyzer(Channel channel_, Category category_) :
  BaseTreeLooper(), channel(channel_), category(category_), allowSSChannel(false), recordCategorizationKDs(false), recordKDVariables(false)
{}
TemplatesEventAnalyzer::TemplatesEventAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_) :
  BaseTreeLooper(inTree), channel(channel_), category(category_), allowSSChannel(false), recordCategorizationKDs(false), recordKDVariables(false)
{}
TemplatesEventAnalyzer::TemplatesEventAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_) :
  BaseTreeLooper(inTreeList), channel(channel_), category(category_), allowSSChannel(false), recordCategorizationKDs(false), recordKDVariables(false)
{}
TemplatesEventAnalyzer::TemplatesEventAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_) :
  BaseTreeLooper(inTreeSet), channel(channel_), category(category_), allowSSChannel(false), recordCategorizationKDs(false), recordKDVariables(false)
{}

bool TemplatesEventAnalyzer::runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product){
  bool validProducts=(tree!=nullptr);
  if (validProducts){
    // Get tree and binning information
    //product.setNamedVal("MH", tree->MHVal);

    // Get main observables
    float& ZZMass = *(valfloats["ZZMass"]);
    product.setNamedVal("ZZMass", ZZMass);
    unordered_map<TString, float*>::const_iterator it_GenHMass;
    if (HelperFunctions::getUnorderedMapIterator("GenHMass", valfloats, it_GenHMass)){
      float const& GenHMass = *(it_GenHMass->second);
      product.setNamedVal("GenHMass", GenHMass);
    }

    if (!mass_boundaries.empty()){
      bool mass_range_found=false;
      for (std::pair<float, float>& mass_boundary:mass_boundaries){
        if (
          (mass_boundary.first<0. || ZZMass>=mass_boundary.first)
          &&
          (mass_boundary.second<0. || ZZMass<mass_boundary.second)
          ){
          mass_range_found=true;
          break;
        }
      }
      validProducts &= mass_range_found;
    }
    if (!validProducts) return validProducts;

    // Construct the weights
    float wgt = externalWgt;
    unordered_map<TString, float*>::const_iterator it_genHEPMCweight;
    if (HelperFunctions::getUnorderedMapIterator("genHEPMCweight", valfloats, it_genHEPMCweight)) // Probably registered all others as well
      wgt *= (*(valfloats["dataMCWeight"]))*(*(valfloats["trigEffWeight"]))*(*(valfloats["PUWeight"]))*(*(it_genHEPMCweight->second));
    for (auto rewgt_it=Rewgtbuilders.cbegin(); rewgt_it!=Rewgtbuilders.cend(); rewgt_it++){
      auto& rewgtBuilder = rewgt_it->second;
      if (rewgt_it->first=="MELARewgt"){
        float mela_wgt_sum = rewgtBuilder->getSumPostThresholdWeights(tree);
        float mela_wgt = (mela_wgt_sum!=0. ? rewgtBuilder->getPostThresholdWeight(tree)/mela_wgt_sum : 0.); // Normalized to unit
        float mela_samplewgt, mela_sumsamplewgts;
        if (ReweightingBuilder::useNeffInNormComponent){
          mela_samplewgt = rewgtBuilder->getSumPostThresholdNeffs(tree);
          mela_sumsamplewgts = rewgtBuilder->getSumAllPostThresholdNeffs(tree);
        }
        else{
          mela_samplewgt = rewgtBuilder->getSumPostThresholdSqWeightInvs(tree);
          mela_sumsamplewgts = rewgtBuilder->getSumAllPostThresholdSqWeightInvs(tree);
        }
        if (mela_sumsamplewgts!=0.) mela_wgt *= mela_samplewgt / mela_sumsamplewgts;
        //unsigned int mela_nevts = rewgtBuilder->getSumNonZeroWgtEvents(tree);
        //unsigned int mela_sumnevts = rewgtBuilder->getSumAllNonZeroWgtEvents(tree);
        //if (mela_sumnevts!=0) mela_wgt *= static_cast<float>(mela_nevts) / static_cast<float>(mela_sumnevts);
        mela_wgt *= rewgtBuilder->getNormComponent(tree);
        wgt *= mela_wgt;
        //product.setNamedVal("MELARewgtWeight", mela_wgt);
        //product.setNamedVal("MELARewgtBin", rewgtBuilder->findBin(tree));
      }
      else wgt *= rewgtBuilder->getPostThresholdWeight(tree);
    }
    for (auto zxfr_it=ZXFakeRateHandlers.cbegin(); zxfr_it!=ZXFakeRateHandlers.cend(); zxfr_it++){
      auto& ZXFRHandle = zxfr_it->second;
      float frwgt = ZXFRHandle->getFakeRateWeight(tree);
      wgt *= frwgt;
    }
    for (auto syst_it=SystVariations.cbegin(); syst_it!=SystVariations.cend(); syst_it++){
      auto& systVar = syst_it->second;
      float systWgt = systVar->eval(tree);
      wgt *= systWgt;
      product.setNamedVal(syst_it->first, systWgt);
    }
    product.setNamedVal("weight", wgt);
    if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.){
      if (wgt!=0.){
        MELAerr << "TemplatesEventAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded at mass " << ZZMass << " for tree " << tree->sampleIdentifier << "." << endl;
        exit(1);
      }
      validProducts=false;
    }
    if (!validProducts) return validProducts;

    // Compute the KDs
    // Reserve the special DjjVBF, DjjZH and DjjWH discriminants
    float DjjVBF[nACHypotheses];
    float DjVBF[nACHypotheses];
    float DjjWH[nACHypotheses];
    float DjjZH[nACHypotheses];
    for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
      DjjVBF[iac]=-1;
      DjjZH[iac]=-1;
      DjjWH[iac]=-1;
    }
    for (auto it=KDbuilders.cbegin(); it!=KDbuilders.cend(); it++){
      auto& KDbuilderpair = it->second;
      auto& KDbuilder = KDbuilderpair.first;
      auto& strKDVarsList = KDbuilderpair.second;
      vector<float> KDBuildVals; KDBuildVals.reserve(strKDVarsList.size());
      for (auto const& s:strKDVarsList){
        KDBuildVals.push_back(*(valfloats[s]));
        if (recordKDVariables) product.setNamedVal(s, KDBuildVals.back());
      }
      float KD = KDbuilder->update(KDBuildVals, ZZMass);
      validProducts &= !(std::isnan(KD) || std::isinf(KD));

      if (it->first.Contains("DjjVBF")){
        ACHypothesisHelpers::ACHypothesis hypo=kSM;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
          TString strComp="DjjVBF"+ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac);
          if (it->first==strComp){ hypo=(ACHypothesisHelpers::ACHypothesis)iac; break; }
        }
        DjjVBF[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else if (it->first.Contains("DjjZH")){
        ACHypothesisHelpers::ACHypothesis hypo=kSM;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
          TString strComp="DjjZH"+ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac);
          if (it->first==strComp){ hypo=(ACHypothesisHelpers::ACHypothesis)iac; break; }
        }
        DjjZH[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else if (it->first.Contains("DjjWH")){
        ACHypothesisHelpers::ACHypothesis hypo=kSM;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
          TString strComp="DjjWH"+ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac);
          if (it->first==strComp){ hypo=(ACHypothesisHelpers::ACHypothesis)iac; break; }
        }
        DjjWH[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else if (it->first.Contains("DjVBF")){
        ACHypothesisHelpers::ACHypothesis hypo=kSM;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
          TString strComp="DjVBF"+ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac);
          if (it->first==strComp){ hypo=(ACHypothesisHelpers::ACHypothesis)iac; break; }
        }
        DjVBF[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else{
        product.setNamedVal(it->first, KD);
        validProducts &= (KD != float(-999.) || it->first.Contains("m4l"));
      }
      //product.setNamedVal(it->first, KD);
    }
    DjjWH[kL1ZGs] = DjjWH[kSM];
    if (!validProducts) return validProducts;

    // Category check
    SimpleEntry catvars;
    bool fitsAtLeastOneCategory=(category==Inclusive);
    if (!fitsAtLeastOneCategory){
      if (globalCategorizationScheme==UntaggedOrJJVBFOrHadVH_WithMultiplicityAndBTag){
        for (unordered_map<TString, short*>::const_iterator it = valshorts.cbegin(); it!=valshorts.cend(); it++){
          if (it->first.Contains("nCleanedJetsPt30BTagged")) catvars.setNamedVal("nJets_BTagged", *(it->second));
          else if (it->first.Contains("nCleanedJetsPt30")) catvars.setNamedVal("nJets", *(it->second));
          else if (it->first.Contains("nExtraLep")) catvars.setNamedVal("nExtraLep", *(it->second));
        }
      }
      bool isRequestedCategory[ACHypothesisHelpers::nACHypotheses]={ false };
      for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
        if (iac!=(int) ACHypothesisHelpers::kSM){
          DjjVBF[iac]=std::max(DjjVBF[iac], DjjVBF[kSM]);
          DjjZH[iac]=std::max(DjjZH[iac], DjjZH[kSM]);
          DjjWH[iac]=std::max(DjjWH[iac], DjjWH[kSM]);
        }
        catvars.setNamedVal("DjjVBF", DjjVBF[iac]);
        catvars.setNamedVal("DjjZH", DjjZH[iac]);
        catvars.setNamedVal("DjjWH", DjjWH[iac]);
        Category catFound = CategorizationHelpers::getCategory(catvars, false);
        isRequestedCategory[iac] = (category==catFound);
        TString catFlagName = TString("is_")
          + CategorizationHelpers::getCategoryName(category)
          + TString("_")
          + ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac);
        product.setNamedVal(catFlagName, isRequestedCategory[iac]);
        fitsAtLeastOneCategory |= isRequestedCategory[iac];
      }
    }
    validProducts &= fitsAtLeastOneCategory;
    if (!validProducts) return validProducts;

    // Channel check
    validProducts &= SampleHelpers::testChannel(channel, *(valshorts["Z1Flav"]), *(valshorts["Z2Flav"]), allowSSChannel);
    if (!validProducts) return validProducts;
  }

  return validProducts;
}

#endif
