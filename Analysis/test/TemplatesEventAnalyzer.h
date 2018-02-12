#ifndef TEMPLATESEVENTANALYZER_H
#define TEMPLATESEVENTANALYZER_H

#include "common_includes.h"

// gg analyzer
class TemplatesEventAnalyzer : public BaseTreeLooper{
protected:
  Channel channel;
  Category category;

  bool recordCategorizationKDs;
  bool recordKDVariables;

  std::vector<std::pair<float, float>> mass_boundaries;

  bool runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  TemplatesEventAnalyzer(Channel channel_, Category category_);
  TemplatesEventAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_);
  TemplatesEventAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_);
  TemplatesEventAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_);

  void setRecordCategorizationKDs(bool flag){ recordCategorizationKDs=flag; }
  void setRecordKDVariables(bool flag){ recordKDVariables=flag; }

  void addMassWindow(std::pair<float, float> const boundaries){ mass_boundaries.push_back(boundaries); }
};

TemplatesEventAnalyzer::TemplatesEventAnalyzer(Channel channel_, Category category_) :
  BaseTreeLooper(), channel(channel_), category(category_), recordCategorizationKDs(false), recordKDVariables(false)
{}
TemplatesEventAnalyzer::TemplatesEventAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_) :
  BaseTreeLooper(inTree), channel(channel_), category(category_), recordCategorizationKDs(false), recordKDVariables(false)
{}
TemplatesEventAnalyzer::TemplatesEventAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_) :
  BaseTreeLooper(inTreeList), channel(channel_), category(category_), recordCategorizationKDs(false), recordKDVariables(false)
{}
TemplatesEventAnalyzer::TemplatesEventAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_) :
  BaseTreeLooper(inTreeSet), channel(channel_), category(category_), recordCategorizationKDs(false), recordKDVariables(false)
{}

bool TemplatesEventAnalyzer::runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product){
  bool validProducts=(tree!=nullptr);
  if (validProducts){
    // Get tree and binning information
    //product.setNamedVal("MH", tree->MHVal);

    // Get main observables
    float& ZZMass = *(valfloats["ZZMass"]);
    float& GenHMass = *(valfloats["GenHMass"]);
    product.setNamedVal("ZZMass", ZZMass);
    product.setNamedVal("GenHMass", GenHMass);

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
    wgt *= (*(valfloats["dataMCWeight"]))*(*(valfloats["trigEffWeight"]))*(*(valfloats["PUWeight"]))*(*(valfloats["genHEPMCweight"]));
    for (auto rewgt_it=Rewgtbuilders.cbegin(); rewgt_it!=Rewgtbuilders.cend(); rewgt_it++){
      auto& rewgtBuilder = rewgt_it->second;
      if (rewgt_it->first=="MELARewgt"){
        float mela_wgt_sum = rewgtBuilder->getSumPostThresholdWeights(tree);
        float mela_wgt = (mela_wgt_sum!=0. ? rewgtBuilder->getPostThresholdWeight(tree)/mela_wgt_sum : 0.); // Normalized to unit
        unsigned int mela_nevts = rewgtBuilder->getSumNonZeroWgtEvents(tree);
        unsigned int mela_sumnevts = rewgtBuilder->getSumAllNonZeroWgtEvents(tree);
        if (mela_sumnevts!=0) mela_wgt *= static_cast<float>(mela_nevts) / static_cast<float>(mela_sumnevts);
        mela_wgt *= rewgtBuilder->getNormComponent(tree);
        wgt *= mela_wgt;
        //product.setNamedVal("MELARewgtWeight", mela_wgt);
        //product.setNamedVal("MELARewgtBin", rewgtBuilder->findBin(tree));
      }
      else wgt *= rewgtBuilder->getPostThresholdWeight(tree);
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
          if (it->first.Contains(ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac))){
            hypo=(ACHypothesisHelpers::ACHypothesis)iac;
            break;
          }
        }
        DjjVBF[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else if (it->first.Contains("DjjZH")){
        ACHypothesisHelpers::ACHypothesis hypo=kSM;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
          if (it->first.Contains(ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac))){
            hypo=(ACHypothesisHelpers::ACHypothesis)iac;
            break;
          }
        }
        DjjZH[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else if (it->first.Contains("DjjWH")){
        ACHypothesisHelpers::ACHypothesis hypo=kSM;
        for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
          if (it->first.Contains(ACHypothesisHelpers::getACHypothesisName((ACHypothesisHelpers::ACHypothesis)iac))){
            hypo=(ACHypothesisHelpers::ACHypothesis)iac;
            break;
          }
        }
        DjjWH[hypo]=KD;
        if (recordCategorizationKDs) product.setNamedVal(it->first, KD);
      }
      else{
        product.setNamedVal(it->first, KD);
        validProducts &= (KD != float(-999.));
      }
      //product.setNamedVal(it->first, KD);
    }
    if (!validProducts) return validProducts;

    // Category check
    bool fitsAtLeastOneCategory=(category==Inclusive);
    if (!fitsAtLeastOneCategory){
      bool isRequestedCategory[ACHypothesisHelpers::nACHypotheses]={ false };
      for (int iac=0; iac<(int) ACHypothesisHelpers::nACHypotheses; iac++){
        if (iac!=(int) ACHypothesisHelpers::kSM){
          DjjVBF[iac]=std::max(DjjVBF[iac], DjjVBF[kSM]);
          DjjZH[iac]=std::max(DjjZH[iac], DjjZH[kSM]);
          DjjWH[iac]=std::max(DjjWH[iac], DjjWH[kSM]);
        }
        Category catFound = CategorizationHelpers::getCategory(DjjVBF[iac], DjjZH[iac], DjjWH[iac], false);
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
    validProducts &= SampleHelpers::testChannel(channel, *(valshorts["Z1Flav"]), *(valshorts["Z2Flav"]));
    if (!validProducts) return validProducts;
  }

  return validProducts;
}

#endif
