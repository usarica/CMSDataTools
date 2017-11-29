#ifndef GGANALYZER_H
#define GGANALYZER_H

#include "common_includes.h"

// gg analyzer
class GGAnalyzer : public BaseTreeLooper{
protected:
  Channel channel;
  Category category;

  bool runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  GGAnalyzer(Channel channel_, Category category_) : BaseTreeLooper(), channel(channel_), category(category_){}
  GGAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_) : BaseTreeLooper(inTree), channel(channel_), category(category_){}
  GGAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_) : BaseTreeLooper(inTreeList), channel(channel_), category(category_){}
  GGAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_) : BaseTreeLooper(inTreeSet), channel(channel_), category(category_){}

};

bool GGAnalyzer::runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product){
  bool validProducts=(tree!=nullptr);
  if (validProducts){
    // Get tree and binning information
    product.setNamedVal("MH", tree->MHVal);

    // Get main observables
    float& ZZMass = *(valfloats["ZZMass"]);
    float& GenHMass = *(valfloats["GenHMass"]);
    product.setNamedVal("ZZMass", ZZMass);
    product.setNamedVal("GenHMass", GenHMass);

    // Construct the weights
    float wgt = externalWgt;
    bool hasPUGenHEPRewgt=false;
    for (auto rewgt_it=Rewgtbuilders.cbegin(); rewgt_it!=Rewgtbuilders.cend(); rewgt_it++){
      auto& rewgtBuilder = rewgt_it->second;
      if (rewgt_it->first=="PUGenHEPRewgt"){
        float pugenhep_wgt_sum = rewgtBuilder->getSumPostThresholdWeights(tree);
        float pugenhep_wgt = (pugenhep_wgt_sum!=0. ? rewgtBuilder->getPostThresholdWeight(tree)/pugenhep_wgt_sum : 0.); // Normalized to unit
        wgt *= pugenhep_wgt;
        hasPUGenHEPRewgt=true;
      }
      else if (rewgt_it->first=="MELARewgt"){
        float mela_wgt_sum = rewgtBuilder->getSumPostThresholdWeights(tree);
        float mela_wgt = (mela_wgt_sum!=0. ? rewgtBuilder->getPostThresholdWeight(tree)/mela_wgt_sum : 0.); // Normalized to unit
        mela_wgt *= rewgtBuilder->getNormComponent(tree);
        wgt *= mela_wgt;
        product.setNamedVal("MELARewgtBin", rewgtBuilder->findBin(tree));
      }
    }
    wgt *= (*(valfloats["dataMCWeight"]))*(*(valfloats["trigEffWeight"]));
    if (!hasPUGenHEPRewgt) wgt *= (*(valfloats["PUWeight"]))*(*(valfloats["genHEPMCweight"]));

    product.setNamedVal("weight", wgt);
    if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.){
      if (wgt!=0.){
        MELAerr << "GGAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded at mass " << ZZMass << " for tree " << tree->sampleIdentifier << "." << endl;
        exit(1);
      }
      validProducts=false;
    }

    // Compute the KDs
    // Reserve the special DjjVBF, DjjZH and DjjWH discriminants
    float DjjVBF=-1;
    float DjjWH=-1;
    float DjjZH=-1;
    for (auto it=KDbuilders.cbegin(); it!=KDbuilders.cend(); it++){
      auto& KDbuilderpair = it->second;
      auto& KDbuilder = KDbuilderpair.first;
      auto& strKDVarsList = KDbuilderpair.second;
      vector<float> KDBuildVals; KDBuildVals.reserve(strKDVarsList.size());
      for (auto const& s : strKDVarsList) KDBuildVals.push_back(*(valfloats[s]));
      float KD = KDbuilder->update(KDBuildVals, ZZMass);
      validProducts &= !(std::isnan(KD) || std::isinf(KD));

      if (it->first=="DjjVBF") DjjVBF=KD;
      else if (it->first=="DjjZH") DjjZH=KD;
      else if (it->first=="DjjWH") DjjWH=KD;
      else{
        product.setNamedVal(it->first, KD);
        validProducts &= (KD != float(-999.));
      }
    }

    // Category check
    Category catFound = CategorizationHelpers::getCategory(DjjVBF, DjjZH, DjjWH, false);
    validProducts &= (category==Inclusive || category==catFound);

    // Channel check
    validProducts &= SampleHelpers::testChannel(channel, *(valshorts["Z1Flav"]), *(valshorts["Z2Flav"]));
  }

  return validProducts;
}


#endif
