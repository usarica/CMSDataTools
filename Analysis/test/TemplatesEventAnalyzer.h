#ifndef TEMPLATESEVENTANALYZER_H
#define TEMPLATESEVENTANALYZER_H

#include "common_includes.h"

// gg analyzer
class TemplatesEventAnalyzer : public BaseTreeLooper{
protected:
  Channel channel;
  Category category;

  bool runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product);

public:
  TemplatesEventAnalyzer(Channel channel_, Category category_) : BaseTreeLooper(), channel(channel_), category(category_){}
  TemplatesEventAnalyzer(CJLSTTree* inTree, Channel channel_, Category category_) : BaseTreeLooper(inTree), channel(channel_), category(category_){}
  TemplatesEventAnalyzer(std::vector<CJLSTTree*> const& inTreeList, Channel channel_, Category category_) : BaseTreeLooper(inTreeList), channel(channel_), category(category_){}
  TemplatesEventAnalyzer(CJLSTSet const* inTreeSet, Channel channel_, Category category_) : BaseTreeLooper(inTreeSet), channel(channel_), category(category_){}

};

bool TemplatesEventAnalyzer::runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product){
  bool validProducts=(tree!=nullptr);
  if (validProducts){
    // Get tree and binning information
    //product.setNamedVal("MH", tree->MHVal);

    // Get main observables
    float& ZZMass = *(valfloats["ZZMass"]);
    float& GenHMass = *(valfloats["GenHMass"]);
    product.setNamedVal("ZZMass", ZZMass);
    //product.setNamedVal("GenHMass", GenHMass);

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
    product.setNamedVal("weight", wgt);
    if (std::isnan(wgt) || std::isinf(wgt) || wgt==0.){
      if (wgt!=0.){
        MELAerr << "TemplatesEventAnalyzer::runEvent: Invalid weight " << wgt << " is being discarded at mass " << ZZMass << " for tree " << tree->sampleIdentifier << "." << endl;
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
      for (auto const& s:strKDVarsList) KDBuildVals.push_back(*(valfloats[s]));
      float KD = KDbuilder->update(KDBuildVals, ZZMass);
      validProducts &= !(std::isnan(KD) || std::isinf(KD));

      if (it->first=="DjjVBF" || it->first=="DjjVBFBSM") DjjVBF=std::max(KD, DjjVBF);
      else if (it->first=="DjjZH" || it->first=="DjjZHBSM") DjjZH=std::max(KD, DjjZH);
      else if (it->first=="DjjWH" || it->first=="DjjWHBSM") DjjWH=std::max(KD, DjjWH);
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
