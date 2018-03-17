#include <algorithm>
#include <utility>
#include <iterator>
#include "BaseTreeLooper.h"
#include "BaseTreeLooper.hpp"


using namespace std;


BaseTreeLooper::BaseTreeLooper() : sampleIdOpt(BaseTreeLooper::kNoStorage), maxNEvents(-1) { setExternalProductList(); setExternalProductTree(); }
BaseTreeLooper::BaseTreeLooper(CJLSTTree* inTree) : sampleIdOpt(BaseTreeLooper::kNoStorage), maxNEvents(-1) { this->addTree(inTree); setExternalProductList(); setExternalProductTree(); }
BaseTreeLooper::BaseTreeLooper(std::vector<CJLSTTree*> const& inTreeList) :
  sampleIdOpt(BaseTreeLooper::kNoStorage),
  treeList(inTreeList),
  maxNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::BaseTreeLooper(CJLSTSet const* inTreeSet) :
  sampleIdOpt(BaseTreeLooper::kNoStorage),
  treeList(inTreeSet->getCJLSTTreeList()),
  maxNEvents(-1)
{
  setExternalProductList();
  setExternalProductTree();
}
BaseTreeLooper::~BaseTreeLooper(){}

void BaseTreeLooper::addTree(CJLSTTree* tree){ this->treeList.push_back(tree); }

void BaseTreeLooper::addDiscriminantBuilder(TString KDname, Discriminant* KDbuilder, std::vector<TString> const& KDvars){
  if (!KDbuilder) return;
  if (KDbuilders.find(KDname)!=KDbuilders.end()) MELAerr << "BaseTreeLooper::addDiscriminantBuilder: " << KDname << " already exists." << endl;
  else{
    KDbuilders[KDname] = std::pair<Discriminant*, std::vector<TString>>(KDbuilder, KDvars);
    for (auto const& v:KDvars) this->addConsumed<float>(v);
  }
}
void BaseTreeLooper::addReweightingBuilder(TString rewgtname, ReweightingBuilder* Rewgtbuilder){
  if (!Rewgtbuilder) return;
  if (Rewgtbuilders.find(rewgtname)!=Rewgtbuilders.end()) MELAerr << "BaseTreeLooper::addReweightingBuilder: " << rewgtname << " already exists but will override it regardless." << endl;
  Rewgtbuilders[rewgtname] = Rewgtbuilder;
}
void BaseTreeLooper::addZXFakeRateHandler(TString frname, ZXFakeRateHandler* ZXFakeRateHandler){
  if (!ZXFakeRateHandler) return;
  if (ZXFakeRateHandlers.find(frname)!=ZXFakeRateHandlers.end()) MELAerr << "BaseTreeLooper::addZXFakeRateHandler: " << frname << " already exists but will override it regardless." << endl;
  ZXFakeRateHandlers[frname] = ZXFakeRateHandler;
}
void BaseTreeLooper::addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&)){
  if (!fcn) return;
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) MELAerr << "BaseTreeLooper::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}
void BaseTreeLooper::addSystematic(TString systname, SystematicsHelpers::SystematicsClass* syst){
  if (!syst) return;
  if (SystVariations.find(systname)!=SystVariations.end()) MELAerr << "BaseTreeLooper::addSystematic: " << systname << " already exists." << endl;
  else SystVariations[systname] = syst;
}

void BaseTreeLooper::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef=extProductListRef;
  else this->productListRef=&(this->productList);
}

void BaseTreeLooper::setExternalProductTree(BaseTree* extTree){
  this->productTree=extTree;
  this->productListRef=&(this->productList); // To make sure product list collects some events before flushing
}

void BaseTreeLooper::setMaximumEvents(int n){ maxNEvents=n; }

void BaseTreeLooper::setSampleIdStorageOption(BaseTreeLooper::SampleIdStorageType opt){ sampleIdOpt=opt; }

void BaseTreeLooper::addProduct(SimpleEntry& product, unsigned int* ev_rec){
  this->productListRef->push_back(product);
  if (ev_rec) (*ev_rec)++;
}

void BaseTreeLooper::recordProductsToTree(){
  if (!this->productTree) return;
  BaseTree::writeSimpleEntries(this->productListRef->begin(), this->productListRef->end(), this->productTree);
  this->clearProducts();
}

bool BaseTreeLooper::linkConsumes(CJLSTTree* tree){
  bool process = tree->isValid();
  if (process){
    process &= this->linkConsumed<short>(tree);
    process &= this->linkConsumed<unsigned int>(tree);
    process &= this->linkConsumed<int>(tree);
    process &= this->linkConsumed<float>(tree);
    process &= this->linkConsumed<double>(tree);
    process &= this->linkConsumed<std::vector<short>>(tree);
    process &= this->linkConsumed<std::vector<unsigned int>>(tree);
    process &= this->linkConsumed<std::vector<int>>(tree);
    process &= this->linkConsumed<std::vector<float>>(tree);
    process &= this->linkConsumed<std::vector<double>>(tree);
    // Silence unused branches
    tree->silenceUnused();
  }
  if (!process) MELAerr << "BaseTreeLooper::linkConsumes: Linking failed for some reason for tree " << tree->sampleIdentifier << endl;
  return process;
}

void BaseTreeLooper::loop(bool loopSelected, bool loopFailed, bool keepProducts){
  // Loop over the trees
  unsigned int ev_acc=0;
  unsigned int ev_rec=0;
  const bool storeSampleIdByMH = (sampleIdOpt==kStoreByMH);
  const bool storeSampleIdByHashVal = (sampleIdOpt==kStoreByHashVal);
  for (CJLSTTree*& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->linkConsumes(tree))) continue;

    float wgtExternal = 1;
    CJLSTSet const* associatedSet = tree->getAssociatedSet();
    if (associatedSet) wgtExternal *= associatedSet->getPermanentWeight(tree);
    if (wgtExternal==0.){
      MELAerr << "BaseTreeLooper::loop: External weights are 0 for the " << tree->sampleIdentifier << " sample. Skipping..." << endl;
      continue;
    }

    size_t sampleId=0;
    if (storeSampleIdByMH) sampleId = tree->MHVal;
    else if (storeSampleIdByHashVal){ std::hash<TString> tmphash; sampleId=tmphash(tree->sampleIdentifier); }

    // Loop over selected events
    if (loopSelected){
      MELAout << "BaseTreeLooper::loop: Looping over " << tree->sampleIdentifier << " selected events" << endl;
      int ev=0;
      const int nevents = tree->getSelectedNEvents();
      while (tree->getSelectedEvent(ev)){
        if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;
        SimpleEntry product;
        if (tree->isValidEvent()){
          if (this->runEvent(tree, wgtExternal, product)){
            if (keepProducts){
              if (storeSampleIdByMH || storeSampleIdByHashVal){
                product.setNamedVal("SampleId", sampleId);
                product.setNamedVal("EventNumber", ev_acc);
              }
              this->addProduct(product, &ev_rec);
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
      }
    }
    // Loop over failed events
    if (loopFailed){
      MELAout << "BaseTreeLooper::loop: Looping over " << tree->sampleIdentifier << " failed events" << endl;
      int ev=0;
      const int nevents = tree->getFailedNEvents();
      while (tree->getFailedEvent(ev)){
        if (maxNEvents>=0 && (int) ev_rec==maxNEvents) break;
        SimpleEntry product;
        if (tree->isValidEvent()){
          if (this->runEvent(tree, wgtExternal, product)){
            if (keepProducts){
              if (storeSampleIdByMH || storeSampleIdByHashVal){
                product.setNamedVal("SampleId", sampleId);
                product.setNamedVal("EventNumber", ev_acc);
              }
              this->addProduct(product, &ev_rec);
            }
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
      }
    }

    // Record products to external tree
    this->recordProductsToTree();

  } // End loop over the trees
  MELAout << "BaseTreeLooper::loop: Total number of products: " << ev_rec << " / " << ev_acc << endl;
}

std::vector<SimpleEntry> const& BaseTreeLooper::getProducts() const{ return *productListRef; }

void BaseTreeLooper::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "BaseTreeLooper::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "BaseTreeLooper::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void BaseTreeLooper::clearProducts(){ std::vector<SimpleEntry> emptyList; std::swap(emptyList, *productListRef); }
