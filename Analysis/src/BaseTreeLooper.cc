#include <algorithm>
#include <utility>
#include <iterator>
#include "BaseTreeLooper.h"
#include "BaseTreeLooper.hpp"


using namespace std;


BaseTreeLooper::BaseTreeLooper(){ setExternalProductList(); }
BaseTreeLooper::BaseTreeLooper(CJLSTTree* inTree){ this->addTree(inTree); setExternalProductList(); }
BaseTreeLooper::BaseTreeLooper(std::vector<CJLSTTree*> const& inTreeList) :
  treeList(inTreeList)
{
  setExternalProductList();
}
BaseTreeLooper::BaseTreeLooper(CJLSTSet const* inTreeSet) :
  treeList(inTreeSet->getCJLSTTreeList())
{
  setExternalProductList();
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
void BaseTreeLooper::addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&)){
  if (!fcn) return;
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) MELAerr << "BaseTreeLooper::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}

void BaseTreeLooper::setExternalProductList(std::vector<SimpleEntry>* extProductListRef){
  if (extProductListRef) this->productListRef=extProductListRef;
  else this->productListRef=&(this->productList);
}

void BaseTreeLooper::addProduct(SimpleEntry& product){ this->productListRef->push_back(product); }

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
  for (CJLSTTree*& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->linkConsumes(tree))) continue;

    // Loop over selected events
    if (loopSelected){
      MELAout << "BaseTreeLooper::loop: Looping over " << tree->sampleIdentifier << " selected events" << endl;
      int ev=0;
      const int nevents = tree->getSelectedNEvents();
      while (tree->getSelectedEvent(ev)){
        SimpleEntry product;
        if (tree->isValidEvent()){
          if (this->runEvent(tree, product)){
            if (keepProducts) this->addProduct(product);
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
        SimpleEntry product;
        if (tree->isValidEvent()){
          if (this->runEvent(tree, product)){
            if (keepProducts) this->addProduct(product);
          }
        }
        HelperFunctions::progressbar(ev, nevents);
        ev++; ev_acc++;
      }
    }

  } // End loop over the trees
  MELAout << "BaseTreeLooper::loop: Total number of products accumulated over " << ev_acc << " events considered: " << productListRef->size() << endl;
}

std::vector<SimpleEntry> const& BaseTreeLooper::getProducts() const{ return *productListRef; }

void BaseTreeLooper::moveProducts(std::vector<SimpleEntry>& targetColl){
  MELAout << "BaseTreeLooper::moveProducts: Moving " << productListRef->size() << " products into a list of initial size " << targetColl.size() << endl;
  std::move(productListRef->begin(), productListRef->end(), std::back_inserter(targetColl));
  clearProducts();
  MELAout << "BaseTreeLooper::moveProducts: Target list final size: " << targetColl.size() << endl;
}

void BaseTreeLooper::clearProducts(){ std::vector<SimpleEntry> emptyList; std::swap(emptyList, *productListRef); }
