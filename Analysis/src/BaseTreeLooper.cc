#include "BaseTreeLooper.h"
#include "BaseTreeLooper.hpp"


using namespace std;


BaseTreeLooper::BaseTreeLooper(){}
BaseTreeLooper::BaseTreeLooper(CJLSTTree* inTree){ this->addTree(inTree); }
BaseTreeLooper::BaseTreeLooper(std::vector<CJLSTTree*> const& inTreeList) : treeList(inTreeList){}
BaseTreeLooper::BaseTreeLooper(CJLSTSet const* inTreeSet) : treeList(inTreeSet->getCJLSTTreeList()){}
void BaseTreeLooper::addTree(CJLSTTree* tree){ this->treeList.push_back(tree); }

BaseTreeLooper::~BaseTreeLooper(){}

void BaseTreeLooper::addDiscriminantBuilder(TString KDname, Discriminant* KDbuilder){
  if (KDbuilders.find(KDname)!=KDbuilders.end()) cerr << "BaseTreeLooper::addDiscriminantBuilder: " << KDname << " already exists but will override it regardless." << endl;
  KDbuilders[KDname] = KDbuilder;
}
void BaseTreeLooper::addReweightingBuilder(TString rewgtname, ReweightingBuilder* Rewgtbuilder){
  if (Rewgtbuilders.find(rewgtname)!=Rewgtbuilders.end()) cerr << "BaseTreeLooper::addReweightingBuilder: " << rewgtname << " already exists but will override it regardless." << endl;
  Rewgtbuilders[rewgtname] = Rewgtbuilder;
}
void BaseTreeLooper::addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&)){
  if (externalFunctions.find(fcnname)!=externalFunctions.end()) cerr << "BaseTreeLooper::addExternalFunction: " << fcnname << " already exists but will override it regardless." << endl;
  externalFunctions[fcnname] = fcn;
}

void BaseTreeLooper::addProduct(SimpleEntry& product){ this->productList.push_back(product); }

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
  return process;
}

void BaseTreeLooper::loop(bool loopSelected, bool loopFailed, bool keepProducts){
  // Loop over the trees
  for (CJLSTTree*& tree:treeList){
    // Skip the tree if it cannot be linked
    if (!(this->linkConsumes(tree))) continue;

    // Loop over selected events
    if (loopSelected){
      int ev=0;
      while (tree->getSelectedEvent(ev)){
        SimpleEntry product;
        if (tree->isValidEvent()){
          this->runEvent(tree, product);
          if (keepProducts) this->addProduct(product);
        }
        ev++;
      }
    }
    // Loop over failed events
    if (loopFailed){
      int ev=0;
      while (tree->getFailedEvent(ev)){
        SimpleEntry product;
        if (tree->isValidEvent()){
          this->runEvent(tree, product);
          if (keepProducts) this->addProduct(product);
        }
        ev++;
      }
    }

  } // End loop over the trees
}

const std::vector<SimpleEntry> BaseTreeLooper::getProducts() const{ return productList; }
