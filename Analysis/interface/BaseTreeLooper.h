#ifndef BASETREELOOPER_H
#define BASETREELOOPER_H

#include "CJLSTSet.h"
#include "ReweightingBuilder.h"
#include "ZXFakeRateHandler.h"
#include "Discriminant.h"
#include "SystematicsHelpers.h"


class BaseTreeLooper{
public:
  enum SampleIdStorageType{
    kNoStorage,
    kStoreByMH,
    kStoreByHashVal
  };

protected:
  SampleIdStorageType sampleIdOpt; // When not kNoStorage, stores a sample identifier and original tree entry index in each output event

  // List of trees to loop over
  std::vector<CJLSTTree*> treeList;

  // Max. events to process
  int maxNEvents;

  // Verbosity flag
  bool verbose;

  // Consumes
  std::unordered_map<TString, short*> valshorts;
  std::unordered_map<TString, unsigned int*> valuints;
  std::unordered_map<TString, int*> valints;
  std::unordered_map<TString, unsigned long*> valulongs;
  std::unordered_map<TString, long*> vallongs;
  std::unordered_map<TString, float*> valfloats;
  std::unordered_map<TString, double*> valdoubles;

  std::unordered_map<TString, std::vector<short>*> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>*> valVuints;
  std::unordered_map<TString, std::vector<int>*> valVints;
  std::unordered_map<TString, std::vector<unsigned long>*> valVulongs;
  std::unordered_map<TString, std::vector<long>*> valVlongs;
  std::unordered_map<TString, std::vector<float>*> valVfloats;
  std::unordered_map<TString, std::vector<double>*> valVdoubles;

  template<typename T> bool linkConsumed(CJLSTTree* tree);
  bool linkConsumes(CJLSTTree* tree);

  // External dependencies
  std::unordered_map<TString, std::pair<Discriminant*, std::vector<TString>>> KDbuilders;
  std::unordered_map<TString, ReweightingBuilder*> Rewgtbuilders;
  std::unordered_map<TString, ZXFakeRateHandler*> ZXFakeRateHandlers;
  std::unordered_map<TString, void(*)(BaseTreeLooper*, SimpleEntry&)> externalFunctions;
  std::unordered_map<TString, SystematicsHelpers::SystematicsClass*> SystVariations;

  // List of products
  std::vector<SimpleEntry> productList;
  std::vector<SimpleEntry>* productListRef;
  BaseTree* productTree;
  void addProduct(SimpleEntry& product, unsigned int* ev_rec=nullptr);

  // Flush product list into tree
  void recordProductsToTree();

  // Abstract function to loop over a single event
  virtual bool runEvent(CJLSTTree* tree, float const& externalWgt, SimpleEntry& product)=0;

public:
  // Constructors
  BaseTreeLooper();
  BaseTreeLooper(CJLSTTree* inTree);
  BaseTreeLooper(std::vector<CJLSTTree*> const& inTreeList);
  BaseTreeLooper(CJLSTSet const* inTreeSet);
  void addTree(CJLSTTree* tree);

  // Destructors
  virtual ~BaseTreeLooper();

  // Set verbosity
  void setVerbosity(bool flag);

  // Add the necessary objects
  template<typename T> void addConsumed(TString name);
  void addDiscriminantBuilder(TString KDname, Discriminant* KDbuilder, std::vector<TString> const& KDvars);
  void addReweightingBuilder(TString rewgtname, ReweightingBuilder* Rewgtbuilder);
  void addZXFakeRateHandler(TString frname, ZXFakeRateHandler* ZXFakeRateHandler);
  void addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&));
  void addSystematic(TString systname, SystematicsHelpers::SystematicsClass* systVar);
  void setExternalProductList(std::vector<SimpleEntry>* extProductListRef=nullptr);
  void setExternalProductTree(BaseTree* extTree=nullptr);

  // Max. events
  void setMaximumEvents(int n);

  // Sample id storage option
  // POWHEG can be stored by mH, but might be better to use hash in others
  void setSampleIdStorageOption(SampleIdStorageType opt);

  // Function to loop over the tree list
  virtual void loop(bool loopSelected, bool loopFailed, bool keepProducts);

  // Get the products
  std::vector<SimpleEntry> const& getProducts() const;
  // Move the products
  void moveProducts(std::vector<SimpleEntry>& targetColl);
  // Clear the products
  void clearProducts();

};


#endif
