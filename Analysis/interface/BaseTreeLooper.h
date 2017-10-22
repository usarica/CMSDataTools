#ifndef BASETREELOOPER_H
#define BASETREELOOPER_H

#include "CJLSTSet.h"
#include "ReweightingBuilder.h"
#include "Discriminant.h"


class BaseTreeLooper{
protected:
  // List of trees to loop over
  std::vector<CJLSTTree*> treeList;

  // Consumes
  std::unordered_map<TString, short*> valshorts;
  std::unordered_map<TString, unsigned int*> valuints;
  std::unordered_map<TString, int*> valints;
  std::unordered_map<TString, float*> valfloats;
  std::unordered_map<TString, double*> valdoubles;

  std::unordered_map<TString, std::vector<short>*> valVshorts;
  std::unordered_map<TString, std::vector<unsigned int>*> valVuints;
  std::unordered_map<TString, std::vector<int>*> valVints;
  std::unordered_map<TString, std::vector<float>*> valVfloats;
  std::unordered_map<TString, std::vector<double>*> valVdoubles;

  template<typename T> bool linkConsumed(CJLSTTree* tree);
  bool linkConsumes(CJLSTTree* tree);

  // External dependencies
  std::unordered_map<TString, Discriminant*> KDbuilders;
  std::unordered_map<TString, ReweightingBuilder*> Rewgtbuilders;
  std::unordered_map<TString, void(*)(BaseTreeLooper*, SimpleEntry&)> externalFunctions;

  // List of products
  std::vector<SimpleEntry> productList;
  void addProduct(SimpleEntry& product);

  // Abstract function to loop over a single event
  virtual void runEvent(CJLSTTree* tree, SimpleEntry& product)=0;

public:
  // Constructors
  BaseTreeLooper();
  BaseTreeLooper(CJLSTTree* inTree);
  BaseTreeLooper(std::vector<CJLSTTree*> const& inTreeList);
  BaseTreeLooper(CJLSTSet const* inTreeSet);
  void addTree(CJLSTTree* tree);

  // Destructors
  virtual ~BaseTreeLooper();

  // Add the necessary objects
  template<typename T> void addConsumed(TString name);
  void addDiscriminantBuilder(TString KDname, Discriminant* KDbuilder);
  void addReweightingBuilder(TString rewgtname, ReweightingBuilder* Rewgtbuilder);
  void addExternalFunction(TString fcnname, void(*fcn)(BaseTreeLooper*, SimpleEntry&));

  // Function to loop over the tree list
  virtual void loop(bool loopSelected, bool loopFailed, bool keepProducts);

  // Get the products
  const std::vector<SimpleEntry> getProducts() const;

};


#endif
