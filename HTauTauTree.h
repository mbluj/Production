#ifndef HTauTauTree_h
#define HTauTauTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

#include "HTTEvent.h"
#include "ScaleFactor.h"
#include <iostream>

// Header of the base class
#include "HTauTauTreeBase.h"

class HTauTauTree : public HTauTauTreeBase {
 public :

  /////////////////////////////////////////////////
  /// MT final state specific
  bool diMuonVeto();
  bool pairSelection(unsigned int index);  
<<<<<<< HEAD
  bool jetSelection(unsigned int index, unsigned int bestPairIndex);
  int getMCMatching(unsigned int index);
  TLorentzVector getGenComponentP4(unsigned int index, unsigned int iAbsCharge);

  template<typename T> T getBranchValue(char *branchAddress, unsigned int index);
  float getProperty(std::string name, unsigned int index);
  std::vector<float> getProperties(const std::vector<std::string> & propertiesList, unsigned int index);
  void writePropertiesHeader(const std::vector<std::string> & propertiesList);
  void writeTriggersHeader(const TH1F*);
=======
  /////////////////////////////////////////////////
>>>>>>> 0c40d85... Changes to accomodate split into a base and derived specialized converter classes
  
  HTauTauTree(TTree *tree=0, std::string prefix="WAWMT");
  virtual ~HTauTauTree();
  
};

#endif

#ifdef HTauTauTree_cxx
HTauTauTree::HTauTauTree(TTree *tree, std::string prefix) : HTauTauTreeBase(tree, prefix) 
{}

HTauTauTree::~HTauTauTree()
{}

#endif // #ifdef HTauTauTree_cxx
