#ifndef HTauhTauhTree_h
#define HTauhTauhTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

#include "HTTEvent.h"
#include <iostream>

// Header of the base class
#include "HTauTauTreeBase.h"


class HTauhTauhTree : public HTauTauTreeBase {
 public :

  /////////////////////////////////////////////////
  /// TT final state specific
  bool pairSelection(unsigned int index);
  unsigned int bestPair(std::vector<unsigned int> &pairIndexes);
  /////////////////////////////////////////////////
  
  HTauhTauhTree(TTree *tree=0, bool doSvFit=false, std::string prefix="WAWTT");
  virtual ~HTauhTauhTree();
  
};

#endif

#ifdef HTauhTauhTree_cxx
HTauhTauhTree::HTauhTauhTree(TTree *tree, bool doSvFit, std::string prefix) : HTauTauTreeBase(tree, doSvFit, prefix) 
{}

HTauhTauhTree::~HTauhTauhTree()
{}

#endif // #ifdef HTauhTauhTree_cxx
