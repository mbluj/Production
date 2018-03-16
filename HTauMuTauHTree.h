#ifndef HTauMuTauHTree_h
#define HTauMuTauHTree_h

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

class HTauMuTauHTree : public HTauTauTreeBase {
 public :

  /////////////////////////////////////////////////
  /// MT final state specific
  bool diMuonVeto();
  bool pairSelection(unsigned int index);
  /////////////////////////////////////////////////
  
  HTauMuTauHTree(TTree *tree=0, bool doSvFit=false, std::string prefix="WAWMT");
  virtual ~HTauMuTauHTree();
  
};

#endif

#ifdef HTauMuTauHTree_cxx
HTauMuTauHTree::HTauMuTauHTree(TTree *tree, bool doSvFit, std::string prefix) : HTauTauTreeBase(tree, doSvFit, prefix) 
{}

HTauMuTauHTree::~HTauMuTauHTree()
{}

#endif // #ifdef HTauMuTauHTree_cxx
