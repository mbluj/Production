#ifndef HTauhTauhTree_h
#define HTauhTauhTree_h

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

<<<<<<< HEAD
/////////////////////////////////////////////////
  ///Added by AK
  void fillEvent();
  void fillPairs(unsigned int bestPairIndex);
  void fillJets(unsigned int bestPairIndex);
  void fillLeptons();
  void fillGenLeptons();
  bool diMuonVeto();
  bool thirdLeptonVeto(unsigned int signalMuonIndex, unsigned int signalTauIndex, int leptonPdg, double dRmin=-1);
  bool extraMuonVeto(unsigned int signalLeptonIndex, unsigned int signalTauIndex, double dRmin=-1);
  bool extraElectronVeto(unsigned int signalLeptonIndex, unsigned int signalTauIndex, double dRmin=-1);
  bool muonSelection(unsigned int index);
  bool electronSelection(unsigned int index);
  bool pairSelection(unsigned int index);  
  bool jetSelection(unsigned int index, unsigned int bestPairIndex);
  int getMCMatching(unsigned int index);
  bool isGoodToMatch(unsigned int ind);
  TLorentzVector getGenComponentP4(unsigned int index, unsigned int iAbsCharge);
=======
class HTauhTauhTree : public HTauTauTreeBase {
 public :
>>>>>>> 0c40d85... Changes to accomodate split into a base and derived specialized converter classes

  /////////////////////////////////////////////////
  /// TT final state specific
  bool pairSelection(unsigned int index);  
  /////////////////////////////////////////////////
  
  HTauhTauhTree(TTree *tree=0, std::string prefix="WAWTT");
  virtual ~HTauhTauhTree();
  
};

#endif

#ifdef HTauhTauhTree_cxx
HTauhTauhTree::HTauhTauhTree(TTree *tree, std::string prefix) : HTauTauTreeBase(tree, prefix) 
{}

HTauhTauhTree::~HTauhTauhTree()
{}

#endif // #ifdef HTauhTauhTree_cxx
