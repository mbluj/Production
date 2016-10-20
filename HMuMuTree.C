#define HMuMuTree_cxx
#include "HMuMuTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HMuMuTree::pairSelection(unsigned int iPair){

  ///Baseline+post sync selection as on
  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015#Baseline_mu_tau_h_AN1
  ///Indexes for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
  ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

  if(!mothers_px->size()) return false;

  int pdgIdLeg1 = PDGIdDaughters->at(indexDau1->at(iPair));
  int pdgIdLeg2 = PDGIdDaughters->at(indexDau2->at(iPair));
  if( std::abs(pdgIdLeg1)!=13 || std::abs(pdgIdLeg2)!=13 ) return 0;
  unsigned int indexLeg1 = indexDau1->at(iPair);
  unsigned int indexLeg2 = indexDau2->at(iPair);

  //Sort muons within the pair
  double pt2_1 = (daughters_px->at(indexLeg1)*daughters_px->at(indexLeg1)+
		  daughters_py->at(indexLeg1)*daughters_py->at(indexLeg1));
  double pt2_2 = (daughters_px->at(indexLeg2)*daughters_px->at(indexLeg2)+
		  daughters_py->at(indexLeg2)*daughters_py->at(indexLeg2));
  if(pt2_2>pt2_1){//Muon with higher-Pt first
    unsigned int indexLegTmp = indexLeg1;
    indexLeg1 = indexLeg2;
    indexLeg2 = indexLegTmp;    
  }
  TLorentzVector mu1P4(daughters_px->at(indexLeg1),
		       daughters_py->at(indexLeg1),
		       daughters_pz->at(indexLeg1),
		       daughters_e->at(indexLeg1));
  
  TLorentzVector mu2P4(daughters_px->at(indexLeg2),
		       daughters_py->at(indexLeg2),
		       daughters_pz->at(indexLeg2),
		       daughters_e->at(indexLeg2));
  
  bool muonBaselineSelection1 = (mu1P4.Pt()>23 && std::abs(mu1P4.Eta())<2.4 &&	
				 std::abs(dz->at(indexLeg1))<0.2 &&
				 std::abs(dxy->at(indexLeg1))<0.045 &&
				 ((daughters_muonID->at(indexLeg1) & (1<<6)) == (1<<6)));//Use Short Term Instructions for ICHEP 2016

  bool muonBaselineSelection2 = (mu2P4.Pt()>10 && std::abs(mu2P4.Eta())<2.4 &&	
				 std::abs(dz->at(indexLeg2))<0.2 &&
				 std::abs(dxy->at(indexLeg2))<0.045 &&
				 ((daughters_muonID->at(indexLeg2) & (1<<6)) == (1<<6)));//Use Short Term Instructions for ICHEP 2016
  
  bool baselinePair = mu1P4.DeltaR(mu2P4) > 0.3;							     
  bool postSynchMuon1 = combreliso->at(indexLeg1)<0.15;
  bool loosePostSynchMuon1 = combreliso->at(indexLeg1)<0.3;
  bool postSynchMuon2 = combreliso->at(indexLeg2)<0.15;
  bool loosePostSynchMuon2 = combreliso->at(indexLeg2)<0.3;
  ///
  bool triggerSelection = (triggerbit & 1<<0) == (1<<0);
  
  httEvent->setSelectionBit(SelectionBitsEnum::muonBaselineSelection,muonBaselineSelection1);
  httEvent->setSelectionBit(SelectionBitsEnum::tauBaselineSelection,muonBaselineSelection2);
  httEvent->setSelectionBit(SelectionBitsEnum::baselinePair,baselinePair);
  httEvent->setSelectionBit(SelectionBitsEnum::postSynchMuon,postSynchMuon1);
  httEvent->setSelectionBit(SelectionBitsEnum::postSynchTau,postSynchMuon2);
  httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,thirdLeptonVeto(indexLeg1, indexLeg2, 13));
  httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto,thirdLeptonVeto(indexLeg1, indexLeg2, 11));
  
  return muonBaselineSelection1 && muonBaselineSelection2 && baselinePair
    && ( (loosePostSynchMuon1 && postSynchMuon2) || (loosePostSynchMuon2 && postSynchMuon1) )
    //&& !thirdLeptonVeto(indexMuonLeg, indexTauLeg, 13)
    && true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
