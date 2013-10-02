/*
 * root -l examples/HWWbb.C\(\"pythia_events.hep.output.root\"\)
 */

//------------------------------------------------------------------------------


#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include "TH1.h"
#include "TLegend.h"
#include "TPaveText.h"

bool isThisJetALepton(TLorentzVector jet, TLorentzVector l1, TLorentzVector l2){
 double DRmax = 0.2;
 bool isLep = false;
 if (jet.DeltaR(l1) < DRmax) isLep = true;
 if (jet.DeltaR(l2) < DRmax) isLep = true;
 return isLep;
}

struct myclassMin {
 bool operator() (std::pair<float, std::pair <int,Int_t> > i, std::pair<float, std::pair <int,Int_t> > j) { 
  return (i.first < j.first);
 }
} myObjMin;


struct myclassMax {
 bool operator() (std::pair<float, std::pair <int,Int_t> > i, std::pair<float, std::pair <int,Int_t> > j) { 
  return (i.first > j.first);
 }
} myObjMax;

struct TestPlots {
 TH1 *hmll;
 TH1 *hdphill;
 
 TH1 *hpt1;
 TH1 *hpt2;
 TH1 *hnleppt10;
 TH1 *hchannel;

 TH1 *hmjj;
 TH1 *hmbb;
 TH1 *hjetpt1;
 TH1 *hjetpt2;
 TH1 *hbjetpt1;
 TH1 *hbjetpt2;
 
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

// # mumu #    channel == 0
// # mue #     channel == 3
// # emu #     channel == 2
// # ee #      channel == 1

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots) {
 TLegend *legend;
 TPaveText *comment;
 
 //---- leptons
 plots->hpt1      = result->AddHist1D( "hpt1",    "p_{T}^{l,max}",   "p_{T}^{l,max}",   "Events / 20 GeV",  50, 0.0, 1000.0);
 plots->hpt2      = result->AddHist1D( "hpt2",    "p_{T}^{l,min}",   "p_{T}^{l,min}",   "Events / 10 GeV",  50, 0.0,  500.0);
 plots->hnleppt10 = result->AddHist1D( "hnleppt10", "num lep pt>10", "num leptons",     "Events"         ,   5, 0.0,    5.0);
 plots->hchannel  = result->AddHist1D( "hchannel", "channel: mm, ee, em, me", "channels",     "Events"   ,   4,-0.5,    3.5);
 plots->hmll      = result->AddHist1D( "hmll",    "m_{ll}",          "m_{ll} [GeV]",    "Events / 10 GeV",  10, 0.0, 100.0);
 plots->hdphill   = result->AddHist1D( "hdphill", "#Delta#phi_{ll}", "#Delta#phi_{ll}", "Events / 0.1 ",    32, 0.0, 3.15);

 //---- jets
 //                                                                                                       bin   min  max  logx logy     
 plots->hmjj      = result->AddHist1D( "hmjj",    "m_{jj}",          "m_{jj} [GeV]",    "Events / 100 GeV", 30, 0.0, 3000.0, 0, 1);
 plots->hmbb      = result->AddHist1D( "hmbb",    "m_{bb}",          "m_{bb} [GeV]",    "Events / 10 GeV",  20, 0.0, 200.0);

 plots->hjetpt1      = result->AddHist1D( "hjetpt1",    "jet-p_{T}^{l,max}",   "jet-p_{T}^{l,max}",   "Events / 20 GeV",  50, 0.0, 1000.0);
 plots->hjetpt2      = result->AddHist1D( "hjetpt2",    "jet-p_{T}^{l,min}",   "jet-p_{T}^{l,min}",   "Events / 10 GeV",  50, 0.0,  500.0);
 
 plots->hbjetpt1      = result->AddHist1D( "hbjetpt1",    "bjet-p_{T}^{l,max}",   "bjet-p_{T}^{l,max}",   "Events / 20 GeV",  50, 0.0, 1000.0);
 plots->hbjetpt2      = result->AddHist1D( "hbjetpt2",    "bjet-p_{T}^{l,min}",   "bjet-p_{T}^{l,min}",   "Events / 10 GeV",  50, 0.0,  500.0);
 
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots) {
 TClonesArray *branchParticle = treeReader->UseBranch("Particle");
 TClonesArray *branchElectron = treeReader->UseBranch("Electron");
 TClonesArray *branchMuon = treeReader->UseBranch("Muon");
 TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
 
 TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
 TClonesArray *branchEFlowTower = treeReader->UseBranch("EFlowTower");
 TClonesArray *branchEFlowMuon = treeReader->UseBranch("EFlowMuon");
 TClonesArray *branchJet = treeReader->UseBranch("Jet");
 
 Long64_t allEntries = treeReader->GetEntries();
 
 cout << "** Chain contains " << allEntries << " events" << endl;
 
 GenParticle *particle;
 Electron *electron;
 Photon *photon;
 Muon *muon;
 
 Track *track;
 Tower *tower;
 
 Jet *jet; // P4 returns a TLorentzVector
 TObject *object;
 
 TLorentzVector momentum;
 
 Float_t Eem, Ehad;
 Bool_t skip;
 
 Long64_t entry;
 
 Int_t i, j, pdgCode;
 
 // Loop over all events
 for(entry = 0; entry < allEntries; entry++) {
  // Load selected branches with data from specified event
  treeReader->ReadEntry(entry);
  
  
  ///---- take the two highest pt leptons in the event (m or e)
  //    maps are ordered in ascending order
  std::map <float, int> m_maxptleptons;
  
  // Loop over all electrons in event
  for(i = 0; i < branchElectron->GetEntriesFast(); i++) {
   electron = (Electron*) branchElectron->At(i);
   double pt = electron->PT;
   m_maxptleptons[-pt] = -(i+1);
  }
  
  // Loop over all muons in event
  for(i = 0; i < branchMuon->GetEntriesFast(); i++) {
   muon = (Muon*) branchMuon->At(i);
   double pt = muon->PT;
   m_maxptleptons[-pt] = (i+1);
  }
  
  //---- at least 2 leptons ----
  if (m_maxptleptons.size() < 2) continue;
  
  // kind = 0/1 if m/e
  
  std::map<float, int>::iterator it_type_m_lepton = m_maxptleptons.begin();
  int flav1 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
  plots->hpt1->Fill(- it_type_m_lepton->first);
  
  it_type_m_lepton++;
  int flav2 = (it_type_m_lepton->second<0);  // m>0, e<0 ---> m=0, e=1
  plots->hpt2->Fill(- it_type_m_lepton->first);
  
  int nlep = 0;
  for(it_type_m_lepton = m_maxptleptons.begin(); it_type_m_lepton != m_maxptleptons.end(); it_type_m_lepton++) {
   if ( -(it_type_m_lepton->first) > 10) nlep++;
  }
  plots->hnleppt10->Fill(nlep);
  
  //                       ee/mm          e   m           m    e
  plots->hchannel->Fill(flav1*flav2+2*(flav1>flav2)+3*(flav1<flav2));
  
  // # mumu #    channel == 0
  // # mue #     channel == 3
  // # emu #     channel == 2
  // # ee #      channel == 1
  
  TLorentzVector l1, l2;
  it_type_m_lepton = m_maxptleptons.begin();
  
  if (it_type_m_lepton->second>0) { l1 = ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->P4();}
  else                            { l1 = ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->P4();}
  
  it_type_m_lepton++;
  if (it_type_m_lepton->second>0) { l2 = ((Muon*)         branchMuon->At(  it_type_m_lepton->second - 1 ))->P4();}
  else                            { l2 = ((Electron*) branchElectron->At(-(it_type_m_lepton->second + 1)))->P4();}
  
  
  plots->hmll->Fill((l1+l2).M());
  plots->hdphill->Fill( l1.DeltaPhi(l2) );
  
  
  
  
  
  // Loop over all jets in event
  //---- two highest mjj are VBF jets
  //---- the other two are "H>bb" jets
  
  //---- at least 4 jets with pt>MINPTJET GeV
  float MINPTJET = 15.;  
  int countJets = 0;
  for(i = 0; i < branchJet->GetEntriesFast(); i++) {
   jet = (Jet*) branchJet->At(i);
   TLorentzVector jetP4 = jet->P4();
   if (jet->PT > MINPTJET && !isThisJetALepton(jetP4, l1, l2)) countJets++;
  }

  if (countJets < 4) {
   continue;
  }
  
  TLorentzVector jet1, jet2, jet3, jet4;
  TLorentzVector Jet1, Jet2, bJet1, bJet2; //---- the VBF jets and the H>bb jets
  int ijet = 0;
  for(i = 0; i < branchJet->GetEntriesFast(); i++) {
   jet = (Jet*) branchJet->At(i);
   TLorentzVector jetP4 = jet->P4();
   
   if (jet->PT > MINPTJET && !isThisJetALepton(jetP4, l1, l2)) {
    if      (ijet == 0) {jet1 = jetP4; ijet++; }
    else if (ijet == 1) {jet2 = jetP4; ijet++; }
    else if (ijet == 2) {jet3 = jetP4; ijet++; }
    else if (ijet == 3) {jet4 = jetP4; ijet++; }
   }
  }
  
  
  if (jet4.Pt() < MINPTJET) {
   std::cout << "We have a problem; countJets = " << countJets << "; ijet = " << ijet << " and jet4.Pt() = " << jet4.Pt() << std::endl; 
  }
  
  float mjj12 = (jet1+jet2).M();
  float mjj13 = (jet1+jet3).M();
  float mjj14 = (jet1+jet4).M();
  float mjj23 = (jet2+jet3).M();
  float mjj24 = (jet2+jet4).M();
  float mjj34 = (jet3+jet4).M();
  
  std::map<float, int> m_maxmjj;
  m_maxmjj[-mjj12] = 1;
  m_maxmjj[-mjj13] = 2;
  m_maxmjj[-mjj14] = 3;
  m_maxmjj[-mjj23] = 4;
  m_maxmjj[-mjj24] = 5;
  m_maxmjj[-mjj34] = 6;
  
  
  std::map<float, int>::iterator it_type_m_maxmjj = m_maxmjj.begin();
  
  plots->hmjj->Fill( - it_type_m_maxmjj->first );
  
  if (it_type_m_maxmjj->second == 1) { Jet1 = jet1; Jet2 = jet2; bJet1 = jet3; bJet2 = jet4; };
  if (it_type_m_maxmjj->second == 2) { Jet1 = jet1; Jet2 = jet3; bJet1 = jet2; bJet2 = jet4; };
  if (it_type_m_maxmjj->second == 3) { Jet1 = jet1; Jet2 = jet4; bJet1 = jet2; bJet2 = jet3; };
  if (it_type_m_maxmjj->second == 4) { Jet1 = jet2; Jet2 = jet3; bJet1 = jet1; bJet2 = jet4; };
  if (it_type_m_maxmjj->second == 5) { Jet1 = jet2; Jet2 = jet4; bJet1 = jet1; bJet2 = jet3; };
  if (it_type_m_maxmjj->second == 6) { Jet1 = jet3; Jet2 = jet4; bJet1 = jet1; bJet2 = jet2; };
  
  
  //---- sub-order in pt
  if (Jet1.Pt() < Jet2.Pt()) {
   TLorentzVector tempjet = Jet1;
   Jet1 = Jet2;
   Jet2 = tempjet;
  }
  
  if (bJet1.Pt() < bJet2.Pt()) {
   TLorentzVector tempjet = bJet1;
   bJet1 = bJet2;
   bJet2 = tempjet;
  }  
  
  plots->hmbb->Fill( (bJet1+bJet2).M() );
  
  plots->hbjetpt1->Fill( bJet1.Pt() );
  plots->hbjetpt2->Fill( bJet2.Pt() );
  
  plots->hjetpt1->Fill( Jet1.Pt() );
  plots->hjetpt2->Fill( Jet2.Pt() );
  
  

  
  
  
  
  //    std::cout << " it_type_m_lepton->first = " << it_type_m_lepton->first << std::endl;
  //    plots->hpt1->Fill(m_maxptleptons.begin()->first);
  
  //    typedef std::map<float, int>::iterator it_type_m_lepton;
  //    for(it_type iterator = m.begin(); iterator != m.end(); iterator++) {
  //       plots->hpt1->Fill(v_maxptleptons.at(0).first);
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  // Loop over all tracks in event
  for(i = 0; i < branchEFlowTrack->GetEntriesFast(); i++) {
   //    track = (Track*) branchEFlowTrack->At(i);
   //    particle = (GenParticle*) track->Particle.GetObject();
   //    
   //    plots->fTrackDeltaPT->Fill((particle->PT - track->PT)/particle->PT);
   //    plots->fTrackDeltaEta->Fill((particle->Eta - track->Eta)/particle->Eta);
  }
  
  // Loop over all towers in event
  for(i = 0; i < branchEFlowTower->GetEntriesFast(); i++) {
   //    tower = (Tower*) branchEFlowTower->At(i);
   //    v_maxptleptons
   //    Eem = 0.0;
   //    Ehad = 0.0;
   //    skip = kFALSE;
   //    for(j = 0; j < tower->Particles.GetEntriesFast(); ++j)
   //    {
   //     particle = (GenParticle*) tower->Particles.At(j);
   //     pdgCode = TMath::Abs(particle->PID);
   //     
   //     // skip muons and neutrinos
   //     if(pdgCode == 12 || pdgCode == 13 || pdgCode == 14 || pdgCode == 16)
   //     {
   //      continue;
   //     }
   //     
   //     // skip K0short and Lambda
   //     if(pdgCode == 310 || pdgCode == 3122)
   //     {
   //      skip = kTRUE;
   //     }
   //     
   //     if(pdgCode == 11 || pdgCode == 22)
   //     {
   // //      Eem += particle->E;
   //     }
   //     else
   //     {
   //      Ehad += particle->E;
   //     }
   //    }
   //    if(skip) continue;
   //    if(Eem > 0.0 && tower->Eem > 0.0) plots->fTowerDeltaEem->Fill((Eem - tower->Eem)/Eem);
   //    if(Ehad > 0.0 && tower->Ehad > 0.0) plots->fTowerDeltaEhad->Fill((Ehad - tower->Ehad)/Ehad);
  }
  
  // Loop over all jets in event
  for(i = 0; i < branchJet->GetEntriesFast(); i++) {
   //    jet = (Jet*) branchJet->At(i);
   //    
   //    momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
   //    
   //    // Loop over all jet's constituents
   //    for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
   //    {
   //     object = jet->Constituents.At(j);
   //     
   //     // Check if the constituent is accessible
   //     if(object == 0) continue;
   //     
   //     if(object->IsA() == GenParticle::Class())
   //     {
   //      momentum += ((GenParticle*) object)->P4();
   //     }
   //     else if(object->IsA() == Track::Class())
   //     {
   //      momentum += ((Track*) object)->P4();
   //     }
   //     else if(object->IsA() == Tower::Class())
   //     {
   //      momentum += ((Tower*) object)->P4();
   //     }
   //     else if(object->IsA() == Muon::Class())
   //     {
   //      momentum += ((Muon*) object)->P4();
   //     }
   //    }
   //    plots->fJetDeltaPT->Fill((jet->PT - momentum.Pt())/jet->PT );
   //   }
   
   
  }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots) {
 result->Print("png");
}

//------------------------------------------------------------------------------

void HHWWbb(const char *inputFile) {
 gSystem->Load("libDelphes");
 
 TChain *chain = new TChain("Delphes");
 chain->Add(inputFile);
 
 ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
 ExRootResult *result = new ExRootResult();
 
 TestPlots *plots = new TestPlots;
 
 BookHistograms(result, plots);
 
 AnalyseEvents(treeReader, plots);
 
 PrintHistograms(result, plots);
 
 result->Write("results_HHWWbb.root");
 
 cout << "** Exiting..." << endl;
 
 delete plots;
 delete result;
 delete treeReader;
 delete chain;
}

//------------------------------------------------------------------------------
