// Dumb selection template
// O. Bondu (June 2013)
// C++ headers
#include <iostream>
#include <string>
// ROOT headers
#include "TROOT.h"
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TClonesArray.h>
// Delphes headers
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "classes/DelphesClasses.h"
// Verbosity
#define DEBUG 0
// namespaces
using namespace std;

int main()
{
	TChain *chain = new TChain("Delphes");
	chain->Add("../delphes_output_test.root");	
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TFile *outfile = new TFile("output.root", "RECREATE");
	TTree *outtree = new TTree("diphoton", "reduced");

	float gg_mass;
	outtree->Branch("gg_mass", &gg_mass, "gg_mass/F");

	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

	for(int ievt=0 ; ievt < treeReader->GetEntries() ; ievt++)
	{
		treeReader->ReadEntry(ievt);
		if(DEBUG) cout << "ievt= " << ievt << endl;
		if(branchPhoton->GetEntries() < 2) continue;
		Photon *photon1 = (Photon*) branchPhoton->At(0);
		Photon *photon2 = (Photon*) branchPhoton->At(1);
		TLorentzVector ph1, ph2;
		ph1.SetPtEtaPhiE(photon1->PT, photon1->Eta, photon1->Phi, photon1->E);
		ph2.SetPtEtaPhiE(photon2->PT, photon2->Eta, photon2->Phi, photon2->E);
		TLorentzVector gg = ph1+ph2;
		
		gg_mass = gg.M();
		outtree->Fill();

	}

	outfile->cd();
	outtree->Write();
	outfile->Close();

	return 0;
}
