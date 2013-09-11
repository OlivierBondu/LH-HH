// Dumb selection template
// O. Bondu (June 2013)
// C++ headers
#include <iostream>
#include <string>
#include <boost/program_options.hpp>
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
namespace po = boost::program_options;

int main(int argc, char* argv[])
{
	// declare arguments
	string inputfile;
	string outputfile;
	string outputtree;
	// print out passed arguments
	copy(argv, argv + argc, ostream_iterator<char*>(cout, " ")); cout << endl;
	// argument parsing
	try
	{
		po::options_description desc("Allowed options");
		desc.add_options()
			("help,h", "produce help message")
			("inputfile,i", po::value<string>(&inputfile)->default_value("../GluGluToHHTo2B2G_M-125_8TeV_madgraph_v2_DEL_v03.root"), "input file")
			("outputfile,o", po::value<string>(&outputfile)->default_value("output.root"), "output file")
			("outputtree,ot", po::value<string>(&outputtree)->default_value("GluGluToHHTo2B2G_8TeV"), "output tree")
		;
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm);
		if (vm.count("help")) {
			cout << desc << "\n";
			return 1;
		}
	} catch(exception& e) {
		cerr << "error: " << e.what() << "\n";
		return 1;
	} catch(...) {
		cerr << "Exception of unknown type!\n";
	}
	// end of argument parsing
  //################################################

	TChain *chain = new TChain("Delphes");
	chain->Add(inputfile.c_str());	
	ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
	TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
	TTree *outtree = new TTree(outputtree.c_str(), "selected events");

// declare output variables
	float pho1_pt, pho1_eta, pho1_phi, pho1_e, pho1_mass;
	float pho2_pt, pho2_eta, pho2_phi, pho2_e, pho2_mass;
	float jet1_pt, jet1_eta, jet1_phi, jet1_e, jet1_mass;
	float jet2_pt, jet2_eta, jet2_phi, jet2_e, jet2_mass;
	float dipho_pt, dipho_eta, dipho_phi, dipho_e, dipho_mass;
	float dijet_pt, dijet_eta, dijet_phi, dijet_e, dijet_mass;
	float tetraphojet_pt, tetraphojet_eta, tetraphojet_phi, tetraphojet_e, tetraphojet_mass;
	outtree->Branch("pho1_pt", pho1_ptpho1_pt, "pho1_pt/F");
	outtree->Branch("pho1_eta", pho1_etapho1_eta, "pho1_eta/F");
	outtree->Branch("pho1_phi", pho1_phipho1_phi, "pho1_phi/F");
	outtree->Branch("pho1_e", pho1_epho1_e, "pho1_e/F");
	outtree->Branch("pho1_mass", pho1_masspho1_mass, "pho1_mass/F");
	outtree->Branch("pho2_pt", pho2_ptpho2_pt, "pho2_pt/F");
	outtree->Branch("pho2_eta", pho2_etapho2_eta, "pho2_eta/F");
	outtree->Branch("pho2_phi", pho2_phipho2_phi, "pho2_phi/F");
	outtree->Branch("pho2_e", pho2_epho2_e, "pho2_e/F");
	outtree->Branch("pho2_mass", pho2_masspho2_mass, "pho2_mass/F");
	outtree->Branch("jet1_pt", jet1_ptjet1_pt, "jet1_pt/F");
	outtree->Branch("jet1_eta", jet1_etajet1_eta, "jet1_eta/F");
	outtree->Branch("jet1_phi", jet1_phijet1_phi, "jet1_phi/F");
	outtree->Branch("jet1_e", jet1_ejet1_e, "jet1_e/F");
	outtree->Branch("jet1_mass", jet1_massjet1_mass, "jet1_mass/F");
	outtree->Branch("jet2_pt", jet2_ptjet2_pt, "jet2_pt/F");
	outtree->Branch("jet2_eta", jet2_etajet2_eta, "jet2_eta/F");
	outtree->Branch("jet2_phi", jet2_phijet2_phi, "jet2_phi/F");
	outtree->Branch("jet2_e", jet2_ejet2_e, "jet2_e/F");
	outtree->Branch("jet2_mass", jet2_massjet2_mass, "jet2_mass/F");
	outtree->Branch("dipho_pt", dipho_ptdipho_pt, "dipho_pt/F");
	outtree->Branch("dipho_eta", dipho_etadipho_eta, "dipho_eta/F");
	outtree->Branch("dipho_phi", dipho_phidipho_phi, "dipho_phi/F");
	outtree->Branch("dipho_e", dipho_edipho_e, "dipho_e/F");
	outtree->Branch("dipho_mass", dipho_massdipho_mass, "dipho_mass/F");
	outtree->Branch("dijet_pt", dijet_ptdijet_pt, "dijet_pt/F");
	outtree->Branch("dijet_eta", dijet_etadijet_eta, "dijet_eta/F");
	outtree->Branch("dijet_phi", dijet_phidijet_phi, "dijet_phi/F");
	outtree->Branch("dijet_e", dijet_edijet_e, "dijet_e/F");
	outtree->Branch("dijet_mass", dijet_massdijet_mass, "dijet_mass/F");
	outtree->Branch("tetraphojet_pt", tetraphojet_pttetraphojet_pt, "tetraphojet_pt/F");
	outtree->Branch("tetraphojet_eta", tetraphojet_etatetraphojet_eta, "tetraphojet_eta/F");
	outtree->Branch("tetraphojet_phi", tetraphojet_phitetraphojet_phi, "tetraphojet_phi/F");
	outtree->Branch("tetraphojet_e", tetraphojet_etetraphojet_e, "tetraphojet_e/F");
	outtree->Branch("tetraphojet_mass", tetraphojet_masstetraphojet_mass, "tetraphojet_mass/F");

	TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

	for(int ievt=0 ; ievt < treeReader->GetEntries() ; ievt++)
	{
		treeReader->ReadEntry(ievt);
		if(DEBUG) cout << "ievt= " << ievt << endl;
    // diphoton selection
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

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
