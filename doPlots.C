{
 //  /////////////////////////////////////////////
  // here we put the options for ploting
  TStyle *defaultStyle = new TStyle("defaultStyle","Default Style");
  defaultStyle->SetOptStat(0000);
  defaultStyle->SetOptFit(000); 
  defaultStyle->SetPalette(1);
  /////// pad ////////////
  defaultStyle->SetPadBorderMode(1);
  defaultStyle->SetPadBorderSize(1);
  defaultStyle->SetPadColor(0);
  defaultStyle->SetPadTopMargin(0.05);
  defaultStyle->SetPadBottomMargin(0.13);
  defaultStyle->SetPadLeftMargin(0.13);
  defaultStyle->SetPadRightMargin(0.02);
  /////// canvas /////////
  defaultStyle->SetCanvasBorderMode(0);
  defaultStyle->SetCanvasColor(0);
  defaultStyle->SetCanvasDefH(600);
  defaultStyle->SetCanvasDefW(600);
  /////// frame //////////
  defaultStyle->SetFrameBorderMode(0);
  defaultStyle->SetFrameBorderSize(1);
  defaultStyle->SetFrameFillColor(0); 
  defaultStyle->SetFrameLineColor(1);
  /////// label //////////
  defaultStyle->SetLabelOffset(0.005,"XY");
  defaultStyle->SetLabelSize(0.05,"XY");
  defaultStyle->SetLabelFont(42,"XY");
  /////// title //////////
  defaultStyle->SetTitleOffset(1.1,"X");
  defaultStyle->SetTitleSize(0.01,"X");
  defaultStyle->SetTitleOffset(1.25,"Y");
  defaultStyle->SetTitleSize(0.05,"Y");
  defaultStyle->SetTitleFont(42, "XYZ");
  /////// various ////////
  defaultStyle->SetNdivisions(505,"Y");
  defaultStyle->SetLegendBorderSize(0);  // For the axis titles:

    defaultStyle->SetTitleColor(1, "XYZ");
    defaultStyle->SetTitleFont(42, "XYZ");
    defaultStyle->SetTitleSize(0.06, "XYZ");
 
    // defaultStyle->SetTitleYSize(Float_t size = 0.02);
    defaultStyle->SetTitleXOffset(0.9);
    defaultStyle->SetTitleYOffset(1.05);
    // defaultStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    defaultStyle->SetLabelColor(1, "XYZ");
    defaultStyle->SetLabelFont(42, "XYZ");
    defaultStyle->SetLabelOffset(0.007, "XYZ");
    defaultStyle->SetLabelSize(0.04, "XYZ");

    // For the axis:
    defaultStyle->SetAxisColor(1, "XYZ");
    defaultStyle->SetStripDecimals(kTRUE);
    defaultStyle->SetTickLength(0.03, "XYZ");
    defaultStyle->SetNdivisions(510, "XYZ");
    defaultStyle->SetPadTickX(1);   // To get tick marks on the opposite side of the frame
    defaultStyle->SetPadTickY(1);
    defaultStyle->cd();
 /////////////////////////////////////////////////////

 int nmass = 1;
 //int mass[nmass]= {300,500,600,700,800,900,1500,2500,3000};

 const char* channel[nmass]={
 "testout.root"
 }; // close channel

 //const char* lege[nmass]={"300 GeV","500 GeV","600 GeV","700 GeV","800 GeV","900 GeV","1500 GeV","2500 GeV","3000 GeV"};

 const char* lege[nmass]={"250 GeV"};
 //for (int k=0; k<nmass; k++) channel[k] = Form("Control_shower_%d.root", k);

 int nplots =85;
 TLegend *leg = new TLegend(0.7,0.80,0.95,0.90);
   leg->SetTextSize(0.03146853);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
 leg->SetHeader("BKG to 8TeV");

 const char* namplots[nplots] = {
 "hbb_pt.png",
 "hbb_eta.png",
 "hbb_phi.png",
 "hbb_e.png",
 "hbb_mass.png",
 "gen_hbb_pt.png",
 "gen_hbb_phi.png",
 "gen_hbb_eta.png",
 "gen_hbb_e.png",
 "gen_hbb_mass.png",
 "gen_hh_mass.png",
 "bjetpt1.png",
 "bjetpt2.png",
 "bjeteta1.png",
 "bjeteta2.png",
 "bjetphi1.png",
 "bjetphi2.png",
 "bjete1.png",
 "bjete2.png",
 "hww_mt.png",
 "hww_pt.png",
 "hww_phi.png",
 "hww_etap.png",
 "hww_etam.png",
 "gen_hww_mt.png",
 "gen_hww_pt.png",
 "gen_hww_phi.png",
 "gen_hww_eta.png",
 "gen_hww_mass.png",
 "hw1_pt.png",
 "hw2_pt.png",
 "hw1_eta.png",
 "hw2_eta.png",
 "hw1_phi.png",
 "hw2_phi.png",
 "hw1_e.png",
 "hw2_e.png",
 "nlep.png",
 "channel.png",
 "mll.png",
 "ptll.png",
 "pzll.png",
 "dphill.png",
 "pfmet.png",
 "gen_pt1.png",
 "gen_pt2.png",
 "gen_nlep.png",
 "gen_channel.png",
 "gen_mll.png",
 "gen_dphill.png",
 "gen_ptll.png",
 "gen_pzll.png",
 "gen_pfmet.png",
 "gen_pfmez.png",
 "gen_mvv.png",
 "xhh_ww_mt.png",
 "xhh_m_ww_pt.png",
 "xhh_m_ww_eta.png",
 "xhh_m_ww_phi.png",
 "xhh_m_ww_m.png",
 "xhh_p_ww_pt.png",
 "xhh_p_ww_eta.png",
 "xhh_p_ww_phi.png",
 "xhh_p_ww_m.png",
 "Njets.png",
 "Ntags.png",
 "Nbtags.png",
 "NgenVBF.png",
 "gen_vbf_m.png",
 "gen_vbf_m1.png",
 "gen_vbf_m2.png",
 "gen_vbf_pt1.png",
 "gen_vbf_pt2.png",
 "gen_vbf_Deta.png",
 "gen_vbf_DR.png",
 "vbf_m.png",
 "vbf_delta_eta.png",
 "vbf_delta_R.png",
 "vbf_pt1.png",
 "vbf_pt2.png",
 "vbf_m1.png",
 "vbf_m2.png",
 "vbf_genB.png",
 "vbf_btagged.png",
 "vbf_fattagged.png"
    }; //25
////////////////////////////////////////////////////////
 const char* plots[nplots] = {
 "hbb_pt",
 "hbb_eta",
 "hbb_phi",
 "hbb_e",
 "hbb_mass",
 "gen_hbb_pt",
 "gen_hbb_phi",
 "gen_hbb_eta",
 "gen_hbb_e",
 "gen_hbb_mass",
 "gen_hh_mass",
 "bjetpt1",
 "bjetpt2",
 "bjeteta1",
 "bjeteta2",
 "bjetphi1",
 "bjetphi2",
 "bjete1",
 "bjete2",
 "hww_mt",
 "hww_pt",
 "hww_phi",
 "hww_etap",
 "hww_etam",
 "gen_hww_mt",
 "gen_hww_pt",
 "gen_hww_phi",
 "gen_hww_eta",
 "gen_hww_mass",
 "hw1_pt",
 "hw2_pt",
 "hw1_eta",
 "hw2_eta",
 "hw1_phi",
 "hw2_phi",
 "hw1_e",
 "hw2_e",
 "nlep",
 "channel",
 "mll",
 "ptll",
 "pzll",
 "dphill",
 "pfmet",
 "gen_pt1",
 "gen_pt2",
 "gen_nlep",
 "gen_channel",
 "gen_mll",
 "gen_dphill",
 "gen_ptll",
 "gen_pzll",
 "gen_pfmet",
 "gen_pfmez",
 "gen_mvv",
 "xhh_ww_mt",
 "xhh_m_ww_pt",
 "xhh_m_ww_eta",
 "xhh_m_ww_phi",
 "xhh_m_ww_m",
 "xhh_p_ww_pt",
 "xhh_p_ww_eta",
 "xhh_p_ww_phi",
 "xhh_p_ww_m",
 "Njets",
 "Ntags",
 "Nbtags",
 "NgenVBF",
 "gen_vbf_m",
 "gen_vbf_m1",
 "gen_vbf_m2",
 "gen_vbf_pt1",
 "gen_vbf_pt2",
 "gen_vbf_Deta",
 "gen_vbf_DR",
 "vbf_m",
 "vbf_delta_eta",
 "vbf_delta_R",
 "vbf_pt1",
 "vbf_pt2",
 "vbf_m1",
 "vbf_m2",
 "vbf_genB",
 "vbf_btagged",
 "vbf_fattagged"
    }; //25
cout<<channel[0]<<endl;
TCanvas* PT_HAT = new TCanvas();
PT_HAT->cd(); 
TFile *file;
file = TFile::Open(channel[0]);
//TTree *MyTree tes; channel[0]->GetObject("test",MyTree);
TTree* teste = (TTree* ) file->Get("test;1");
cout<<teste->GetEntriesFast()<<endl;
//TBranch *branch = teste.GetBranch("vbf_m");
//TGraph *gr = new TGraph(branch->GetSelectedRows(),branch->GetV2(), branch->GetV1());
// Float_t gen_vbf_DR;

//teste->SetBranchAddress("gen_vbf_DR",&gen_vbf_DR);
//TBranch *brep = teste->GetBranch("gen_vbf_DR");

//Int_t nevent = teste->GetEntries();
//Int_t nselected = 0;
//Float_t repold = -1;
//Float_t *gX  = new Float_t[nevent];
//Float_t *grep = new Float_t[nevent];
cout<<"GetSelectedRows: "<<teste->GetSelectedRows()<<endl;

for(int i=1;i<85;i++) {
	leg->SetHeader("260 GeV graviton - no cuts");
//	for (Int_t j=0; j<nevent; j++) {
//	  brep->GetEvent(j);    // read branch 'rep' only
//	  teste->GetEvent(j);  //read complete accepted event in memory. This
//	  gX[j]  = j;
//	  grep[j] = gen_vbf_DR;	  
//	}
//	TGraph *graph = new TGraph(nevent, gX, grep);
//	graph->Draw("ac");

//	branch->Draw();
	teste->Draw(plots[i]);
        //teste->Draw(plots[i-1], "","same");
//        TGraph *gr = new TGraph(teste->GetSelectedRows(),teste->GetV1(),teste->GetV1());
//        gr->Draw("ap"); //draw graph in current pad
	leg->Draw("same");
//	teste->Draw("gen_vbf_DR,same");
	PT_HAT->SaveAs(namplots[i]);
	PT_HAT->Clear();
	leg->Clear();

} // for nplots

PT_HAT->Close(); 

}
