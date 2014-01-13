{
  /////////////////////////////////////////////
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

 int nmass = 8;
 //int mass[nmass]= {300,500,600,700,800,900,1500,2500,3000};

 const char* channel[nmass]={
 "bulk_graviton/CMS_MG260.root",
 "bulk_graviton/kgraviton_cg0_0137.root",
 "bulk_graviton/kgraviton_cg1.root",
 "bulk_graviton/kgraviton_cg0.root",
 "bulk_graviton/CMS_MG450.root",
 "bulk_graviton/CMS_MG550.root",
 "bulk_graviton/CMS_MG600.root",
 "bulk_graviton/CMS_MG650.root"};/*,
 "bulk_graviton/CMS_MG700.root",
 "bulk_graviton/CMS_MG750.root",
 "bulk_graviton/CMS_MG800.root"
 };*/ // close channel

 //const char* lege[nmass]={"300 GeV","500 GeV","600 GeV","700 GeV","800 GeV","900 GeV","1500 GeV","2500 GeV","3000 GeV"};

 const char* lege[nmass]={"Calc cg =0","Mad cg = 0.0137","Mad cg = 1","Mad cg =0","450 GeV","550 GeV","600 GeV","650 GeV"};/*,
	"700 GeV","750 GeV","800 GeV"};*/
 const int color[nmass]={91,8,99,42,106,223,66,3};//,227,205,32};
 //for (int k=0; k<nmass; k++) channel[k] = Form("Control_shower_%d.root", k);

 int nplots =75;
 TLegend *leg = new TLegend(0.65,0.49,0.99,0.99);
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
 "gen_total_mass.png",
 "hw1_pt.png",
 "hw2_pt.png",
 "hw1_eta.png",
 "hw2_eta.png",
 "hw1_phi.png",
 "hw2_phi.png",
 "hw1_e.png",
 "hw2_e.png",
 "nlep.png",
 "Njets.png",
 "Ntags.png",
 "Nbtags.png",
 "NgenVBF.png",
 "gen_vbf_m.png",
 "gen_vbf_m1.png",
 "gen_vbf_m2.png",
 "gen_vbf_pt1.png",
 "gen_vbf_pt2.png",
 "gen_vbf_eta1.png",
 "gen_vbf_eta2.png",
 "gen_vbf_Deta.png",
 "gen_vbf_DR.png",
 "gen_vbf1H1.png",
 "gen_vbf1H2.png",
 "gen_vbf2H1.png",
 "gen_vbf2H2.png",
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
 "hww_mt",
 "hww_pt",
 "hww_phi",
 "hww_etap",
 "hww_etam",
 "gen_hww_mt",
 "gen_hww_pt",
 "gen_hww_phi",
 "gen_hww_phi",
 "gen_hww_mass",
 "gen_total_mass",
 "hw1_pt",
 "hw2_pt",
 "hw1_eta",
 "hw2_eta",
 "hw1_phi",
 "hw2_phi",
 "hw1_e",
 "hw2_e",
 "nlep",
 "Njets",
 "Ntags",
 "Nbtags",
 "NgenVBF",
 "gen_vbf_m",
 "gen_vbf_m1",
 "gen_vbf_m2",
 "gen_vbf_pt1",
 "gen_vbf_pt2",
 "gen_vbf_eta1",
 "gen_vbf_eta2",
 "gen_vbf_Deta",
 "gen_vbf_DR",
 "gen_vbf1H1",
 "gen_vbf1H2",
 "gen_vbf2H1",
 "gen_vbf2H2",
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
TCanvas* PT_HAT = new TCanvas();
PT_HAT->cd();
for(int i=0;i<22;i++) { // 0-22 /// 22 -46 /// 46 - 57
	leg->SetHeader("M = 260 GeV");
//	for (Int_t j=0; j<nevent; j++) {
//	  brep->GetEvent(j);    // read branch 'rep' only
//	  teste->GetEvent(j);  //read complete accepted event in memory. This
//	  gX[j]  = j;
//	  grep[j] = gen_vbf_DR;	  
//	}
//	TGraph *graph = new TGraph(nevent, gX, grep);
//	graph->Draw("ac");
/*  float y;
  cout<<plots[i]<<endl;
  teste->SetBranchAddress(plots[i] , &y);
  cout<<teste.GetMinimum()<<endl;
  //teste->SetBranchAddress(plots[i-1] , &x);
  float ymin=-4, ymax= 4, ybin = (ymax - ymin)/100;
  TH1F* scan = new TH1F("scan", "", ybin, ymin, ymax);
  int nevent = teste->GetEntries();
  for(int i=0; i<nevent; ++i){
    teste->GetEvent(i);
    //if(scan->GetBinContent(scan->FindBin(y))==0){
    scan->Fill(y);
    //}
    cout<<"GetSelectedRows: "<<y<<" "<<endl;
  }
  scan->Draw();
*/
//	branch->Draw();
   TFile *file[8];
   for(int j=0;j<4;j++) { // files
	file[j] = TFile::Open(channel[j]);
	//TTree *MyTree tes; channel[0]->GetObject("test",MyTree);
	TTree* teste = (TTree* ) file[j]->Get("test;1");	
	teste->SetLineColor(color[j]);
	teste->SetLineWidth(3);
	if(j==0) teste->Draw(plots[i]); else teste->Draw(plots[i], "","same");
        leg->AddEntry(teste,lege[j]);
/*        TGraph *gr = new TGraph(teste->GetSelectedRows(),teste->GetV1(),teste->GetV1());
        gr->Draw("ap"); //draw graph in current pad
*/
//	teste->Draw("gen_vbf_DR,same");
        cout<<teste->GetEntriesFast()<<j<<" "<<i<<" "<<namplots[i]<<" "<<plots[i]<<endl;
  }// for channel
	leg->Draw("same");
	PT_HAT->SetLogy(0);
	PT_HAT->SaveAs(namplots[i]);
	PT_HAT->Clear();
	leg->Clear();
} // for nplots
PT_HAT->Close();

}
