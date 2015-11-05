#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

TCanvas* example_plot( int iPeriod, int iPos );

//###### EDIT THIS PART #######

int histoFromFile = 0; //=1 histograms available in root file, =0 create histo from tree

double xmin = 156;
double xmax = 4010;
double ymin = 0;
double ymax = 1.3;

//double xminZoom = 944;
//double xmaxZoom = 1530;
double xminZoom = 1181;
double xmaxZoom = 2438;
double yminZoom = 0.6;
double ymaxZoom = 1.2;


double threshold = 1213;
//double threshold = 1700;

//TString myinputFile = "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/data/Run2015B_plus_Run2015C_goldenJSON_JEC-Summer15_50nsV4_29Aug2015/rootfile_JetHT__Run2015B_plus_Run2015C__MINIAOD_Run2015B_goldenJSON_JEC-Summer15_50nsV4_29Aug2015_reduced_skim.root";
//TString myinputFile = "/cmshome/santanas/CMS/Releases/CMSSW_7_4_3/src/CMSDIJET/DijetRootTreeAnalyzer/test/santanas/output/output_test_reduced_skim.root";
//--------
//TString myinputFile =  "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134/merged/rootfile_SingleMuon__Run2015B-17Jul2015-v1__MINIAOD_santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134_reduced_skim.root";
//TString myinputFile =  "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134/merged/rootfile_SingleMuon__Run2015B-PromptReco-v1__MINIAOD_santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134_reduced_skim.root";
//TString myinputFile =  "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134/merged/rootfile_SingleMuon__Run2015C-PromptReco-v1__MINIAOD_santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134_reduced_skim.root";
//TString myinputFile =  "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134/merged/rootfile_SingleMuon__Run2015BandC-All-v1__MINIAOD_santanas__SingleMu__65pb-1_50ns_11_09_2015_20150911_162134_reduced_skim.root";
//TString myinputFile =  "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/santanas__SingleMu__15p5pb-1_25ns_22_09_2015_20150922_155905/merged/rootfile_SingleMuon__Run2015C-PromptReco-v1__MINIAOD_santanas__SingleMu__15p5pb-1_25ns_22_09_2015_20150922_155905_reduced_skim.root";
//TString myinputFile =  "/t3/users/santanas/Dijet/reducedRootTrees/rootfile_SingleMuon__Run2015D-PromptReco-v4__MINIAOD_santanas__SingleMu__721pb-1_25ns_28_10_2015_20151029_003601_reduced_skim.root";
TString myinputFile =  "/t3/users/santanas/Dijet/reducedRootTrees/rootfile_SingleMuon__Run2015D-PromptReco-v4__MINIAOD_santanas__SingleMu__728pb-1_25ns_29_10_2015_20151030_125632_reduced_skim.root";

//TString mynumerator = "h_mjj_HLTpass_PFHT800"; // only used if histoFromFile = 1
TString mynumerator = "h_mjj_HLTpass_PFHT800AndMu45Eta2p1"; // only used if histoFromFile = 1
//TString mytitlelegendNum = "PFHT800 AND Mu45Eta2p1";
TString mytitlelegendNum = "(PFHT800 OR PFJET500) AND Mu45Eta2p1";
//TString mytitlelegendNum = "(All JetHT triggers) AND Mu45Eta2p1";

//TString mydenominator = "h_mjj_HLTpass_PFHT475"; // only used if histoFromFile = 1
TString mydenominator = "h_mjj_HLTpass_Mu45Eta2p1"; // only used if histoFromFile = 1
TString mytitlelegendDen = "Mu45Eta2p1";

//TString mytitle = ";Dijet Mass [GeV];Relative Efficiency";
//TString mytitle = ";Dijet Mass [GeV];Trigger Efficiency";
TString mytitle = ";Dijet Mass [GeV];Trigger Efficiency (control region)";
//TString mytitlelegend = "PF H_{T} > 800 GeV";
TString mytitlelegend = "PF H_{T} > 800 GeV OR p_{T}^{PF jet} > 500 GeV ";
//TString mytitlelegend = "All JetHT triggers";
//TString mytitlelegend = "Selection of JetHT triggers";

TString xAxisTitle = "Dijet Mass [GeV]";
TString yAxisTitle = "Number of events";

//TString myoutputfilename = "triggerEfficiency";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015B-17Jul2015";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015B-PromptReco";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015C-PromptReco";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015BAndC-All";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015C-25ns-PromptReco";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015Dv4-25ns-PromptReco_HT_DetaJJLess1p3";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015Dv4-25ns-PromptReco_HTOrJetPt_DetaJJLess1p3";
TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015Dv4-25ns-PromptReco_HTOrJetPt_DetaJJFrom1p3To2p6";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015Dv4-25ns-PromptReco_AllJetHT_DetaJJFrom1p3To2p6";
//TString myoutputfilename = "triggerEfficiency_SingleMu_Run2015Dv4-25ns-PromptReco_JetHTBestGuess_DetaJJFrom1p3To2p6";
//TString myoutputfilename = "pippo";

//####### NOTE: #######
//Change the style settings under triggerEfficiency()
//#####################


void triggerEfficiency()
{
  //=== General CMS Style ===

  gROOT->ForceStyle();

  //gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //gROOT->LoadMacro("CMS_lumi.C");
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  writeExtraText = true;       // remove or keep "Preliminary"

  lumi_13TeV  = "2015";  // for trigger
  //lumi_13TeV  = "65 pb^{-1}, 50ns";  // for trigger
  //lumi_13TeV  = "15.5 pb^{-1}, 25ns";  // for trigger

  //lumi_13TeV  = "65 pb^{-1}";  // default is ""
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 11;     // 0=out , 11=left, 22=center, 33=right


  //====================================================================================
  // Style

  int W = 600;
  int H = 600;

  // 
  // Simple example of macro: plot with CMS name and lumi text
  //  (this script does not pretend to work in all configurations)
  // iPeriod = 1*(0/1 7 TeV) + 2*(0/1 8 TeV)  + 4*(0/1 13 TeV) 
  // For instance: 
  //               iPeriod = 3 means: 7 TeV + 8 TeV
  //               iPeriod = 7 means: 7 TeV + 8 TeV + 13 TeV 
  // Initiated by: Gautier Hamel de Monchenault (Saclay)
  // Updated by:   Dinko Ferencek (Rutgers)
  //
  int H_ref = 600; 
  int W_ref = 600; 

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.15*W_ref;
  float R = 0.04*W_ref;

  TString canvName = "trigger";
  canvName += W;
  canvName += "-";
  canvName += H;
  canvName += "_";  
  canvName += iPeriod;
  if( writeExtraText ) canvName += "-prelim";
  if( iPos%10==0 ) canvName += "-out";
  else if( iPos%10==1 ) canvName += "-left";
  else if( iPos%10==2 )  canvName += "-center";
  else if( iPos%10==3 )  canvName += "-right";

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);

  canv->SetGridx(true);
  canv->SetGridy(true);

  //====================================================================================
  // Efficiency

  TFile *fileInput = TFile::Open(myinputFile);

  fileInput->ls();

  TEfficiency* h_efficiency = 0;

  TH1F *h_denominator;
  TH1F *h_numerator;
  if(histoFromFile==1)
    {
      //== taking histo from file       
      h_denominator = (TH1F*)fileInput->Get(mydenominator);
      h_numerator   = (TH1F*)fileInput->Get(mynumerator);
    }
  else
    {
      //== creating histo from tree   
      TTree *thistree = (TTree*)fileInput->Get("rootTupleTree/tree");
      thistree->Print();
      TH1F *h_denominator_tmp = (TH1F*)fileInput->Get(mydenominator);
      h_denominator = (TH1F*)h_denominator_tmp->Clone();
      h_numerator = (TH1F*)h_denominator_tmp->Clone();
      h_denominator->Reset();
      h_numerator->Reset();
      h_denominator->SetName("h_denominator");
      h_numerator->SetName("h_numerator");
      //fill histograms
      //--
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_PFHT475==1"); //signal region
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1"); //signal region
      thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)>1.3 && fabs(deltaETAjj)<2.6 && passHLT_Mu45==1"); //control region
      //--
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_PFHT475==1 && passHLT_PFHT800==1");
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1 && passHLT_PFHT800==1");
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1 && (passHLT_PFHT800==1 || passHLT_PFJET500==1)");
      thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)>1.3 && fabs(deltaETAjj)<2.6 && passHLT_Mu45==1 && (passHLT_PFHT800==1 || passHLT_PFJET500==1)");
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)>1.3 && fabs(deltaETAjj)<2.6 && passHLT_Mu45==1 && (passHLT_PFHT800==1 || passHLT_PFJET500==1 || passHLT_PFHT650MJJ950==1 || passHLT_PFHT650MJJ900==1 || passHLT_AK8DiPFJet280200TrimMass30Btag==1 || passHLT_AK8PFHT600TriMass50Btag==1 || passHLT_AK8PFHT700TriMass50==1 || passHLT_AK8PFJet360TrimMass50==1 || passHLT_CaloJet500NoJetID==1 || passHLT_DiPFJetAve300HFJEC==1 || passHLT_DiPFJetAve500==1 || passHLT_PFHT400SixJet30Btag==1 || passHLT_PFHT450SixJet40Btag==1 || passHLT_PFHT750FourJetPt50==1 || passHLT_QuadPFJetVBF==1 || passHLT_PFHT650==1 || passHLT_PFHT475==1 || passHLT_PFHT200==1 || passHLT_PFJET450==1)");
      //-- option placeholder 
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1","",10000);
      //--  
    }

  //==========================

  if(TEfficiency::CheckConsistency(*h_numerator,*h_denominator))
    {
      h_efficiency = new TEfficiency(*h_numerator,*h_denominator);    
      //stat option, see https://root.cern.ch/root/html/TEfficiency.html#TEfficiency:SetStatisticOption
      h_efficiency->SetStatisticOption(TEfficiency::kFWilson);  
      //h_efficiency->SetStatisticOption(TEfficiency::kFCP); //default  
      h_efficiency->SetTitle(mytitle);
      h_efficiency->Draw();
      gPad->Update();
      h_efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(xmin,xmax);
      h_efficiency->GetPaintedGraph()->GetXaxis()->SetNdivisions(505);
      h_efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(ymin,ymax);
      //h_efficiency->GetPaintedGraph()->GetYaxis()->SetTitleOffset(0.9);
      // h_efficiency->GetPaintedGraph()->GetYaxis()->SetLabelSize(0.04);

      for (int bin=0;bin<h_efficiency->GetPaintedGraph()->GetN();bin++)
	{
	  double x=-1; 
	  double y=-1;
	  double eyh=-1;
	  double eyl=-1;

	  h_efficiency->GetPaintedGraph()->GetPoint(bin,x,y);
	  eyh = h_efficiency->GetPaintedGraph()->GetErrorYhigh(bin);
	  eyl = h_efficiency->GetPaintedGraph()->GetErrorYlow(bin);
	  cout << "bin = " << bin << ": x= " << x << " , y = " << y << " + " << eyh << " - " << eyl << endl;       
	}

      // draw the legend
      TLegend *legend=new TLegend(0.35,0.22,0.89,0.32);
      //legend->SetTextFont(72);
      //legend->SetTextSize(0.04);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetShadowColor(0);
      legend->AddEntry(h_efficiency,mytitlelegend,"lpe");
      legend->Draw();
    }

  //====================================================================================
  //Draw

  //## Trigger Efficiency plot ##
  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos ); 
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  canv->Print(myoutputfilename+".pdf",".pdf");
  canv->Print(myoutputfilename+".png",".png");
  canv->Print(myoutputfilename+".root",".root");
  
  //## Trigger Efficiency plot (zoom) ## 
  h_efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(xminZoom,xmaxZoom);
  h_efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(yminZoom,ymaxZoom);

  CMS_lumi( canv, iPeriod, iPos ); 
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  canv->Print(myoutputfilename+"_zoom.pdf",".pdf");
  canv->Print(myoutputfilename+"_zoom.png",".png");
  
  //Integral above threshold
  int totalNev =  h_denominator->Integral(h_denominator->FindFixBin(threshold),h_denominator->GetNbinsX());
  int passedNev =  h_numerator->Integral(h_numerator->FindFixBin(threshold),h_numerator->GetNbinsX());
  //int totalNev =  h_denominator->Integral(h_denominator->FindFixBin(threshold),h_denominator->FindFixBin(threshold));
  //int passedNev =  h_numerator->Integral(h_numerator->FindFixBin(threshold),h_numerator->FindFixBin(threshold));
  float effIntegrated = float(passedNev)/float(totalNev);
  cout << "totalNev = " << totalNev <<  " , passedNev=" << passedNev << " , efficiency=" << effIntegrated << endl;       
  TEfficiency* pEff = 0;
  float effIntegrated_errUp = pEff->Wilson(totalNev,passedNev,0.683,true) - effIntegrated;
  float effIntegrated_errDown = pEff->Wilson(totalNev,passedNev,0.683,false) - effIntegrated;
  cout << "efficiency integrated above threshold of "<< threshold <<" = " << effIntegrated << " + " << effIntegrated_errUp << " - " << effIntegrated_errDown << endl;


  //## Mjj Spectra ##
  canv->SetGridx(false);
  canv->SetGridy(false);
  canv->SetLogy(true);
  h_denominator->UseCurrentStyle();  
  h_denominator->SetLineColor(2);
  h_numerator->SetLineColor(1);
  h_denominator->Draw();
  h_numerator->Draw("same");
  h_denominator->GetXaxis()->SetRangeUser(xmin,xmax);
  h_denominator->GetXaxis()->SetTitle(xAxisTitle);
  h_denominator->GetYaxis()->SetTitle(yAxisTitle);
  h_denominator->GetYaxis()->SetTitleOffset(1.3);

  CMS_lumi( canv, iPeriod, iPos ); 
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  gPad->SetTickx(1);
  gPad->SetTicky(1);

  // draw the legend
  TLegend *legend1=new TLegend(0.4,0.65,0.91,0.83);
  //legend->SetTextFont(72);
  //legend->SetTextSize(0.06);
  legend1->SetFillStyle(0);
  legend1->SetLineColor(0);
  legend1->SetShadowColor(0);
  legend1->AddEntry(h_denominator,mytitlelegendDen,"l");
  legend1->AddEntry(h_numerator,mytitlelegendNum,"l");
  legend1->Draw();

  canv->Print(myoutputfilename+"_histo.pdf",".pdf");
  canv->Print(myoutputfilename+"_histo.png",".png");
  



  //-----------------------------------------------------------------------------

}
