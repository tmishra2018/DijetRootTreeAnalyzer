#include "tdrstyle.C"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

TCanvas* example_plot( int iPeriod, int iPos );

//###### EDIT THIS PART #######

double xmin = 526;
double xmax = 2231;
double ymin = 0;
double ymax = 1.3;

double xminZoom = 944;
double xmaxZoom = 1530;
double yminZoom = 0.6;
double ymaxZoom = 1.2;

TString myinputFile = "dcap://cmsrm-se01.roma1.infn.it/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/data/Run2015B_plus_Run2015C_goldenJSON_JEC-Summer15_50nsV4_29Aug2015/rootfile_JetHT__Run2015B_plus_Run2015C__MINIAOD_Run2015B_goldenJSON_JEC-Summer15_50nsV4_29Aug2015_reduced_skim.root";

TString mynumerator = "h_mjj_HLTpass_PFHT800";
TString mydenominator = "h_mjj_HLTpass_PFHT475";
TString mytitle = ";Dijet Mass [GeV];Relative Efficiency";
TString mytitlelegend = "PF H_{T} > 800 GeV";

TString myoutputfilename = "triggerEfficiency";
TString myoutputfilename_zoom = "triggerEfficiencyZoom";

//####### NOTE: #######
//Change the style settings under triggerEfficiency()
//#####################


void triggerEfficiency()
{
  //=== General CMS Style ===

  //gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //gROOT->LoadMacro("CMS_lumi.C");
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  writeExtraText = true;       // remove or keep "Preliminary"

  lumi_13TeV  = "2015";  // for trigger
  //lumi_13TeV  = "65 pb^{-1}";  // default is ""
  //lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  //lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)

  //=== Draw ===

  //  example_plot( iPeriod, 0 );   // out of frame (in exceptional cases)
  example_plot( iPeriod, 11 );  // default: left-aligned
  //  example_plot( iPeriod, 22 );  // centered
  //  example_plot( iPeriod, 33 );  // right-aligned  
  // mode generally : 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
}

TCanvas* example_plot( int iPeriod, int iPos )
{ 
  //  if( iPos==0 ) relPosX = 0.12;

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
  float L = 0.12*W_ref;
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

  //----------------------------------------------------------------------------

  TFile *fileInput = TFile::Open(myinputFile);

  fileInput->ls();

  TEfficiency* h_efficiency = 0;

  TH1F *h_denominator = (TH1F*)fileInput->Get(mydenominator);
  TH1F *h_numerator   = (TH1F*)fileInput->Get(mynumerator);

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
      // h_efficiency->GetPaintedGraph()->GetYaxis()->SetTitleOffset(0.9);
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
	  cout << "x= " << x << " , y = " << y << " + " << eyh << " - " << eyl << endl;       
	}

      // draw the legend
      TLegend *legend=new TLegend(0.5,0.5,0.85,0.55);
      //legend->SetTextFont(72);
      //legend->SetTextSize(0.04);
      legend->SetFillStyle(0);
      legend->SetLineColor(0);
      legend->SetShadowColor(0);
      legend->AddEntry(h_efficiency,mytitlelegend,"lpe");
      legend->Draw();
    }

  //-----------------------------------------------------------------------------

  //## Trigger Efficiency plot ##
  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos ); 
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  //canv->Print(canvName+"Effieciency"+".pdf",".pdf");
  canv->Print(myoutputfilename+".pdf",".pdf");
  canv->Print(myoutputfilename+".png",".png");

  //## Trigger Efficiency plot (zoom) ## 
  h_efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(xminZoom,xmaxZoom);
  h_efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(yminZoom,ymaxZoom);

  CMS_lumi( canv, iPeriod, iPos ); 
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  canv->Print(myoutputfilename_zoom+".pdf",".pdf");
  canv->Print(myoutputfilename_zoom+".png",".png");
  
  return canv;
}
