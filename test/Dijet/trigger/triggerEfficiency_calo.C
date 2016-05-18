#include "test/santanas/DijetScouting/style/tdrstyle.C"
#include "test/santanas/DijetScouting/style/CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"

TCanvas* example_plot( int iPeriod, int iPos );

//###### EDIT THIS PART #######

int histoFromFile = 0; //=1 histograms available in root file, =0 create histo from tree

//high-mass
//double xmin = 156;
//double xmax = 3704;
//scouting HT450
double xmin = 156;
double xmax = 2037;
//scouting L1HTT
//double xmin = 23;
//double xmax = 526;
//
double ymin = 0;
double ymax = 1.3;

//high-mass
//double xminZoom = 944;
//double xmaxZoom = 1530;
//scouting HT450
double xminZoom = 156;
double xmaxZoom = 490;
//double xmaxZoom = 1455;
//scouting L1HTT
//double xminZoom = 88;
//double xmaxZoom = 325;
//
double yminZoom = 0.1;
double ymaxZoom = 1.2;

//high-mass
//double threshold = 1213;
//scouting
double threshold = 119;

//
int doFit = 1;
double xminFit = 250;
//double xmaxFit = 2037;
double xmaxFit = 890;
//double xmaxFit = 1455;


//high-mass
// TString myinputFile =  "rootfile_CaloScoutingCommissioning_JEC_CaloL1L2L3_PFL2L3Residual_20160505_183413_59_reduced_skim.root";
TString myinputFile =  "eos/cms/store/group/phys_exotica/dijet/Dijet13TeV/apresyan/reduced_skims/JEC_CaloL1L2L3_PFL2L3Residual_20160505_201702/CaloScoutingCommissioning_JEC_CaloL1L2L3_PFL2L3Residual_reduced_skim.root";
TString mybaselinehisto = "h_mjj_NoTrigger_1GeVbin"; // needed to define the x-axis range 

TString mynumerator = "h_mjj_HLTpass_CaloScoutingHT250"; // only used if histoFromFile = 1
TString mytitlelegendNum = "Calo HT250";//scouting HT450
TString mydenominator = "h_mjj_HLTpass_CaloJet40_CaloScouting_PFScouting"; // only used if histoFromFile = 1
TString mytitlelegendDen = "Calo Jet40"; //scouting HT450

TString mytitle = ";Dijet Mass [GeV];Trigger efficiency";
//TString mytitle = ";Dijet Mass [GeV];Trigger efficiency (control region)";
//TString mytitlelegend = "PF H_{T} > 800 GeV";
//TString mytitlelegend = "PF H_{T} > 800 GeV OR p_{T}^{PF jet} > 500 GeV ";
//TString mytitlelegend = "PF H_{T} > 0.8 TeV OR p_{T}^{PF jet} > 0.5 TeV ";//paper
TString mytitlelegend = "H_{T} > 250 GeV";//scouting HT450
//TString mytitlelegend = "L1HTT";//scouting L1HTT
//TString mytitlelegend = "All JetHT triggers";
//TString mytitlelegend = "Selection of JetHT triggers";

TString xAxisTitle = "Dijet Mass [GeV]";
TString yAxisTitle = "Number of events";

TString myoutputfilename = "triggerEfficiency";//scouting HT450


//####### NOTE: #######
//Change the style settings under triggerEfficiency()
//#####################


void triggerEfficiency_calo()
{
  //=== General CMS Style ===

  gROOT->ForceStyle();

  //gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();

  //gROOT->LoadMacro("CMS_lumi.C");
  //extraText  = "Preliminary";  // default extra text is "Preliminary"
  extraText  = "Supplementary";  // default extra text is "Preliminary"
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
      //thistree->Print();
      //TH1F *h_denominator_tmp = (TH1F*)fileInput->Get(mybaselinehisto);
      TH1F *h_denominator_tmp = new TH1F("h_denominator_tmp","h_denominator_tmp",10000,0,10000) ;
      h_denominator = (TH1F*)h_denominator_tmp->Clone();
      h_numerator = (TH1F*)h_denominator_tmp->Clone();
      h_denominator->Reset();
      h_numerator->Reset();
      h_denominator->SetName("h_denominator");
      h_numerator->SetName("h_numerator");
      //fill histograms
      //--
      cout << "filling denominator" << endl;
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_PFHT475==1"); //signal region
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1"); //signal region //paper      
      thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && PassJSON==1 && (passHLT_L1HTT_CaloScouting_PFScouting)"); //signal region //scouting HT450
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_ZeroBias==1"); //signal region //scouting L1HTT
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)>1.3 && fabs(deltaETAjj)<2.6 && passHLT_Mu45==1"); //control region
      cout << "filled denominator" << endl;
      //--
      cout << "filling numerator" << endl;
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_PFHT475==1 && passHLT_PFHT800==1");
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1 && passHLT_PFHT800==1");
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1 && (passHLT_PFHT800==1 || passHLT_PFJET500==1)");//paper
      thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && PassJSON==1 && passHLT_CaloScoutingHT250 && (passHLT_L1HTT_CaloScouting_PFScouting) ");//scouting HT450
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)<1.3 && passHLT_ZeroBias==1 && passHLT_L1HTT==1");//scouting L1HTT
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)>1.3 && fabs(deltaETAjj)<2.6 && passHLT_Mu45==1 && (passHLT_PFHT800==1 || passHLT_PFJET500==1)");
      //thistree->Draw("mjj >> h_numerator","fabs(deltaETAjj)>1.3 && fabs(deltaETAjj)<2.6 && passHLT_Mu45==1 && (passHLT_PFHT800==1 || passHLT_PFJET500==1 || passHLT_PFHT650MJJ950==1 || passHLT_PFHT650MJJ900==1 || passHLT_AK8DiPFJet280200TrimMass30Btag==1 || passHLT_AK8PFHT600TriMass50Btag==1 || passHLT_AK8PFHT700TriMass50==1 || passHLT_AK8PFJet360TrimMass50==1 || passHLT_CaloJet500NoJetID==1 || passHLT_DiPFJetAve300HFJEC==1 || passHLT_DiPFJetAve500==1 || passHLT_PFHT400SixJet30Btag==1 || passHLT_PFHT450SixJet40Btag==1 || passHLT_PFHT750FourJetPt50==1 || passHLT_QuadPFJetVBF==1 || passHLT_PFHT650==1 || passHLT_PFHT475==1 || passHLT_PFHT200==1 || passHLT_PFJET450==1)");
      cout << "filled numerator" << endl;
      //-- option placeholder 
      //thistree->Draw("mjj >> h_denominator","fabs(deltaETAjj)<1.3 && passHLT_Mu45==1","",10000);
      //--  
      
    }

  //==========================

  if(!TEfficiency::CheckConsistency(*h_numerator,*h_denominator))
    {
      cout << "Numerator and denominator are not consistent! Exit" << endl;
      gApplication->Terminate(0);
    }

  //  if(TEfficiency::CheckConsistency(*h_numerator,*h_denominator))
  //    {
  h_efficiency = new TEfficiency(*h_numerator,*h_denominator);    
  h_efficiency->SetName("efficiency"); 
  //stat option, see https://root.cern.ch/root/html/TEfficiency.html#TEfficiency:SetStatisticOption
  h_efficiency->SetStatisticOption(TEfficiency::kFWilson);  
  //h_efficiency->SetStatisticOption(TEfficiency::kFCP); //default  
  h_efficiency->SetTitle(mytitle);
  
  // //fit efficiency
  TF1* f1 = new TF1("f1","([0]/2)* ( 1 + TMath::Erf((x-[1])/[2]))",xminFit,xmaxFit);      
  //TF1* f1 = new TF1("f1","(1/2)* ( 1 + TMath::Erf((x-[0])/[1]))",xminFit,xmaxFit);      
  TGraphAsymmErrors* graph_efficiency = h_efficiency->CreateGraph();
  TFitResultPtr fitResult;
  f1->SetLineWidth(2);
  if(doFit==1)
    {
      f1->SetParameters(1,300,95);
      //f1->SetParLimits(0,0.99,1);//unconstrained
      //f1->SetParLimits(0,0.999999,1);//constrained
      f1->FixParameter(0,1);//fixed
      f1->SetParLimits(1,200,600);
      f1->SetParLimits(2,30,120);
      //fitResult = h_efficiency->Fit(f1,"S V R I");
      fitResult = graph_efficiency->Fit(f1,"VLRS");            
    }

  //fitResult->Print("V");


  //int numberOfParameters = f1->GetNumberFreeParameters();  
  int numberOfParameters = fitResult->NFreeParameters();  

  //Correlation matrix
  TMatrixDSym corrMatrix = fitResult->GetCorrelationMatrix();
  // cout << corrMatrix[1][1] << endl;
  // cout << corrMatrix[1][2] << endl;
  // cout << corrMatrix[2][1] << endl;
  // cout << corrMatrix[2][2] << endl;
  TMatrixDSym covMatrix = fitResult->GetCovarianceMatrix();
  //cout << covMatrix[1][1] << endl;
  //cout << covMatrix[1][2] << endl;
  //cout << covMatrix[2][1] << endl;
  //cout << covMatrix[2][2] << endl;

  graph_efficiency->Draw("AP");
  // gStyle->SetOptFit(0);
  gPad->Update();
  graph_efficiency->GetXaxis()->SetRangeUser(xmin,xmax);
  graph_efficiency->GetXaxis()->SetNdivisions(505);
  graph_efficiency->GetYaxis()->SetRangeUser(ymin,ymax);

  /*
  h_efficiency->Draw();
  gPad->Update();
  h_efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(xmin,xmax);
  h_efficiency->GetPaintedGraph()->GetXaxis()->SetNdivisions(505);
  h_efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(ymin,ymax);
  */
  
  double chi2 = 0;
  int npoints = 0;
  int ndf = 0;
  double chi2Norm = 0;
  TGraph *g_residual = new TGraph(0);
  g_residual->SetName("");
  TGraphErrors *g_errorBand = new TGraphErrors(0);
  g_errorBand->SetName("errorBand_graph");
  TGraphErrors *g_fitValue = new TGraphErrors(0);
  g_fitValue->SetName("fitValue_graph");

  //  for (int bin=0;bin<h_efficiency->GetPaintedGraph()->GetN();bin++)
  for (int bin=0;bin<graph_efficiency->GetN();bin++)
    {
      double x=-1; 
      double exh=-1;
      double exl=-1;
      double y=-1;
      double eyh=-1;
      double eyl=-1;
      double fitValue = -1;
      double residual = -1;
      double fitError = -1;
      double fitErrorRel = -1;
      
      /*
      h_efficiency->GetPaintedGraph()->GetPoint(bin,x,y);
      eyh = h_efficiency->GetPaintedGraph()->GetErrorYhigh(bin);
      eyl = h_efficiency->GetPaintedGraph()->GetErrorYlow(bin);
      exh = h_efficiency->GetPaintedGraph()->GetErrorXhigh(bin);
      exl = h_efficiency->GetPaintedGraph()->GetErrorXlow(bin);      
      //cout << exl << " , " << exh << endl;
      */

      graph_efficiency->GetPoint(bin,x,y);
      eyh = graph_efficiency->GetErrorYhigh(bin);
      eyl = graph_efficiency->GetErrorYlow(bin);
      exh = graph_efficiency->GetErrorXhigh(bin);
      exl = graph_efficiency->GetErrorXlow(bin);      
      //cout << exl << " , " << exh << endl;

      g_residual->SetPoint(bin,x,0);
      g_residual->RemovePoint(bin);

      g_errorBand->SetPoint(bin,x,0);
      g_errorBand->SetPointError(bin,0,0);
      //g_errorBand->RemovePoint(bin);

      g_fitValue->SetPoint(bin,x,0);
      g_fitValue->SetPointError(bin,0,0);

      if(doFit==1)
	{
	  fitValue = (f1->Integral(x-exl , x+exl))/(exl+exl);

	  if(y<fitValue)
	    residual = (y - fitValue) / eyh;
	  else
	    residual = (y - fitValue) / eyl;

	  double dfdP1   = f1->GradientPar(1,&x);
	  double dfdP2   = f1->GradientPar(2,&x);
	  double VarP1   = covMatrix[1][1] ; 
	  double VarP2   = covMatrix[2][2] ;
	  double CovP1P2 = covMatrix[1][2];

	  //fitError = 0.05;
	  fitError = sqrt(dfdP1*dfdP1*VarP1 + dfdP2*dfdP2*VarP2 + 2*dfdP1*dfdP2*CovP1P2);
	  fitErrorRel = sqrt(dfdP1*dfdP1*VarP1 + dfdP2*dfdP2*VarP2 + 2*dfdP1*dfdP2*CovP1P2)/fitValue;
	}
      
      //cout << "bin = " << bin << ": x= " << x << " , y = " << y << " + " << eyh << " - " << eyl << endl;   	  
      if(x>xminFit && x<xmaxFit && doFit==1)
	{
	  cout << "bin = " << bin << ": x bin = (" << x-exl << "-" << x+exh 
	       << ") , y = " << y << " + " << eyh << " - " << eyl 
	       << " , fit = " << fitValue 
	       << " , fit error = " << fitError 
	       << " , (y-fit) / err = " << residual 
	       << endl ;  
	  chi2 += pow(residual,2); 
	  npoints++;

	  g_residual->SetPoint(bin,x,residual);
	  //g_residual->SetPointError(bin,,0);

	  //cout << x << " - " << fitValue << " +/- " << fitError << endl;

	  /*
	  g_errorBand->SetPoint(bin,x,fitValue);
	  g_errorBand->SetPointError(bin,0,fitError);
	  */

	  g_errorBand->SetPoint(bin,x,0);
	  //g_errorBand->SetPointError(bin,0,fitErrorRel);
	  g_errorBand->SetPointError(bin,0,fitError);

	  g_fitValue->SetPoint(bin,x,fitValue);
	  g_fitValue->SetPointError(bin,0,fitError);
	}
      else
	{
	  /*
	  cout << "bin = " << bin << ": x bin = (" << x-exl << "-" << x+exh 
	       << ") , y = " << y << " + " << eyh << " - " << eyl 
	       << endl ;       	      
	  */
	}
    }
  
  if(doFit==1)
    {
      ndf = npoints - numberOfParameters;
      chi2Norm=chi2/ndf;
      cout << "=== chi2 = " << chi2 << endl;
      cout << "=== npoints = " << npoints << endl;
      cout << "=== numberOfParameters = " << numberOfParameters << endl;
      cout << "=== normalized chi2 = " << chi2Norm << endl;
    }
  
  // draw the legend
  TLegend *legend=new TLegend(0.35,0.22,0.95,0.32);
  //legend->SetTextFont(72);
  //legend->SetTextSize(0.04);
  legend->SetFillStyle(0);
  legend->SetLineColor(0);
  legend->SetShadowColor(0);
  //  legend->AddEntry(h_efficiency,mytitlelegend,"lpe");
  legend->AddEntry(graph_efficiency,mytitlelegend,"lpe");

  legend->Draw();
  //    }
  
  //====================================================================================
  //Draw
  
  //## Trigger Efficiency plot ##
  // writing the lumi information and the CMS "logo"
  CMS_lumi( canv, iPeriod, iPos ); 
  TText *pave = new TText(2346.12,1.19462,"arxiv:1512.01224");
  pave->SetTextColor(4);
  pave->SetTextSize(0.04);
  pave->SetTextFont(42);
  pave->Draw();
  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  canv->Print(myoutputfilename+".eps",".eps");
  canv->Print(myoutputfilename+".pdf",".pdf");
  canv->Print(myoutputfilename+".png",".png");
  canv->Print(myoutputfilename+".root",".root");

  //## Trigger Efficiency plot (zoom) ## 
  /*
  h_efficiency->GetPaintedGraph()->GetXaxis()->SetRangeUser(xminZoom,xmaxZoom);
  h_efficiency->GetPaintedGraph()->GetYaxis()->SetRangeUser(yminZoom,ymaxZoom);
  */
  graph_efficiency->GetXaxis()->SetRangeUser(xminZoom,xmaxZoom);
  graph_efficiency->GetYaxis()->SetRangeUser(yminZoom,ymaxZoom);


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
  

  //## Fit result ##
  TCanvas* canv2 = new TCanvas("canv2","canv2",50,50,W,H);

  // //pad 1
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0.15); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->SetGridy();         // Vertical grid
  pad1->Draw();
  pad1->cd(); 

  //TGraphAsymmErrors *g_efficiency = (TGraphAsymmErrors*) h_efficiency->GetPaintedGraph();
  graph_efficiency->Draw("ap");
  f1->Draw("same");
  gPad->Update();

  TPaveStats *st = (TPaveStats*)graph_efficiency->FindObject("stats");
  st->Delete();
  TPaveText *fit_stat = new TPaveText(0.52,0.22,0.91,0.58,"NDC");    
  char chi2text[100]; 
  sprintf (chi2text, "chi2/ndf = %.1f/%d = %.1f", chi2, ndf, chi2Norm);
  fit_stat->AddText(chi2text);
  char par0text[100]; 
  sprintf (par0text, "p0 = %.6f #pm %.6f", f1->GetParameter(0), f1->GetParError(0));
  fit_stat->AddText(par0text);
  char par1text[100]; 
  sprintf (par1text, "p1 = %.1f #pm %.1f", f1->GetParameter(1), f1->GetParError(1));
  fit_stat->AddText(par1text);
  char par2text[100]; 
  sprintf (par2text, "p2 = %.1f #pm %.1f", f1->GetParameter(2), f1->GetParError(2));
  fit_stat->AddText(par2text);

  fit_stat->SetFillColor(0);
  fit_stat->SetLineColor(1);
  fit_stat->SetFillStyle(1001);
  fit_stat->SetBorderSize(1);
  fit_stat->SetTextFont(42);
  fit_stat->SetTextSize(0.040);
  fit_stat->SetTextAlign(12); 
  fit_stat->SetTextColor(1); 
  fit_stat->Draw("SAME");

  // //pad 2
  canv2->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0.); // Upper and lower plot are joined
  pad2->SetBottomMargin(0.15); // Upper and lower plot are joined
  pad2->SetGridx();         // Vertical grid
  pad2->SetGridy();         // Vertical grid
  pad2->Draw();
  pad2->cd(); 
  g_residual->SetFillColor(2);
  g_residual->SetLineColor(0);
  g_residual->GetXaxis()->SetRangeUser(pad1->GetUxmin(),pad1->GetUxmax()); 
  g_residual->GetXaxis()->SetNdivisions(505);
  g_residual->GetYaxis()->SetRangeUser(-3.5,3.5); 
  g_residual->GetYaxis()->SetLabelSize(0.1); 
  g_residual->GetYaxis()->SetTitleSize(0.1); 
  g_residual->GetYaxis()->SetTitleOffset(0.5); 
  g_residual->GetYaxis()->SetTitle("#frac{(Data-Fit)}{#sigma_{stat}}"); 
  g_residual->Draw("AP");

  canv2->Update();
  canv2->Modified();
  gPad->Update();
  canv2->Print(myoutputfilename+"_fit.root",".root");
  canv2->Print(myoutputfilename+"_fit.png",".png");

  //## Store TEfficiency ##
  TFile outputFile(myoutputfilename+"_output.root","recreate");
  h_efficiency->Write();
  g_fitValue->Write();
  g_errorBand->Write();

  //## Fit error ##
  TCanvas* canv3 = new TCanvas("canv3","canv3",50,50,W,H);
  canv3->SetGridx();
  canv3->SetGridy();
  g_errorBand->SetLineColor(4);
  g_errorBand->SetFillColor(4);
  g_errorBand->SetFillStyle(3001);
  g_errorBand->Draw("aE3");
  g_errorBand->GetXaxis()->SetRangeUser(xminFit,xmaxFit);  
  g_errorBand->GetXaxis()->SetTitle("Dijet mass [GeV]");  
  g_errorBand->GetYaxis()->SetTitle("Tot. fit uncertainty");  
  g_errorBand->GetYaxis()->SetTitleOffset(1.3);  
  TF1* ftot = new TF1("ftot","(1/2)* ( 1 + TMath::Erf((x-495.7)/95.0))",xminFit,xmaxFit);      
  TF1* period1 = new TF1("period1","ftot - (1/2)* ( 1 + TMath::Erf((x-493.6)/92.1))",xminFit,xmaxFit);      
  TF1* period2 = new TF1("period2","ftot - (1/2)* ( 1 + TMath::Erf((x-496.6)/93.9))",xminFit,xmaxFit);      
  TF1* period3 = new TF1("period3","ftot - (1/2)* ( 1 + TMath::Erf((x-497.1)/96.1))",xminFit,xmaxFit);      
  TF1* period4 = new TF1("period4","ftot - (1/2)* ( 1 + TMath::Erf((x-494.5)/92.4))",xminFit,xmaxFit);      
  period1->SetLineColor(1);
  period2->SetLineColor(2);
  period3->SetLineColor(3);
  period4->SetLineColor(4);
  period1->Draw("lsame");
  period2->Draw("lsame");
  period3->Draw("lsame");
  period4->Draw("lsame");
  leg = new TLegend(0.55,0.71,0.93,0.91);
  leg->AddEntry(period1,"(tot - period1): run 257968-258432 (415/pb)","l");
  leg->AddEntry(period2,"(tot - period2): run 258434-258745 (461/pb)","l");
  leg->AddEntry(period3,"(tot - period3): run 258749-260425 (457/pb)","l");
  leg->AddEntry(period4,"(tot - period4): run 260426-260627 (471/pb)","l");
  leg->Draw();


  //-----------------------------------------------------------------------------

}
