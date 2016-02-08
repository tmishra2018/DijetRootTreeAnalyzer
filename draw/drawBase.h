#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TFrame.h>
#include <TLine.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TStyle.h>
#include <TMath.h>

#define maxIterations 300

using namespace std;

void DrawTH1(TFile* dataFile, const char *nameHisto, const char *XTitle, const char *YTitle, bool log, double xmin, double xmax){
  
  TH1F *h1 = (TH1F*)dataFile->Get(nameHisto );

  TString HistoName =TString::Format(nameHisto);

  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1    -> SetTitle(" ");
  //  h1    -> SetLineColor(kColor);
  if(log) gPad->SetLogy();
  h1    -> SetStats(0);
  h1    -> SetXTitle(XTitle);
  h1    -> SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.4);
  h1    -> GetXaxis()->SetRangeUser(xmin, xmax);
  h1    -> Draw();
  canvas ->SaveAs("Plot/"+HistoName+".png");
  canvas ->Destructor();
}

void FitTH1(TFile* dataFile, const char *nameHisto, const char *XTitle, const char *YTitle, bool log, double xmin, double xmax){
  
  TH1F *h1 = (TH1F*)dataFile->Get(nameHisto );

  TString HistoName =TString::Format(nameHisto);

  //  h1->Rebin(5);

  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1    -> SetTitle(" ");
  if(log) gPad->SetLogy();
  h1    -> SetStats(111); // with statistics
  h1    -> SetXTitle(XTitle);
  h1    -> SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.4);
  h1    -> GetXaxis()->SetRangeUser(xmin, xmax);

  
  double nSigma;
  nSigma=1.5;

  TF1* gaussian = new TF1( "gaussian", "gaus" );
  
  Float_t histMean = h1 -> GetMean();
  Float_t histRMS  = h1 -> GetRMS();

  gaussian->SetParameter(0, h1 ->GetMaximum());
  gaussian->SetParameter(1, histMean);
  gaussian->SetParameter(2, histRMS);
 
  gaussian->SetParLimits(1, 0., 2.*histMean);

  Float_t lowerBound = histMean - nSigma * histRMS;
  Float_t upperBound = histMean + nSigma * histRMS; 

  gaussian->SetRange(lowerBound, upperBound);

  h1 -> Fit(gaussian, "R");

  int n_iter = 4;
  for (int i = 0; i < n_iter; ++i) {
    
    Float_t lowerBound = gaussian->GetParameter(1) - nSigma * gaussian->GetParameter(2);
    Float_t upperBound = gaussian->GetParameter(1) + nSigma * gaussian->GetParameter(2);
    
    gaussian->SetRange(lowerBound, upperBound);
    
    h1 ->Fit(gaussian, "R");
  }

  double meanFit       =  gaussian->GetParameter(1);
  double meanFitErr   = gaussian->GetParError(1);
  double sigmaFit       =  gaussian->GetParameter(2);
  double sigmaFitErr  = gaussian->GetParError(2);
  
  std::cout<<"Mean "<<meanFit<<" +- "<<meanFitErr<<std::endl;
  std::cout<<"Sigma "<<sigmaFit<<" +- "<<sigmaFitErr<<std::endl;
  std::cout<<"ChiSquare/NDF "<<gaussian->GetChisquare()<<" / "<<gaussian->GetNDF()<<std::endl;
  
  TPaveText* fitlabel = new TPaveText(0.75, 0.75, 0.60, 0.60, "brNDC");
  fitlabel->SetTextSize(0.03);
  fitlabel->SetFillColor(0);
  TString Text_Mean = TString::Format("Mean: %.4f #pm %.4f", meanFit, meanFitErr);
  TString Text_Sigma = TString::Format("Sigma: %.4f #pm %.4f", sigmaFit, sigmaFitErr);
  fitlabel->AddText(Text_Mean);
  fitlabel->AddText(Text_Sigma);
  
  h1    -> Draw();
  //  fitlabel->Draw("same");
  
  canvas ->SaveAs("Plot/Fit_"+HistoName+".png");
  canvas ->Destructor();
}


void DrawTH2(TFile* dataFile, const char *nameHisto, const char *XTitle, const char *YTitle, bool log, bool colz, double xmin, double xmax, double ymin, double ymax){
  
  TH2F *h1 = (TH2F*)dataFile->Get(nameHisto );

  TString HistoName =TString::Format(nameHisto);

  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1    -> SetTitle(" ");
  //  h1    -> SetLineColor(kColor);
  if(log) gPad->SetLogz();
  h1    -> SetStats(0);
  h1    -> SetXTitle(XTitle);
  h1    -> SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.4);
  h1    -> GetYaxis()->SetRangeUser(ymin, ymax);
  h1    -> GetXaxis()->SetRangeUser(xmin, xmax);
  if(colz){ 
    gStyle->SetPalette(55);
    h1->SetContour(100);
    h1-> Draw("colz"); 
  }else{ 
    h1->Draw(); 
  }
  canvas ->SaveAs("Plot/"+HistoName+".png");
  canvas ->Destructor();
}


void DrawComparison(TFile* dataFile, const char *nameHisto, const char *XTitle, const char *YTitle, bool log, bool colz, double xmin, double xmax){

  TString HistoName =TString::Format(nameHisto);
  double ymax;

  TH1F *h1 = (TH1F*)dataFile->Get(HistoName+"HLT" );
  TH1F *h2 = (TH1F*)dataFile->Get(HistoName+"Reco" );
  
  TLegend* legend = new TLegend(0.78, 0.78, 0.90, 0.90);
  legend->SetTextFont(42);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextSize( 0.036);
  legend->AddEntry(h1, "HLT", "LP");
  legend->AddEntry(h2, "Reco", "LP");
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  if(log) gPad->SetLogy();
  h1->SetStats(0);
  h1->SetTitle(" ");
  h1->SetXTitle(XTitle);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetYTitle(YTitle);
  if (h1->GetMaximum() >= h2->GetMaximum() ) ymax = h1->GetMaximum() ; 
  if (h1->GetMaximum() < h2->GetMaximum() ) ymax = h2->GetMaximum() ; 
  h1 -> SetMaximum(ymax+ 0.1*ymax);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1->SetLineColor(kBlue);
  h1->SetMarkerColor(kBlue);
  h1->SetMarkerStyle(2);
  h2->SetLineColor(kRed);
  h2->SetMarkerColor(kRed);
  h2->SetMarkerStyle(3);
  h1->Draw();
  h2->Draw("same");
  legend->Draw();
  canvas ->SaveAs("Plot/Compare_"+HistoName+".png");
  canvas ->Destructor();
}

void DrawTProfile(TFile* dataFile, const char *nameProfile, const char *XTitle, const char *YTitle, bool log, double xmin, double xmax){
  
  TProfile *p1 = (TProfile*)dataFile->Get(nameProfile );

  TString ProfileName =TString::Format(nameProfile);

  TCanvas *canvas = new TCanvas("canvas","",800,800);
  p1    -> SetTitle(" ");
  //  p1    -> SetLineColor(kColor);
  if(log) gPad->SetLogy();
  p1    -> SetStats(0);
  p1    -> SetMarkerStyle(20);
  p1    -> SetXTitle(XTitle);
  p1    -> SetYTitle(YTitle);
  p1    -> GetYaxis()->SetTitleOffset(1.4);
  p1    -> GetXaxis()->SetRangeUser(xmin, xmax);
  p1    -> Draw();
  canvas ->SaveAs("Plot/"+ProfileName+".png");
  canvas ->Destructor();
}

void DrawCompareProfile(TFile* dataFile, const char *nameProfile, const char *XTitle, const char *YTitle, bool log, bool colz, double xmin, double xmax){

  TString ProfileName =TString::Format(nameProfile);
  double ymax;

  TProfile *p1 = (TProfile*)dataFile->Get(ProfileName+"HLT" );
  TProfile *p2 = (TProfile*)dataFile->Get(ProfileName+"Reco" );
  
  TLegend* legend = new TLegend(0.78, 0.78, 0.90, 0.90);
  legend->SetTextFont(42);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->SetFillStyle(0);
  legend->SetTextSize( 0.036);
  legend->AddEntry(p1, "HLT", "LP");
  legend->AddEntry(p2, "Reco", "LP");
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  //  if(log) gPad->SetLogy();
  p1->SetStats(0);
  p1->SetTitle(" ");
  p1->SetXTitle(XTitle);
  p1->GetXaxis()->SetRangeUser(xmin,xmax);
  p1->SetYTitle(YTitle);
  if (p1->GetMaximum() >= p2->GetMaximum() ) ymax = p1->GetMaximum() ; 
  if (p1->GetMaximum() < p2->GetMaximum() ) ymax = p2->GetMaximum() ; 
  p1 -> SetMaximum(ymax+ 0.1*ymax);
  p1    -> GetYaxis()->SetTitleOffset(1.55);
  p1->SetLineColor(kBlue);
  p1->SetMarkerColor(kBlue);
  p1->SetMarkerStyle(2);
  p2->SetLineColor(kRed);
  p2->SetMarkerColor(kRed);
  p2->SetMarkerStyle(3);
  p1->Draw();
  p2->Draw("same");
  legend->Draw();
  canvas ->SaveAs("Plot/Compare_"+ProfileName+".png");
  canvas ->Destructor();
}




//////////////////////////////// old stuff
void DrawComparison(const char *output_dir, const char *nameFile ,TH1D* h1, TH1D* h2, int xmin, int xmax, const char *XTitle, const char *YTitle, int h1Color, int h2Color){
   
  char * nameSaved ;
  nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, output_dir);
  //unisco szStringa3 con szStinga2
  strcat(nameSaved, nameFile);
  
  TCanvas *canvas = new TCanvas("canvas","",800,800);
  h1->SetStats(0);
  h1->SetTitle(" ");
  h1->SetXTitle(XTitle);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->SetYTitle(YTitle);
  h1    -> GetYaxis()->SetTitleOffset(1.55);
  h1->SetLineColor(h1Color);
  h1->SetMarkerColor(h1Color);
  h1->SetMarkerStyle(2);
  h2->SetLineColor(h2Color);
  h2->SetMarkerColor(h2Color);
  h2->SetMarkerStyle(3);
  h1->Draw();
  h2->Draw("same");
  //leg->Draw();
  canvas ->SaveAs(nameSaved);
  canvas ->Destructor();
  }

void DrawRatio(const char *output_dir, const char *nameFile,  TH1F* h1, TH1F* h2, const char *XTitle, const char *YTitle,  TLegend *legend ){

    char * nameSaved ;
    nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
    //copio stringa 1 in stringa nuova  
    strcpy(nameSaved, output_dir);
    //unisco szStringa3 con szStinga2
    strcat(nameSaved, nameFile);
    
  int size_h1 = h1 -> GetSize();
  //  int size_h2 = h2 -> GetSize();
  double binWidth = h1 -> GetBinWidth(1);
  double SemibinWidth = binWidth / 2.;
  double x_min = h1->GetXaxis()->GetBinCenter(1) - SemibinWidth;
  double x_max = h1->GetXaxis()->GetBinCenter(size_h1 - 1) - SemibinWidth;

  
  TCanvas *c = new TCanvas("c","",800,800);
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    h1->SetStats(0);          // No statistics on upper plot
    h1 -> SetTitle(" ");
    h1 -> SetYTitle(YTitle);
    h1 -> GetYaxis()->SetTitleOffset(1.55);
    h1 -> SetLineColor(kBlue);
    h1 -> SetMinimum(-0.001);
    double ymax;
    if (h1->GetMaximum() >= h2->GetMaximum() ) ymax = h1->GetMaximum() ; 
    if (h2->GetMaximum() > h1->GetMaximum() ) ymax = h2->GetMaximum() ; 
    h1 -> SetMaximum(ymax+ 0.1*ymax);
    //    h1 -> SetMaximum(ymax+200);
    //    h1 -> GetXaxis()->SetRangeUser(xmin,xmax);
    //    h1 -> SetMinimum(-0.001);
    h2 ->SetLineColor(kRed);
    h1 ->Draw();               // Draw h1
    h2 ->Draw("same");         // Draw h2 on top of h1
    legend -> Draw();   
 // lower plot will be in pad
    c->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad

    int npoints = 15;
    float x [npoints];
    float ex [npoints];
    float y [npoints];
    float yu95[npoints];
    float yd95[npoints];
    for(Int_t ii = 0; ii < npoints; ii++){
      x[ii] = ii * 1000 - 1010;
      ex[ii] = 0;
      y[ii] = 1;  
      yd95[ii] = 0.05;
      yu95[ii] = 0.05;
    }
    TGraph* Line = new TGraph (npoints, x, y);
    Line->SetLineColor(kRed);
    TGraphAsymmErrors *Band95 = new TGraphAsymmErrors(npoints, x, y, ex, ex, yd95, yu95);
    Band95->GetXaxis()->SetLimits(x_min, x_max);
    Band95->GetYaxis()->SetRangeUser(0.8, 1.2);
    Band95->SetTitle("");
    Band95->GetXaxis()->SetTitle(XTitle);
    Band95->GetXaxis()->SetTitleSize(20);
    Band95->GetXaxis()->SetTitleFont(43);
    Band95->GetXaxis()->SetTitleOffset(4.);
    Band95->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Band95->GetXaxis()->SetLabelSize(15);
    Band95->GetYaxis()->SetTitle("Ratio");
    Band95->GetYaxis()->SetTitleOffset(1.0);
    Band95->GetYaxis()->SetNdivisions(505);
    Band95->GetYaxis()->SetTitleSize(20);
    Band95->GetYaxis()->SetTitleFont(43);
    Band95->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Band95->GetYaxis()->SetLabelSize(15);
    Band95->SetFillColor(8);
    Band95->SetLineColor(10);
    Band95->Draw("a3");
    Line -> Draw("lSAME");
    // Define the ratio plot
    TH1D *h3 = (TH1D*)h2->Clone("h3");
    h3->Divide(h1);
    h3->SetLineColor(kBlack);
    h3->SetMarkerStyle(21);
    h3->Draw("same ep");       // Draw the ratio plot
    c -> SaveAs(nameSaved);
    c->Destructor();
  }

void DrawPull(const char *output_dir, const char *nameFile, TH1F* h1, TH1F *h2, const char *XTitle, const char *YTitle, TLegend *legend){

  bool verbose = false;  

  char * nameSaved ;
  nameSaved = new char[strlen(output_dir) + strlen(nameFile) +20] ;
  //copio stringa 1 in stringa nuova  
  strcpy(nameSaved, output_dir);
  //unisco szStringa3 con szStinga2
  strcat(nameSaved, nameFile);
  
  //calcolo il pull
  int size_h1 = h1 -> GetSize();
  int size_h2 = h2 -> GetSize();

  if(verbose)  std::cout<<"size_h1  "<<size_h1 << std::endl;

  double binWidth = h1 -> GetBinWidth(1);
  double SemibinWidth = binWidth / 2.;

  double x_min = h1->GetXaxis()->GetBinCenter(1) - SemibinWidth;
  double x_max = h1->GetXaxis()->GetBinCenter(size_h1 - 1) - SemibinWidth;

  if(verbose)  std::cout<<"x_min  "<<x_min << std::endl;
  if(verbose)  std::cout<<"x_max  "<<x_max << std::endl;
  
    if(size_h1 == size_h2){
    TH1D *H_Pull = new TH1D("H_Pull"," ", 1000, -5, 5);
    TH2D *H_Mass_Pull = new TH2D("H_Mass_Pull"," ",1000, x_min, x_max,100, -5, 5);
    
    for(int nbin = 0; nbin<size_h1; nbin++){
      double binCenter          = h1->GetXaxis()->GetBinCenter(nbin);
      double binContent_h1  = h1 -> GetBinContent(nbin);
      double binError_h1       = h1 -> GetBinError(nbin);
      double binContent_h2  = h2 -> GetBinContent(nbin);
      double binError_h2       = h2 -> GetBinError(nbin);
      if(binError_h1 ==0 && binError_h2 == 0) continue;
      double error_diff = sqrt( binError_h1 * binError_h1 + binError_h2* binError_h2 );
      double pull = (binContent_h2 - binContent_h1 ) / error_diff ;

      if(verbose){
	cout<<"nbin "<<nbin<<endl;
	cout<<"binCenter " <<binCenter<<endl;
	cout<<"binContent_h1 " <<binContent_h1<<endl;
	cout<<"binError_h1 " <<binError_h1<<endl;
	cout<<"binContent_h2 " <<binContent_h2<<endl;
	cout<<"binError_h2 " <<binError_h2<<endl;
	cout<<"error_diff " <<error_diff<<endl;
	cout<<"pull " <<pull<<endl;
      }

      H_Pull -> Fill(pull);
      H_Mass_Pull -> Fill(binCenter, pull);
    }
    
      //      TCanvas *canvas1 = new TCanvas("canvas1","",800,800);
      //H_Pull ->SetMarkerSize(1);
      //H_Pull ->SetMarkerStyle(21);
      //H_Pull->Draw("ep");
      //canvas -> SaveAs("H_Pull.png");
           
    TCanvas *Canv = new TCanvas("Canv","",800,800);
    // Upper plot will be in pad1
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetGridx();         // Vertical grid
    pad1->Draw();             // Draw the upper pad: pad1
    pad1->cd();               // pad1 becomes the current pad
    h1->SetStats(0);          // No statistics on upper plot
    h1 -> SetTitle(" ");
    h1 -> SetYTitle(YTitle);
    h1 -> GetYaxis()->SetTitleOffset(1.55);
    h1 -> SetLineColor(kBlue);
    h1 -> SetMinimum(-0.01);
    double ymax;
    if(verbose) cout<< "h1->GetMaximum() : "<< h1->GetMaximum() << endl;
    if(verbose) cout<< "h2->GetMaximum() : "<< h2->GetMaximum() << endl;
    if (h1->GetMaximum() >= h2->GetMaximum() ) ymax = h1->GetMaximum() ; 
    if (h2->GetMaximum() > h1->GetMaximum() ) ymax = h2->GetMaximum() ; 
    h1 -> SetMaximum(ymax+ 0.1*ymax);
    if(verbose) cout<< "ymax:  "<< ymax+0.1*ymax << endl;
    h2 ->SetLineColor(kRed);
    h1 ->Draw();             
    h2 ->Draw("same");   
    legend->Draw();
    // lower plot will be in pad
    Canv->cd();          // Go back to the main canvas before defining pad2
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();       // pad2 becomes the current pad
    //////////////////////
    int npoints = 15;
    float x [npoints];
    float ex [npoints];
    float y [npoints];
    float yu95[npoints];
    float yd95[npoints];
    float yu68[npoints];
    float yd68[npoints];
    for(Int_t ii = 0; ii < npoints; ii++){
      x[ii] = ii *1000 -1010;
      ex[ii] = 0;
      y[ii] = 0;  
      yd95[ii] = 1;
      yu95[ii] = 1;
      yd68[ii] = 2;
      yu68[ii] = 2;
    }
    TGraph* Line = new TGraph (npoints, x, y);
    Line->SetLineColor(kRed);
    TGraphAsymmErrors *Band95 = new TGraphAsymmErrors(npoints, x, y, ex, ex, yd95, yu95);
    TGraphAsymmErrors *Band68 = new TGraphAsymmErrors(npoints, x, y, ex, ex, yd68, yu68);
    //    Band68->GetXaxis()->SetRangeUser(x_min, x_max);
    Band68->GetXaxis()->SetLimits(x_min, x_max);

    Band68->GetYaxis()->SetRangeUser(-4, 4);
    //    Band68->SetStats(0);
    Band68->SetTitle(""); // Remove the ratio title  
    Band68->GetXaxis()->SetTitle(XTitle);
    Band68->GetXaxis()->SetTitleSize(20);
    Band68->GetXaxis()->SetTitleFont(43);
    Band68->GetXaxis()->SetTitleOffset(4.);
    Band68->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Band68->GetXaxis()->SetLabelSize(15);
    // Y axis ratio plot settings
    Band68->GetYaxis()->SetTitle("Pull");
    Band68->GetYaxis()->SetNdivisions(512);
    Band68->GetYaxis()->SetTitleSize(20);
    Band68->GetYaxis()->SetTitleFont(43);
    Band68->GetYaxis()->SetTitleOffset(1.55);
    Band68->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    Band68->GetYaxis()->SetLabelSize(15);
    Band68->SetFillColor(8);
    Band68->SetLineColor(10);
    Band68->Draw("a3");
    Band95->SetFillColor(5);
    Band95->SetLineColor(5);
    Band95->Draw("3SAME");

    H_Mass_Pull->SetLineColor(kBlack);
    H_Mass_Pull ->SetMarkerSize(1);
    H_Mass_Pull ->SetMarkerStyle(21);
    H_Mass_Pull->Draw("pSAME");
    Line -> Draw("lSAME");
    Canv -> SaveAs(nameSaved);
    delete H_Pull;
    delete H_Mass_Pull;
    delete Canv;
  }
}

///////////////////////////////////////////////////////////

 void Normalizer(TH1F* h1){
   //  char HistoName[200] = GetTitle(h1);
  int factor=1;
  h1 -> Scale(factor/h1->Integral());
  //  cout<<" Integral normalized to: "<<h1->Integral() << endl;

}

////////////////////////////////////////////////////////////
















