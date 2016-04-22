#include <stdlib.h>
#include <sys/stat.h>  
#include <time.h>

#include "TParameter.h"
#include "TPaveText.h"
#include "TError.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include <TMap.h>
#include <TMath.h>
#include "ptBinning.h"
#include "etaBinning.h"
#include <TColor.h>
#include <TStyle.h>

using namespace std;

//Crystal ball function for signal, parameters are:  0:normalizatio, 1:mean, 2:sigma, 3:alpha, 4:n;
Double_t CrystalBall(Double_t *x,Double_t *par) {

  Double_t t = (x[0]-par[1])/par[2];
  //  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);

  if (x[0] >= par[1] - absAlpha*par[2]) {  
    //  if (t >= -absAlpha) {
    return par[0]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[4]/absAlpha, par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[4]/absAlpha - absAlpha;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
}
////////////////////////////////
// double Crystal ball
Double_t doubleCrystalBall(Double_t *x,Double_t *par) {
  // par[0]:  N, normalization
  // par[1]:  mean
  // par[2]:  sigma
  // par[3]:  alphaR, crossover point
  // par[4]:  nR, length of tail
  // par[5]:  alphaL, crossover point
  // par[6]:  nL, length of tail

  double cb = 0.0;
  Double_t t = (x[0]-par[1])/par[2] ;
  Double_t absAlphaR = fabs(par[3]);
  Double_t absAlphaL = fabs(par[5]);

  if (x[0] > par[1] + absAlphaR*par[2]) {  
    Double_t a = pow(par[4]/absAlphaR, par[4]) * exp(-0.5*absAlphaR*absAlphaR);
    Double_t b = par[4]/absAlphaR - absAlphaR;    
    cb= ( a/pow(b - t, par[4]) );
  }
  else if (x[0] < par[1] - absAlphaL*par[2]) {  
    Double_t a = pow(par[6]/absAlphaL, par[6]) * exp(-0.5*absAlphaL*absAlphaL);
    Double_t b = par[6]/absAlphaL - absAlphaL;
    cb = (a/pow(b - t, par[6]) );
  }
  else{
    cb =exp(-0.5*t*t);
  }
  
  if (par[0] > 0.) {
    cb *= par[0];
  }
  return cb;
}

//////////////////////////////////////////

map<std::string , TGraphErrors*> userTGraphs_;
void CreateAndFillUserTGraph(const char* nameAndTitle, Int_t npoint, Double_t xvalue, Double_t yvalue, Double_t xErr, Double_t yErr)
{
  map<std::string , TGraphErrors*>::iterator nh_h = userTGraphs_.find(std::string(nameAndTitle));
  TGraphErrors * h;
  if( nh_h == userTGraphs_.end() )
    {
      h = new TGraphErrors(0);
      userTGraphs_[std::string(nameAndTitle)] = h;
      h->SetPoint(npoint, xvalue, yvalue);
      h->SetPointError(npoint, xErr, yErr);
    }
  else
    {
      nh_h->second->SetPoint(npoint, xvalue, yvalue);
      nh_h->second->SetPointError(npoint, xErr, yErr);
    }
}


double pTMeanArray[] = { 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500};
int binnum = sizeof(pTMeanArray)/sizeof(double) - 1;
double etaBins[] = {-2.50, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.50};  
int  binnumEta = sizeof(etaBins)/sizeof(double) -1;                                                                                                      

TH2D *Histo_Corrections = new TH2D("Histo_Corrections","Corrections", binnumEta, etaBins, binnum, pTMeanArray);

/////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

  std::string outputDir = "WideJets_Plot";                 
  mkdir(outputDir.c_str(), 0777); 

  //  TString dataFileName_DoubleMu;
  //  dataFileName_DoubleMu = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_DoubleMu.root");
  // TFile* dataFile_DoubleMu = TFile::Open(dataFileName_DoubleMu);
  // if (dataFile_DoubleMu) {
  //  std::cout << "Opened data file '" << dataFileName_DoubleMu << "'." << std::endl;
  // }

  TString dataFileName_HT450;
  dataFileName_HT450 = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_HT450.root");
  TFile* dataFile_HT450 = TFile::Open(dataFileName_HT450);
  if (dataFile_HT450) {
    std::cout << "Opened data file '" << dataFileName_HT450 << "'." << std::endl;
  }

  TH1D *pT_WideJetHLT = (TH1D*)dataFile_HT450->Get("pT_WideJetHLT");

  if( !pT_WideJetHLT) cout<<"Not found" <<endl;

  gErrorIgnoreLevel = kWarning;

  // procedura
  // prendo gli istogramma 1D del bias
  // per ogni bin di pT e eta
  // li fitto con una cristal ball
  // poi grafico il picco in funzione del pT in ogni bin di eta
  // fitto con una funzione logaritmica
  // mi creo l'istogramma 2D (eta,pT) con correzione sull'asse z

  //////////////////////// pTBias 1D
  EtaBinning mEtaBinning;
  size_t etaBinningSize = mEtaBinning.size();
  PtBinning mPtBinning;
  size_t pTBinningSize = mPtBinning.size();

  cout<< etaBinningSize << endl;

  for (size_t ii = 0; ii < etaBinningSize; ii++) { 
    // for (size_t ii = 9; ii < 10; ii++) {
    for (size_t jj = 0; jj < pTBinningSize; jj++) {

      cout<< "ii = "<< ii <<endl;
      cout<< "jj = "<< jj <<endl;

      if(jj < 2) continue; // pT<300
      if(jj > 7) continue; // pT>900

      std::string etaName = mEtaBinning.getBinName(ii);
      std::pair<float, float> ptBins = mPtBinning.getBinValue(jj);	       
      
      std::string HistoName = TString::Format("pTBias_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
      
      cout<< HistoName.c_str() << endl;
      
      TH1D *h;
      h = (TH1D*)dataFile_HT450->Get( HistoName.c_str() );
            
      double hN = h->GetMaximum();
      double hMean = h->GetMean();
      double hRMS = h->GetRMS();

      //histo->GetMean() - sigma * histo->GetRMS(), histo->GetMean() + sigma * histo->GetRMS() 

      //      TF1 *crystal = new TF1("crystal", CrystalBall, hMean-1*hRMS , hMean+1*hRMS, 5);
      TF1 *crystal = new TF1("crystal", CrystalBall, -0.1, 0.15, 5);
      //      TF1 *crystal = new TF1("crystal", CrystalBall, -0.05, 0.1, 5);
      crystal->SetParameters(hN, hMean, hRMS, 1.0, 2.0);
      crystal->SetParNames("N", "Mean", "sigma","#alpha","n");
      
      //      TF1 *doublecrystal = new TF1("doublecrystal", doubleCrystalBall, hMean-1*hRMS , hMean+1*hRMS, 7);
      TF1 *doublecrystal = new TF1("doublecrystal", doubleCrystalBall, -0.2 , 0.2, 7);
      doublecrystal->SetParameters(hN, hMean, hRMS, 3, 3, 3, 3);
      doublecrystal->SetParNames("N", "Mean", "sigma","#alphaR","nR","#alphaL", "nL");
      
      TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
      h -> SetMarkerStyle(20);
      h -> SetMarkerSize(1.5);
      h -> GetXaxis()->SetRangeUser(-0.4, 0.4);
      // h -> Fit(crystal, "R");
      h -> Fit(doublecrystal, "R");
      h -> Draw();
      c1 -> SaveAs( Form("WideJets_Plot/%s.png", HistoName.c_str() ) );
      c1-> Destructor();

      /////////////////////////////////////////////////
      
      pT_WideJetHLT ->GetXaxis()->SetRangeUser(ptBins.first, ptBins.second);
      double pTMean = pT_WideJetHLT->GetMean();
      std::cout<< "Bin " << ptBins.first<< "-"<<ptBins.second<<" -> Mean  "<< pTMean << std::endl; 
      
      double pTBias = doublecrystal->GetParameter(1);
      double pTBiasErr = doublecrystal->GetParError(1);
      
      cout<< "pT mean = "<< pTMean <<endl;
      cout<< "pT Bias = "<< pTBias <<endl;
      cout<< "sigma(pT Bias) = "<< pTBiasErr <<endl;

      std::string TGraphName = TString::Format("pTBias_vs_pT_etabin_%i", ii).Data();
      cout<< TGraphName.c_str() << endl;
      CreateAndFillUserTGraph(TGraphName.c_str(), jj, pTMean, pTBias, 0, pTBiasErr);
      
    }
  }
  

  for(int ii = 0; ii < etaBinningSize ; ii++){ 
   
    std::string TGraphName = TString::Format("pTBias_vs_pT_etabin_%i", ii).Data();
    cout<< TGraphName.c_str() << endl;
    map<std::string , TGraphErrors*>::iterator nh_h = userTGraphs_.find( TGraphName );
    if( nh_h == userTGraphs_.end()) continue;

    TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)*log(x)) + [1]*(log(x)*log(x)) + [2]*log(x) + [3]", 300, 1000);
    //TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)) + [1]*log(x) + [2]", 300, 1000);

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    nh_h->second->SetMarkerStyle(20);
    nh_h->second->SetMarkerSize(1.5);
    nh_h->second->GetXaxis()->SetRangeUser(300, 1000);
    nh_h->second->SetMaximum(0.1);
    nh_h->second->SetMinimum(-0.1);
    nh_h->second->Fit(functionFit, "R");
    nh_h->second->Draw("ZAP");
    c1 -> SaveAs( Form("WideJets_Plot/%s.png", TGraphName.c_str() ) );
    c1 -> Destructor();

    double Par0       =  functionFit->GetParameter(0);
    double Par1       =  functionFit->GetParameter(1);
    double Par2       =  functionFit->GetParameter(2);                                                                                                              
    double Par3       =  functionFit->GetParameter(3); 
    double NDF       =  functionFit->GetNDF();
    double Chisquare      =  functionFit->GetChisquare();
    cout<< "Par 0 = " << Par0 << endl;
    cout<< "Par 1 = " << Par1 << endl;
    cout<< "Par 2 = " << Par2 << endl;
    cout<< "Par 3 = " << Par3 << endl;                                                                                                                                                                 
    cout<< "Chisquare / NDF = " << Chisquare<<"/"<<NDF << endl;

    for(int jj=0; jj < binnum; jj++){
      
      int pTMin = pTMeanArray[jj];
      int pTMax = pTMeanArray[jj+1];
      double pTMean= (pTMin + pTMax) /2. ;
      cout<<"pT Mean = " << pTMean <<endl;

      //      double pTBiasFit = Par0*(log(pTMean)*log(pTMean)) + Par1*log(pTMean) + Par2;
      double pTBiasFit = Par0*(log(pTMean)*log(pTMean)*log(pTMean)) + Par1*log(pTMean)*log(pTMean) + Par2*log(pTMean) + Par3;
      cout<<"pT Bias from fit = " << pTBiasFit <<endl;

      Histo_Corrections->SetBinContent(ii+1, jj+1, pTBiasFit);
    }
  } 

    TCanvas *c2 = new TCanvas("c2", "c2", 900, 800);
    Histo_Corrections->SetStats(0);
    gStyle->SetPalette(55);
    Histo_Corrections->SetContour(100);
    Histo_Corrections->Draw("colz");
    c2->SetRightMargin(0.13);
    c2->SaveAs("WideJets_Plot/Corrections.png");

    TFile* outputFile= new TFile("WideJetsCorrections_JEC_HLT_v7.root","RECREATE");
    outputFile->cd();
    Histo_Corrections -> Write();
    outputFile->Close();

   
  return 0;
}

      /*      
	      if( jj < 2){ // pT<300
	      h = (TH1D*)dataFile_DoubleMu->Get( HistoName.c_str() );
	      cout << " Lo prendo dal DoubleMu" << endl;
	      }else if(jj < 6){ // pT<700 due to low statistics
	      h = (TH1D*)dataFile_HT450->Get( HistoName.c_str() );
	      cout << " Lo prendo dal HT450" << endl;
	      }else{ continue; }
      */
