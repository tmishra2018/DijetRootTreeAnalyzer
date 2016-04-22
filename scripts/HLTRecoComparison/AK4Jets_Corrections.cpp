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
  if (par[3] < 0) t = -t;

  Double_t absAlpha = fabs((Double_t)par[3]);

  if (t >= -absAlpha) {
    return par[0]*exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(par[4]/absAlpha, par[4])*exp(-0.5*absAlpha*absAlpha);
    Double_t b= par[4]/absAlpha - absAlpha;

    return par[0]*(a/TMath::Power(b - t, par[4]));
  }
}
////////////////////////////////
Double_t doubleCrystalBall(Double_t *x,Double_t *par) {
  // par[0]:  N, normalization
  // par[1]:  mean
  // par[2]:  sigma
  // par[3]:  alphaR, crossover point
  // par[4]:  nR, length of tail
  // par[5]:  alphaL, crossover point
  // par[6]:  nL, length of tail

  Double_t cb = 0.0;
  Double_t t = (x[0]-par[1])/par[2] ;
  Double_t absAlphaR = fabs((Double_t)par[3]);
  Double_t absAlphaL  = fabs((Double_t)par[5]);
 
  if(x[0] >= par[1] + absAlphaR*par[2]) {  
    //    cout<< "Sto qua 1"<<endl;
    Double_t a = pow(par[4]/absAlphaR, par[4]) * exp(-0.5*absAlphaR*absAlphaR);
    Double_t b = (par[4]/absAlphaR) - absAlphaR;    
    cb= ( a/pow(b - t, par[4]) );
  }
  else if(x[0] <= par[1] - absAlphaL*par[2]) {  
    //    cout<< "Sto qua 2"<<endl;
    Double_t a = pow(par[6]/absAlphaL, par[6]) * exp(-0.5*absAlphaL*absAlphaL);
    Double_t b = (par[6]/absAlphaL) - absAlphaL;
    cb = ( a/pow(b - t, par[6]) );
  }
  else{
    //    cout<< "Sto qua Gauss"<<endl;
    cb =exp(-0.5*t*t);
  }
  
  if (par[0] > 0.) {
    cb *= par[0];
  }
  return cb;
}
/////////////////////////////////////////
Double_t ExpGaussExp(Double_t *x,Double_t *par) {
  // par[0]:  N, normalization
  // par[1]:  mean
  // par[2]:  sigma
  // par[3]:  kL
  // par[4]:  kR

  Double_t function = 0.0;
  Double_t exponent = 0.0;
  Double_t t = (x[0]-par[1])/par[2] ;
  Double_t abskL  = fabs((Double_t)par[3]);
  Double_t abskR = fabs((Double_t)par[4]);
  
 if (x[0] <= par[1] - abskL*par[2]) {  
   exponent = 0.5*abskL*abskL + abskL*t ;
    function = exp(exponent);
 }else if( x[0] >= par[1] + abskR*par[2] ) {  
   exponent = 0.5*abskR*abskR - abskR*t ;
    function = exp(exponent);
 }else{
   function =exp(-0.5*t*t);
 }
 
 if (par[0] > 0.) {
   function *= par[0];
 }
 return function;
}
///////////////////////////////////////

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


double pTMeanArray[] = {30, 100, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000};
int binnum = sizeof(pTMeanArray)/sizeof(double) - 1;
//double etaBins[] = {-2.50, -1.50, -1.25, -1.00, -0.75, -0.50, -0.25, 0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 2.50};  
double etaBins[] = {0.00, 0.50, 1.00, 1.50, 2.50};  
int  binnumEta = sizeof(etaBins)/sizeof(double) -1;                                                                                                      

TH2D *Histo_Corrections = new TH2D("Histo_Corrections","Corrections", binnumEta, etaBins, binnum, pTMeanArray);

/////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

  std::string outputDir = "AK4Jets_Plot";                 
  mkdir(outputDir.c_str(), 0777); 

  TString dataFileName_DoubleMu;
  //  dataFileName_DoubleMu = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_DoubleMu.root");
  //  dataFileName_DoubleMu = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_WideEtaBin_DoubleMu.root");
  dataFileName_DoubleMu = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_WideEtaBin_DoubleMu.root");
  TFile* dataFile_DoubleMu = TFile::Open(dataFileName_DoubleMu);
  if (dataFile_DoubleMu) {
    std::cout << "Opened data file '" << dataFileName_DoubleMu << "'." << std::endl;
  }
  TString dataFileName_HT450;
  //  dataFileName_HT450 = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_HT450.root");
  //  dataFileName_HT450 = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_WideEtaBin_HT450.root");
  dataFileName_HT450 = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts10GeV_WideEtaBin_HT450.root");
  TFile* dataFile_HT450 = TFile::Open(dataFileName_HT450);
  if (dataFile_HT450) {
    std::cout << "Opened data file '" << dataFileName_HT450 << "'." << std::endl;
  }

  TH1D *pT_AK4JetHLT = (TH1D*)dataFile_HT450->Get("AK4_pT_JetHLT");

  gErrorIgnoreLevel = kWarning;

  TF1* lineup1      = new TF1("lineup1", " 0.01", 0, 3000);
  TF1* linedown1 = new TF1("linedown1", "-0.01", 0, 3000);
  lineup1->SetLineColor(kBlack);                                                                                                                                        
  lineup1->SetLineStyle(2);                                                                                                                                             
  linedown1->SetLineColor(kBlack);                                                                                                                                      
  linedown1->SetLineStyle(2);

  TF1* lineup2      = new TF1("lineup2", " 0.02", 0, 3000);
  TF1* linedown2 = new TF1("linedown2", "-0.02", 0, 3000);
  lineup2->SetLineColor(kBlack);                                                                                                                                        
  lineup2->SetLineStyle(2);                                                                                                                                             
  linedown2->SetLineColor(kBlack);                                                                                                                                      
  linedown2->SetLineStyle(2);


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
    // for (size_t ii = 0; ii < 1; ii++) {
    for (size_t jj = 0; jj < pTBinningSize; jj++) {
      // for (size_t jj = 8; jj < 9; jj++) {
      
      std::string etaName = mEtaBinning.getBinName(ii);
      std::pair<float, float> ptBins = mPtBinning.getBinValue(jj);	       
      
      std::string HistoName = TString::Format("AK4_pTBias_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();

      cout<< "("<<ii << ";"<< jj<<")" << endl;
      cout<< HistoName.c_str() << endl;


      TH1D *h;

      double sigma;
      double sigma2;
      
      if( jj < 2){ // pT<300
	sigma = 3;
	sigma2 = 4;
	h = (TH1D*)dataFile_DoubleMu->Get( HistoName.c_str() );
	cout << " Lo prendo dal DoubleMu" << endl;
      }else if( jj<8 ){
	sigma = 4; // 4
	sigma2 = 5; //5
	h = (TH1D*)dataFile_HT450->Get( HistoName.c_str() );
	cout << " Lo prendo dal HT450" << endl;
      }else{
	sigma = 4;
	sigma2 = 4; //3
	h = (TH1D*)dataFile_HT450->Get( HistoName.c_str() );
	cout << " Lo prendo dal HT450" << endl;
      }
      
      if( jj > 8) continue; // >8 prima
      if( ii==3 && jj > 6) continue; // ultimo fit in endcap non funziona
      if(!h) continue;
      int Nentries = h->GetEntries();
      if(Nentries<100) continue; //200

      h->Rebin(2);

      double hN = h->GetMaximum();
      double hMean = h->GetMean();
      double hRMS = h->GetRMS();

      /*
	TF1 *expon = new TF1("expon", "[0] + [1]*exp(-[2]*x)", 0.1, 0.4);
	h -> Fit(expon, "R");
      */
      
      /*
      TF1 *crystal = new TF1("crystal", CrystalBall, hMean-2*hRMS , hMean+1*hRMS, 5);
      crystal->SetParameters(hN, hMean, hRMS, 1.0, 2.0);
      crystal->SetParNames("N", "Mean", "sigma","#alpha","n");
      crystal->SetParLimits(0, hN-0.05*hN, hN+0.05*hN);
      h -> Fit(crystal, "R");
      */

      /*
	TF1 *doublecrystal = new TF1("doublecrystal", doubleCrystalBall, hMean-3*hRMS , hMean+4*hRMS, 7);
	doublecrystal->SetParNames("N", "Mean", "sigma","#alphaR","nR","#alphaL", "nL");
	doublecrystal->SetParameters(hN, 0.01, 0.03, 1.2, 30., 15, 40.);
	doublecrystal->SetParLimits(1, 0.005, 0.015);
	doublecrystal->SetParLimits(2, 0.015, 0.035);
	doublecrystal->SetParLimits(3, 0 , 1.4);
	doublecrystal->SetParLimits(4, 1, 10);
	doublecrystal->SetParLimits(5, 0, 10);
	doublecrystal->SetParLimits(6, 0, 15);
	h -> Fit(doublecrystal, "R");
      */

      /*
      // secondo giro
      TF1 *doublecrystal2 = new TF1("doublecrystal2", doubleCrystalBall, doublecrystal->GetParameter(1)-2* doublecrystal->GetParameter(2),doublecrystal->GetParameter(1)+1*doublecrystal->GetParameter(2), 7);      
      doublecrystal2->SetLineColor(kBlue);
      doublecrystal2->SetParNames("N", "Mean", "sigma","#alphaR","nR","#alphaL", "nL");
      doublecrystal2->SetParameters(doublecrystal->GetParameter(0), doublecrystal->GetParameter(1), doublecrystal->GetParameter(2),  doublecrystal->GetParameter(3),  doublecrystal->GetParameter(4),  doublecrystal->GetParameter(5), doublecrystal->GetParameter(6));  
      h -> Fit(doublecrystal2, "+R");

      // terzo giro
      TF1 *doublecrystal3 = new TF1("doublecrystal3", doubleCrystalBall, doublecrystal2->GetParameter(1)-2* doublecrystal2->GetParameter(2),doublecrystal2->GetParameter(1)+1*doublecrystal2->GetParameter(2), 7);      
      doublecrystal3->SetLineColor(kGreen);
      doublecrystal3->SetParNames("N", "Mean", "sigma","#alphaR","nR","#alphaL", "nL");
      doublecrystal3->SetParameters(doublecrystal2->GetParameter(0), doublecrystal2->GetParameter(1), doublecrystal2->GetParameter(2),  doublecrystal2->GetParameter(3),  doublecrystal2->GetParameter(4),  doublecrystal2->GetParameter(5), doublecrystal2->GetParameter(6));  
      h -> Fit(doublecrystal3, "+R");
      */      
      /*
      // provo a fittare solo la coda --> non viene manco cosi
      TF1 *Function = new TF1("Function", "[0]*(  (pow([4]/fabs([3]),[4])*exp(-0.5*[3]*[3])) / pow( ([4]/fabs([3])- fabs([3]) - (x - [1])/[2]),[4]))", 0.03 , 0.4);
      Function->SetParNames("N", "Mean", "sigma","aR","nR");
      Function->SetParameters(hN, hMean, hRMS, 1 , 3.);
      Function->SetParLimits(0, hN-0.05*hN, hN+0.05*hN);
      h -> Fit(Function, "R");
      */

      TF1 *ExponentialGauss = new TF1("ExponentialGauss", ExpGaussExp, hMean-sigma*hRMS , hMean+sigma*hRMS, 5);
      ExponentialGauss->SetParNames("N", "Mean", "sigma","kL","kR");
      ExponentialGauss->SetParameters(hN, hMean, hRMS, 1., 1.);
      ExponentialGauss->SetParLimits(0, hN-0.2*hN, hN+0.2*hN); // se non funziona rimettilo
      if( jj>= 8 ) ExponentialGauss->SetParLimits(2, 0.005, 0.025);
      // if( jj>= 9 ) ExponentialGauss->SetParLimits(1, -0.01, 0.01);
      // ExponentialGauss->FixParameter(0, hN);
      h -> Fit(ExponentialGauss, "R");

      // secondo giro
      TF1 *ExponentialGauss2 = new TF1("ExponentialGauss2", ExpGaussExp, ExponentialGauss->GetParameter(1)-sigma2*ExponentialGauss->GetParameter(2) , ExponentialGauss->GetParameter(1)+sigma2*ExponentialGauss->GetParameter(2), 5);
      ExponentialGauss2 -> SetParNames("N", "Mean", "sigma","kL","kR");
      ExponentialGauss2 -> SetParameters( ExponentialGauss->GetParameter(0), ExponentialGauss->GetParameter(1), ExponentialGauss->GetParameter(2), ExponentialGauss->GetParameter(3), ExponentialGauss->GetParameter(4));
      h -> Fit(ExponentialGauss2, "R");

      TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
      //      gPad->SetLogy();
      // h -> SetStats(0);
      gStyle->SetOptFit(1111);
      h -> SetMarkerStyle(20);
      h -> SetMarkerSize(1.5);
      h -> GetXaxis()->SetRangeUser(-0.4, 0.4);
      h  -> Draw();
      c1 -> SaveAs( Form("AK4Jets_Plot/%s.png", HistoName.c_str() ) );
      c1 -> Destructor();

      /////////////////////////////////////////////////

      //      int pTMin = ptBins.first;
      //      int pTMax = ptBins.second;
      //      double pTMean= (pTMin + pTMax) /2. ;
           
      pT_AK4JetHLT ->GetXaxis()->SetRangeUser(ptBins.first, ptBins.second);                                                                                                
      double pTMean = pT_AK4JetHLT->GetMean();                                                                                                                             
      std::cout<< "Bin " << ptBins.first<< "-"<<ptBins.second<<" -> Mean  "<< pTMean << std::endl;

      double pTBias     = ExponentialGauss2 -> GetParameter(1);
      double pTBiasErr = ExponentialGauss2 -> GetParError(1);

      cout<< "ii = "<< ii <<endl;
      cout<< "pT mean = "<< pTMean <<endl;
      cout<< "pT Bias = "<< pTBias <<endl;
      cout<< "sigma(pT Bias) = "<< pTBiasErr <<endl;

      std::string TGraphName = TString::Format("AK4_pTBias_vs_pT_etabin_%s", etaName.c_str()).Data();
      cout<< TGraphName.c_str() << endl;
      CreateAndFillUserTGraph(TGraphName.c_str(), jj, pTMean, pTBias, 0, pTBiasErr);

    }
  }

  cout << endl;
  cout<< "Fitto i TGraph" << endl;
  
  for(int ii = 0; ii < etaBinningSize ; ii++){ 
    // for(int ii = 0; ii < 1 ; ii++){ 

    std::string etaName = mEtaBinning.getBinName(ii);
      
    std::string TGraphName = TString::Format("AK4_pTBias_vs_pT_etabin_%s", etaName.c_str()).Data();
    map<std::string , TGraphErrors*>::iterator nh_h = userTGraphs_.find( TGraphName );
    if( nh_h == userTGraphs_.end()) continue;

    cout<< TGraphName.c_str() << endl;
    
    // TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)*log(x)) + [1]*(log(x)*log(x)) + [2]*log(x) + [3]", 65, 700);
    // TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)) + [1]*(log(x)) + [2]", 60, 700);
    // TF1* functionFit = new TF1( "functionFit", "[0]*log(x) + [1]", 65, 700);

    TF1* functionFit = new TF1("functionFit"," ([0]/2)*( [3] + TMath::Erf( (x-[1])/[2]) )", 0, 1000);
    //    TF1* functionFit = new TF1("functionFit"," [0]*(TMath::Erf( (x-[1])/[2])) ", 0, 1000);
    functionFit -> SetParameters(0.02, 100, 200, -0.01);
    functionFit -> SetParLimits(0, 0.005, 0.05); //0.04
    functionFit -> SetParLimits(1, 150, 300);
    functionFit -> SetParLimits(2, 100, 350); //125

    if(ii==3){
      functionFit -> SetParLimits(0, 0.005, 0.03); //0.04
      functionFit -> SetParLimits(1, 150, 300);
      functionFit -> SetParLimits(2, 20, 350); //125
    }
 
   // functionFit -> SetParLimits(3, 0, 10000);
    

    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    //    gPad -> SetLogx();
    gStyle->SetOptFit(0);
    nh_h->second->SetMarkerStyle(20);
    nh_h->second->SetMarkerSize(1.5);
    nh_h->second->SetMaximum(0.05);
    nh_h->second->SetMinimum(-0.05);
    nh_h->second->Fit(functionFit, "R");
    nh_h->second->Draw("ZAP");
    lineup1 -> Draw("same");
    linedown1 -> Draw("same");
    lineup2 -> Draw("same");
    linedown2 -> Draw("same");
    c1 -> SaveAs( Form("AK4Jets_Plot/%s.png", TGraphName.c_str() ) );
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
      //double pTBiasFit = Par0*(log(pTMean)*log(pTMean)*log(pTMean)) + Par1*log(pTMean)*log(pTMean) + Par2*log(pTMean) + Par3;

      //    functionFit  = ([0]/2)   *( [3]  + TMath::Erf( (x-[1])/[2]) )", 0, 1000);
      double pTBiasFit = (Par0/2)*( Par3 + TMath::Erf( (pTMean -Par1)/Par2));
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
    c2->SaveAs("AK4Jets_Plot/Corrections.png");

    TFile* outputFile= new TFile("AK4JetsCorrections_JEC_HLT_v7.root","RECREATE");
    outputFile->cd();
    Histo_Corrections -> Write();
    outputFile->Close();

 
  return 0;
}
