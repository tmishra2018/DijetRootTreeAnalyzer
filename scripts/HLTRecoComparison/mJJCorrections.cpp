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
#include "mJJBinning.h"
#include <TColor.h>
#include <TStyle.h>

using namespace std;

bool ClosureTest = true;
bool ComputeCorrections = false;

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

//TH2D *Histo_Corrections = new TH2D("Histo_Corrections","Corrections", binnumEta, etaBins, binnum, pTMeanArray);
TH1D *Histo_Corrections = new TH1D("Histo_Corrections","Corrections", 26, 400, 3000);

/////////////////////////////////////////////////////////////////
int main(int argc, char* argv[]) {

  std::string outputDir = "Plot_mJJ";                 
  mkdir(outputDir.c_str(), 0777); 
  /*
    TString dataFileName_DoubleMu;
    //  dataFileName_DoubleMu = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts_pTCut10GeV_FinerPtBin_DoubleMu.root");
    dataFileName_DoubleMu = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts_pTCut10GeV_ChangedEtaPtBin_DoubleMu.root");
    TFile* dataFile_DoubleMu = TFile::Open(dataFileName_DoubleMu);
    if (dataFile_DoubleMu) {
    std::cout << "Opened data file '" << dataFileName_DoubleMu << "'." << std::endl;
    }
  */

  TString dataFileName_HT450;
  //  dataFileName_HT450 = TString::Format("../../output/HLTReco_Comparison/rootFile_JEC_HLT_v7_RecoLowerCuts_pTCut10GeV_FinerPtBin_HT450.root");
  //  dataFileName_HT450 = TString::Format("../../output/HLTReco_mJJComparison/rootFile_JEC_HLT_v7_SamePtCuts_HT450.root");

  if(ComputeCorrections)  dataFileName_HT450 = TString::Format("../../output/HLTReco_mJJComparison/rootFile_JEC_HLT_v7_SamePtCuts_mJJCut400_L1HTT150_HT450.root");
  if(ClosureTest)  dataFileName_HT450 = TString::Format("../../output/HLTReco_mJJClosureTest/rootFile_mJJClosureTest_mJJCut400_L1HTT150_HT450_Check.root");

  TFile* dataFile_HT450 = TFile::Open(dataFileName_HT450);
  if (dataFile_HT450) {
    std::cout << "Opened data file '" << dataFileName_HT450 << "'." << std::endl;
  }

  gErrorIgnoreLevel = kWarning;

  TF1* lineup      = new TF1("lineup", " 0.02", -2.5, 3000);
  TF1* linedown = new TF1("linedown", "-0.02", -2.5, 3000);
  lineup->SetLineColor(kBlack);
  lineup->SetLineStyle(2);
  linedown->SetLineColor(kBlack);
  linedown->SetLineStyle(2);


  if(ComputeCorrections){
    mJJBinning mMjjBinning;
    size_t mJJBinningSize = mMjjBinning.size();

    // disegno bias 1D 
    for(size_t jj = 0; jj < mJJBinningSize ; jj++){ 
    
      std::pair<float, float> mJJBins = mMjjBinning.getBinValue(jj);	       
      std::string Histo = TString::Format("mJJBias_mJJ_%i_%i", (int) mJJBins.first, (int) mJJBins.second ).Data();
      cout<< Histo.c_str() << endl;
    
      TH1D *h  = (TH1D*)dataFile_HT450->Get( Histo.c_str() );
      if( !h) continue;
    
      TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
      h->GetXaxis()->SetRangeUser(-0.5, 0.5);
      h->SetMarkerStyle(20);
      h->SetMarkerSize(1.5);
      h->Draw();
      c1->SaveAs( Form("Plot_mJJ/%s.png", Histo.c_str() ) );
    }

    ///////
    TH1D *h2  = (TH1D*)dataFile_HT450->Get( "mJJBias_vs_Eta" );
    TCanvas *c3 = new TCanvas("c3", "c3", 800, 800);
    h2->SetStats(0);
    h2->SetTitle("");
    h2->SetXTitle("#eta");
    h2->SetYTitle("mJJ Bias");
    h2-> GetYaxis()->SetTitleOffset(1.4);
    h2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    h2->SetMaximum(0.1);
    h2->SetMinimum(-0.1);
    h2 -> SetMarkerStyle(20);
    h2 -> SetMarkerSize(1.5);
    c3->SetLeftMargin(0.13);
    h2 -> Draw();
    lineup->Draw("same");
    linedown->Draw("same");
    c3-> SaveAs( "Plot_mJJ/mJJBias_vs_Eta.png" );   

    //Histo da fittare per le correzioni
    TH1D *h  = (TH1D*)dataFile_HT450->Get( "mJJBias_vs_Mjj" );
  
    //  TF1* functionFit = new TF1( "functionFit", "-[0]*(log(x)*log(x)) + [1]*log(x)+[2]", 500, 2500); // quadratica log
    TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)*log(x)) + [1]*(log(x)*log(x)) + [2]*log(x) + [3]", 400, 2100); // quadratica log
    // functionFit->SetParLimits(0, 0., 10000);
    // functionFit->SetParLimits(1, 0., 10000);
  
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    h->SetStats(0);
    h->SetTitle("");
    h->SetXTitle("mJJ [GeV]");
    h->SetYTitle("mJJ Bias");
    h-> GetYaxis()->SetTitleOffset(1.4);
    h->GetXaxis()->SetRangeUser(400,2100);
    h->SetMaximum(0.1);
    h->SetMinimum(-0.1);
    h -> SetMarkerStyle(20);
    h -> SetMarkerSize(1.5);
    h -> Fit(functionFit, "R");
    c1->SetLeftMargin(0.13);
    h -> Draw();
    lineup->Draw("same");
    linedown->Draw("same");
    c1-> SaveAs( "Plot_mJJ/mJJBias_vs_Mjj.png" );
  
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
  
    for(size_t jj = 0; jj < mJJBinningSize ; jj++){ 
      std::pair<float, float> mJJBins = mMjjBinning.getBinValue(jj);	       
      int mJJMin = mJJBins.first;
      int mJJMax = mJJBins.second;
      double mJJMean= (mJJMin + mJJMax) /2. ;
    
      //    double mJJBiasFit = -Par0*(log(mJJMean)*log(mJJMean)) + Par1*log(mJJMean) + Par2;                                                                                                      
      double mJJBiasFit = Par0*(log(mJJMean)*log(mJJMean)*log(mJJMean)) + Par1*(log(mJJMean)*log(mJJMean)) + Par2*log(mJJMean) + Par3; 
      cout<<"mJJ Mean = " << mJJMean <<endl;
      cout<<"mJJ Bias from fit = " << mJJBiasFit <<endl;    
      Histo_Corrections -> SetBinContent(jj+1, mJJBiasFit);
    }
  
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);                                                                                       
    Histo_Corrections->SetStats(0);                                                                          
    Histo_Corrections->SetTitle("");    
    Histo_Corrections->SetXTitle("mJJ [GeV]");              
    Histo_Corrections->SetYTitle("Correction");                                                         
    Histo_Corrections->SetMaximum(0.1);
    Histo_Corrections->SetMinimum(-0.1);
    Histo_Corrections-> GetYaxis()->SetTitleOffset(1.4);                                                
    Histo_Corrections-> GetXaxis()->SetRangeUser(400, 3000);                                 
    Histo_Corrections-> SetMarkerStyle(20);                                                                          
    Histo_Corrections-> SetMarkerSize(1.5);                                                     
    c2->SetLeftMargin(0.13);                                                                                               
    Histo_Corrections-> Draw();                                                                                                           
    lineup->Draw("same");
    linedown->Draw("same");
    c2-> SaveAs( "Plot_mJJ/Histo_Corrections.png" );   

    TFile* outputFile= new TFile("mJJCorrections_JEC_HLT_v7.root","RECREATE");                                                                                                      
    outputFile->cd();                                                                                                                                                            
    Histo_Corrections -> Write();                                                                                                                                                
    outputFile->Close();   
  }

  if(ClosureTest){

    TH1D *h  = (TH1D*)dataFile_HT450->Get( "mJJBias_vs_Mjj" );    
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
    h->SetStats(0);
    h->SetTitle("");
    h->SetXTitle("mJJ [GeV]");
    h->SetYTitle("mJJ Bias");
    h-> GetYaxis()->SetTitleOffset(1.4);
    h->GetXaxis()->SetRangeUser(400,2100);
    h->SetMaximum(0.1);
    h->SetMinimum(-0.1);
    h -> SetMarkerStyle(20);
    h -> SetMarkerSize(1.5);
    c1->SetLeftMargin(0.13);
    h -> Draw();
    lineup->Draw("same");
    linedown->Draw("same");
    c1-> SaveAs( "Plot_mJJ/ClosureTest_mJJBias_vs_Mjj.png" );   

    TH1D *h2  = (TH1D*)dataFile_HT450->Get( "mJJBias_vs_Eta" );
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 800);
    h2->SetStats(0);
    h2->SetTitle("");
    h2->SetXTitle("#eta");
    h2->SetYTitle("mJJ Bias");
    h2-> GetYaxis()->SetTitleOffset(1.4);
    h2->GetXaxis()->SetRangeUser(-2.5, 2.5);
    h2->SetMaximum(0.1);
    h2->SetMinimum(-0.1);
    h2 -> SetMarkerStyle(20);
    h2 -> SetMarkerSize(1.5);
    c2->SetLeftMargin(0.13);
    h2 -> Draw();
    lineup->Draw("same");
    linedown->Draw("same");
    c2-> SaveAs( "Plot_mJJ/ClosureTest_mJJBias_vs_Eta.png" );   


  }

  
  return 0;
  
}




/*
  for(int ii = 1; ii < (etaBinningSize+1) ; ii++){ 
   
  std::string TGraphName = TString::Format("Test_AK4_pTBias_vs_pT_etabin_%i", ii).Data();
  map<std::string , TGraphErrors*>::iterator nh_h = userTGraphs_.find( TGraphName );
  if( nh_h == userTGraphs_.end()) continue;

  TF1* functionFit = new TF1( "functionFit", "-[0]*(log(x)*log(x)) + [1]*log(x)", 30,700); // quadratica log
  // functionFit->SetParLimits(0, 0., 10000);
  // functionFit->SetParLimits(1, 0., 10000);

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  nh_h->second->SetMarkerStyle(20);
  nh_h->second->SetMarkerSize(1.5);
  //    nh_h->second->Fit(functionFit, "R");
  nh_h->second->Draw("ZAP");
  c1->SaveAs( Form("Plot/%s.png", TGraphName.c_str() ) );
  }




  /////////////////////////////////////////////////////////////////////
  // creating TGraph
  // PtBinning mPtBinning;
  //  size_t pTBinningSize = mPtBinning.size();

  for(int ii = 1; ii < (etaBinningSize+1) ; ii++){ //# bins in eta
  cout<< "Bin in eta # "<< ii << endl;
  for(size_t jj=0; jj < pTBinningSize; jj++){

  std::pair<float, float> ptBins = mPtBinning.getBinValue(jj);	       
  std::string Histo = TString::Format("AK4_pTBias_vs_Eta_VariableBin_pT_%i_%i", (int) ptBins.first, (int) ptBins.second ).Data();
  cout<< Histo.c_str() << endl;

  int pTMin = ptBins.first;
  int pTMax = ptBins.second;
  double pTMean= (pTMin + pTMax) /2. ;

  TH1D *h;
  if( jj < 2){ // pT<300
  h = (TH1D*)dataFile_DoubleMu->Get( Histo.c_str() );
  cout << " Lo prendo dal DoubleMu" << endl;
  }else if(jj < 6){ // pT<700 due to low statistics
  h = (TH1D*)dataFile_HT450->Get( Histo.c_str() );
  cout << " Lo prendo dal HT450" << endl;
  }else{ continue; }
      
  double pTBias      = h->GetBinContent(ii);
  double pTBiasErr = h->GetBinError(ii);

  cout<< "pT mean = "<< pTMean <<endl;
  cout<< "pT Bias = "<< pTBias <<endl;
  cout<< "sigma(pT Bias) = "<< pTBiasErr <<endl;

  std::string TGraphName = TString::Format("AK4_pTBias_vs_pT_etabin_%i", ii).Data();
  CreateAndFillUserTGraph(TGraphName.c_str(), jj, pTMean, pTBias, 0, pTBiasErr);

  }// for sui pT
  }// for bin in eta
  /////////////////////////////// costruiti i TGraph da fittare

  for(int ii = 1; ii < (etaBinningSize+1) ; ii++){ 
   
  std::string TGraphName = TString::Format("AK4_pTBias_vs_pT_etabin_%i", ii).Data();
  map<std::string , TGraphErrors*>::iterator nh_h = userTGraphs_.find( TGraphName );
  if( nh_h == userTGraphs_.end()) continue;

  //    TF1* functionFit = new TF1( "functionFit", "[0] + [1]* TMath::Log(x)", 30,700); // linearLog
  TF1* functionFit = new TF1( "functionFit", "-[0]*(log(x)*log(x)) + [1]*log(x)", 30,700); // quadratica log
  //    TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)) + [1]*log(x)", 30,700); // quadratica log
  // TF1* functionFit = new TF1( "functionFit", "[0]*(log(x)*log(x)) + [1]*log(x) +[2]", 30,700); // quadratica log

  // functionFit->SetParLimits(0, 0., 10000);
  // functionFit->SetParLimits(1, 0., 10000);

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 800);
  nh_h->second->SetMarkerStyle(20);
  nh_h->second->SetMarkerSize(1.5);
  nh_h->second->Fit(functionFit, "R");
  nh_h->second->Draw("ZAP");
  c1->SaveAs( Form("Plot/%s.png", TGraphName.c_str() ) );
    
  double Par0       =  functionFit->GetParameter(0);
  double Par1       =  functionFit->GetParameter(1);
  //    double Par2       =  functionFit->GetParameter(2);
  cout<< "Par 0 = " << Par0 << endl;
  cout<< "Par 1 = " << Par1 << endl;
  //    cout<< "Par 2 = " << Par2 << endl;

  for(int jj=0; jj < binnum; jj++){
      
  int pTMin = pTMeanArray[jj];
  int pTMax = pTMeanArray[jj+1];
  double pTMean= (pTMin + pTMax) /2. ;
  cout<<"pT Mean = " << pTMean <<endl;

  //   double pTBiasFit = Par0 + Par1* TMath::Log(pTMean);
  //      double pTBiasFit = -Par0*(log(pTMean)*log(pTMean)) + Par1*log(pTMean) +Par2;
  double pTBiasFit = Par0*(log(pTMean)*log(pTMean)) + Par1*log(pTMean);
  cout<<"pT Bias from fit = " << pTBiasFit <<endl;

  Histo_Corrections->SetBinContent(ii, jj+1, pTBiasFit);

  } //pT bin
  } // eta bin


  TCanvas *c2 = new TCanvas("c2", "c2", 900, 800);
  Histo_Corrections->SetStats(0);
  gStyle->SetPalette(55);
  Histo_Corrections->SetContour(100);
  Histo_Corrections->Draw("colz");
  c2->SetRightMargin(0.13);
  c2->SaveAs("Plot/Corrections.png");


  TFile* outputFile= new TFile("Corrections_JEC_HLT_v7.root","RECREATE");
  outputFile->cd();
  Histo_Corrections -> Write();
  outputFile->Close();
*/





/*
  EtaBinning mEtaBinning;
  size_t etaBinningSize = mEtaBinning.size();
  PtBinning mPtBinning;
  size_t pTBinningSize = mPtBinning.size();

  for (size_t ii = 0; ii < etaBinningSize; ii++) {
  for (size_t jj = 0; jj < pTBinningSize; jj++) {

  std::string etaName = mEtaBinning.getBinName(ii);
  std::pair<float, float> ptBins = mPtBinning.getBinValue(jj);	       
      
  std::string HistoName = TString::Format("pTResolution_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
  std::string HistoName2 = TString::Format("pt_WideJet_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();

  std::string HistoName3 = TString::Format("AK4_pTResolution_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();
  std::string HistoName4 = TString::Format("AK4_pt_Jets_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();

  std::cout << HistoName.c_str()<< std::endl;

  // se facciamo un binning piu stretto e non esistono tutti gli istogramma possiamo fare cosi
  //      TH1D *h1 = (TH1D*)dataFile1->Get(HistoName.c_str());
  //      if( !h1) continue;

  }// pt binning
  }// eta binning

*/
