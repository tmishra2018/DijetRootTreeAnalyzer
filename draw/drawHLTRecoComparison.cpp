#include <stdlib.h>
#include <sys/stat.h>  

#include "TParameter.h"
#include "TPaveText.h"
#include "TError.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "drawBase.h"

//#include "etaBinning.h"
//#include "ptBinning.h"

#include <TColor.h>

int main(int argc, char* argv[]) {

  if (argc != 2 && argc != 3) {

    std::cout << "USAGE: ./drawPhotonJet [file.root]" << std::endl;
    exit(23);
  }

  std::string outputDir = "Plot";                 
  mkdir(outputDir.c_str(), 0777); 

  std::string data1_dataset(argv[1]);

  TString dataFileName1;
  dataFileName1 = TString::Format("%s", data1_dataset.c_str() );

  TFile* dataFile1 = TFile::Open(dataFileName1);

  if (dataFile1) {
    std::cout << "Opened data file '" << dataFileName1 << "'." << std::endl;
  }


  gErrorIgnoreLevel = kWarning;

  /*  
  vector<TString> HistoNameHLT;
  HistoNameHLT.push_back(TString::Format("cutHisto_noCuts_________________Dijet_MassAK4"));
 
  vector<TString> HistoNameReco;
  HistoNameReco.push_back(TString::Format("cutHisto_noCuts_________________Dijet_MassAK4reco"));

  vector<TString> XAxis;
  XAxis.push_back(TString::Format("M(jj) [GeV]"));
  
  size_t sizeHLT = HistoNameHLT.size();
  std::cout<< sizeHLT << std::endl;

  size_t sizeReco = HistoNameReco.size();
  std::cout<< sizeReco << std::endl;

  for(size_t ii = 0; ii < sizeHLT; ii++){
    std::cout<< HistoNameHLT.at(ii) << std::endl;
    std::cout<< HistoNameReco.at(ii) << std::endl;
    
    TH1F *h1 = (TH1F*)dataFile1->Get(HistoNameHLT.at(ii) );
    TH1F *h2 = (TH1F*)dataFile1->Get(HistoNameReco.at(ii) );
    double mean_h1 = h1->GetMean();
    double mean_h2 = h2->GetMean();

    std::cout<< mean_h1 << std::endl; 
    std::cout<< mean_h2 << std::endl;

    h1->Rebin(10);
    h2->Rebin(10);
    Normalizer(h1);
    Normalizer(h2);

    TLegend* legend = new TLegend(0.70, 0.80, 0.90, 0.90);
    legend->SetTextFont(42);
    legend->SetBorderSize(0);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.036);
    legend->AddEntry(h1, "HLT", "PL");
    legend->AddEntry(h2, "Reco", "PL");
   	 
    DrawPull("Plot/", HistoNameHLT.at(ii)+"_Pull.png", h1, h2, XAxis.at(ii) , "Events", legend);
    //    DrawRatio("Plot/", HistoName.at(ii)+"_Ratio.png", h1, h2, XAxis.at(ii) , "Events", legend);
  }
  */

  bool log = true;
  bool colz = true;

  DrawTH1(dataFile1, "deltaR_minimum",  "dR",  "Events", log,  0., 1.);
  //  DrawTH1(dataFile1, "pTResolution",  "pT resolution",  "Events", false,  -2., 2.);
  //  DrawTH1(dataFile1, "pTResolutionCut",  "pT resolution",  "Events", false,  -2., 2.);
  FitTH1(dataFile1, "pTResolution",  "pT resolution",  "Events", false,  -0.5, 0.5);
  DrawTH2(dataFile1, "ptWJ_Reco_vs_HLT",  "pT(Reco)",  "pT(HLT)", log, colz,  0., 2000., 0., 2000.);
  FitTH1(dataFile1, "pTResolution_Cut1",  "pT resolution",  "Events", false,  -0.5, 0.5);
  DrawTH2(dataFile1, "ptWJ_Reco_vs_HLT_Cut1",  "pT(Reco)",  "pT(HLT)", log, colz,  0., 2000., 0., 2000.);
  FitTH1(dataFile1, "pTResolution_Cut2",  "pT resolution",  "Events", false,  -0.5, 0.5);
  DrawTH2(dataFile1, "ptWJ_Reco_vs_HLT_Cut2",  "pT(Reco)",  "pT(HLT)", log, colz,  0., 2000., 0., 2000.);


  return 0;
  
}

