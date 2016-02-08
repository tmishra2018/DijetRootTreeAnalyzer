#include <stdlib.h>
#include <sys/stat.h>  
#include <time.h>

#include "TParameter.h"
#include "TPaveText.h"
#include "TError.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TAxis.h"
#include "drawBase.h"
#include "ptBinning.h"
#include "etaBinning.h"

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

  bool colz = true;
  bool log= true;

  DrawTH1(dataFile1, "deltaR_minimum",  "#Delta R",  "Events", log,  0, 1);

  DrawComparison(dataFile1, "NVertex", "n vertex", "Events", false, false, 0, 35);
  DrawComparison(dataFile1, "met", "MET [GeV]", "Events", log, false, 0, 1500);
  DrawComparison(dataFile1, "metSig", "MET Significance", "Events", log, false, 0, 1);
  DrawComparison(dataFile1, "deltaETAjj", "#Delta #eta (jj)", "Events", log, false, 0, 3);
  DrawComparison(dataFile1, "AK4_deltaETAjj", "#Delta #eta (jj)", "Events", log, false, 0, 3);

  DrawComparison(dataFile1, "chargedHadEnFrac", "CHF", "Events", log, false, 0, 1);
  DrawComparison(dataFile1, "neutrHadEnFrac", "NHF", "Events", log, false, 0, 1);
  DrawComparison(dataFile1, "chargedElectromFrac", "CEMF", "Events", log, false, 0, 1);
  DrawComparison(dataFile1, "neutrElectromFrac", "NEMF", "Events", log, false, 0, 1);
  DrawComparison(dataFile1, "muEnFrac", "MUF", "Events", log, false, 0, 1);
  DrawComparison(dataFile1, "chargedMult", "CM", "Events", log, false, 0, 100);
  DrawComparison(dataFile1, "neutrMult", "NM", "Events", log, false, 0, 100);
  DrawComparison(dataFile1, "photonMult", "PM", "Events", log, false, 0, 100);
  DrawComparison(dataFile1, "jetCSVAK4", "CSV", "Events", log, false, 0, 1);

  DrawComparison(dataFile1, "pt_WideJet", "pT [GeV]", "Events", false, false, 0, 1500);
  DrawComparison(dataFile1, "eta_WideJet", "#eta", "Events", false, false, -2.5, 2.5);
  DrawComparison(dataFile1, "phi_WideJet", "#phi", "Events", false, false, -3.5, 3.5);
  DrawComparison(dataFile1, "MJJ_WideJet", "m(jj) [GeV]", "Events", false, false, 0, 3000);

  DrawTH2(dataFile1, "NVertex",  "nVtx(Reco)",  "nVtx(HLT)", log, colz,  0, 35, 0, 35);
  DrawTH2(dataFile1, "met",  "MET(Reco) [GeV]",  "MET(HLT) [GeV]", log, colz,  0., 1500., 0., 1500.);
  DrawTH2(dataFile1, "metSig",  "MET Sig (Reco)",  "MET Sig (HLT)", log, colz,  0., 1., 0., 1.);
  DrawTH2(dataFile1, "deltaETAjj",  "#Delta #eta (jj) (Reco)",  "#Delta #eta (jj) (HLT)", log, colz,  0., 3., 0., 3.);
  DrawTH2(dataFile1, "AK4_deltaETAjj",  "#Delta #eta (jj) (Reco)",  "#Delta #eta (jj) (HLT)", log, colz,  0., 3., 0., 3.);

  DrawTH2(dataFile1, "chargedHadEnFrac",  "CHF Reco",  "CHF HLT", log, colz,  0., 1., 0., 1.);
  DrawTH2(dataFile1, "neutrHadEnFrac",  "NHF Reco",  "NHF HLT", log, colz,  0., 1., 0., 1.);
  DrawTH2(dataFile1, "chargedElectromFrac",  "CEMF Reco",  "CEMF HLT", log, colz,  0., 1., 0., 1.);
  DrawTH2(dataFile1, "neutrElectromFrac",  "NEMF Reco",  "NEMF HLT", log, colz,  0., 1., 0., 1.);
  DrawTH2(dataFile1, "muEnFrac",  "MUF Reco",  "MUF HLT", log, colz,  0., 1., 0., 1.);
  DrawTH2(dataFile1, "chargedMult",  "CM Reco",  "CM HLT", log, colz,  0., 100., 0., 100.);
  DrawTH2(dataFile1, "neutrMult",  "NM Reco",  "NM HLT", log, colz,  0., 100., 0., 100.);
  DrawTH2(dataFile1, "photonMult",  "PM Reco",  "PM HLT", log, colz,  0., 100., 0., 100.);
  DrawTH2(dataFile1, "jetCSVAK4",  "CSV Reco",  "CSV HLT", log, colz,  0., 1., 0., 1.);

  DrawTH2(dataFile1, "pt_WideJet",  "pT(Reco) [GeV]",  "pT(HLT) [GeV]", log, colz,  0., 1500., 0., 1500.);
  DrawTH2(dataFile1, "eta_WideJet",  "#eta (Reco)",  "#eta (HLT)", log, colz,  -2.5, 2.5, -2.5, 2.5);
  DrawTH2(dataFile1, "phi_WideJet",  "#phi (Reco)",  "#phi (HLT)", log, colz,  -3.5, 3.5, -3.5, 3.5);
  DrawTH2(dataFile1, "MJJ_WideJet",  "mjj(Reco) [GeV]",  "mjj(HLT) [GeV]", log, colz,  0, 3000, 0, 3000);

  DrawTH1(dataFile1, "pTResolution",  "pT bias",  "Events", false,  -0.5, 0.5);
  FitTH1(dataFile1, "pTResolution",  "pT bias",  "Events", false,  -0.5, 0.5);

  DrawTH2(dataFile1, "AK4_pt_Jets",  "pT(Reco) [GeV]",  "pT(HLT) [GeV]", log, colz,  0., 1500., 0., 1500.);
  DrawTH1(dataFile1, "AK4_pTResolution",  "pT bias",  "Events", false,  -0.5, 0.5);
  FitTH1(dataFile1, "AK4_pTResolution",  "pT bias",  "Events", false,  -0.5, 0.5);

  // vs N Vertex
  DrawCompareProfile(dataFile1, "AK4_pT_vs_Nvtx",  "nVtx",  "pT [GeV]", false, false, 0, 35);
  DrawTProfile(dataFile1, "AK4_pTResolution_vs_Nvtx",  "nVtx",  "pT bias", false, 0, 35);
  DrawCompareProfile(dataFile1, "pT_vs_Nvtx",  "nVtx",  "pT [GeV]", false, false,  0, 35);
  DrawTProfile(dataFile1, "pTResolution_vs_Nvtx",  "nVtx",  "pT bias", false,  0, 35);

  DrawCompareProfile(dataFile1, "pT_vs_Nvtx_Barrel",  "nVtx",  "pT [GeV]", false, false,  0, 35);
  DrawTProfile(dataFile1, "pTResolution_vs_Nvtx_Barrel",  "nVtx",  "pT bias", false,  0, 35);
  FitTH1(dataFile1, "pTResolution_Barrel",  "pT bias",  "Events", false,  -0.5, 0.5);

  DrawCompareProfile(dataFile1, "pT_vs_Nvtx_Endcap",  "nVtx",  "pT [GeV]", false, false,  0, 35);
  DrawTProfile(dataFile1, "pTResolution_vs_Nvtx_Endcap",  "nVtx",  "pT bias", false,  0, 35);
  FitTH1(dataFile1, "pTResolution_Endcap",  "pT bias",  "Events", false,  -0.5, 0.5);

  //////////////////// drawing in in eta and pT bins
  EtaBinning mEtaBinning;
  PtBinning mPtBinning;
  
  size_t etaBinningSize = mEtaBinning.size();
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

      //anche piu stretto va bene
      DrawTH2(dataFile1, HistoName2.c_str(),  "pT(Reco) [GeV]",  "pT(HLT) [GeV]", log, colz,  0., 1500., 0., 1500.);
      DrawTH1(dataFile1, HistoName.c_str(),  "pT bias",  "Events", false,  -0.5, 0.5);
      FitTH1(dataFile1, HistoName.c_str(),  "pT bias",  "Events", false,  -0.5, 0.5);

      DrawTH2(dataFile1, HistoName4.c_str(),  "pT(Reco) [GeV]",  "pT(HLT) [GeV]", log, colz,  0., 1500., 0., 1500.);
      DrawTH1(dataFile1, HistoName3.c_str(),  "pT bias",  "Events", false,  -0.5, 0.5);
      FitTH1(dataFile1, HistoName3.c_str(),  "pT bias",  "Events", false,  -0.5, 0.5);

    }
  }

  /*

  EtaBinning mEtaBinning;
  PtmakBinning mPtBinning;
  
  size_t etaBinningSize = mEtaBinning.size();
  size_t pTBinningSize = mPtBinning.size();

  double nEntriesTot=0;

  for (size_t ii = 0; ii < etaBinningSize; ii++) {
    for (size_t jj = 0; jj < pTBinningSize; jj++) {

      std::string etaName = mEtaBinning.getBinName(ii);
      std::pair<float, float> ptBins = mPtBinning.getBinValue(jj);	       
      
      std::string HistoName = TString::Format("pTResolution_%s_pT_%i_%i", etaName.c_str(), (int) ptBins.first, (int) ptBins.second ).Data();

      std::cout << HistoName.c_str()<< std::endl;

      TH1D *h1 = (TH1D*)dataFile1->Get(HistoName.c_str());
   
      if( !h1) continue;

      double nEntries = h1->GetEntries();
      std::cout<< "nEntries = "<< nEntries << std::endl;

      nEntriesTot +=nEntries;

    }
  }

      std::cout<< "nEntries TOTALI = "<< nEntriesTot << std::endl;

  */



  return 0;
  
}

