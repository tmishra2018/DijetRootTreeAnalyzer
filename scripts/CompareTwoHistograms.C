// %%%%%%%%%%%%%%%%%%  EXAMPLE OF USAGE %%%%%%%%%%%%%%%%%%%
// ##RSqq
// # root -b -l << EOF
// # .L CompareTwoHistograms.C
// # CompareTwoHistograms("../rootfile_RSGravitonToQuarkQuark_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8__RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2__MINIAODSIM.root","../rootfile_RSGravitonToQuarkQuark_kMpl01_M_1000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM.root","cutHisto_allCuts________________Diparton_M","resonance mass (gen level) [GeV]","RSqq1000.png","RSqq1000Log.png","RSqq - Spring15","RSqq - Phys14")
// ....
// ....
// # CompareTwoHistograms("../rootfile_RSGravitonToQuarkQuark_kMpl01_M_9000_TuneCUETP8M1_13TeV_pythia8__RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v1__MINIAODSIM.root","../rootfile_RSGravitonToQuarkQuark_kMpl01_M_9000_Tune4C_13TeV_pythia8__Phys14DR-PU20bx25_PHYS14_25_V1-v1__MINIAODSIM.root","cutHisto_allCuts________________Diparton_M","resonance mass (gen level) [GeV]","RSqq9000.png","RSqq9000Log.png","RSqq - Spring15","RSqq - Phys14")
// # EOF
// %%%%%%%%%%%%%%%%%  END EXAMPLE OF USAGE %%%%%%%%%%%%%%%%%%%

void CompareTwoHistograms(const char* inputFile1, const char* inputFile2, const char* histoname, const char* XaxisName, const char *plotName, const char *plotNameLog, const char* histo1Label, const char* histo2Label)
{
  gROOT->Reset();

  gStyle->SetOptStat(0);

  TFile *f1 = TFile::Open(inputFile1);
  TFile *f2 = TFile::Open(inputFile2);

  //TFile *f1 = TFile::Open("../rootFile_Qstar7000_Spring15.root");
  //TFile *f2 = TFile::Open("../rootFile_Qstar7000_Phys14.root");
    
  //char histoname[200];
  //sprintf(histoname,"cutHisto_allCuts________________Diparton_M");

  //char plotName[200];
  //sprintf(plotName,"Qstar7000.png");
  //char plotNameLog[200];
  //sprintf(plotNameLog,"Qstar7000Log.png");

  //char XaxisName[200];
  //sprintf(XaxisName,"resonance mass (gen level) [GeV]");

  //char histo1Label[200];
  //sprintf(histo1Label,"Qstar - Spring15");
  //sprintf(histo1Label,&pippo);
  //char histo2Label[200];
  //sprintf(histo2Label,"Qstar - Phys14");

  //==============================================================

  TCanvas c1;
  //  c1.Divide(1,2);
  //c1.SetLogy();

  //  c1.cd(1);

  f1->cd();
  TH1F *histo1 = (TH1F*)f1->Get(histoname);
  histo1->SetLineColor(1);
  histo1->GetXaxis()->SetRangeUser(histo1->GetMean()-5*histo1->GetRMS(),histo1->GetMean()+5*histo1->GetRMS());
  histo1->GetXaxis()->SetTitle(XaxisName);
  histo1->DrawNormalized();

  f2->cd();
  TH1F *histo2 = (TH1F*)f2->Get(histoname);
  histo2->SetLineColor(2);
  histo2->DrawNormalized("sames");

  TLegend *legend=new TLegend(0.,0.85,0.35,0.99);
  legend->SetTextFont(72);
  legend->SetTextSize(0.035);
  legend->AddEntry(histo1,histo1Label,"l");
  legend->AddEntry(histo2,histo2Label,"l");
  legend->Draw();

  c1.SaveAs(plotName);

  c1.SetLogy();

  c1.SaveAs(plotNameLog);

  // c1.cd(2);
  // TH1F *histoRatio = new TH1F("1 / 2","histoRatio",histo1->GetNbinsX(),histo1->GetXaxis()->GetXmin(),histo1->GetXaxis()->GetXmax());
  // histoRatio->Sumw2();
  // histoRatio->Divide(histo1,histo2,1,1);
  // histoRatio->GetYaxis()->SetRangeUser(0.,2.0);
  // histoRatio->Draw();

}
