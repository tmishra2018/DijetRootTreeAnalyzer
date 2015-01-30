// Default fit
Double_t fitQCD1( Double_t *m, Double_t *p)
{
  double x=m[0]/13000.;
  return p[0]*pow(1.-x,p[1])/pow(x,p[2]+p[3]*log(x));
}


int makePlot_coolStyle(string inputfile){

  // Set TDR Style	
  gROOT->ProcessLine(".L tdrStyle.C");
  gROOT->ProcessLine("setTDRStyle()");

  TFile* input = new TFile(inputfile.c_str(),"read");
  TH1F* h_QCD_mc = (TH1F*)input->Get("hist_allCutsQCD");
  TH1F* h_QCD_1GeV = new TH1F("h_QCD_1GeV","",10000,1.,10000.);
  h_QCD_1GeV->Sumw2();
  h_QCD_1GeV->FillRandom(h_QCD_mc,h_QCD_mc->Integral());
 
  
  const int nMassBins = 94;
  
  double massBoundaries[nMassBins+1] = {0, 1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325,	 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10000};
  
  TH1F* h_QCD_rebin = (TH1F*)h_QCD_1GeV->Rebin(nMassBins,"h_QCD_rebin",massBoundaries);
  cout << "n bins = " << h_QCD_rebin->GetNbinsX() << endl; 
  h_QCD_rebin->Draw();
 

  TH1F* h_QCD = (TH1F*)h_QCD_rebin->Clone("h_QCD");
  for (int i=1; i <=h_QCD_rebin->GetNbinsX(); i++){
    h_QCD->SetBinContent(i,h_QCD_rebin->GetBinContent(i)/h_QCD_rebin->GetBinWidth(i) );
    h_QCD->SetBinError(i,h_QCD_rebin->GetBinError(i)/h_QCD_rebin->GetBinWidth(i));
  }
  //QCD Fit -- fit to qcd

  TF1 *f_qcd = new TF1("fit_qcd",fitQCD1, 0.,10000.0,4); 
  gStyle->SetOptFit(0000); 
  f_qcd->SetParameter(0,0.0001456);
  f_qcd->SetParameter(1,3.413);
  f_qcd->SetParameter(2,7.772);
  f_qcd->SetParameter(3,0.4801);
  //h_QCD->Fit("fit_qcd", "R");  
  
  TH1F* hPulls = new TH1F("hPulls", "", h_QCD->GetNbinsX(), h_QCD->GetBinLowEdge(1), h_QCD->GetBinLowEdge(h_QCD->GetNbinsX()));
  double pull = 0;
  for(int i=1; i <= h_QCD->GetNbinsX(); i++){
    pull = (h_QCD->GetBinContent(i)-f_qcd->Eval(h_QCD->GetBinCenter(i)))/h_QCD->GetBinError(i);
    hPulls->SetBinContent(i,pull);

  }

  
  TCanvas* c = new TCanvas("c","DijetMass cross section with Fit and QCD MC",600,650);
  c->GetWindowHeight();
  c->GetWindowWidth();
  c->SetLogy();
  c->Divide(1,2,0,0,0);
  c->cd(1);
  p11_1 = (TPad*)c->GetPad(1);
  p11_1->SetPad(0.01,0.23,0.99,0.98);
  p11_1->SetLogy();
  p11_1->SetRightMargin(0.05);
  p11_1->SetTopMargin(0.05);

  // Pave text
  TPaveText *pave_fit = new TPaveText(0.18,0.15,0.40,0.27,"NDC");
  pave_fit->AddText(" #sqrt{s} = 8 TeV");
  pave_fit->AddText("|#eta| < 2.5, |#Delta#eta| < 1.3");
  pave_fit->AddText("M_{jj} > 890 GeV");
  pave_fit->AddText("Wide Jets");
  pave_fit->SetFillColor(0);
  pave_fit->SetLineColor(0);
  pave_fit->SetFillStyle(0);
  pave_fit->SetBorderSize(0);
  pave_fit->SetTextFont(42);
  pave_fit->SetTextSize(0.03);
  pave_fit->SetTextAlign(12); 


  TPaveText *pt1 = new TPaveText(0.1284756,0.9602144,0.3887139,0.9902251,"brNDC");
  pt1->SetBorderSize(0);
  pt1->SetFillColor(0);
  pt1->SetFillStyle(0);
  pt1->SetLineColor(0);
  pt1->SetTextAlign(12);
  pt1->SetTextSize(0.035);
  TText *text = pt1->AddText("CMS");

  TPaveText *pt2 = new TPaveText(0.45,0.96,0.65,0.99,"brNDC");
  pt2->SetBorderSize(0);
  pt2->SetFillColor(0);
  pt2->SetFillStyle(0);
  pt2->SetLineColor(0);
  pt2->SetTextAlign(12);
  pt2->SetTextSize(0.035);
  TText *text2 = pt2->AddText("#sqrt{s} = 13 TeV");

  TPaveText *pt3 = new TPaveText(0.7687988,0.9602144,0.9297357,0.9902251,"brNDC");
  pt3->SetBorderSize(0);
  pt3->SetFillColor(0);
  pt3->SetFillStyle(0);
  pt3->SetLineColor(0);
  pt3->SetTextAlign(12);
  pt3->SetTextSize(0.035);
  TText *text3 = pt3->AddText("L= 1 fb^{-1}");

  TH1F *vFrame = p11_1->DrawFrame(890.0,0.000000001,5500.0,10.0);
  vFrame->SetTitle("");
  vFrame->SetXTitle("Dijet Mass (GeV)");
  vFrame->GetXaxis()->SetTitleSize(0.06);
  vFrame->SetYTitle("events / bin width");
  //h_QCD->GetYaxis()->SetRangeUser(0.000000002,10);
  h_QCD->GetXaxis()->SetRangeUser(1000,8000.0);
  h_QCD->SetTitle("");
  h_QCD->SetLineColor(1);
  h_QCD->SetFillColor(1);
  h_QCD->SetMarkerColor(1);
  h_QCD->SetMarkerStyle(20);
  h_QCD->Draw("AP");
  f_qcd->SetLineWidth(2);
  f_qcd->SetLineStyle(2);
  f_qcd->SetLineColor(2);
  f_qcd->Draw("same");

  TLegend *leg = new TLegend(0.60,0.638,0.934,0.915);
  leg->SetTextSize(0.03146853);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(0);
  leg->SetMargin(0.35);
  leg->AddEntry(h_QCD,"pseudo data" ,"PL");
  leg->AddEntry(f_qcd,"fit to pseudodata","L");
  leg->Draw("same");

  pt1->Draw("same"); 
  pt2->Draw("same");
  pt3->Draw("same");


//redraw axis
  p11_1->RedrawAxis();
  p11_1->Update();
  cout << "MIN: " << p11_1->GetUxmin() << endl;
  cout << "MAX: " << p11_1->GetUxmax() << endl;

  //---- Next PAD

  c->cd(2);
  p11_2 = (TPad*)c->GetPad(2);
  p11_2->SetPad(0.01,0.02,0.99,0.24);
  p11_2->SetBottomMargin(0.35);
  p11_2->SetRightMargin(0.05);
  p11_2->SetGridx();
  p11_2->SetGridy();
  //c3_2->SetTickx(50);

  TH1F *vFrame2 = p11_2->DrawFrame(p11_1->GetUxmin(), -4.5, p11_1->GetUxmax(), 4.5);

  vFrame2->SetTitle("");
  vFrame2->SetXTitle("Dijet Mass (GeV)");
  vFrame2->GetXaxis()->SetTitleSize(0.06);
  vFrame2->SetYTitle("(pseudoData-Fit)/#sigma_{pseudoData}");
  vFrame2->GetYaxis()->SetTitleSize(0.12);
  vFrame2->GetYaxis()->SetLabelSize(0.07);
  vFrame2->GetYaxis()->SetTitleOffset(0.50);
  vFrame2->GetXaxis()->SetTitleOffset(0.90);
  vFrame2->GetXaxis()->SetTitleSize(0.18);
  vFrame2->GetXaxis()->SetLabelSize(0.1);

  hPulls->GetXaxis()->SetRangeUser(1000.,8000.);
  hPulls->GetYaxis()->SetRangeUser(-4.,4.);
  hPulls->SetLineWidth(0);
  hPulls->SetFillColor(2);
  hPulls->SetLineColor(1);
  hPulls->Draw("SAMEHIST");

  TLine *line = new TLine(890.,0,5500,0);
  line->Draw("");

  return 0;

}


