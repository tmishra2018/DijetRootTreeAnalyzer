{

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TFile *_file0 = TFile::Open("triggerEfficiency_L1HTT150seed_HT450_DetaJJLess1p3_HLTv7_output.root");
  _file0->ls();
  TEfficiency *h_efficiency = (TEfficiency*)_file0->Get("h_mjj_HLTpass_L1HTT150_1GeVbin_clone");

  h_efficiency->Draw();

  //-- sigmoid function --
  //TF1* f1 = new TF1("f1","1/( 1 + exp( -[0]*(x-[1]) ) )",386,1118);
  //TF1* f1 = new TF1("f1","1/( 1 + exp( -[0]*(x-[1]) ) )",220,2037);
  //TF1* f1 = new TF1("f1","1/( 1 + exp( -[0]*(x-[1]) ) )",606,2037);
  //f1->SetParameters(0.02,500);

  gPad->Update();

  TGraphAsymmErrors *g_eff = (TGraphAsymmErrors*) h_efficiency->GetPaintedGraph();

  // //-- err function --
  // TF1* f1 = new TF1("f1","([0]/2)* ( 1 + TMath::Erf((x-[1])/[2]))",453,890);
  // f1->SetParameters(1,500,95);
  // //f1->SetParLimits(0,0.99,1);//unconstrained
  // f1->SetParLimits(0,0.999999,1);//constrained
  // f1->SetParLimits(1,450,550);
  // f1->SetParLimits(2,80,120);

  // //-----------------
  //g_eff->Fit(f1,"VLRI");
  //g_eff->Draw();

}
