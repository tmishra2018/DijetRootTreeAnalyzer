

// Write background pdfs and datasets into the workspace 
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName, double mass) {

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  
  //RooWorkspace *wAll = w;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all"); 
     
  // Retrieve the datasets and PDFs 
  //RooDataSet* data[NCAT]; 
  //RooDataHist* dataBinned[NCAT]; 
  RooDataSet* data;
  RooDataHist* dataBinned;
  RooArgSet* argset;

  data = (RooDataSet*) w->data(TString::Format("Data_cat0"));
  dataBinned = (RooDataHist*) w->data(TString::Format("data_obs_cat0"));
  

  //// Create binning object with range (-10,10)
  //const int nMassBins = 103;

  //double massBoundaries[nMassBins+1] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000};
  //RooBinning abins(nMassBins,massBoundaries,"mjj_binning") ;
  //data = (RooDataSet*) w->data(TString::Format("Data_cat0"));
  //RooRealVar* mjj = (RooRealVar*) data->get()->find("mjj");

  //mjj->setBinning(abins) ; 
  //RooDataHist dataBinned(TString::Format("data_obs_cat0"),"",RooArgSet(*mjj),*data) ;

  //for (int c=0; c<ncat; ++c) { 
  //  cout << "load dataset" << endl;
  //  
  //  data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c)); 
  //  if(!data[c]) cout << "WARNING: unbinned dataset not found" << endl ;
  //  data[c]->Print();
  //  
  //  cout << "create binned dataset" << endl;
  //  
  //  argset = (RooArgSet*)data[c]->get();
  //  cout << "debug " << endl;
  //  //dataBinned[c]->Print("v") ;
      
    // //giulia: gives "bad alloc error.... 
    //dataBinned[c]= (RooDataHist*)data[c]->binnedClone(); 
    //if(!dataBinned[c]) cout << "WARNING: binned dataset not ncreated" << endl ;
    //else cout << "assigned dataBinned" << endl;
    //dataBinned[c]->Print();

    ////giulia
    ////wAll->import(*w->pdf(TString::Format("DijetMassBkg_cat%d",c))); 
    //cout << "import binned dataset" << endl;
    //wAll->import(*dataBinned[c]); 
  //  cout << "import unbinned dataset" << endl;
  //  wAll->import(*w->data(TString::Format("Data_cat%d",c)),RooFit::Rename(TString::Format("data_unbinned_obs_cat%d",c))); 
  //  cout << "done with WriteWorkspace" << endl;
  // }
  wAll->import(*data); 
  wAll->import(*dataBinned); 
  std::cout << "test: done with importing background dataset only" << std::endl;     
  //std::cout << "done with importing background pdfs" << std::endl;


  TString filename;
  filename = (wsDir+TString(fileBaseName)+TString::Format("_m%.2f.root",mass));


  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;
  
//  std::cout << std::endl; 
//  std::cout << "observation:" << std::endl;
//  for (int c=0; c<ncat; ++c) {
//    std::cout << "cat " << c << ", " << wAll->data(TString::Format("data_unbinned_obs_cat%d",c))->sumEntries() << endl;
//    wAll->data(TString::Format("data_unbinned_obs_cat%d",c))->Print();
//  }
//  std::cout << std::endl;
//  
  return;
}

////////////////////////////////////////////////////////////////////////
//NOT YET ADAPTED
/////////////////////////////////////////////////////////////////////////
/*
// Write signal pdfs and datasets into the workspace 
void MakeSigWS(RooWorkspace* w, const char* fileBaseName, Float_t width, std::string model){
  
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");  


  // Retrieve P.D.F.s
   //w->Print("V");
  for (int c=0; c<ncat; ++c) {
    std::cout<<"flag"<<std::endl;
      wAll->import(*w->pdf("PhotonsMassSig"+TString::Format("_cat%d",c)));//*w->pdf("PhotonsMassSigCBCExt"+TString::Format("_cat%d",c))
     
    
      wAll->import(*w->data(TString::Format("SigWeight_cat%d",c)));
      wAll->import(*w->function("PhotonsMassSig"+TString::Format("_cat%d_norm",c)));
                                                 
  }
  std::cout << "done with importing signal pdfs" << std::endl;
  wAll->import(*w->var("massReduced"));
  wAll->import(*w->var("PhotonsMassTrue"));
  // (2) Systematics on energy scale and resolution // chiara: per ora tutte le sistematiche non hanno senso
  // wAll->factory("CMS_hgg_sig_m0_absShift[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat0(massggnewvtx_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat1(massggnewvtx_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");

  // (3) Systematics on resolution: create new sigmas
  // wAll->factory("CMS_hgg_sig_sigmaScale[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat0(massggnewvtx_sig_sigma0_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(massggnewvtx_sig_sigma1_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat1(massggnewvtx_sig_sigma0_cat1, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(massggnewvtx_sig_sigma1_cat1, CMS_hgg_sig_sigmaScale)")

  TString filename(wsDir+TString(fileBaseName)+TString::Format(("_m%.2f_w%.2f.inputsig_"+model+".root").c_str(),w->var("MH")->getVal(),width));
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  
  return;
}


*/
