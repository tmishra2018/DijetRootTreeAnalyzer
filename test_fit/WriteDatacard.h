// preparing datacards
void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, Float_t width,int iChan, std::string model) {

  TString cardDir = "/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/datacardWithAllSyst/"+filePOSTfix;
  Int_t ncat = NCAT;
  TString wsDir   = "/afs/cern.ch/work/s/soffi/CMSSW611-Analysis/src/h2gglobe/ChiaraFitLimits/workspaces/"+filePOSTfix;

  // **********************
  // Retrieve the datasets
  cout << "Start retrieving dataset" << endl;
  
  RooDataSet* data[9];
  RooDataSet* signal[9];
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
  }

  RooRealVar*  lumi = w->var("lumi");
  /*
  // *****************************
  // Print Expected event yields
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data: " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("SigWeight")->sumEntries()  << endl;
  Float_t siglikeErr[6];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << signal[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*signal[c]->sumEntries();
  }
  cout << "====================================================" << endl;  

 */
  // *****************************
  // Printdata Data Card int file
  TString filename(cardDir+TString(fileBaseName)+"_"+"8TeV"+Form(("_m%.2f_w%.2f_channel%d_"+model+".txt").c_str(),w->var("MH")->getVal(),width,iChan));
  ofstream outFile(filename);

  outFile << "#CMS-HGG HighMass DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d datacardName.txt -U -m *mass* -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax *" << endl;
  outFile << "jmax *" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << "shapes data_obs * " << wsDir+TString(fileBkgName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()) << " w_all:data_obs_$CHANNEL" << endl;
  outFile << "shapes sig * "      << wsDir+TString(fileBaseName)+"_8TeV"+TString::Format(("_m%.2f_w%.2f.inputsig_"+model+".root").c_str(),w->var("MH")->getVal(),width) << " w_all:PhotonsMassSig_$CHANNEL" << endl;
  outFile << "shapes bkg * "      << wsDir+TString(fileBkgName)+TString::Format("_m%.2f.root",w->var("MH")->getVal()) << " w_all:PhotonsMassBkg_$CHANNEL" << endl;

  outFile << "---------------" << endl;
  outFile << Form("bin          cat%d", iChan) << endl;
  outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                 " << Form("cat%d      cat%d", iChan, iChan) << endl;
  outFile << "process                 sig      bkg" << endl;
  outFile << "process                   0        1" << endl;
  // if(signalScaler==1.)
  // signalScaler=1./signal[2]->sumEntries()*20;
  outFile << "rate                   " 
    //	  << signal[iChan]->sumEntries()*signalScaler << " " << data[iChan]->sumEntries() << endl;
    // << 1 << " " << data[iChan]->sumEntries() << endl;
	  <<19620. << "  "<<data[iChan]->sumEntries() << endl;
  outFile << "--------------------------------" << endl;
  outFile << "# signal scaled by " << signalScaler << endl;

 
  outFile << "lumi_8TeV     lnN     1.026000  - " << endl;
  outFile << "eff_trig     lnN     1.010000  - " << endl;
  outFile << "global_syst     lnN     1.050000  - " << endl;
  if(width==0.1) outFile << "bw_syst     lnN     1.00100  - " << endl;
  if(width==2) outFile << "bw_syst     lnN     1.0200  - " << endl;
  if(width==5) outFile << "bw_syst     lnN     1.0500  - " << endl;
  if(width==7) outFile << "bw_syst     lnN     1.0700  - " << endl;
  if(width==10) outFile << "bw_syst     lnN     1.1000  - " << endl;
  if(iChan==0){    
    outFile << "id_eff_eb     lnN     1.005000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.000000  - " << endl;    
    outFile << "r9Eff   lnN   1.0145/0.9915   - " << endl;
    outFile << "vtxEff   lnN   0.996/1.008   - " << endl; 
  }else if(iChan==1){    
    outFile << "id_eff_eb     lnN     1.005000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.000000  - " << endl;    
    outFile << "r9Eff   lnN   0.985/1.0085   - " << endl;
    outFile << "vtxEff   lnN   0.998/1.005   - " << endl; 
  }else if(iChan==2){    
    outFile << "id_eff_eb     lnN     1.00000 - " << endl;    
    outFile << "id_eff_ee     lnN     1.02600  - " << endl;    
    outFile << "r9Eff   lnN   1.0193/0.964   - " << endl;
    outFile << "vtxEff   lnN   0.996/1.007   - " << endl; 
  }else if(iChan==3){    
    outFile << "id_eff_eb     lnN     1.000000  - " << endl;    
    outFile << "id_eff_ee     lnN     1.026000  - " << endl;    
    outFile << "r9Eff   lnN   0.981/1.0357   - " << endl;
    outFile << "vtxEff   lnN   0.998/1.003   - " << endl; 
  }
    
    outFile << Form("mShift_cat%d    param   0 1 ",iChan) << endl;
    outFile << Form("mSmear_cat%d    param   0 1 ",iChan) << endl;



  // outFile << "CMS_VV_eff_g         lnN  0.8/1.20      - # Signal Efficiency" << endl;
  // outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  // outFile << Form("CMS_hgg_sig_m0_absShift    param   1   0.0125   # displacement of the mean w.r.t. nominal in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_sig_sigmaScale     param   1   0.1   # multiplicative correction to sigmas in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("rooHistFunc_cat%d_norm       ",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope2_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope3_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // if (iChan != 2 )  outFile << Form("CMS_hgg_bkg_8TeV_slope1_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
}


