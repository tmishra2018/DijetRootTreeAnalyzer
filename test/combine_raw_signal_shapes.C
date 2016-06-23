void combine_raw_signal_shapes()
{
    TFile *f_gg_500 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_500_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    //TFile *f_gg_750 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_750_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_1000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_2000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_2000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_3000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_3000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_4000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_4000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_5000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_5000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_6000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_6000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_7000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_7000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_8000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_8000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_gg_9000 = new TFile("output/RSGravToGluGlu/RSGravitonToGluonGluon_kMpl01_M_9000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_500 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_500_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_750 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_750_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_1000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_1000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    //TFile *f_qq_2000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_2000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_3000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_3000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    //TFile *f_qq_4000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_4000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_5000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_5000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_6000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_6000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_7000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_7000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_8000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_8000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qq_9000 = new TFile("output/RSGravToQQ/RSGravitonToQuarkQuark_kMpl01_M_9000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_500 = new TFile("output/QstarToJJ/QstarToJJ_M_500_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    //TFile *f_qg_750 = new TFile("output/QstarToJJ/QstarToJJ_M_750_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_1000 = new TFile("output/QstarToJJ/QstarToJJ_M_1000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_2000 = new TFile("output/QstarToJJ/QstarToJJ_M_2000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_3000 = new TFile("output/QstarToJJ/QstarToJJ_M_3000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_4000 = new TFile("output/QstarToJJ/QstarToJJ_M_4000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_5000 = new TFile("output/QstarToJJ/QstarToJJ_M_5000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_6000 = new TFile("output/QstarToJJ/QstarToJJ_M_6000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_7000 = new TFile("output/QstarToJJ/QstarToJJ_M_7000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_8000 = new TFile("output/QstarToJJ/QstarToJJ_M_8000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");
    TFile *f_qg_9000 = new TFile("output/QstarToJJ/QstarToJJ_M_9000_TuneCUETP8M1_13TeV_pythia8__RunIISpring16MiniAODv2-PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1__MINIAODSIM_JEC_Spring16_25nsV3.root");

    TH1F *h_RSGgg_Spring16_M500_WJ, *h_RSGgg_Spring16_M750_WJ,
        *h_RSGgg_Spring16_M1000_WJ, *h_RSGgg_Spring16_M2000_WJ,
	*h_RSGgg_Spring16_M3000_WJ, *h_RSGgg_Spring16_M4000_WJ,
	*h_RSGgg_Spring16_M5000_WJ, *h_RSGgg_Spring16_M6000_WJ,
	*h_RSGgg_Spring16_M7000_WJ, *h_RSGgg_Spring16_M8000_WJ,
	*h_RSGgg_Spring16_M9000_WJ;
    TH1F *h_RSGqq_Spring16_M500_WJ, *h_RSGqq_Spring16_M750_WJ,
        *h_RSGqq_Spring16_M1000_WJ, *h_RSGqq_Spring16_M2000_WJ,
	*h_RSGqq_Spring16_M3000_WJ, *h_RSGqq_Spring16_M4000_WJ,
	*h_RSGqq_Spring16_M5000_WJ, *h_RSGqq_Spring16_M6000_WJ,
	*h_RSGqq_Spring16_M7000_WJ, *h_RSGqq_Spring16_M8000_WJ,
	*h_RSGqq_Spring16_M9000_WJ;
    TH1F *h_QstarToJJ_Spring16_M500_WJ, *h_QstarToJJ_Spring16_M750_WJ,
        *h_QstarToJJ_Spring16_M1000_WJ, *h_QstarToJJ_Spring16_M2000_WJ,
	*h_QstarToJJ_Spring16_M3000_WJ, *h_QstarToJJ_Spring16_M4000_WJ,
	*h_QstarToJJ_Spring16_M5000_WJ, *h_QstarToJJ_Spring16_M6000_WJ,
	*h_QstarToJJ_Spring16_M7000_WJ, *h_QstarToJJ_Spring16_M8000_WJ,
	*h_QstarToJJ_Spring16_M9000_WJ;
    //TString histogram_name = "h_mjj_ratio";
    TString histogram_name = "h_mjj_ratio_nom";
    f_gg_500->GetObject(histogram_name, h_RSGgg_Spring16_M500_WJ);
    //f_gg_750->GetObject(histogram_name, h_RSGgg_Spring16_M750_WJ);
    f_gg_1000->GetObject(histogram_name, h_RSGgg_Spring16_M1000_WJ);
    f_gg_2000->GetObject(histogram_name, h_RSGgg_Spring16_M2000_WJ);
    f_gg_3000->GetObject(histogram_name, h_RSGgg_Spring16_M3000_WJ);
    f_gg_4000->GetObject(histogram_name, h_RSGgg_Spring16_M4000_WJ);
    f_gg_5000->GetObject(histogram_name, h_RSGgg_Spring16_M5000_WJ);
    f_gg_6000->GetObject(histogram_name, h_RSGgg_Spring16_M6000_WJ);
    f_gg_7000->GetObject(histogram_name, h_RSGgg_Spring16_M7000_WJ);
    f_gg_8000->GetObject(histogram_name, h_RSGgg_Spring16_M8000_WJ);
    f_gg_9000->GetObject(histogram_name, h_RSGgg_Spring16_M9000_WJ);
    f_qq_500->GetObject(histogram_name, h_RSGqq_Spring16_M500_WJ);
    f_qq_750->GetObject(histogram_name, h_RSGqq_Spring16_M750_WJ);
    f_qq_1000->GetObject(histogram_name, h_RSGqq_Spring16_M1000_WJ);
    //f_qq_2000->GetObject(histogram_name, h_RSGqq_Spring16_M2000_WJ);
    f_qq_3000->GetObject(histogram_name, h_RSGqq_Spring16_M3000_WJ);
    //f_qq_4000->GetObject(histogram_name, h_RSGqq_Spring16_M4000_WJ);
    f_qq_5000->GetObject(histogram_name, h_RSGqq_Spring16_M5000_WJ);
    f_qq_6000->GetObject(histogram_name, h_RSGqq_Spring16_M6000_WJ);
    f_qq_7000->GetObject(histogram_name, h_RSGqq_Spring16_M7000_WJ);
    f_qq_8000->GetObject(histogram_name, h_RSGqq_Spring16_M8000_WJ);
    f_qq_9000->GetObject(histogram_name, h_RSGqq_Spring16_M9000_WJ);
    f_qg_500->GetObject(histogram_name, h_QstarToJJ_Spring16_M500_WJ);
    //f_qg_750->GetObject(histogram_name, h_QstarToJJ_Spring16_M750_WJ);
    f_qg_1000->GetObject(histogram_name, h_QstarToJJ_Spring16_M1000_WJ);
    f_qg_2000->GetObject(histogram_name, h_QstarToJJ_Spring16_M2000_WJ);
    f_qg_3000->GetObject(histogram_name, h_QstarToJJ_Spring16_M3000_WJ);
    f_qg_4000->GetObject(histogram_name, h_QstarToJJ_Spring16_M4000_WJ);
    f_qg_5000->GetObject(histogram_name, h_QstarToJJ_Spring16_M5000_WJ);
    f_qg_6000->GetObject(histogram_name, h_QstarToJJ_Spring16_M6000_WJ);
    f_qg_7000->GetObject(histogram_name, h_QstarToJJ_Spring16_M7000_WJ);
    f_qg_8000->GetObject(histogram_name, h_QstarToJJ_Spring16_M8000_WJ);
    f_qg_9000->GetObject(histogram_name, h_QstarToJJ_Spring16_M9000_WJ);

    h_RSGgg_Spring16_M500_WJ->SetTitle("Mjj_WJ/M");
    //h_RSGgg_Spring16_M750_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M1000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M2000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M3000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M4000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M5000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M6000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M7000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M8000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGgg_Spring16_M9000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M500_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M750_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M1000_WJ->SetTitle("Mjj_WJ/M");
    //h_RSGqq_Spring16_M2000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M3000_WJ->SetTitle("Mjj_WJ/M");
    //h_RSGqq_Spring16_M4000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M5000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M6000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M7000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M8000_WJ->SetTitle("Mjj_WJ/M");
    h_RSGqq_Spring16_M9000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M500_WJ->SetTitle("Mjj_WJ/M");
    //h_QstarToJJ_Spring16_M750_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M1000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M2000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M3000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M4000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M5000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M6000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M7000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M8000_WJ->SetTitle("Mjj_WJ/M");
    h_QstarToJJ_Spring16_M9000_WJ->SetTitle("Mjj_WJ/M");

    //TFile *f_gg_out = new TFile("InputShapes_RSGgg_Spring16.root", "RECREATE");
    TFile *f_gg_out = new TFile("InputShapes_RSGgg_CaloScouting_Spring16.root", "RECREATE");
    f_gg_out->cd();
    h_RSGgg_Spring16_M500_WJ->Write("h_RSGgg_Spring16_M500_WJ");
    //h_RSGgg_Spring16_M750_WJ->Write("h_RSGgg_Spring16_M750_WJ");
    h_RSGgg_Spring16_M1000_WJ->Write("h_RSGgg_Spring16_M1000_WJ");
    h_RSGgg_Spring16_M2000_WJ->Write("h_RSGgg_Spring16_M2000_WJ");
    h_RSGgg_Spring16_M3000_WJ->Write("h_RSGgg_Spring16_M3000_WJ");
    h_RSGgg_Spring16_M4000_WJ->Write("h_RSGgg_Spring16_M4000_WJ");
    h_RSGgg_Spring16_M5000_WJ->Write("h_RSGgg_Spring16_M5000_WJ");
    h_RSGgg_Spring16_M6000_WJ->Write("h_RSGgg_Spring16_M6000_WJ");
    h_RSGgg_Spring16_M7000_WJ->Write("h_RSGgg_Spring16_M7000_WJ");
    h_RSGgg_Spring16_M8000_WJ->Write("h_RSGgg_Spring16_M8000_WJ");
    h_RSGgg_Spring16_M9000_WJ->Write("h_RSGgg_Spring16_M9000_WJ");
    f_gg_out->Close();

    //TFile *f_qq_out = new TFile("InputShapes_RSGqq_Spring16.root", "RECREATE");
    TFile *f_qq_out = new TFile("InputShapes_RSGqq_CaloScouting_Spring16.root", "RECREATE");
    f_qq_out->cd();
    h_RSGqq_Spring16_M500_WJ->Write("h_RSGqq_Spring16_M500_WJ");
    h_RSGqq_Spring16_M750_WJ->Write("h_RSGqq_Spring16_M750_WJ");
    h_RSGqq_Spring16_M1000_WJ->Write("h_RSGqq_Spring16_M1000_WJ");
    //h_RSGqq_Spring16_M2000_WJ->Write("h_RSGqq_Spring16_M2000_WJ");
    h_RSGqq_Spring16_M3000_WJ->Write("h_RSGqq_Spring16_M3000_WJ");
    //h_RSGqq_Spring16_M4000_WJ->Write("h_RSGqq_Spring16_M4000_WJ");
    h_RSGqq_Spring16_M5000_WJ->Write("h_RSGqq_Spring16_M5000_WJ");
    h_RSGqq_Spring16_M6000_WJ->Write("h_RSGqq_Spring16_M6000_WJ");
    h_RSGqq_Spring16_M7000_WJ->Write("h_RSGqq_Spring16_M7000_WJ");
    h_RSGqq_Spring16_M8000_WJ->Write("h_RSGqq_Spring16_M8000_WJ");
    h_RSGqq_Spring16_M9000_WJ->Write("h_RSGqq_Spring16_M9000_WJ");
    f_qq_out->Close();

    //TFile *f_qg_out = new TFile("InputShapes_QstarToJJ_Spring16.root", "RECREATE");
    TFile *f_qg_out = new TFile("InputShapes_QstarToJJ_CaloScouting_Spring16.root", "RECREATE");
    f_qg_out->cd();
    h_QstarToJJ_Spring16_M500_WJ->Write("h_QstarToJJ_Spring16_M500_WJ");
    //h_QstarToJJ_Spring16_M750_WJ->Write("h_QstarToJJ_Spring16_M750_WJ");
    h_QstarToJJ_Spring16_M1000_WJ->Write("h_QstarToJJ_Spring16_M1000_WJ");
    h_QstarToJJ_Spring16_M2000_WJ->Write("h_QstarToJJ_Spring16_M2000_WJ");
    h_QstarToJJ_Spring16_M3000_WJ->Write("h_QstarToJJ_Spring16_M3000_WJ");
    h_QstarToJJ_Spring16_M4000_WJ->Write("h_QstarToJJ_Spring16_M4000_WJ");
    h_QstarToJJ_Spring16_M5000_WJ->Write("h_QstarToJJ_Spring16_M5000_WJ");
    h_QstarToJJ_Spring16_M6000_WJ->Write("h_QstarToJJ_Spring16_M6000_WJ");
    h_QstarToJJ_Spring16_M7000_WJ->Write("h_QstarToJJ_Spring16_M7000_WJ");
    h_QstarToJJ_Spring16_M8000_WJ->Write("h_QstarToJJ_Spring16_M8000_WJ");
    h_QstarToJJ_Spring16_M9000_WJ->Write("h_QstarToJJ_Spring16_M9000_WJ");
    f_qg_out->Close();

    return;
}
