#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string *inputList, string *cutFile,
                             string *treeName, string *outputFileName,
                             string *cutEfficFile)
    : baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
    std::cout << "analysisClass::analysisClass(): begins " << std::endl;

    std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
    std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

    std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
    std::cout << "analysisClass::Loop() begins" << std::endl;

    if (fChain == 0)
        return;

    //////////book histos here

    // variable binning for mjj plots
    const int nMassBins = 103;

    double massBoundaries[nMassBins+1]
        = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176,
           197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606,
           649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246,
           1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132,
           2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416,
           3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253,
           5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866,
           8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179,
           11571, 11977, 12395, 12827, 13272, 13732, 14000};

    TH1F* h_mjj_fullSel_varBin = new TH1F("h_mjj_fullSel_varBin", "",
                                          103, massBoundaries);
    TH1F* h_mjj_fullSel_fixBin = new TH1F("h_mjj_fullSel_fixBin", "",
                                          getHistoNBins("mjj"),
                                          getHistoMin("mjj"), getHistoMax("mjj"));

    TH1F* h_mjj_ratio = new TH1F("h_mjj_ratio", "", getHistoNBins("mjj_ratio"),
                                 getHistoMin("mjj_ratio"), getHistoMax("mjj_ratio"));

    /////////initialize variables

    Long64_t nentries = fChain->GetEntriesFast();
    std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;

    ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
    ////// If the root version is updated and rootNtupleClass regenerated,     /////
    ////// these lines may need to be updated.                                 /////
    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
            break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if(jentry < 10 || jentry%100000 == 0)
            std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;

        ////////////////////// User's code starts here ///////////////////////

        ///Stuff to be done for every event

        resetCuts();

        double mjj_ratio = mjj/getPreCutValue1("resonanceMass");

        //== Fill Variables ==
        fillVariableWithValue("deltaETAjj", deltaETAjj);
        fillVariableWithValue("mjj", mjj);
        fillVariableWithValue("mjj_ratio", mjj_ratio);

        // Evaluate cuts (but do not apply them)
        evaluateCuts();

        if (passedCut("all")) {
            h_mjj_fullSel_varBin->Fill(getVariableValue("mjj"));
            h_mjj_fullSel_fixBin->Fill(getVariableValue("mjj"));
            h_mjj_ratio->Fill(getVariableValue("mjj_ratio"));
        }

        ////////////////////// User's code ends here ///////////////////////

    } // End loop over events

    //////////write histos

    h_mjj_fullSel_varBin->Write();
    h_mjj_fullSel_fixBin->Write();
    h_mjj_ratio->Write();

    std::cout << "analysisClass::Loop() ends" << std::endl;
}
