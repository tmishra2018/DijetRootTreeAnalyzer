#!/usr/bin/env python
import sys, os, copy, re
from array import array
from argparse import ArgumentParser
import math
from ROOT import *
import CMS_lumi, setTDRStyle

# Global variables
readExistingHisto = 1
inputHistoFileName = "rawhistV7_Run2015D_scoutingPFHT_UNBLINDED_649_838_JEC_HLTplusV7_Mjj_cor_smooth.root"
inputEfficiencyFileName = "triggerEfficiency_L1HTT150seed_HT450_DetaJJLess1p3_HLTv7Corr_output.root"
inputDataFileName = "/t3/users/santanas/Dijet13TeVScouting/rootTrees_reduced/ScoutingPFHT__15_01_2016_20160115_192039/merged/rootfile_ScoutingPFHT__Run2015D-v1__RAW_ScoutingPFHT__15_01_2016_20160115_192039_reduced_skim.root"
treeName = "rootTupleTree/tree"
outputLabel = "output"
var = "mjj"
sqrtS = 13000.0
lumiValue = 1866.0 #[pb]
showCrossSection = 1 #1=cross section [pb] , 0=number of events/GeV
fixedRange = 1 #1=YES , 0=NO  (the option works only if showCrossSection=1; otherwise=0)
minY = 0.00000003
maxY = 20
if showCrossSection==1:
    lumi = lumiValue
else:
    lumi = 1
massMin = 565
massMax = 2037
#massMax = 6328
efficiency_massMin = 453
efficiency_massMax = 890
blindRegionMassMin = 0
#blindRegionMassMin = 649
blindRegionMassMax = 0
#blindRegionMassMax = 838
doBlind = False
if blindRegionMassMin != blindRegionMassMax:
    doBlind = True
xaxisTitle = "Dijet Mass [GeV]"
doSimultaneousFit = True

input_sig_gg_file_name = "ResonanceShapes_gg_13TeV_Scouting_Spring15.root"
input_sig_RSGgg750 = "h_gg_750"
xsec_RSGgg750 = 15.0 # pb

if showCrossSection==1:
    yaxisTitle_main = "d#sigma / dm_{jj}   [pb / GeV]"
else:
    yaxisTitle_main = "Number of events / GeV"
yaxisTitle_efficiency = "Trigger Efficiency"
yaxisTitle_secondary = "#frac{(Data-Fit)}{#sigma_{Data}}   "
range_residual = 3.5
eff_range_residual = 3.5
MinNumEvents = 10.
nParFit = 4
massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
################### DO NOT MODIFY ######################
massBins_list_actual = []
firstBin = -1
lastBin = -1
index_bin_boundary=-1
for bin_boundary in massBins_list:
    index_bin_boundary = index_bin_boundary+1
    if (bin_boundary>=massMin and firstBin==-1 ):
        massMin = bin_boundary
        firstBin=1
        print "FIRST BIN is "+str(massMin)
    if (bin_boundary>(massMax+0.0000000001) and lastBin==-1 ):
        massMax = massBins_list[index_bin_boundary-1]
        lastBin=1
        print "LAST BIN is "+str(massMax)
    if (bin_boundary>=massMin and bin_boundary<=massMax):
        massBins_list_actual.append(bin_boundary)
print massBins_list_actual
massBins = array("d",massBins_list_actual)
N_massBins = len(massBins)-1
eff_massBins_list_actual = []
eff_firstBin = -1
eff_lastBin = -1
index_bin_boundary=-1
for bin_boundary in massBins_list:
    index_bin_boundary = index_bin_boundary+1
    if (bin_boundary>=efficiency_massMin and eff_firstBin==-1):
        efficiency_massMin = bin_boundary
        eff_firstBin=1
        print "FIRST EFFICIENCY BIN is "+str(efficiency_massMin)
    if (bin_boundary>(efficiency_massMax+0.0000000001) and eff_lastBin==-1 ):
        efficiency_massMax = massBins_list[index_bin_boundary-1]
        eff_lastBin=1
        print "LAST EFFICIENCY BIN is "+str(efficiency_massMax)
    if (bin_boundary>=efficiency_massMin and bin_boundary<=efficiency_massMax):
        eff_massBins_list_actual.append(bin_boundary)
print eff_massBins_list_actual
eff_massBins = array("d", eff_massBins_list_actual)
N_eff_massBins = len(eff_massBins)-1
#############################################
#sel = "mjj>"+str(massMin)+" && (mjj<693 || mjj>838) && fabs(deltaETAjj)<1.3 && PassJSON==1 && (passHLT_HT450_BtagSeq==1 || passHLT_HT450==1)"
#sel = "((mjj>1181 && mjj<2037) || (mjj>4686 && mjj<6328)) && fabs(deltaETAjj)<1.3 && PassJSON==1 && (passHLT_HT450_BtagSeq==1 || passHLT_HT450==1)"
sel = "((mjj>"+str(massMin)+" && mjj<"+str(blindRegionMassMin)+") || (mjj>"+str(blindRegionMassMax)+" && mjj<"+str(massMax)+")) && fabs(deltaETAjj)<1.3 && PassJSON==1 && (passHLT_HT450_BtagSeq==1 || passHLT_HT450==1)"

#==================
# CMS style and lumi
#==================

#set tdr style
setTDRStyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "%.1f fb^{-1}" % (lumiValue/1000.0)
CMS_lumi.lumi_8TeV = "%.1f fb^{-1}" % (lumiValue/1000.0)
CMS_lumi.lumi_13TeV = "%.1f fb^{-1}" % (lumiValue/1000.0)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.cmsTextSize = 1.0

iPos = 11
if iPos==0:
    CMS_lumi.relPosX = 0.12
iPeriod = 4

W = 600
H = 700
H_ref = 700
W_ref = 600
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

#######################################################################

def main():

    # import ROOT stuff
    #from ROOT import TFile, TH1F, TH1D, TGraph, kTRUE, kFALSE, TF1, TCanvas
    #from ROOT import RooRealVar, RooDataHist, RooArgList, RooArgSet, RooAddPdf, RooFit, RooGenericPdf, RooWorkspace, RooMsgService, RooHistPdf, RooExtendPdf

    # input data
    if readExistingHisto==1:
        inputDataFile = TFile(inputHistoFileName)
    else:
        inputDataFile = TFile(inputDataFileName)
        treeData = inputDataFile.Get(treeName)
        print "number of entries in the data tree: ", treeData.GetEntries()
    inputEfficiencyFile = TFile(inputEfficiencyFileName)
    input_sig_RSGgg750_file = TFile(input_sig_gg_file_name)

    # mass variable
    mjj = RooRealVar('mjj', 'mjj', min(float(massMin), float(efficiency_massMin)),
                     max(float(massMax), float(efficiency_massMax)))
    mjj.setRange("MassFit", massMin, massMax)
    mjj.setRange("EfficiencyFit", efficiency_massMin, efficiency_massMax)
    mjj.setRange("RangeLow", massMin, blindRegionMassMin)
    mjj.setRange("RangeHigh", blindRegionMassMax, massMax)

    # data
    dataInt = 1
    if (readExistingHisto==1):
        h_data = (inputDataFile.Get("mjj_mjjcor_gev")).Clone("h_data")
        if doBlind:
            dataInt = h_data.Integral(massMin + 1, blindRegionMassMin) \
                      + h_data.Integral(blindRegionMassMax + 1, massMax)
        else:
            dataInt = h_data.Integral(massMin + 1, massMax)
    else:
        h_data = TH1F("h_data", "", 13000, 0.0, 13000.0);
        treeData.Project("h_data", var, sel);
        dataInt = h_data.GetEntries()
    print "number of events in the fit range:", int(dataInt)
    efficiency = inputEfficiencyFile.Get("efficiency")

    h_sig_RSGgg750 = input_sig_RSGgg750_file.Get(input_sig_RSGgg750)
    integral_sig_RSGgg750 = h_sig_RSGgg750.Integral()
    h_sig_RSGgg750.Scale(xsec_RSGgg750/integral_sig_RSGgg750)
    h_sig_RSGgg750.SetTitle("G_{RS} 750 GeV 15 pb")

    h_data_roo = RooDataHist('h_data_roo', 'h_data_roo', RooArgList(mjj), h_data)
    h_data_roo.Print()

    # Fill trigger efficiency dataset
    h_passed = efficiency.GetPassedHistogram()
    h_total = efficiency.GetTotalHistogram()
    cut = RooCategory("cut", "cutr")
    cut.defineType("pass", 1)
    cut.defineType("fail", 0)
    trigger_data = RooDataSet("trigger_data", "trigger data", RooArgSet(mjj, cut))
    for i in range(efficiency_massMin+1, efficiency_massMax+1):
        mjj.setVal(h_total.GetXaxis().GetBinCenter(i))
        if h_passed.GetBinContent(i) > 0:
            cut.setLabel("pass")
            weight = h_passed.GetBinContent(i)
            for j in range(int(weight)):
                trigger_data.add(RooArgSet(mjj, cut))
        if h_total.GetBinContent(i) > h_passed.GetBinContent(i):
            cut.setLabel("fail")
            weight = h_total.GetBinContent(i) - h_passed.GetBinContent(i)
            for j in range(int(weight)):
                trigger_data.add(RooArgSet(mjj, cut))
    trigger_data.Print("v")

    # Trigger efficiency model
    m_eff = RooRealVar('m_eff', 'm_eff', 494.404950525, 450.0, 550.0)
    sigma_eff = RooRealVar('sigma_eff', 'sigma_eff', 96.2722550919, 80.0, 120.0)
    effFunc = RooFormulaVar("effFunc",
                            "0.5*(1.0 + TMath::Erf((mjj - m_eff)/sigma_eff))",
                            RooArgList(mjj, m_eff, sigma_eff))
    effPdf = RooEfficiency("effPdf", "effPdf", effFunc, cut, "pass")
    if not doSimultaneousFit:
        effPdf.fitTo(trigger_data)
        effPdf.Print()
        m_eff.setConstant(kTRUE)
        sigma_eff.setConstant(kTRUE)

    # Background model
    norm = RooRealVar('norm', 'norm', dataInt, 0.0, 1.0e9)
    p1 = RooRealVar('p1', 'p1', 5.57031294347, 0.0, 100.0)
    p2 = RooRealVar('p2', 'p2', 7.21504246896, 0.0, 100.0)
    p3 = RooRealVar('p3', 'p3', 0.370141801135, -10.0, 10.0)
    background = RooGenericPdf('background',
                               '(pow(1-@0/%.1f,@1)/pow(@0/%.1f,@2+@3*log(@0/%.1f)))*(0.5*(1.0 + TMath::Erf((@0 - @4)/@5)))'%(sqrtS,sqrtS,sqrtS),
                               RooArgList(mjj, p1, p2, p3, m_eff, sigma_eff))
    background.Print()
    background_ext = RooExtendPdf("background_ext", "", background, norm)

    # fit
    num_regions = 1.0
    if doBlind:
        num_regions = 2.0
    if doBlind:
        nll_mass = background_ext.createNLL(h_data_roo,
                                            RooFit.Range("RangeLow,RangeHigh"),
                                            RooFit.Extended(kTRUE),
                                            RooFit.Save(kTRUE),
                                            RooFit.Strategy(1))
    else:
        nll_mass = background_ext.createNLL(h_data_roo, RooFit.Range("MassFit"),
                                            RooFit.Extended(kTRUE),
                                            RooFit.Save(kTRUE),
                                            RooFit.Strategy(1))

    if doSimultaneousFit:
        nll_efficiency = effPdf.createNLL(trigger_data,
                                          RooFit.Range("EfficiencyFit"),
                                          RooFit.Save(kTRUE), RooFit.Strategy(1))
        nll_simultaneous = RooAddition("nll_simultaneous", "simultaneous nll",
                                       RooArgList(nll_mass, nll_efficiency))
        res_b = RooMinuit(nll_simultaneous)
        res_b.migrad()
        res_b.simplex()
    else:
        res_b = RooMinuit(nll_mass)
        res_b.migrad()
        res_b.simplex()
    res_b.Print()

    # fit results
    norm_b = norm
    p1_b = p1
    p2_b = p2
    p3_b = p3

    background_noNorm = TF1("background_noNorm",
                            "(TMath::Power(1-x/%.1f,[0]))/(TMath::Power(x/%.1f,[1]+[2]*log(x/%.1f)))*(0.5*(1.0 + TMath::Erf((x - [3])/[4])))"%(sqrtS,sqrtS,sqrtS),
                            float(massMin), float(massMax))
    background_noNorm.SetParameter(0, p1_b.getVal())
    background_noNorm.SetParameter(1, p2_b.getVal())
    background_noNorm.SetParameter(2, p3_b.getVal())
    background_noNorm.SetParameter(3, m_eff.getVal())
    background_noNorm.SetParameter(4, sigma_eff.getVal())
    int_b = 1.0
    if doBlind:
        int_b = background_noNorm.Integral(float(massMin), float(blindRegionMassMin)) \
            + background_noNorm.Integral(float(blindRegionMassMax), float(massMax))
    else:
        int_b = background_noNorm.Integral(float(massMin), float(massMax))
    print "rescale for background function:", float(int_b)
    p0_b = norm_b.getVal() / (int_b*lumi) * num_regions
    print "p0_b = {0} + {1} - {2}".format(p0_b,
                                          norm_b.getErrorHi()/int_b*math.sqrt(num_regions),
                                          norm_b.getErrorLo()/int_b*math.sqrt(num_regions))
    print "p1_b = {0} + {1} - {2}".format(p1_b.getVal(), p1_b.getErrorHi(),
                                          p1_b.getErrorLo())
    print "p2_b = {0} + {1} - {2}".format(p2_b.getVal(), p2_b.getErrorHi(),
                                          p2_b.getErrorLo())
    print "p3_b = {0} + {1} - {2}".format(p3_b.getVal(), p3_b.getErrorHi(),
                                          p3_b.getErrorLo())
    print "m_eff = {0} + {1} - {2}".format(m_eff.getVal(), m_eff.getErrorHi(),
                                           m_eff.getErrorLo())
    print "sigma_eff = {0} + {1} - {2}".format(sigma_eff.getVal(),
                                               sigma_eff.getErrorHi(),
                                               sigma_eff.getErrorLo())

    background = TF1("background",
                     "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)))*(0.5*(1.0 + TMath::Erf((x - [4])/[5])))"%(sqrtS,sqrtS,sqrtS),
                     float(massMin), float(massMax))
    background.SetParameter(0, p0_b)
    background.SetParameter(1, p1_b.getVal())
    background.SetParameter(2, p2_b.getVal())
    background.SetParameter(3, p3_b.getVal())
    background.SetParameter(4, m_eff.getVal())
    background.SetParameter(5, sigma_eff.getVal())

    eff_fit = TF1("eff_fit", "0.5*(1.0 + TMath::Erf((x - [0])/[1]))",
                  float(efficiency_massMin), float(efficiency_massMax))
    eff_fit.SetParameter(0, m_eff.getVal())
    eff_fit.SetParameter(1, sigma_eff.getVal())

    f_passed = TF1("f_passed",
                   "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)))*(0.5*(1.0 + TMath::Erf((x - [4])/[5])))"%(sqrtS,sqrtS,sqrtS),
                   float(efficiency_massMin), float(efficiency_massMax))
    f_passed.SetParameter(0, p0_b)
    f_passed.SetParameter(1, p1_b.getVal())
    f_passed.SetParameter(2, p2_b.getVal())
    f_passed.SetParameter(3, p3_b.getVal())
    f_passed.SetParameter(4, m_eff.getVal())
    f_passed.SetParameter(5, sigma_eff.getVal())

    f_total = TF1("f_total",
                  "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)))"%(sqrtS,sqrtS,sqrtS),
                  float(efficiency_massMin), float(efficiency_massMax))
    f_total.SetParameter(0, p0_b)
    f_total.SetParameter(1, p1_b.getVal())
    f_total.SetParameter(2, p2_b.getVal())
    f_total.SetParameter(3, p3_b.getVal())

    # data graph
    h_data.Rebin(N_massBins, "h_data_varBin", massBins)
    g_data = TGraphAsymmErrors(h_data_varBin)

    alpha = 1-0.6827
    for i in range(0, g_data.GetN()):
        N = g_data.GetY()[i]
        binWidth = g_data.GetEXlow()[i] + g_data.GetEXhigh()[i]
        #print str(g_data.GetX()[i]-g_data.GetEXlow()[i])+" "+str(g_data.GetX()[i]+g_data.GetEXhigh()[i])+" "+str(N) 

        L = 0
        if N != 0:
            L = ROOT.Math.gamma_quantile(alpha/2, N, 1.0)
        U = ROOT.Math.gamma_quantile_c(alpha/2, N+1, 1)

        g_data.SetPointEYlow(i, (N-L)/(binWidth*lumi));
        g_data.SetPointEYhigh(i, (U-N)/(binWidth*lumi));
        g_data.SetPoint(i, g_data.GetX()[i], N/(binWidth*lumi))

        if (g_data.GetX()[i]>blindRegionMassMin
            and g_data.GetX()[i]<blindRegionMassMax):
            g_data.SetPointEYlow(i, 0);
            g_data.SetPointEYhigh(i, 0);
            g_data.SetPoint(i, g_data.GetX()[i], 0)

    # output
    output = TFile("{0}.root".format(outputLabel), "RECREATE");
    h_data.Write()

    list_parameter = [p0_b, p1_b.getVal(), p2_b.getVal(), p3_b.getVal(),
                      m_eff.getVal(), sigma_eff.getVal(),
                      (norm_b.getErrorHi() - norm_b.getErrorLo())/(2.0*int_b)*math.sqrt(num_regions),
                      (p1_b.getErrorHi() - p1_b.getErrorLo())/2.0,
                      (p2_b.getErrorHi() - p2_b.getErrorLo())/2.0,
                      (p3_b.getErrorHi() - p3_b.getErrorLo())/2.0,
                      (m_eff.getErrorHi() - m_eff.getErrorLo())/2.0,
                      (sigma_eff.getErrorHi() - sigma_eff.getErrorLo())/2.0]

    mass_1GeV_list = []
    for bin in range(efficiency_massMin, efficiency_massMax+1):
        mass_1GeV_list.append(bin)
    eff_1GeV_massBins = array("d", mass_1GeV_list)
    # plot
    h_background = convertFunctionToHisto(background,"h_background",N_massBins,massBins)
    h_efficiency = convertEffFunctionToHisto(eff_fit, "h_efficiency",
                                             len(eff_1GeV_massBins)-1, eff_1GeV_massBins)
    g_efficiency_binned = convertEfficiencyToGraph(h_total, h_passed, "g_efficiency_binned",
                                                   N_eff_massBins, eff_massBins)
    g_efficiency = convertEfficiencyToGraph(h_total, h_passed, "g_efficiency",
                                            len(eff_1GeV_massBins)-1, eff_1GeV_massBins)
    h_fit_residual_vs_mass = TH1D("h_fit_residual_vs_mass","h_fit_residual_vs_mass",N_massBins,massBins)
    h_efficiency_residual_vs_mass = TH1D("h_efficiency_residual_vs_mass",
                                         "h_efficiency_residual_vs_mass",
                                         efficiency_massMax - efficiency_massMin,
                                         efficiency_massMin, efficiency_massMax)
    list_chi2AndNdf_background \
        = calculateChi2AndFillResiduals(g_data, h_background,
                                        h_fit_residual_vs_mass, 6, 0)
    list_chi2AndNdf_efficiency \
        = calculateChi2AndFillResiduals(g_efficiency, h_efficiency,
                                        h_efficiency_residual_vs_mass, 2, 0)
    h_efficiency_residual_vs_mass_binned = TH1D("h_efficiency_residual_vs_mass_binned",
                                                "h_efficiency_residual_vs_mass_binned",
                                                N_eff_massBins, eff_massBins)
    h_efficiency_binned = TH1D("h_efficiency_binned", "h_efficiency_binned",
                               N_eff_massBins, eff_massBins)
    list_chi2AndNdf_efficiency_binned \
        = bin_efficiency(g_efficiency, h_efficiency_residual_vs_mass, eff_fit,
                         h_efficiency_residual_vs_mass_binned,
                         h_efficiency_binned, g_efficiency_binned, f_passed,
                         f_total)
    list_likelihood = []
    if doSimultaneousFit:
        list_likelihood = [nll_simultaneous.getVal(), nll_mass.getVal(),
                           nll_efficiency.getVal()]
    else:
        nll_efficiency = effPdf.createNLL(trigger_data,
                                          RooFit.Range("EfficiencyFit"),
                                          RooFit.Save(kTRUE), RooFit.Strategy(1))
        list_likelihood = [nll_mass.getVal() + nll_efficiency.getVal(),
                           nll_mass.getVal(), nll_efficiency.getVal()]
    signal_hist_list = [h_sig_RSGgg750]
    drawAndSavePlot_background(g_data, h_background, h_fit_residual_vs_mass,
                               outputLabel, list_chi2AndNdf_background,
                               list_parameter, g_efficiency, eff_fit,
                               h_efficiency_residual_vs_mass,
                               list_chi2AndNdf_efficiency_binned,
                               h_efficiency_residual_vs_mass_binned,
                               g_efficiency_binned, list_likelihood,
                               signal_hist_list)

    output.Close()

    return


#==============================================================================


def convertFunctionToHisto(background_,name_,N_massBins_,massBins_):

    background_hist_ = TH1D(name_,name_,N_massBins_,massBins_)

    for bin in range (0,N_massBins_):
        xbinLow = massBins_[bin]
        xbinHigh = massBins_[bin+1]
        binWidth_current = xbinHigh - xbinLow
        value = background_.Integral(xbinLow , xbinHigh) / binWidth_current
        # print "{0}: {1} {2} {3} {4}".format(massBins_[bin], xbinLow, xbinHigh,
        #                                     background_.Integral(xbinLow, xbinHigh),
        #                                     value)
        background_hist_.SetBinContent(bin+1,value)

    return background_hist_

def convertEffFunctionToHisto(efficiency_, name_, N_massBins_, massBins_):

    efficiency_hist_ = TH1D(name_, name_, N_massBins_, massBins_)

    for bin in range (0,N_massBins_):
        xbinLow = massBins_[bin]
        xbinHigh = massBins_[bin+1]
        xbinCenter = (xbinLow + xbinHigh)/2.0
        value = efficiency_.Eval(xbinCenter)
        efficiency_hist_.SetBinContent(bin+1,value)

    return efficiency_hist_

def convertEfficiencyToGraph(h_total_, h_passed_, name_, N_massBins_, massBins_):
    h_total_rebinned = h_total_.Rebin(N_massBins_, "h_total_rebinned", massBins_)
    h_passed_rebinned = h_passed_.Rebin(N_massBins_, "h_passed_rebinned", massBins_)
    efficiency_rebinned = TEfficiency(h_passed_rebinned, h_total_rebinned)
    return efficiency_rebinned.CreateGraph()


def calculateChi2AndFillResiduals(data_obs_TGraph_,background_hist_,hist_fit_residual_vsMass_, params_, prinToScreen_=0):
    print "-- "+str(background_hist_.GetName())

    N_massBins_ = data_obs_TGraph_.GetN()

    chi2_FullRangeAll = 0
    chi2_PlotRangeAll = 0
    chi2_PlotRangeNonZero = 0
    chi2_PlotRangeMinNumEvents = 0

    N_FullRangeAll = 0
    N_PlotRangeAll = 0
    N_PlotRangeNonZero = 0
    N_PlotRangeMinNumEvents = 0

    if(prinToScreen_):
        print ""
        print ""
        print "======== Number of events / GeV (data, errors, fit, residuals) ========"
        print ""

    for bin in range (0, N_massBins_):

        ## Values and errors

        value_data = data_obs_TGraph_.GetY()[bin]
        err_low_data = data_obs_TGraph_.GetEYlow()[bin]
        err_high_data = data_obs_TGraph_.GetEYhigh()[bin]
        xbinCenter = data_obs_TGraph_.GetX()[bin]
        xbinLow = data_obs_TGraph_.GetX()[bin]-data_obs_TGraph_.GetEXlow()[bin]
        xbinHigh = data_obs_TGraph_.GetX()[bin]+data_obs_TGraph_.GetEXhigh()[bin]
        binWidth_current = xbinHigh - xbinLow
        value_fit = background_hist_.GetBinContent(bin+1)

        ## Fit residuals

        err_tot_data = 0
        if (value_fit >= value_data):
            err_tot_data = err_high_data
        else:
            err_tot_data = err_low_data
        if (xbinCenter<blindRegionMassMin or xbinCenter>blindRegionMassMax):
            fit_residual = (value_data - value_fit) / err_tot_data
            err_fit_residual = 1
        else:
            fit_residual = 0
            err_fit_residual = 0

        ## Fill histo with residuals

        hist_fit_residual_vsMass_.SetBinContent(bin+1,fit_residual)
        hist_fit_residual_vsMass_.SetBinError(bin+1,err_fit_residual)

        ## Chi2

        chi2_FullRangeAll += pow(fit_residual,2)
        N_FullRangeAll += 1
        if (xbinLow >= massMin and xbinHigh<=massMax):
            chi2_PlotRangeAll += pow(fit_residual,2)
            N_PlotRangeAll += 1
            if (value_data > 0):
                chi2_PlotRangeNonZero += pow(fit_residual,2)
                N_PlotRangeNonZero += 1
                if(value_data * binWidth_current * lumi > MinNumEvents):
                    chi2_PlotRangeMinNumEvents += pow(fit_residual,2)
                    N_PlotRangeMinNumEvents += 1

        if(prinToScreen_):
            print str(xbinLow)+" "+str(xbinHigh)+" "+str(binWidth_current)+" : "+str(value_data)+" "+str(value_data * binWidth_current * lumi)+" - "+str(err_low_data)+" + "+str(err_high_data)+" fit: "+str(value_fit)+" fit residual: "+str(fit_residual) 

    #==================
    # Calculate chi2/ndf
    #==================

    # ndf
    ndf_FullRangeAll = N_FullRangeAll - params_
    ndf_PlotRangeAll = N_PlotRangeAll - params_
    ndf_PlotRangeNonZero = N_PlotRangeNonZero - params_
    ndf_PlotRangeMinNumEvents = N_PlotRangeMinNumEvents - params_

    chi2_ndf_FullRangeAll = chi2_FullRangeAll / ndf_FullRangeAll
    chi2_ndf_PlotRangeAll = chi2_PlotRangeAll / ndf_PlotRangeAll
    chi2_ndf_PlotRangeNonZero = chi2_PlotRangeNonZero / ndf_PlotRangeNonZero
    chi2_ndf_PlotRangeMinNumEvents = chi2_PlotRangeMinNumEvents / ndf_PlotRangeMinNumEvents

    print "chi2/ndf FullRangeAll : %.8f / %d = %.2f" % ( chi2_FullRangeAll , ndf_FullRangeAll , chi2_ndf_FullRangeAll ) 
    print "chi2/ndf PlotRangeAll : %.1f / %d = %.2f" % ( chi2_PlotRangeAll , ndf_PlotRangeAll , chi2_ndf_PlotRangeAll ) 
    print "chi2/ndf PlotRangeNonZero : %.1f / %d = %.2f" % ( chi2_PlotRangeNonZero , ndf_PlotRangeNonZero , chi2_ndf_PlotRangeNonZero ) 
    print "chi2/ndf PlotRangeMinNumEvents : %.1f / %d = %.2f" % ( chi2_PlotRangeMinNumEvents , ndf_PlotRangeMinNumEvents , chi2_ndf_PlotRangeMinNumEvents ) 

    return [chi2_FullRangeAll, ndf_FullRangeAll, chi2_PlotRangeAll, ndf_PlotRangeAll, chi2_PlotRangeNonZero, ndf_PlotRangeNonZero, chi2_PlotRangeMinNumEvents, ndf_PlotRangeMinNumEvents]


def bin_efficiency(g_efficiency_, h_efficiency_residual_vs_mass_, eff_fit_,
                   h_efficiency_residual_vs_mass_binned_, h_efficiency_binned_,
                   g_efficiency_binned_, f_passed_, f_total_):
    chi2 = 0
    ndf = 0

    bin = 0
    for i in range(0, N_eff_massBins):
        sum_value = 0.0
        sum_residual = 0.0
        sum_error2 = 0.0
        for j in range(int(eff_massBins[i]), int(eff_massBins[i+1])):
            bin += 1
            sum_value += g_efficiency_.GetY()[bin]
            sum_residual += h_efficiency_residual_vs_mass_.GetBinContent(bin)
            if g_efficiency_.GetEYhigh()[bin] > 0:
                sum_error2 += g_efficiency_.GetEYlow()[bin]*g_efficiency_.GetEYhigh()[bin]
            else:
                sum_error2 += g_efficiency_.GetEYlow()[bin]*g_efficiency_.GetEYlow()[bin]
            if bin == 437:
                sum_value += 1.0
        h_efficiency_binned_.SetBinContent(i+1, sum_value/(eff_massBins[i+1] - eff_massBins[i]))
        # p = sum_value/(eff_massBins[i+1] - eff_massBins[i])
        # err = math.sqrt(sum_error2)/(eff_massBins[i+1] - eff_massBins[i])
        # if err == 0:
        #     err = 0.00001
        # fit_residual = (p - eff_fit_.Eval((eff_massBins[i] + eff_massBins[i+1])/2.0))/err
        fit_residual = 0.0
        difference = g_efficiency_binned_.GetY()[i] - f_passed_.Integral(int(eff_massBins[i]), int(eff_massBins[i+1]))/f_total_.Integral(int(eff_massBins[i]), int(eff_massBins[i+1]))
        if difference > 0:
            err = g_efficiency_binned_.GetEYlow()[i]
            if err != 0:
                fit_residual = difference/err
        else:
            err = g_efficiency_binned_.GetEYhigh()[i]
            if err != 0:
                fit_residual = difference/err
        h_efficiency_residual_vs_mass_binned_.SetBinContent(i+1, fit_residual)

        chi2 += pow(fit_residual, 2)
        ndf += 1

    ndf -= 2

    return [chi2, ndf]


def fillSignalPulls(h_signal_, g_data_, h_pulls_):
    h_signal_binned = h_signal_.Rebin(N_massBins, "h_data_varBin", massBins)

    for i in range(0, h_pulls_.GetXaxis().GetNbins()):
        value = h_signal_binned.GetBinContent(i+1)
        h_signal_binned.SetBinContent(i+1, value/h_signal_binned.GetXaxis().GetBinWidth(i+1))
        if g_data_.GetEYlow()[i] != 0.0:
            h_pulls_.SetBinContent(i+1, h_signal_binned.GetBinContent(i+1)/g_data_.GetEYlow()[i])
        else:
            h_pulls_.SetBinContent(i+1, 0.0)
        h_pulls_.SetBinError(i+1, 0.0)
    return


def drawAndSavePlot_background(data_obs_TGraph_, background_TH1_,
                               hist_fit_residual_vsMass_, outputLabel_,
                               list_chi2AndNdf_, list_parameter_, g_efficiency_,
                               eff_fit_, hist_eff_residual_vsMass_,
                               list_chi2AndNdf_eff_,
                               h_efficiency_residual_vs_mass_binned_,
                               h_efficiency_binned_, list_likelihood_,
                               signal_hist_list_):

    global minY, maxY

    eff_canvas = TCanvas("eff_canvas","canvas",W,H)
    eff_canvas.GetWindowHeight()
    eff_canvas.GetWindowWidth()
    eff_canvas.SetTitle("")
    eff_canvas.Divide(1,2,0,0,0)

    eff_canvas.cd(1)
    eff_pad_1 = eff_canvas.GetPad(1)
    eff_pad_1.SetPad(0.01,0.36,0.99,0.98)
    eff_pad_1.SetRightMargin(0.05)
    eff_pad_1.SetTopMargin(0.05)
    eff_pad_1.SetFillColor(0)
    eff_pad_1.SetBorderMode(0)
    eff_pad_1.SetFrameFillStyle(0)
    eff_pad_1.SetFrameBorderMode(0)

    effFrame = eff_pad_1.DrawFrame(efficiency_massMin, 0.0, efficiency_massMax, 1.04)
    effFrame.SetTitle("")
    effFrame.GetXaxis().SetTitleSize(0.06)
    effFrame.GetXaxis().SetTitleOffset(0.95)
    effFrame.GetXaxis().SetLabelSize(0.05)
    effFrame.GetXaxis().SetTitle(xaxisTitle)
    effFrame.GetYaxis().SetTitleSize(0.06)
    effFrame.GetYaxis().SetLabelSize(0.05)
    effFrame.GetYaxis().SetTitle(yaxisTitle_efficiency)

    gStyle.SetErrorX(1)

    eff_fit_.SetLineColor(2)

    g_efficiency_.SetMarkerStyle(20)
    g_efficiency_.SetMarkerColor(17)
    g_efficiency_.SetMarkerSize(0.5)
    g_efficiency_.SetLineColor(17)
    g_efficiency_.SetFillColor(18)
    g_efficiency_.SetTitle("")
    g_efficiency_.GetXaxis().SetTitle(xaxisTitle)
    g_efficiency_.GetXaxis().SetTitle("Trigger Efficiency")
    g_efficiency_.GetXaxis().SetLimits(efficiency_massMin, efficiency_massMax)

    h_efficiency_binned_.SetLineColor(1)
    h_efficiency_binned_.SetMarkerStyle(20)

    #draw objects
    g_efficiency_.Draw("P E2")
    eff_fit_.Draw("C same")
    h_efficiency_binned_.Draw("P E1 same")

    #draw text
    eff_pave_chi2 = TPaveText(0.475, 0.3, 0.825, 0.35, "NDC")
    eff_pave_chi2.SetFillColor(0)
    eff_pave_chi2.SetBorderSize(0)
    eff_pave_chi2.SetFillStyle(0)
    eff_pave_chi2.AddText(0.5, 0.0,
                          "#chi^{{2}} / ndf = {0:.2f} / {1:d} = {2:.2f}".format(
                          list_chi2AndNdf_eff_[0], list_chi2AndNdf_eff_[1],
                          list_chi2AndNdf_eff_[0]/list_chi2AndNdf_eff_[1]))
    eff_pave_chi2.Draw()

    eff_pave_fit = TPaveText(0.45, 0.1, 0.85, 0.3, "NDC")
    eff_pave_fit.SetFillColor(0)
    eff_pave_fit.SetBorderSize(0)
    eff_pave_fit.SetFillStyle(0)
    eff_pave_fit.SetTextColor(13)
    eff_pave_fit.AddText(0.5, 0.8,
                         "p0 = {0:.5g} #pm {1:.5g}".format(list_parameter_[0],
                                                           list_parameter_[6]))
    eff_pave_fit.AddText(0.5, 0.6,
                         "p1 = {0:.4f} #pm {1:.4f}".format(list_parameter_[1],
                                                           list_parameter_[7]))
    eff_pave_fit.AddText(0.5, 0.4,
                         "p2 = {0:.4f} #pm {1:.4f}".format(list_parameter_[2],
                                                           list_parameter_[8]))
    eff_pave_fit.AddText(0.5, 0.2,
                         "p3 = {0:.4f} #pm {1:.4f}".format(list_parameter_[3],
                                                           list_parameter_[9]))
    eff_pave_fit.AddText(0.5, 0.05,
                         "m_{{eff}} = {0:.1f} #pm {1:.1f}".format(list_parameter_[4],
                                                                  list_parameter_[10]))
    eff_pave_fit.AddText(0.5, 0.0,
                         "#sigma_{{eff}} = {0:.2f} #pm {1:.2f}".format(list_parameter_[5],
                                                                       list_parameter_[11]))
    eff_pave_fit.Draw()

    #draw legend
    eff_leg = TLegend(0.5564991,0.58,0.9203575,0.835812)
    eff_leg.SetTextSize(0.03546853)
    eff_leg.SetLineColor(0)
    eff_leg.SetLineStyle(1)
    eff_leg.SetLineWidth(1)
    eff_leg.SetFillColor(0)
    eff_leg.SetFillStyle(0)
    eff_leg.SetMargin(0.35)
    eff_leg.AddEntry(h_efficiency_binned_,"Data" ,"EPL")
    eff_leg.AddEntry(eff_fit_,"Fit","L")
    eff_leg.Draw("SAME")

    #draw pad
    eff_pad_1.RedrawAxis()
    eff_pad_1.Update()
    eff_pad_1.GetFrame().Draw()
    CMS_lumi.CMS_lumi(eff_pad_1, iPeriod, iPos)


    #pad2 - residuals
    eff_canvas.cd(2)
    eff_pad_2 = eff_canvas.GetPad(2)
    eff_pad_2.SetPad(0.01,0.02,0.99,0.37)
    eff_pad_2.SetBottomMargin(0.35)
    eff_pad_2.SetRightMargin(0.05)
    eff_pad_2.SetGridx()
    eff_pad_2.SetGridy()

    effFrame2 = eff_pad_2.DrawFrame(eff_pad_1.GetUxmin(), -eff_range_residual,
                                    eff_pad_1.GetUxmax(), +eff_range_residual)
    effFrame2.SetTitle("")
    effFrame2.SetXTitle(xaxisTitle)
    effFrame2.GetXaxis().SetTitleSize(0.06)
    effFrame2.SetYTitle(yaxisTitle_secondary)
    effFrame2.GetYaxis().SetTitleSize(0.12)
    effFrame2.GetYaxis().SetTitleOffset(0.60)
    effFrame2.GetYaxis().SetLabelSize(0.09)
    effFrame2.GetXaxis().SetTitleSize(0.15)
    effFrame2.GetXaxis().SetTitleOffset(0.90)
    effFrame2.GetXaxis().SetLabelSize(0.12)
    effFrame2.GetXaxis().SetTitle(xaxisTitle)
    effFrame2.GetYaxis().SetTitle(yaxisTitle_secondary)

    #style residuals
    h_efficiency_residual_vs_mass_binned_.GetXaxis().SetRangeUser(efficiency_massMin,efficiency_massMax)
    h_efficiency_residual_vs_mass_binned_.GetYaxis().SetRangeUser(-eff_range_residual,+eff_range_residual)
    h_efficiency_residual_vs_mass_binned_.GetYaxis().SetNdivisions(206,kFALSE)
    h_efficiency_residual_vs_mass_binned_.SetLineWidth(0)
    h_efficiency_residual_vs_mass_binned_.SetFillColor(2)
    h_efficiency_residual_vs_mass_binned_.SetLineColor(1)
    h_efficiency_residual_vs_mass_binned_.GetXaxis().SetTitle(xaxisTitle)
    h_efficiency_residual_vs_mass_binned_.GetYaxis().SetTitle(yaxisTitle_secondary)
    h_efficiency_residual_vs_mass_binned_.GetYaxis().SetTitleSize(0.12)
    h_efficiency_residual_vs_mass_binned_.GetYaxis().SetTitleOffset(0.60)
    h_efficiency_residual_vs_mass_binned_.GetXaxis().SetTitleSize(0.15)
    h_efficiency_residual_vs_mass_binned_.GetXaxis().SetTitleOffset(0.90)
    h_efficiency_residual_vs_mass_binned_.GetYaxis().SetLabelSize(0.09)
    h_efficiency_residual_vs_mass_binned_.GetXaxis().SetLabelSize(0.12)
    h_efficiency_residual_vs_mass_binned_.Draw("HIST")
    hist_eff_residual_vsMass_.GetYaxis().SetRangeUser(-eff_range_residual,
                                                       +eff_range_residual)
    hist_eff_residual_vsMass_.GetYaxis().SetNdivisions(206, kFALSE)
    hist_eff_residual_vsMass_.SetLineWidth(0)
    hist_eff_residual_vsMass_.SetFillColor(2)
    hist_eff_residual_vsMass_.SetLineColor(0)
    hist_eff_residual_vsMass_.SetMarkerStyle(20)
    hist_eff_residual_vsMass_.SetMarkerColor(17)
    hist_eff_residual_vsMass_.SetMarkerSize(0.5)
    hist_eff_residual_vsMass_.Draw("P HIST same")

    eff_line = TLine(efficiency_massMin, 0.0, efficiency_massMax, 0.0)
    eff_line.Draw("")

    #draw pad
    eff_pad_2.RedrawAxis()

    eff_canvas.SaveAs(outputLabel_+"_eff"+".root")
    eff_canvas.SaveAs(outputLabel_+"_eff"+".png")
    #eff_canvas.SaveAs(outputLabel_+"_eff"+".pdf")

    ######################
    ######################
    ######################

    canvas = TCanvas("canvas","canvas",W,H)
    canvas.GetWindowHeight()
    canvas.GetWindowWidth()
    #canvas.SetLogy()
    canvas.SetTitle("")
    canvas.Divide(1,2,0,0,0)

    #pad1 - data spectrum
    canvas.cd(1)
    pad_1 = canvas.GetPad(1)
    #pad_1.SetPad(0.01,0.26,0.99,0.98) #FIXME
    pad_1.SetPad(0.01,0.36,0.99,0.98)
    pad_1.SetLogy()
    pad_1.SetRightMargin(0.05)
    pad_1.SetTopMargin(0.05)
    pad_1.SetFillColor(0)
    pad_1.SetBorderMode(0)
    pad_1.SetFrameFillStyle(0)
    pad_1.SetFrameBorderMode(0)

    if fixedRange==1 and showCrossSection==1:
        #minY = 0.0001/lumi
        minY = 1.0/lumi
        maxY = data_obs_TGraph_.GetY()[0]*10

    #vFrame = pad_1.DrawFrame(minX_mass_plot,0.0001/lumi,maxX_mass_plot,data_obs_TGraph_.GetY()[0]*10)
    vFrame = pad_1.DrawFrame(massMin, minY, massMax, maxY)
    vFrame.SetTitle("")
    #vFrame.SetXTitle(xaxisTitle)
    #vFrame.SetYTitle(yaxisTitle_main)
    vFrame.GetXaxis().SetTitleSize(0.06)
    vFrame.GetXaxis().SetTitleOffset(0.95)
    vFrame.GetXaxis().SetLabelSize(0.05)
    vFrame.GetYaxis().SetTitleSize(0.06)
    #vFrame.GetYaxis().SetTitleOffset(1.0)
    vFrame.GetYaxis().SetLabelSize(0.05)

    #style data spectrum
    gStyle.SetErrorX(1)

    data_obs_TGraph_.SetMarkerStyle(20)
    data_obs_TGraph_.SetMarkerColor(1)
    data_obs_TGraph_.SetLineColor(1)
    data_obs_TGraph_.SetTitle("")
    data_obs_TGraph_.GetXaxis().SetTitle(xaxisTitle)
    data_obs_TGraph_.GetYaxis().SetTitle(yaxisTitle_main)
    data_obs_TGraph_.GetXaxis().SetLimits(massMin,massMax)
    #data_obs_TGraph_.GetYaxis().SetRangeUser(0.0001/lumi,data_obs_TGraph_.GetY()[0]*10)
    data_obs_TGraph_.GetYaxis().SetRangeUser(minY,maxY)

    #style background function
    background_TH1_.SetLineColor(2)
    background_TH1_.SetLineWidth(2)
    background_TH1_.SetTitle("")
    background_TH1_.GetXaxis().SetTitle(xaxisTitle)
    background_TH1_.GetYaxis().SetTitle(yaxisTitle_main)
    background_TH1_.GetXaxis().SetLimits(massMin,massMax)
    background_TH1_.GetYaxis().SetRangeUser(minY,maxY)

    signal_colors = [4, 51, 50, 209]
    signal_styles = [2, 4, 6, 8]

    #draw objects
    #data_obs_TGraph_.Draw("A P E0")
    #background_TH1_.Draw("C SAME")
    background_TH1_.Draw("C")
    data_obs_TGraph_.Draw("P E0 SAME")

    for i, h_sig in enumerate(signal_hist_list_):
        h_sig.SetLineColor(signal_colors[i])
        h_sig.SetLineStyle(signal_styles[i])
        h_sig.SetLineWidth(2)
        h_sig.Draw("HIST ][ SAME")

    #draw text
    pave_sel = TPaveText(0.55, 0.0817972, 0.8, 0.254608, "NDC")
    pave_sel.SetFillColor(0)
    pave_sel.SetBorderSize(0)
    pave_sel.SetFillStyle(0)
    pave_sel.AddText(0.5,1.2,"Wide Jets")
    pave_sel.AddText(0.5,0.5,"m_{jj} > 565 GeV")
    pave_sel.AddText(0.5,0.,"|#eta| < 2.5, |#Delta#eta| < 1.3")
    pave_sel.Draw("SAME")

    pave_chi2 = TPaveText(0.5, 0.25, 0.85, 0.3, "NDC")
    pave_chi2.SetFillColor(0)
    pave_chi2.SetBorderSize(0)
    pave_chi2.SetFillStyle(0)
    pave_chi2.AddText(0.5, 0.0,
                      "#chi^{{2}} / ndf = {0:.2f} / {1:d} = {2:.2f}".format(
                          list_chi2AndNdf_[4], list_chi2AndNdf_[5],
                          list_chi2AndNdf_[4]/list_chi2AndNdf_[5]))
    pave_chi2.Draw("SAME")

    pave_simultaneous_chi2 = TPaveText(0.529489, 0.575, 0.85, 0.675, "NDC")
    pave_simultaneous_chi2.SetFillColor(0)
    pave_simultaneous_chi2.SetBorderSize(0)
    pave_simultaneous_chi2.SetFillStyle(0)
    pave_simultaneous_chi2.SetTextColor(209)
    pave_simultaneous_chi2.AddText(0.5, 0.5, "Combined")
    pave_simultaneous_chi2.AddText(0.5, 0.0,
                                   "#chi^{{2}} / ndf = {0:.2f} / {1:d} = {2:.2f}".format(
                                   list_chi2AndNdf_[4] + list_chi2AndNdf_eff_[0],
                                   list_chi2AndNdf_[5] + list_chi2AndNdf_eff_[1] + 2,
                                   (list_chi2AndNdf_[4] + list_chi2AndNdf_eff_[0])/(list_chi2AndNdf_[5] + list_chi2AndNdf_eff_[1] + 2.0)))
    pave_simultaneous_chi2.Draw()

    pave_simultaneous_nll = TPaveText(0.65, 0.465, 0.925, 0.565, "NDC")
    pave_simultaneous_nll.SetFillColor(0)
    pave_simultaneous_nll.SetBorderSize(0)
    pave_simultaneous_nll.SetFillStyle(0)
    pave_simultaneous_nll.SetTextColor(209)
    pave_simultaneous_nll.AddText(0.5, 0.67,
                                  "NLL_{{eff}} = {0}".format(list_likelihood_[2]))
    pave_simultaneous_nll.AddText(0.5, 0.33,
                                  "NLL_{{mjj}} = {0}".format(list_likelihood_[1]))
    pave_simultaneous_nll.AddText(0.5, 0.0,
                                  "NLL_{{tot}} = {0}".format(list_likelihood_[0]))
    pave_simultaneous_nll.Draw()

    pave_fit = TPaveText(0.15, 0.34, 0.55, 0.54, "NDC")
    pave_fit.SetFillColor(0)
    pave_fit.SetBorderSize(0)
    pave_fit.SetFillStyle(0)
    pave_fit.SetTextColor(13)
    pave_fit.AddText(0.5, 0.8,
                     "p0 = {0:.5g} #pm {1:.5g}".format(list_parameter_[0],
                                                       list_parameter_[6]))
    pave_fit.AddText(0.5, 0.6,
                     "p1 = {0:.4f} #pm {1:.4f}".format(list_parameter_[1],
                                                       list_parameter_[7]))
    pave_fit.AddText(0.5, 0.4,
                     "p2 = {0:.4f} #pm {1:.4f}".format(list_parameter_[2],
                                                       list_parameter_[8]))
    pave_fit.AddText(0.5, 0.2,
                     "p3 = {0:.4f} #pm {1:.4f}".format(list_parameter_[3],
                                                       list_parameter_[9]))
    pave_fit.AddText(0.5, 0.05,
                     "m_{{eff}} = {0:.1f} #pm {1:.1f}".format(list_parameter_[4],
                                                              list_parameter_[10]))
    pave_fit.AddText(0.5, 0.0,
                     "#sigma_{{eff}} = {0:.2f} #pm {1:.2f}".format(list_parameter_[5],
                                                                   list_parameter_[11]))
    pave_fit.Draw("SAME")

    # pave_toy = TPaveText(0.6,0.4,0.9,0.6,"NDC")
    # pave_toy.SetFillColor(0)
    # pave_toy.SetBorderSize(0)
    # pave_toy.SetFillStyle(0)
    # pave_toy.SetTextColor(2)
    # pave_toy.AddText(0.5, 0.0, "Toy MC (sbtoy1)")
    # pave_toy.Draw("SAME")

    #draw legend
    leg = TLegend(0.5564991,0.65,0.9203575,0.89)
    leg.SetTextSize(0.03546853)
    leg.SetLineColor(0)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetMargin(0.35)
    leg.AddEntry(data_obs_TGraph_,"Data" ,"EPL")
    leg.AddEntry(background_TH1_,"Fit","L")
    for h_sig in signal_hist_list_:
        leg.AddEntry(h_sig, h_sig.GetTitle(), "L")
    leg.Draw("SAME")

    signal_pulls_list = []
    for i, h_sig in enumerate(signal_hist_list_):
        h_sig_pulls = hist_fit_residual_vsMass_.Clone("{0}_pulls".format(h_sig.GetName()))
        fillSignalPulls(h_sig, data_obs_TGraph_, h_sig_pulls)
        h_sig_pulls.SetLineColor(signal_colors[i])
        h_sig_pulls.SetLineStyle(signal_styles[i])
        h_sig_pulls.SetLineWidth(2)
        signal_pulls_list.append(h_sig_pulls)

    #draw pad
    pad_1.RedrawAxis()
    pad_1.Update()
    pad_1.GetFrame().Draw()
    CMS_lumi.CMS_lumi(pad_1, iPeriod, iPos)


    #pad2 - residuals
    canvas.cd(2)
    pad_2 = canvas.GetPad(2)
    #pad_2.SetPad(0.01,0.02,0.99,0.27) #FIXME
    pad_2.SetPad(0.01,0.02,0.99,0.37)
    pad_2.SetBottomMargin(0.35)
    pad_2.SetRightMargin(0.05)
    pad_2.SetGridx()
    pad_2.SetGridy()

    vFrame2 = pad_2.DrawFrame(pad_1.GetUxmin(), -range_residual, pad_1.GetUxmax(), +range_residual)    
    vFrame2.SetTitle("")
    vFrame2.SetXTitle(xaxisTitle)
    vFrame2.GetXaxis().SetTitleSize(0.06)
    vFrame2.SetYTitle(yaxisTitle_secondary)
    vFrame2.GetYaxis().SetTitleSize(0.12)
    vFrame2.GetYaxis().SetTitleOffset(0.60)
    vFrame2.GetYaxis().SetLabelSize(0.09)
    vFrame2.GetXaxis().SetTitleSize(0.15)
    vFrame2.GetXaxis().SetTitleOffset(0.90)
    vFrame2.GetXaxis().SetLabelSize(0.12)

    #style residuals
    hist_fit_residual_vsMass_.GetXaxis().SetRangeUser(massMin, massMax)
    hist_fit_residual_vsMass_.GetYaxis().SetRangeUser(-range_residual,
                                                      +range_residual)
    hist_fit_residual_vsMass_.GetYaxis().SetNdivisions(206, kFALSE)
    hist_fit_residual_vsMass_.SetLineWidth(0)
    hist_fit_residual_vsMass_.SetFillColor(2)
    hist_fit_residual_vsMass_.SetLineColor(1)
    hist_fit_residual_vsMass_.Draw("SAME HIST")
    for h_sig_pulls in signal_pulls_list:
        h_sig_pulls.Draw("SAME HIST")

    line = TLine(massMin,0,massMax,0)
    line.Draw("")

    #draw pad
    pad_2.RedrawAxis()

    #============

    #write canvas
    canvas.SaveAs(outputLabel_+"_B"+".root")
    canvas.SaveAs(outputLabel_+"_B"+".png")
    #canvas.SaveAs(outputLabel_+"_B"+".pdf")

    raw_input("Press Enter to exit...")


#==============================================================================



if __name__ == '__main__':
    main()
