#!/usr/bin/env python
import sys, os, copy, re
from array import array
from argparse import ArgumentParser
import math
from ROOT import *
import CMS_lumi, setTDRStyle

# Global variables
dataset = 0
if dataset == 0:
    inputHistoFileName = "rawhistV10_RECO_2016BCD_rounds45678_Spring16_25nsV6_ICHEP_12.root"
elif dataset == 1:
    inputHistoFileName = "data_CaloScoutingHT_Run2016BCD_NewBiasCorrectedFlat_Golden12910pb_CaloDijet2016.root"
inputEfficiencyFileName = "triggerEfficiency_L1HTT150seed_HT450_DetaJJLess1p3_HLTv7Corr_output.root"
inputDataFileName = "/t3/users/santanas/Dijet13TeVScouting/rootTrees_reduced/ScoutingPFHT__15_01_2016_20160115_192039/merged/rootfile_ScoutingPFHT__Run2015D-v1__RAW_ScoutingPFHT__15_01_2016_20160115_192039_reduced_skim.root"
treeName = "rootTupleTree/tree"
outputLabel = "output"
var = "mjj"
sqrtS = 13000.0
lumiValue = 12910.0 #[pb]
showCrossSection = 1 #1=cross section [pb] , 0=number of events/GeV
fixedRange = 1 #1=YES , 0=NO  (the option works only if showCrossSection=1; otherwise=0)
minY = 0.00000003
maxY = 20
if showCrossSection==1:
    lumi = lumiValue
else:
    lumi = 1
if dataset == 0:
    massMin = 1058
    massMax = 7866
elif dataset == 1:
    massMin = 453
    massMax = 2037
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
doSimultaneousFit = False

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

#iPos = 11
iPos = 0
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
    inputDataFile = TFile(inputHistoFileName)
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
    if dataset == 0:
        h_data = (inputDataFile.Get("mjj_gev")).Clone("h_data")
    elif dataset == 1 or dataset == 2:
        h_data = (inputDataFile.Get("h_mjj_HLTpass_HT250_1GeVbin")).Clone("h_data")
    dataInt = h_data.Integral(massMin + 1, massMax)
    print "number of events in the fit range:", int(dataInt)

    background2 = TF1("background2",
                      "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)+[4]*log(x/%.1f)*log(x/%.1f)))"%(sqrtS,sqrtS,sqrtS,sqrtS,sqrtS),
                      float(massMin), float(massMax))
    if dataset == 0:
        background2.SetParameter(0, 3.16766240835e-07)
        background2.SetParameter(2, 6.27895318762)
    elif dataset == 1:
        background2.SetParameter(0, 3.3110576299e-06)
        background2.SetParameter(2, 5.40474642534)
    background2.SetParameter(1, 0.0)
    background2.SetParameter(3, 0.0)
    background2.SetParameter(4, 0.0)

    background3 = TF1("background3",
                      "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)+[4]*log(x/%.1f)*log(x/%.1f)))"%(sqrtS,sqrtS,sqrtS,sqrtS,sqrtS),
                      float(massMin), float(massMax))
    if dataset == 0:
        background3.SetParameter(0, 1.51335086414e-05)
        background3.SetParameter(1, 9.05904632837)
        background3.SetParameter(2, 5.02571873249)
    elif dataset == 1:
        background3.SetParameter(0, 8.58951893126e-05)
        background3.SetParameter(1, 14.3144967628)
        background3.SetParameter(2, 4.57510842619)
    background3.SetParameter(3, 0.0)
    background3.SetParameter(4, 0.0)

    background4 = TF1("background4",
                      "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)+[4]*log(x/%.1f)*log(x/%.1f)))"%(sqrtS,sqrtS,sqrtS,sqrtS,sqrtS),
                      float(massMin), float(massMax))
    if dataset == 0:
        background4.SetParameter(0, 3.46397990934e-06)
        background4.SetParameter(1, 7.34961109006)
        background4.SetParameter(2, 5.96523987147)
        background4.SetParameter(3, 0.163614755034)
    elif dataset == 1:
        background4.SetParameter(0, 9.49994522633e-07)
        background4.SetParameter(1, 5.84090837867)
        background4.SetParameter(2, 6.83225395106)
        background4.SetParameter(3, 0.299718539449)
    background4.SetParameter(4, 0.0)

    background5 = TF1("background5",
                      "([0]*TMath::Power(1-x/%.1f,[1]))/(TMath::Power(x/%.1f,[2]+[3]*log(x/%.1f)+[4]*log(x/%.1f)*log(x/%.1f)))"%(sqrtS,sqrtS,sqrtS,sqrtS,sqrtS),
                      float(massMin), float(massMax))
    if dataset == 0:
        background5.SetParameter(0, 2.38002753645e-05)
        background5.SetParameter(1, 8.87638444646)
        background5.SetParameter(2, 4.18327834544)
        background5.SetParameter(3, -0.457885218595)
        background5.SetParameter(4, -0.0784871491669)
    elif dataset == 1:
        background5.SetParameter(0, 9.52739808228e-07)
        background5.SetParameter(1, 5.84779556315)
        background5.SetParameter(2, 6.83125819802)
        background5.SetParameter(3, 0.299773743681)
        background5.SetParameter(4, 3.51085801569e-05)

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

    list_parameter = [9.5006023175e-07, 6.23424107417, 6.82675441147, 0.312944823058, 0.0,
                      0.0, 0.0,
                      0.0, 0.0, 0.0, 0.0, 0.0,
                      0.0, 0.0]

    mass_1GeV_list = []
    for bin in range(efficiency_massMin, efficiency_massMax+1):
        mass_1GeV_list.append(bin)
    eff_1GeV_massBins = array("d", mass_1GeV_list)
    # plot
    h_background2 = convertFunctionToHisto(background2,"h_background2",N_massBins,massBins)
    h_background3 = convertFunctionToHisto(background3,"h_background3",N_massBins,massBins)
    h_background4 = convertFunctionToHisto(background4,"h_background4",N_massBins,massBins)
    h_background5 = convertFunctionToHisto(background5,"h_background5",N_massBins,massBins)
    h_fit_residual_vs_mass2 = TH1D("h_fit_residual_vs_mass2","h_fit_residual_vs_mass2",N_massBins,massBins)
    h_fit_residual_vs_mass3 = TH1D("h_fit_residual_vs_mass3","h_fit_residual_vs_mass3",N_massBins,massBins)
    h_fit_residual_vs_mass4 = TH1D("h_fit_residual_vs_mass4","h_fit_residual_vs_mass4",N_massBins,massBins)
    h_fit_residual_vs_mass5 = TH1D("h_fit_residual_vs_mass5","h_fit_residual_vs_mass5",N_massBins,massBins)
    list_chi2AndNdf_background2 = calculateChi2AndFillResiduals(g_data,
                                                               h_background2,
                                                               h_fit_residual_vs_mass2,
                                                               2, 0)
    list_chi2AndNdf_background3 = calculateChi2AndFillResiduals(g_data,
                                                               h_background3,
                                                               h_fit_residual_vs_mass3,
                                                               3, 0)
    list_chi2AndNdf_background4 = calculateChi2AndFillResiduals(g_data,
                                                               h_background4,
                                                               h_fit_residual_vs_mass4,
                                                               4, 0)
    list_chi2AndNdf_background5 = calculateChi2AndFillResiduals(g_data,
                                                               h_background5,
                                                               h_fit_residual_vs_mass5,
                                                               5, 0)
    backgrounds = [h_background2, h_background3, h_background4, h_background5]
    residuals = [h_fit_residual_vs_mass2, h_fit_residual_vs_mass3,
                 h_fit_residual_vs_mass4, h_fit_residual_vs_mass5]
    drawAndSavePlot_background(g_data, backgrounds, residuals, outputLabel)

    #output.Close()

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
    print "chi2/ndf PlotRangeNonZero : %.8f / %d = %.2f" % ( chi2_PlotRangeNonZero , ndf_PlotRangeNonZero , chi2_ndf_PlotRangeNonZero ) 
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
                               hist_fit_residual_vsMass_, outputLabel_):

    global minY, maxY

    canvas = TCanvas("canvas","canvas",W,H)
    canvas.GetWindowHeight()
    canvas.GetWindowWidth()
    canvas.SetTitle("")
    canvas.Divide(1,2,0,0,0)

    canvas.cd(1)
    pad_1 = canvas.GetPad(1)
    pad_1.SetPad(0.01,0.36,0.99,0.98)
    pad_1.SetLogy()
    pad_1.SetRightMargin(0.05)
    pad_1.SetTopMargin(0.05)
    pad_1.SetFillColor(0)
    pad_1.SetBorderMode(0)
    pad_1.SetFrameFillStyle(0)
    pad_1.SetFrameBorderMode(0)

    if fixedRange==1 and showCrossSection==1:
        if dataset == 0:
            minY = 0.0002/lumi
            maxY = data_obs_TGraph_.GetY()[0]*3
        elif dataset == 1 or dataset == 2:
            minY = 200.0/lumi
            maxY = data_obs_TGraph_.GetY()[0]*2

    vFrame = pad_1.DrawFrame(massMin, minY, massMax, maxY)
    vFrame.SetTitle("")
    vFrame.GetXaxis().SetTitleSize(0.06)
    vFrame.GetXaxis().SetTitleOffset(0.95)
    vFrame.GetXaxis().SetLabelSize(0.05)
    vFrame.GetYaxis().SetTitleSize(0.06)
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
    data_obs_TGraph_.GetYaxis().SetRangeUser(minY,maxY)

    #style background function
    for background_TH1 in background_TH1_:
        background_TH1.SetLineWidth(2)
        background_TH1.SetTitle("")
        background_TH1.GetXaxis().SetTitle(xaxisTitle)
        background_TH1.GetYaxis().SetTitle(yaxisTitle_main)
        background_TH1.GetXaxis().SetLimits(massMin,massMax)
        background_TH1.GetYaxis().SetRangeUser(minY,maxY)
    background_TH1_[0].SetLineColor(209)
    background_TH1_[1].SetLineColor(4)
    background_TH1_[2].SetLineColor(2)
    background_TH1_[3].SetLineColor(90)
    background_TH1_[0].SetLineStyle(5)
    background_TH1_[1].SetLineStyle(9)
    background_TH1_[2].SetLineStyle(1)
    background_TH1_[3].SetLineStyle(2)

    signal_colors = [4, 51, 50, 209]
    signal_styles = [2, 4, 6, 8]

    #draw objects
    background_TH1_[0].Draw("C")
    background_TH1_[1].Draw("C SAME")
    background_TH1_[2].Draw("C SAME")
    background_TH1_[3].Draw("C SAME")
    data_obs_TGraph_.Draw("P E0 Z SAME")

    #draw text
    pave_sel = TPaveText(0.02,0.015,0.42,0.215,"NDC")
    pave_sel.SetFillColor(0)
    pave_sel.SetBorderSize(0)
    pave_sel.SetFillStyle(0)
    pave_sel.SetTextAlign(11);
    pave_sel.SetTextFont(42);
    pave_sel.SetTextSize(0.045);
    pave_sel.AddText(0.5,1.2,"Wide Jets")
    if dataset == 0:
        pave_sel.AddText(0.5,0.5,"1.1 < m_{jj} < 7.9 TeV")
    elif dataset == 1 or dataset == 2:
        pave_sel.AddText(0.5,0.5,"453 < m_{jj} < 2037 GeV")
    pave_sel.AddText(0.5,0.,"|#eta| < 2.5, |#Delta#eta| < 1.3")
    pave_sel.Draw("SAME")

    leg = TLegend(0.6,0.52,0.89,0.88);
    leg.SetBorderSize(1);
    leg.SetLineColor(0);
    leg.SetLineStyle(1);
    leg.SetLineWidth(0);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    entry_data = leg.AddEntry(data_obs_TGraph_,"Data","pe");
    entry_data.SetLineColor(1);
    entry_data.SetLineStyle(1);
    entry_data.SetLineWidth(1);
    entry_data.SetMarkerColor(1);
    entry_data.SetMarkerStyle(20);
    entry_data.SetMarkerSize(0.9);
    entry_data.SetTextFont(42);
    entry_fit2 = leg.AddEntry(background_TH1_[0],"2 Parameter Fit","l");
    entry_fit2.SetLineColor(2);
    entry_fit2.SetLineStyle(1);
    entry_fit2.SetLineWidth(2);
    entry_fit2.SetMarkerColor(1);
    entry_fit2.SetMarkerStyle(21);
    entry_fit2.SetMarkerSize(1);
    entry_fit2.SetTextFont(42);
    entry_fit3 = leg.AddEntry(background_TH1_[1],"3 Parameter Fit","l");
    entry_fit3.SetLineColor(2);
    entry_fit3.SetLineStyle(1);
    entry_fit3.SetLineWidth(2);
    entry_fit3.SetMarkerColor(1);
    entry_fit3.SetMarkerStyle(21);
    entry_fit3.SetMarkerSize(1);
    entry_fit3.SetTextFont(42);
    entry_fit4 = leg.AddEntry(background_TH1_[2],"4 Parameter Fit","l");
    entry_fit4.SetLineColor(2);
    entry_fit4.SetLineStyle(1);
    entry_fit4.SetLineWidth(2);
    entry_fit4.SetMarkerColor(1);
    entry_fit4.SetMarkerStyle(21);
    entry_fit4.SetMarkerSize(1);
    entry_fit4.SetTextFont(42);
    entry_fit5 = leg.AddEntry(background_TH1_[3],"5 Parameter Fit","l");
    entry_fit5.SetLineColor(2);
    entry_fit5.SetLineStyle(1);
    entry_fit5.SetLineWidth(2);
    entry_fit5.SetMarkerColor(1);
    entry_fit5.SetMarkerStyle(21);
    entry_fit5.SetMarkerSize(1);
    entry_fit5.SetTextFont(42);
    leg.Draw();

    signal_pulls_list = []

    #draw pad
    pad_1.RedrawAxis()
    pad_1.Update()
    pad_1.GetFrame().Draw()
    CMS_lumi.CMS_lumi(pad_1, iPeriod, iPos)


    #pad2 - residuals
    canvas.cd(2)
    pad_2 = canvas.GetPad(2)
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
    for hist_fit_residual_vsMass in hist_fit_residual_vsMass_:
        hist_fit_residual_vsMass.GetXaxis().SetRangeUser(massMin, massMax)
        hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-range_residual,
                                                         +range_residual)
        hist_fit_residual_vsMass.GetYaxis().SetNdivisions(206, kFALSE)
        hist_fit_residual_vsMass.SetLineWidth(2)
    hist_fit_residual_vsMass_[0].SetLineColor(209)
    hist_fit_residual_vsMass_[1].SetLineColor(4)
    hist_fit_residual_vsMass_[2].SetLineColor(2)
    hist_fit_residual_vsMass_[3].SetLineColor(90)
    hist_fit_residual_vsMass_[0].SetLineStyle(5)
    hist_fit_residual_vsMass_[1].SetLineStyle(9)
    hist_fit_residual_vsMass_[2].SetLineStyle(1)
    hist_fit_residual_vsMass_[3].SetLineStyle(2)
    hist_fit_residual_vsMass_[0].SetFillColor(209)
    hist_fit_residual_vsMass_[1].SetFillColor(4)
    hist_fit_residual_vsMass_[2].SetFillColor(2)
    hist_fit_residual_vsMass_[3].SetFillColor(90)
    hist_fit_residual_vsMass_[0].SetFillStyle(3004)
    hist_fit_residual_vsMass_[1].SetFillStyle(3005)
    hist_fit_residual_vsMass_[2].SetFillStyle(3006)
    hist_fit_residual_vsMass_[3].SetFillStyle(3007)
    for hist_fit_residual_vsMass in hist_fit_residual_vsMass_:
        hist_fit_residual_vsMass.Draw("SAME HIST")
    for h_sig_pulls in signal_pulls_list:
        h_sig_pulls.Draw("SAME HIST")

    line = TLine(massMin,0,massMax,0)
    line.Draw("")

    #draw pad
    pad_2.RedrawAxis()

    #============

    #write canvas
    #canvas.SaveAs(outputLabel_+"_B"+".root")
    #canvas.SaveAs(outputLabel_+"_B"+".png")
    #canvas.SaveAs(outputLabel_+"_B"+".pdf")

    canvas.Update()

    raw_input("Press Enter to exit...")


#==============================================================================



if __name__ == '__main__':
    main()
