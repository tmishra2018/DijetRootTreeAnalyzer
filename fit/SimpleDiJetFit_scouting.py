#!/usr/bin/python

from ROOT import *
from array import array
import CMS_lumi, tdrstyle
import subprocess
import os
import imp
import multiprocessing
from itertools import repeat
import math
import optparse

gStyle.SetOptFit(1111)
#set the tdr style
tdrstyle.setTDRStyle()

lumi=1806
#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_13TeV = str(int(lumi))+" pb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = ""
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--plotSig",action="store_true",dest="plotSig",default=False)

(options, args) = parser.parse_args()
plotSig = options.plotSig


#Fit function
# 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#
#1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#
#2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#   --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
#
#3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
#   --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
#
#4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )"
#   --> "log" extension wrt to DEFAULT
#
#5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )"
#   --> "log" extension wrt to DEFAULT
#
#6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )"
#   --> "log" extension wrt to DEFAULT
#
number_of_variableWidth_bins = 103

massBins =[1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430,10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000];

v_massBins = array("d",massBins)

minX_mass = 606
maxX_mass = 2037

FunctionType = 10
fileNameSuffix = "RunD_4param_JEC_Summer15_25nsV6_Silver_1806pb-1_plotSig"

####### INPUT #############
# data
input_root_file="histograms.root"

### input file and 1D histo
file0 = TFile.Open( input_root_file )
input_1Dhistogram = "h_mjj"

hist_mass_original = file0.Get(input_1Dhistogram)
hist_binned = hist_mass_original.Rebin(number_of_variableWidth_bins,"hist_binned",v_massBins)
hist_mass = TH1F("hist_mass","",number_of_variableWidth_bins,v_massBins)

##########OUTPUT########
outputDir="fit_directory/"
os.system("mkdir -p "+outputDir)

#================================================================================================================

def main():
    for  i in range (1, number_of_variableWidth_bins):
        #data
        bincontent = hist_binned.GetBinContent(i)
        binwidth = hist_binned.GetBinWidth(i)
        binerror = hist_binned.GetBinError(i)
        hist_mass.SetBinContent(i,bincontent/(binwidth*lumi))

    #######################################################
    #data in TGraph format (hist binned)
    alpha = 1 - 0.6827;

    x=[]
    y=[]
    exl=[]
    exh=[]
    eyl=[]
    eyh=[]
    x_mc=[]
    y_mc=[]

    for i in range(0,number_of_variableWidth_bins):
        n    = hist_binned.GetBinContent(i+1)
        dm   = hist_binned.GetBinWidth(i+1)
        mass = hist_binned.GetBinCenter(i+1)
        xl   = hist_binned.GetBinLowEdge(i+1)
        xh   = xl+dm
        x.append( (xl+xh)/2.)
        exl.append( dm/2.)
        exh.append( dm/2.)
        y.append( n / (dm*lumi))
        l = 0.5*TMath.ChisquareQuantile(alpha/2,2*n)
        h = 0.5*TMath.ChisquareQuantile(1-alpha/2,2*(n+1))
        eyl.append( (n-l)/(lumi*dm) )
        eyh.append( (h-n)/(lumi*dm) )

    vx = array("f",x)
    vy = array("f",y)
    vexl = array("f",exl)
    vexh = array("f",exh)
    veyl = array("f",eyl)
    veyh = array("f",eyh)
    vx_mc = array("f",x_mc)
    vy_mc = array("f",y_mc)

    #data in TGraph format
    g = TGraphAsymmErrors(number_of_variableWidth_bins,vx,vy,vexl,vexh,veyl,veyh)
    g.SetName("g_data")
    #g.Print()


    nBins_fit = hist_mass.FindBin(maxX_mass)- hist_mass.FindBin(minX_mass)
    ##count zero bins
    zeroBins = 0
    for i in range(0,nBins_fit):
        if hist_mass.GetBinContent(hist_mass.FindBin(minX_mass)+i)==0:
            zeroBins +=1

    nBins_fit = nBins_fit-zeroBins
    fitresult = doFitAndChi2(FunctionType,hist_mass,g,hist_mass_original)
    M1Bkg = fitresult[3]
    hist_fit_residual_vsMass = fitresult[4]
    nPar = nBins_fit - fitresult[1]
    chi2 = fitresult[2]
    ndof = fitresult[1]
    DrawFit(g, M1Bkg, hist_fit_residual_vsMass, FunctionType, nPar,
            fileNameSuffix, chi2, ndof)


def blinded4p(x, par):
    if x[0] >= 693.0 and x[0] < 838.0:
        TF1.RejectPoint()
        return 0.0;
    return par[0]*TMath.Power(1.0 - x[0]/13000.0, par[1])/TMath.Power(x[0]/13000.0, par[2] + par[3]*TMath.Log(x[0]/13000.0))*(0.5*(1.0 + TMath.Erf((x[0] - 507.1)/94.2)))
    #return par[0]*TMath.Power(1.0 - x[0]/13000.0, par[1])/TMath.Power(x[0]/13000.0, par[2] + par[3]*TMath.Log(x[0]/13000.0))

def doFitAndChi2(FunctionType,hist_mass,g,hist_mass_original):
    ### fit mass histogram with background function
    # -2: VARIATION-1 (3 par.) - " [0] / ( TMath::Power(x/13000,[2]) )"
    if( FunctionType==-2 ):
        nPar=2
        M1Bkg = TF1("M1Bkg"," [0] / ( TMath::Power(x/13000,[1]) )",minX_mass,maxX_mass)
        M1Bkg.SetParameter(0,0.08)
        M1Bkg.SetParameter(1,2)
        M1Bkg.SetParLimits(0,0.,10);
        M1Bkg.SetParLimits(1,0.,20);


    # -1: VARIATION-1 (3 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]) )"
    if( FunctionType==-1 ):
        nPar=3
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]) )",minX_mass,maxX_mass)
        M1Bkg.SetParameter(0,0.00001)
        M1Bkg.SetParameter(1,10)
        M1Bkg.SetParameter(2,5)

        #1 fb-1
        M1Bkg.SetParLimits(0,0.,1000);
        #10 fb-1
        #M1Bkg->SetParLimits(0,0,1.);
        #M1Bkg.SetParLimits(1,0.,100)
        #M1Bkg.SetParLimits(2,0.,20)

    # 0: DEFAULT (4 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    if( FunctionType==0 ):
        nPar=4
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
        #M1Bkg.SetParameter(0,0.0014)
        M1Bkg.SetParameter(0, 0.00274215) # 1.4
        M1Bkg.SetParameter(1, 6.96048) # 12
        M1Bkg.SetParameter(2, 6.6761) # 2
        M1Bkg.SetParameter(3, 0.270368) # -0.5

        #M1Bkg.SetParLimits(0,0.,1000);
        #1 fb-1
        #M1Bkg.SetParLimits(0,0.,10);
        #10 fb-1
        #M1Bkg->SetParLimits(0,0,1.);
        # M1Bkg.SetParLimits(1,0.,10000)
        # M1Bkg.SetParLimits(2,1.,10000)
        # M1Bkg.SetParLimits(1,0,100)
        # M1Bkg.SetParLimits(2,0.,20)
        # M1Bkg.SetParLimits(3,-10.,10)
        # M1Bkg.FixParameter(3,0.)


    # 1: VARIATION 1 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    if( FunctionType==1 ):
        nPar=5
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass)
        M1Bkg.SetParameter(0,0.005)
        M1Bkg.SetParameter(1,9.3)
        M1Bkg.SetParameter(2,7.2)
        M1Bkg.SetParameter(3,0.4)
        M1Bkg.SetParameter(4,3.1)
    #   M1Bkg.SetParLimits(4,-5,5)

    # 2: VARIATION 2 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    #  --> 2nd order poly extension : inspired by HERA PDF 1.0 [http://arxiv.org/abs/arXiv:0911.0884 , Eq. 4.1]
    if( FunctionType==2 ):
        nPar=6
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*(1+[4]*x/13000+[5]*pow(x/13000,2)) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )",minX_mass,maxX_mass);
        M1Bkg.SetParameter(0,0.005)
        M1Bkg.SetParameter(1,9.3)
        M1Bkg.SetParameter(2,7.2)
        M1Bkg.SetParameter(3,0.4)
        M1Bkg.SetParameter(4,3.1)
        M1Bkg.SetParameter(5,25.6)
        #M1Bkg.SetParLimits(5,10,50)


    # 3: VARIATION 3 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    #    --> "exponential" extension wrt to DEFAULT - inspired by CTEQ 2008 [http://arxiv.org/pdf/hep-ph/0201195v3.pdf , Eq. 4]
    if( FunctionType==3 ) :
        nPar=7
        M1Bkg =  TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1])*exp([4]*x/13000)*TMath::Power(1+exp([5])*x/13000,[6]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )" ,minX_mass,maxX_mass)
        M1Bkg.SetParameter(0,0.005)
        M1Bkg.SetParameter(1,15.1)
        M1Bkg.SetParameter(2,7.2)
        M1Bkg.SetParameter(3,0.4)
        M1Bkg.SetParameter(4,13.0)
        M1Bkg.SetParameter(5,-4.0)
        M1Bkg.SetParameter(6,70.0)
        # M1Bkg.SetParLimits(4,-1,1)


    # 4: VARIATION 4 (5 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )"
    #    --> "log" extension wrt to DEFAULT
    if( FunctionType==4 ) :
        nPar=5
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)) )",minX_mass,maxX_mass);
        M1Bkg.SetParameter(0,0.08)
        M1Bkg.SetParameter(1,12)
        M1Bkg.SetParameter(2,2)
        M1Bkg.SetParameter(3,-0.5)
        M1Bkg.SetParameter(4,0.)
        #M1Bkg.SetParLimits(3,0,0.4)
        #M1Bkg.FixParameter(3,0.)
        #M1Bkg.FixParameter(4,0.)
        #M1Bkg.SetParLimits(1,0.,10000)
        #M1Bkg.SetParLimits(2,1.,10000)
        #M1Bkg.SetParLimits(0,0.,10);
        #M1Bkg.SetParLimits(1,0.,100)
        #M1Bkg.SetParLimits(2,0.,20)
        #M1Bkg.SetParLimits(3,-10,10)
        #M1Bkg.SetParLimits(4,-10,10)


    # 5: VARIATION 5 (6 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )"
    #    --> "log" extension wrt to DEFAULT
    if( FunctionType==5 ):
        nPar=6
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)) )",minX_mass,maxX_mass)
        M1Bkg.SetParameter(0,0.08)
        M1Bkg.SetParameter(1,12)
        M1Bkg.SetParameter(2,2)
        M1Bkg.SetParameter(3,-0.5)
        M1Bkg.SetParameter(4,0.)
        M1Bkg.SetParameter(5,0.)
        #M1Bkg.SetParLimits(3,0,0.4)
        #M1Bkg.FixParameter(3,0.)
        #M1Bkg.FixParameter(4,0.)
        #M1Bkg.FixParameter(5,0.)
        #M1Bkg.SetParLimits(1,0.,10000)
        #M1Bkg.SetParLimits(2,1.,10000)
        #M1Bkg.SetParLimits(3,0.,10000)
        #M1Bkg.SetParLimits(0,0.,10);
        #M1Bkg.SetParLimits(1,1.,100)
        #M1Bkg.SetParLimits(2,1.,20)
        #M1Bkg.SetParLimits(3,-10,10)
        #M1Bkg.SetParLimits(4,-10,10)
        #M1Bkg.SetParLimits(5,-10,10)


    # 6: VARIATION 6 (7 par.) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )"
    #  //    --> "log" extension wrt to DEFAULT
    if( FunctionType==6 ) :
        nPar=7
        M1Bkg = TF1("M1Bkg","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)+[4]*TMath::Power(log(x/13000),2)+[5]*TMath::Power(log(x/13000),3)+[6]*TMath::Power(log(x/13000),4)) )",minX_mass,maxX_mass)
        M1Bkg.SetParameter(0,0.08)
        M1Bkg.SetParameter(1,12)
        M1Bkg.SetParameter(2,2)
        M1Bkg.SetParameter(3,-0.5)
        M1Bkg.SetParameter(4,0.)
        M1Bkg.SetParameter(5,0.)
        M1Bkg.SetParameter(6,0.)
        #M1Bkg.SetParLimits(3,0,0.4)
        #M1Bkg.FixParameter(3,0.)
        #M1Bkg.FixParameter(4,0.)
        #M1Bkg.FixParameter(5,0.)
        #M1Bkg.FixParameter(6,0.)
        #M1Bkg.SetParLimits(1,0.,10000)
        #M1Bkg.SetParLimits(2,1.,10000)
        #M1Bkg.SetParLimits(3,0.,10000)
        #M1Bkg.SetParLimits(0,0.,10);
        #M1Bkg.SetParLimits(1,0.,100)
        #M1Bkg.SetParLimits(2,0.,20)
        #M1Bkg.SetParLimits(3,-10,10)
        #M1Bkg.SetParLimits(4,-10,10)
        #M1Bkg.SetParLimits(5,-10,10)
        #M1Bkg.SetParLimits(6,-10,10)

    # 10: (4 par. blind) - "( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )"
    if( FunctionType==10 ):
        nPar=4
        M1Bkg = TF1("M1Bkg", blinded4p, minX_mass, maxX_mass, 4)
        # M1Bkg.SetParameter(0, 0.00113356)
        # M1Bkg.SetParameter(1, 6.3042)
        # M1Bkg.SetParameter(2, 6.95734)
        # M1Bkg.SetParameter(3, 0.310343)
        M1Bkg.SetParameter(0, 0.000425498)
        M1Bkg.SetParameter(1, 4.91305)
        M1Bkg.SetParameter(2, 7.52418)
        M1Bkg.SetParameter(3, 0.399293)



    #TFitResultPtr r;
    stopProgram=1;
    for loop in range (0,10):
        r = hist_mass_original.Fit("M1Bkg", "ELSR", "", minX_mass, maxX_mass)
        #r = hist_mass_original.Fit("M1Bkg", "ELLSR", "", minX_mass, maxX_mass)
        #r = hist_mass_original.Fit("M1Bkg", "MSR", "", minX_mass, maxX_mass)
        fitStatus = int(r)
        print "fit status : %d" % fitStatus
        if(fitStatus==0):
            stopProgram=0
            r.Print("V")
            break


    if(stopProgram==1):
        print "######################"
        print"######################"
        print "ERROR : Fit failed!!!!"
        print "######################"
        print "######################"


    # fit residuals and chi2
    hist_fit_residual_vsMass =  TH1D("hist_fit_residual_vsMass","hist_fit_residual_vsMass",number_of_variableWidth_bins,v_massBins)
    hist_fit_residual = TH1D("hist_fit_residual","hist_fit_residual",10,-5,5)
    NumberOfVarBins = 0
    NumberOfObservations_VarBin = 0
    NumberOfObservations_VarBin_5entries = 0
    chi2_VarBin = 0.
    chi2_VarBin_notNorm = 0.
    chi2_VarBin_5entries = 0.
    chi2_VarBin_zeroes = 0.

    for bin in range (1,number_of_variableWidth_bins):

        if( hist_mass.GetXaxis().GetBinLowEdge(bin)>=minX_mass and hist_mass.GetXaxis().GetBinUpEdge(bin)<=maxX_mass ):
            NumberOfVarBins += 1
            #print "bin content = " + str(hist_mass.GetBinContent(bin)) + "   graph y = " + str(vy[bin-1]) + "  error y low = " + str(g.GetErrorYlow(bin-1))
            data = hist_mass.GetBinContent(bin)
            err_data_low = g.GetErrorYlow(bin-1)
            err_data_high= g.GetErrorYhigh(bin-1)
            fit = M1Bkg.Integral(hist_mass.GetXaxis().GetBinLowEdge(bin) , hist_mass.GetXaxis().GetBinUpEdge(bin) )
            fit = fit / ( hist_mass.GetBinWidth(bin) * lumi )
            if(fit > data): err_tot = err_data_high
            else: err_tot = err_data_low
            fit_residual = 0.0
            if err_tot != 0.0:
                fit_residual = (data - fit) / err_tot
            err_fit_residual = 1
            ##skip bin with zero entries
            #print "res = "+str(pow( (data - fit) , 2 ) / pow( err_tot , 2 ))
            if err_tot != 0.0:
                chi2_VarBin_zeroes += pow( (data - fit) , 2 ) / pow( err_tot , 2 )
            if (hist_mass.GetBinContent(bin)>0):
                NumberOfObservations_VarBin+=1
                if err_tot != 0.0:
                    chi2_VarBin += pow( (data - fit) , 2 ) / pow( err_tot , 2 )
                chi2_VarBin_notNorm += pow( (data - fit) , 2 )

            ##skip bin with fewer than 5 entries
            if (hist_mass.GetBinContent(bin)*hist_mass.GetBinWidth(bin)*lumi>=5):
                NumberOfObservations_VarBin_5entries+=1
                if err_tot != 0.0:
                    chi2_VarBin_5entries += pow( (data - fit) , 2 ) / pow( err_tot , 2 )

            hist_fit_residual_vsMass.SetBinContent(bin,fit_residual)
            hist_fit_residual_vsMass.SetBinError(bin,err_fit_residual)
            hist_fit_residual.Fill(fit_residual)

    ndf_VarBin_5entries = NumberOfObservations_VarBin_5entries - nPar# -1
    ndf_VarBin = NumberOfObservations_VarBin - nPar# -1
    ndf_VarBin_withzeroes = NumberOfVarBins - nPar# -1
    #expected = M1Bkg.Integral(5058,10072) * lumi
    expected = M1Bkg.Integral(5058,10072)
    observed = hist_binned.Integral(hist_mass.FindBin(5058),hist_mass.FindBin(10072))
    print "============================"
    print "NumberOfObservations_VarBin: %d" %  NumberOfObservations_VarBin
    print "chi2 / dof: %f / %d" % (chi2_VarBin,ndf_VarBin)
    print "chi2 / dof (with zeroes) : %f / %d " % (chi2_VarBin_zeroes,ndf_VarBin_withzeroes)
    print "chi2 / dof (> 5 evt) : %f / %d " % (chi2_VarBin_5entries,ndf_VarBin_5entries)
    #print "ndf_VarBin with 5entries: %d" % ndf_VarBin_5entries
    #print "chi2_VarBin_5entries: %f" % chi2_VarBin_5entries
    #print "chi2_VarBin_notNorm: %f" % chi2_VarBin_notNorm
    print "expected events > 5 TeV (fit) : %f" % expected
    print "observed events > 5 TeV  : %f" % observed
    print "============================"
    return [chi2_VarBin_notNorm,ndf_VarBin,chi2_VarBin,M1Bkg,hist_fit_residual_vsMass]




def DrawFit(g, M1Bkg, hist_fit_residual_vsMass, FunctionType, nPar,
            fileNameSuffix, chi2, ndof):
    ##rescale fit function
    M1Bkg_xsec = TF1("M1Bkg_xsec","( [0]*TMath::Power(1-x/13000,[1]) ) / ( TMath::Power(x/13000,[2]+[3]*log(x/13000)) )*(0.5*(1.0 + TMath::Erf((x - 507.1)/94.2)))", minX_mass, maxX_mass)
    M1Bkg_xsec.SetParameter(0,M1Bkg.GetParameter(0)/lumi)
    M1Bkg_xsec.SetParameter(1,M1Bkg.GetParameter(1))
    M1Bkg_xsec.SetParameter(2,M1Bkg.GetParameter(2))
    M1Bkg_xsec.SetParameter(3,M1Bkg.GetParameter(3))
    M1Bkg_xsec.SetLineColor(2)
    M1Bkg_xsec.SetLineWidth(2)

    #  //### Draw plots
    W = 600
    H = 700
    H_ref = 700
    W_ref = 600
    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.04*W_ref

    c = TCanvas("c","DijetMass cross section with Fit and QCD MC",W,H)
    c.GetWindowHeight()
    c.GetWindowWidth()
    c.SetLogy()
    c.Divide(1,2,0,0,0)


    #------------ pad 1  ----------------
    c.cd(1)
    p11_1 = c.GetPad(1)
    p11_1.SetPad(0.01,0.26,0.99,0.98)
    p11_1.SetLogy()
    p11_1.SetRightMargin(0.05)
    p11_1.SetTopMargin(0.05)
    p11_1.SetFillColor(0)
    p11_1.SetBorderMode(0)
    p11_1.SetFrameFillStyle(0)
    p11_1.SetFrameBorderMode(0)

    #Pave text
    #pave_fit = TPaveText(0.1558691,0.30735043,0.3750171,0.4070085,"NDC")
    pave_fit = TPaveText(0.2358691,0.04035043,0.5050171,0.1870085,"NDC")

    pave_fit.AddText("#chi^{{2}}/ndof = {0:.4g}/{1}".format(chi2, ndof))
    pave_fit.AddText("Wide Jets")
    pave_fit.AddText("|#eta| < 2.5, |#Delta#eta_{jj}| < 1.3")
    #pave_fit.AddText("M_{jj} > 1.2 TeV")
    pave_fit.SetFillColor(0)
    pave_fit.SetLineColor(0)
    pave_fit.SetFillStyle(0)
    pave_fit.SetBorderSize(0)
    pave_fit.SetTextFont(42)
    pave_fit.SetTextSize(0.040)
    pave_fit.SetTextAlign(12)


    vFrame = p11_1.DrawFrame(minX_mass, 0.005, maxX_mass, 100.0)

    vFrame.SetTitle("")
    vFrame.SetXTitle("Dijet Mass [GeV]")
    vFrame.SetYTitle("d#sigma / dm_{jj}   [pb / GeV]")
    vFrame.GetXaxis().SetTitleSize(0.06)
    vFrame.GetXaxis().SetTitleOffset(0.95)
    vFrame.GetXaxis().SetLabelSize(0.05)
    vFrame.GetYaxis().SetTitleSize(0.06)
    vFrame.GetYaxis().SetLabelSize(0.05)

    g.GetXaxis().SetNdivisions(405)
    g.SetMarkerSize(0.9)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.9)
    g.SetMarkerStyle(20)
    g.Draw("pe0 same")
    M1Bkg_xsec.Draw("same")
    g.Draw("pe0 same")

    leg =  TLegend(0.4564991,0.62,0.9303575,0.90)
    leg.SetTextSize(0.038)
    leg.SetLineColor(0)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetMargin(0.35)
    leg.AddEntry(hist_mass,"Data" ,"PLE");
    leg.AddEntry(M1Bkg_xsec,"Fit","L");
    leg.Draw("same")
    pave_fit.Draw("same")

    # writing the lumi information and the CMS "logo"
    #  CMS_lumi( p11_1, iPeriod, iPos );
    #redraw axis
    p11_1.RedrawAxis()
    p11_1.Update()
    p11_1.GetFrame().Draw()
    #draw the lumi text on the canvas
    CMS_lumi.CMS_lumi(p11_1, iPeriod, iPos)

    #--- Next PAD

    c.cd(2)
    p11_2 = c.GetPad(2)
    p11_2.SetPad(0.01,0.02,0.99,0.27)
    p11_2.SetBottomMargin(0.35)
    p11_2.SetRightMargin(0.05)
    p11_2.SetGridx()
    p11_2.SetGridy()

    vFrame2 = p11_2.DrawFrame(p11_1.GetUxmin(), -2.2, p11_1.GetUxmax(), 2.2)

    vFrame2.SetTitle("")
    vFrame2.SetXTitle("Dijet Mass [GeV]")
    vFrame2.GetXaxis().SetTitleSize(0.06)
    vFrame2.SetYTitle("#frac{(Data-Fit)}{#sigma_{Data}}")
    vFrame2.GetYaxis().SetTitleSize(0.15)
    vFrame2.GetYaxis().SetTitleOffset(0.40)
    vFrame2.GetYaxis().SetLabelSize(0.09)
    vFrame2.GetXaxis().SetTitleSize(0.15)
    vFrame2.GetXaxis().SetTitleOffset(0.90)
    vFrame2.GetXaxis().SetLabelSize(0.12)

    hist_fit_residual_vsMass.GetXaxis().SetNdivisions(405)
    hist_fit_residual_vsMass.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
    hist_fit_residual_vsMass.GetYaxis().SetRangeUser(-2.5,2.5)
    hist_fit_residual_vsMass.SetLineWidth(0)
    hist_fit_residual_vsMass.SetFillColor(2)
    hist_fit_residual_vsMass.SetLineColor(1)
    hist_fit_residual_vsMass.Draw("SAMEHIST")

    if plotSig:
        hist_sig_significance.Draw("c same")
        hist_sig_significance2.Draw("c same")
        hist_sig_significance3.Draw("c same")

    line = TLine(minX_mass,0,maxX_mass,0)
    line.Draw("")
    p11_2.RedrawAxis()
    line2=TLine()
    line2.DrawLine(p11_2.GetUxmin(), p11_2.GetUymax(), p11_2.GetUxmax(), p11_2.GetUymax())
    line2.DrawLine(p11_2.GetUxmax(), p11_2.GetUymin(), p11_2.GetUxmax(), p11_2.GetUymax())

    ### Output files

    output_root_file = "dijetFitResults_FuncType%d_nParFit%d_%s.root" % (FunctionType,nPar,fileNameSuffix)

    f_output = TFile(output_root_file,"RECREATE")
    f_output.cd()
    g.Write()
    hist_mass.Write()
    c.Write()
    f_output.Close()
    c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.png" %(FunctionType,nPar,fileNameSuffix)
    c.SaveAs(outputDir+"/"+c_fileName)
    c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.pdf" %(FunctionType,nPar,fileNameSuffix)
    c.SaveAs(outputDir+"/"+c_fileName)
    c_fileName = "fitAndResiduals_FuncType%d_nParFit%d_%s.root" %(FunctionType,nPar,fileNameSuffix)
    c.SaveAs(outputDir+"/"+c_fileName)


def histSignificance(M1Bkg,hist_mass,h_sig):
    nbins_sig = h_sig.GetNbinsX()
    massBinsSig_list=[]
    for i in range(0,nbins_sig+1):
        massBinsSig_list.append(h_sig.GetXaxis().GetBinLowEdge(i+1))

    massBinsSig = array("d",massBinsSig_list)
    ##histogram signal significance
    hist_sig_significance = TH1D("hist_sig_significance","hist_sig_significance",nbins_sig,massBinsSig)

    for  i in range(1,number_of_variableWidth_bins+1):
        if( hist_mass.GetXaxis().GetBinLowEdge(i)>=minX_mass and hist_mass.GetXaxis().GetBinUpEdge(i)<=maxX_mass ):
            fit = M1Bkg.Integral(hist_mass.GetXaxis().GetBinLowEdge(i), hist_mass.GetXaxis().GetBinUpEdge(i) )
            significance = 0

            hist_sig_significance.SetBinContent(i,significance);

    return hist_sig_significance



#----- keep the GUI alive ------------
if __name__ == '__main__':
    main()
    # rep = ''
    # while not rep in ['q','Q']:
    #   rep = raw_input('enter "q" to quit: ')
    #   if 1 < len(rep):
    #     rep = rep[0]
