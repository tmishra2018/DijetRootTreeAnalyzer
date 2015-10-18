#!usr/bin/python

from setTDRStyle import setTDRStyle
import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi

import optparse


gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)
gStyle.SetTitleFont(42, "XYZ")
gStyle.SetTitleSize(0.06, "XYZ")
gStyle.SetLabelFont(42, "XYZ")
gStyle.SetLabelSize(0.05, "XYZ")
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetPadLeftMargin(0.15)
gStyle.SetPadRightMargin(0.05)
gStyle.SetPadTopMargin(0.05)
gStyle.SetPadBottomMargin(0.15)
gROOT.ForceStyle()

#######################################################

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--var",action="store",type="string",dest="var",default='ptHat')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--bins",action="store",type="int",dest="bins",default=100)
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList_1",action="store",type="string",default="list1.txt",dest="inputList_1")
parser.add_option("--inputList_2",action="store",type="string",default="list2.txt",dest="inputList_2")
parser.add_option("--name1",action="store",type="string",default="name1",dest="name1")
parser.add_option("--name2",action="store",type="string",default="name2",dest="name2")
parser.add_option("--lumi",action="store",type="float",default="1000.",dest="lumi")

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
bins = options.bins
xtitle = options.xtitle
rebin = options.rebin
logy = options.logy
outputDir = options.outputDir
inputList_1 = options.inputList_1
inputList_2 = options.inputList_2
name1 = options.name1
name2 = options.name2
lumi = options.lumi

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')
#####################
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "(13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

minX_mass = 1181
maxX_mass = 5877 
#maxX_mass = 8447 
#minX_mass = 119
#maxX_mass = 4010 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
    
massBins = array("d",massBins_list)
hist_allCuts      = []
hist_allCuts_2 = []

LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines1 = [line1.strip() for line1 in open(inputList_1)]
lines2 = [line2.strip() for line2 in open(inputList_2)]

#---- split sample name and xsec
fileNames = []
fileNames2 = []
xsecs = []
ii = 0
for line in lines1:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset 1 : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  ii+=1

ii = 0
for line in lines2:
  parts2 = line.split()
  fileNames2.append(parts2[0])
  xsecs.append(parts[1])
  print ("dataset 2 : %s    xsec : %s" % (fileNames2[ii], xsecs[ii]))


  ii+=1

dataset1 = []
dataset2 = []

i_f = 0
for f in fileNames:
  inf = TFile.Open(f)
  print inf.GetName()
  dataset1.append ( os.path.basename(fileNames[i_f]) )
  dataset1[i_f] = dataset1[i_f].split("rootfile_")[1]
  dataset1[i_f] = dataset1[i_f].split("_reduced_skim.root")[0]
  dataset1[i_f] = dataset1[i_f].split("__")[0]

  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events list1 : %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts = TH1F("h_allCuts", "", bins, xmin, xmax)
  h_allCuts.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  ##standard
  tree.Project(h_allCuts.GetName(), var,' deltaETAjj<1.3 && mjj>1181')
  ##blind m > 4TeV and look high mass
  #tree.Project(h_allCuts.GetName(), var,' deltaETAjj<1.3 && mjj<4010 && mjj>1181 && pTWJ_j1>1000 && pTWJ_j2>1000')
  ##look at low mass requiring prescaled trigger
  #tree.Project(h_allCuts.GetName(), var,' deltaETAjj<1.3 && passHLT_PFHT475')

  h_allCuts.SetDirectory(0)
  #h_allCuts.SetFillStyle(3005)
  #h_allCuts.SetFillColor(kRed)
  h_allCuts.SetLineColor(kRed)
  h_allCuts.SetMarkerColor(kRed)
  h_allCuts.SetMarkerSize(0)
  hist_allCuts.append(h_allCuts)
   
  i_f += 1

#------- open the files in 2nd list -------------

i_f = 0
for f in fileNames2:
  inf = TFile.Open(f)
  print inf.GetName()
  dataset2.append ( os.path.basename(fileNames2[i_f]) )
  dataset2[i_f] = dataset2[i_f].split("rootfile_")[1]
  dataset2[i_f] = dataset2[i_f].split("_reduced_skim.root")[0]
  dataset2[i_f] = dataset2[i_f].split("__")[0]
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events list 2: %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts_2 = TH1F("h_allCuts_2", "", bins, xmin, xmax)
  h_allCuts_2.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  ##standard
  tree.Project(h_allCuts_2.GetName(), var,' deltaETAjj<1.3 && mjj>1181 && PassJSON==1')
  ##blind m > 4TeV and look high mass
  #tree.Project(h_allCuts_2.GetName(), var,' deltaETAjj<1.3 && mjj<4010 && mjj>1181 && pTWJ_j1>1000 && pTWJ_j2>1000')
  ##look at low mass requiring prescaled trigger
  #tree.Project(h_allCuts_2.GetName(), var,' deltaETAjj<1.3 && passHLT_PFHT475')
  h_allCuts_2.SetDirectory(0)
  #h_allCuts_2.SetFillStyle(3004)
  #h_allCuts_2.SetFillColor(kBlue-9)
  h_allCuts_2.SetLineColor(kBlue-9)
  h_allCuts_2.SetMarkerColor(kBlue-9)
  h_allCuts_2.SetMarkerSize(0)
  hist_allCuts_2.append(h_allCuts_2)
   
  i_f += 1

hist_allCuts_tot = hist_allCuts[0].Clone(name1)
hist_allCuts_tot_2 = hist_allCuts_2[0].Clone(name2)

for i in range(1,len(fileNames)):
  hist_allCuts_tot.Add(hist_allCuts[i])
for i in range(1,len(fileNames2)):
  hist_allCuts_tot_2.Add(hist_allCuts_2[i])



#----- Drawing  and save on file -----------------------

outFile = TFile(outputDir+"histo_signal_bkg_"+var+".root", "recreate")
outFile.cd()
hist_allCuts_tot.Write()
hist_allCuts_tot_2.Write()
for i in range(9, len(fileNames)):
  hist_allCuts[i].Write()
  hist_allCuts_2[i].Write()

can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,600,650)

leg = TLegend(0.4, 0.75, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
#leg.AddEntry(hist_allCuts_tot, name1+" JEC Summer15_50ns_V4" , "le")
#leg.AddEntry(hist_allCuts_tot_2, name2+"JEC Summer15_25ns_V3", "le")
leg.AddEntry(hist_allCuts_tot, name1 , "le")
leg.AddEntry(hist_allCuts_tot_2, name2, "le")
#leg.AddEntry(hist_allCuts_tot_2, name2+" scaled to "+name1, "le")

can_allCuts.cd()


#----- pad 1 -----------
pad1 = TPad("pad1", "pad1",0.01,0.13,1.,1.)  
pad1.SetRightMargin(0.1)

pad1.Draw()
pad1.cd()
pad1.Clear()      

#title_y = "Events (L=150 pb^{-1})"
#title_y = "Events (L=15.5 pb^{-1})"
title_y = "Events"

#hist_allCuts_tot.Reset()
hist_allCuts_tot.GetXaxis().SetRangeUser(xmin,xmax)
hist_allCuts_tot.GetXaxis().SetTitle(xtitle)
hist_allCuts_tot_2.GetXaxis().SetRangeUser(xmin,xmax)
hist_allCuts_tot_2.GetXaxis().SetTitle(xtitle)
size = hist_allCuts_tot.GetBinWidth(1)

hist_allCuts_tot.GetYaxis().SetTitle(title_y)
hist_allCuts_tot_2.GetYaxis().SetTitle(title_y)


#maximumBin = array('f',  [hist_allCuts_tot.GetBinContent(hist_allCuts_tot.GetMaximumBin()), hist_allCutsSig_1000.GetBinContent(hist_allCutsSig_1000.GetMaximumBin()), hist_allCutsSig_5000.GetBinContent(hist_allCutsSig_5000.GetMaximumBin()), hist_allCutsSig_8000.GetBinContent(hist_allCutsSig_8000.GetMaximumBin())])
#max = TMath.MaxElement(4, maximumBin)
#hist_allCuts_tot.SetMaximum(1.2*max)
#hist_allCuts_tot.SetMinimum(0.00001)
#hist_allCuts_tot.Rebin()
#hist_allCuts_tot_PHYS14_feb15.Rebin()

hist_allCuts_rebin =[]
hist_allCuts_2_rebin =[]

if (var=="mjj" or var=="Dijet_MassAK4"):
  hist_allCuts_tot_rebin = hist_allCuts_tot.Rebin(len(massBins_list)-1,name1+"_rebin",massBins)
  hist_allCuts_tot_2_rebin = hist_allCuts_tot_2.Rebin(len(massBins_list)-1,name2+"_rebin",massBins)
  hist_allCuts_tot_rebin.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
  hist_allCuts_tot_2_rebin.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
  #hist_allCuts_tot_rebin.GetYaxis().SetRangeUser(0.1,1000000)
  #hist_allCuts_tot_2_rebin.GetYaxis().SetRangeUser(0.1,1000000)
else :
  hist_allCuts_tot.Rebin(rebin)
  hist_allCuts_tot_2.Rebin(rebin)
  hist_allCuts_tot_rebin = hist_allCuts_tot.Clone(name1+"_rebin")
  hist_allCuts_tot_2_rebin = hist_allCuts_tot_2.Clone(name2+"_rebin")

#hist_allCuts_tot_rebin.Scale(hist_allCuts_tot_2_rebin.Integral()/hist_allCuts_tot_rebin.Integral())
#hist_allCuts_tot_2_rebin.Scale((23.186+22.33+19.50)/16.089)
#hist_allCuts_tot_2_rebin.Scale(68./6.5)
#hist_allCuts_tot_rebin.Scale(150./68.)
#hist_allCuts_tot_rebin.Scale(150./16.)


if ((var=="mjj" or var=="Dijet_MassAK4") and logy):
  min=0.2
  hist_allCuts_tot_2_rebin.SetMinimum(min)
max=hist_allCuts_tot_2_rebin.GetBinContent(hist_allCuts_tot_2_rebin.GetMaximumBin())
if logy:
  pad1.SetLogy(1)
  hist_allCuts_tot_2_rebin.SetMaximum(50*max)
else:
  pad1.SetLogy(0)
  hist_allCuts_tot_2_rebin.SetMaximum(2*max)

hist_allCuts_tot_2_rebin.Draw("pe")
hist_allCuts_tot_rebin.Draw("pe same")
leg.Draw()
pad1.RedrawAxis()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

can_allCuts.cd()
#-------pad 2------
pad2 = TPad("pad2", "pad2",0.01,0.001,1,0.25)
pad2.SetGrid()
	      
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.4)
pad2.SetRightMargin(0.1)
pad2.Draw()	       
pad2.cd()

ratio = hist_allCuts_tot_2_rebin.Clone("ratio")
ratio.Divide(hist_allCuts_tot_rebin)

#for i in range(ratio.GetNbinsX()):
#  if(hist_allCuts_tot_2_rebin.GetBinContent(i+1) and hist_allCuts_tot_rebin.GetBinContent(i+1)):
#    err=TMath.Sqrt(TMath.Power(((hist_allCuts_tot_2_rebin.GetBinError(i+1)/((23.186+22.33+19.50)/8.103))/hist_allCuts_tot_2_rebin.GetBinContent(i+1)),2) + TMath.Power((hist_allCuts_tot_rebin.GetBinError(i+1)/hist_allCuts_tot_rebin.GetBinContent(i+1)),2) * ratio.GetBinContent(i+1))
#    print "error ratio: "+str(err) + "   err old: "+str(ratio.GetBinError(i+1))
#    #ratio.SetBinError(i+1,err)

ratio.SetFillColor(0)
ratio.SetLineColor(kBlack)
ratio.SetMarkerSize(0.9)
ratio.SetMarkerColor(kBlack)
ratio.Draw("p")
ratio.GetYaxis().SetRangeUser(0., 2.)
ratio.GetYaxis().SetNdivisions(405, kTRUE)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitle("#frac{"+name2+"}{"+name1+"}")
ratio.GetXaxis().SetTitleSize(0.2)
ratio.GetXaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetTitleSize(0.15)
ratio.GetYaxis().SetTitleOffset(0.5)
#ratio.GetXaxis().SetTitleOffset(0.8)
ratio.Draw()

can_allCuts.Write()   
if logy:
  can_allCuts.SaveAs(outputDir+var+'_allCuts_logy.png')
  can_allCuts.SaveAs(outputDir+var+'_allCuts_logy.pdf')

else:
  
  can_allCuts.SaveAs(outputDir+var+'_allCuts.pdf')

     

can_allCuts.Close()

outFile.Close()

