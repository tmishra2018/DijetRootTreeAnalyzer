#!usr/bin/python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi
import optparse
from setTDRStyle import setTDRStyle
import math

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

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')


#######################################################

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList_mc",action="store",type="string",default="list_mc.txt",dest="inputList_mc")
parser.add_option("--inputList_data",action="store",type="string",default="list_data.txt",dest="inputList_data")
parser.add_option("--lumi",action="store",type="float",default="1000.",dest="lumi")
parser.add_option("--tag",action="store",type="string",default="DCS",dest="tag")

(options, args) = parser.parse_args()

outputDir = options.outputDir
inputList_mc = options.inputList_mc
inputList_data = options.inputList_data
lumi = options.lumi
tag = options.tag
var = 'mjj'
os.system("mkdir -p "+outputDir)
#############################

CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = str(options.lumi)+" pb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0
#######################
bins = 13999
xmin = 1
xmax= 14000
#minX_mass = 526.
minX_mass = 1181.
#maxX_mass = 5877. 
maxX_mass = 7320. 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
    
massBins = array("d",massBins_list)
hist_allCuts      = []
hist_allCuts_2    = []
#LUMI      = lumi
#PATH      = inputDir

#---- read the list -----------------
lines = [line.strip() for line in open(inputList_mc)]

#---- split sample name and xsec (mc list)
fileNames = []
xsecs = []
ii = 0
for line in lines:
  parts = line.split()
  fileNames.append(parts[0])
  xsecs.append(parts[1])
  print ("dataset : %s    xsec : %s" % (fileNames[ii], xsecs[ii]))
  ii+=1

fileNames_data = []
lines_data = [line.strip() for line in open(inputList_data)]
for line in lines_data:
  fileNames_data.append(line)

#------------  MC  ---------------
i_f = 0
for f in fileNames:
  inf = TFile.Open(f)
  print inf.GetName()
  
  Nev = inf.Get('DijetFilter/EventCount/EventCounter').GetBinContent(1)
  print ('processed events: %s' % Nev)
  wt = 1.0
  #if i_f < 3:

  h_allCuts = TH1F("h_allCuts", "", bins, xmin, xmax)
  h_allCuts.Sumw2()
  h_allCuts_2 = TH1F("h_allCuts_2", "", bins, xmin, xmax)
  h_allCuts_2.Sumw2()
  tree = inf.Get('rootTupleTree/tree')
  #standard
  tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > '+str(minX_mass))
  tree.Project(h_allCuts_2.GetName(), var,'deltaETAjj<0.8 && mjj > '+str(minX_mass))
  
  Npassed = h_allCuts.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  print('(not using efficiency in the weight)')
  #if not (i_f == 9):
  wt = options.lumi*float(xsecs[i_f])/Nev
  print('weight : %f' % wt)
  h_allCuts.Scale(wt)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetLineWidth(2)
  h_allCuts.SetLineColor(kBlue-9)
  h_allCuts.SetMarkerColor(kBlue-9)
  hist_allCuts.append(h_allCuts)
  h_allCuts_2.Scale(wt)
  h_allCuts_2.SetDirectory(0)
  h_allCuts_2.SetLineWidth(2)
  h_allCuts_2.SetLineColor(kBlue-9)
  h_allCuts_2.SetMarkerColor(kBlue-9)
  hist_allCuts_2.append(h_allCuts_2)
  h_allCuts.Print()
  h_allCuts_2.Print()
  print "entries (deta < 1.3): %d   (deta < 1.8): %d" % (h_allCuts.GetEntries(),h_allCuts_2.GetEntries())
  print "integral (deta <1.3) : %f   (deta < 1.8): %f" % (h_allCuts.Integral(),h_allCuts_2.Integral())
  i_f += 1
 
print hist_allCuts  
print hist_allCuts_2  

#--------------------------		   
hist_allCutsQCD=hist_allCuts[0].Clone("hist_allCutsQCD") 
hist_allCutsQCD_2=hist_allCuts_2[0].Clone("hist_allCutsQCD_2") 
for i in range(1,9):
  hist_allCutsQCD.Add(hist_allCuts[i]) 
  hist_allCutsQCD_2.Add(hist_allCuts_2[i]) 
print "------------------------------------------------"
#print "tot entries (deta < 1.3): %d   (deta < 1.8): %d" % (hist_allCutsQCD.Add.GetEntries(),hist_allCutsQCD.Add_2.GetEntries())
print "QCD tot integral (deta <1.3) : %f   (deta < 1.8): %f" % (hist_allCutsQCD.Integral(),hist_allCutsQCD_2.Integral())

#-------  signal -------
hist_allCutsSig = hist_allCuts[9]
hist_allCutsSig_2 = hist_allCuts_2[9]

hist_s_plus_b = hist_allCutsQCD.Clone("hist_s_plus_b")
hist_s_plus_b_2 = hist_allCutsQCD_2.Clone("hist_s_plus_b_2")
hist_s_plus_b.Add(hist_allCutsSig)
hist_s_plus_b_2.Add(hist_allCutsSig_2)
hist_s_plus_b.SetLineColor(2)
hist_s_plus_b.SetLineWidth(2)
hist_s_plus_b_2.SetLineColor(2)
hist_s_plus_b_2.SetLineWidth(2)

#--------- data  --------------
chain = TChain("rootTupleTree/tree")
for i in range(0,len(fileNames_data)):
  chain.Add(fileNames_data[i])
  print fileNames_data[i]

h_dat = TH1F("h_dat", "", bins, xmin, xmax)
h_dat.Sumw2()
h_dat_2 = TH1F("h_dat_2", "", bins, xmin, xmax)
h_dat_2.Sumw2()
#passJSON
#chain.Project("h_dat",var,'deltaETAjj<1.3 && mjj > '+str(minX_mass)+' && PassJSON==1')
#chain.Project("h_dat_2",var,'deltaETAjj<0.8 && mjj > '+str(minX_mass)+' && PassJSON==1')
#standard
chain.Project("h_dat",var,'deltaETAjj<1.3 && mjj > '+str(minX_mass))
chain.Project("h_dat_2",var,'deltaETAjj<0.8 && mjj > '+str(minX_mass))


#h_dat = hist_allCuts[9].Clone()
#h_dat.SetName("h_dat")
#h_dat.SetLineStyle(2)
#h_dat.SetLineWidth(2)
h_dat.SetMarkerColor(kBlack)
h_dat.SetMarkerSize(0.9)
h_dat.SetLineColor(kBlack)
h_dat_2.SetMarkerColor(kBlack)
h_dat_2.SetLineColor(kBlack)

                   
#--- rebin ---------
hist_allCutsQCD_rebin = hist_allCutsQCD.Rebin(len(massBins_list)-1,var+"_rebin",massBins)
hist_allCutsQCD_2_rebin = hist_allCutsQCD_2.Rebin(len(massBins_list)-1,var+"_rebin",massBins)
hist_s_plus_b_rebin = hist_s_plus_b.Rebin(len(massBins_list)-1,var+"_rebin",massBins)
hist_s_plus_b_2_rebin = hist_s_plus_b_2.Rebin(len(massBins_list)-1,var+"_rebin",massBins)

h_dat_rebin = h_dat.Rebin(len(massBins_list)-1,var+"_data_rebin",massBins)
h_dat_2_rebin = h_dat_2.Rebin(len(massBins_list)-1,var+"_data_rebin",massBins)


#---- divide  ------
h_mjjQCD = hist_allCutsQCD_2_rebin.Clone("h_mjjQCD") 
h_mjjQCD.Divide(hist_allCutsQCD_rebin)
h_mjjSB = hist_s_plus_b_2_rebin.Clone("h_mjjSB") 
h_mjjSB.Divide(hist_s_plus_b_rebin)

h_mjjDat = h_dat_2_rebin.Clone("h_mjjDat")
h_mjjDat.Divide(h_dat_rebin)

h_mjjQCD.Print()
h_mjjDat.Print()
##---- set bin error ----
#for i in range(1,h_mjjQCD.GetNbinsX()):
#  f_QCD = h_mjjQCD.GetBinContent(i)
#  N_QCD = hist_allCutsQCD_rebin.GetBinContent(i)
#  h_mjjQCD.SetBinError(i,TMath.Sqrt( f_QCD*(1-f_QCD)/ N_QCD ) 
#  f_Dat = h_mjjDat.GetBinContent(i)
#  N_Dat = hist_allCutsDat_rebin.GetBinContent(i)
#  h_mjjDat.SetBinError(i,TMath.Sqrt( f_Dat*(1-f_Dat)/ N_Dat )


h_mjj_vector_1 = [hist_allCutsQCD_rebin, h_dat_rebin]
h_mjj_vector_2 = [hist_allCutsQCD_2_rebin, h_dat_2_rebin]

#- TGraph -----
g_eff = []
MODE=1
##### efficiency and error for trigger curves
for ii in range(0,2):
  x = []
  y = []
  exl = []
  exh = []
  eyl = []
  eyh = []
  scale = 1.0
  for i in range(0,len(massBins)-1):
    N1 = h_mjj_vector_1[ii].GetBinContent(i+1)
    N2 = h_mjj_vector_2[ii].GetBinContent(i+1)
    p  = 0
    eU = 0
    eL = 0
    if (N1 > 0):
      p = N2/N1
      n = N1+N2
      w = N2/n
      if (MODE==1): ## Wilson for binomial
	scale = 1.0 # makes sense only for the unprescaled trigger 
        d = math.sqrt(p*(1-p)/N1+0.25/(N1*N1))
        eU = (p+0.5/N1+d)/(1+1/N1)-p
        eL = p-(p+0.5/N1-d)/(1+1/N1)
         
      else: ## Wilson for Poisson ratio
        d  = math.sqrt(w*(1-w)/n+0.25/(n*n))
        UB = (w+0.5/n+d)/(1+1/n)
        LB = (w+0.5/n-d)/(1+1/n)
        eU = UB/(1-UB)-p
        eL = p-LB/(1-LB)     

    x.append( h_mjj_vector_1[0].GetBinCenter(i+1))
    y.append( p*scale)
    exl.append( h_mjj_vector_1[0].GetBinWidth(i+1)/2)
    exh.append( h_mjj_vector_1[0].GetBinWidth(i+1)/2)      
    eyl.append( eL*scale)
    eyh.append( eU*scale)
       
  vx = array("f",x)
  vy = array("f",y)
  vexl = array("f",exl)
  vexh = array("f",exh)
  veyl = array("f",eyl)
  veyh = array("f",eyh)

  g_eff.append(TGraphAsymmErrors(len(massBins)-1,vx,vy,vexl,vexh,veyl,veyh))  


ytitle = "d#sigma_{#eta<0.8} / d#sigma_{#eta<1.3}"
#ytitle = "d#sigma_{#eta<1.3} / d#sigma_{#eta>1.3}"
h_mjjDat.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
h_mjjDat.GetXaxis().SetTitle("Dijet Mass [GeV]")
h_mjjQCD.SetLineColor(kBlue-9)
h_mjjQCD.SetLineWidth(2)
h_mjjQCD.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
h_mjjQCD.GetYaxis().SetRangeUser(0,1.5)
h_mjjQCD.GetXaxis().SetTitle("Dijet Mass [GeV]")
h_mjjQCD.GetYaxis().SetTitle(ytitle)
g_eff[0].SetName("g_QCD")
g_eff[0].SetLineColor(kBlue-9)
g_eff[0].SetLineWidth(2)
g_eff[0].GetXaxis().SetRangeUser(minX_mass,maxX_mass)
g_eff[0].GetXaxis().SetTitle("Dijet Mass [GeV]")
g_eff[0].GetYaxis().SetTitle(ytitle)
g_eff[1].SetName("g_data")
g_eff[1].SetMarkerColor(1)
g_eff[1].SetMarkerSize(0.9)
g_eff[1].GetXaxis().SetRangeUser(minX_mass,maxX_mass)
g_eff[1].GetXaxis().SetTitle("Dijet Mass [GeV]")
g_eff[1].GetYaxis().SetTitle(ytitle)

#----- Draw -------
file_out = TFile(outputDir+"/deltaeta_check_"+tag+".root","recreate")
c = TCanvas('c','c',600,600)
c.cd()
leg = TLegend(0.45, 0.7, 0.85, 0.9)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.AddEntry(h_mjjQCD, "QCD", "l")
leg.AddEntry(h_mjjSB,  "QCD + q* (4 TeV, 0.06 pb^{-1})", "l")
leg.AddEntry(g_eff[1], "data", "p")

h_mjjQCD.Draw("hist")
#g_eff[0].Draw("al")
h_mjjSB.Draw("hist same")
g_eff[1].Draw("pe0 same")
leg.Draw()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c, iPeriod, iPos)

c.SaveAs(outputDir+"/deltaeta_check_"+tag+".png")
c.SaveAs(outputDir+"/deltaeta_check_"+tag+".pdf")
c.Write()
h_mjjQCD.Write()
g_eff[0].Write()
g_eff[1].Write()

#-----------  with ratio -----------
c2 = TCanvas('c2','c2',600,650)
c2.cd()
#----- pad 1 -----------
pad1 = TPad("pad1", "pad1",0,0.15,1,1)  
#pad1.SetRightMargin(0.1)

pad1.Draw()
pad1.Clear()
pad1.cd()
h_mjjQCD.Draw("hist")
h_mjjSB.Draw("hist same")
g_eff[1].Draw("pe0 same")
leg.Draw()
#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)
gPad.RedrawAxis()

#-------pad 2------
pad2 = TPad("pad2", "pad2",0,0,1,0.15)
pad2.SetGrid()
	      
pad2.SetTopMargin(0)
#pad2.SetBottomMargin(0.4)
#pad2.SetRightMargin(0.1)
pad2.Draw()	       
pad2.cd()

ratio = h_mjjDat.Clone("ratio")
ratio.GetXaxis().SetRangeUser(minX_mass,maxX_mass)
ratio.Divide(h_mjjQCD)
for i in range(1,len(massBins)):
  error = 0 
  data = h_mjjDat.GetBinContent(i)
  mc = h_mjjQCD.GetBinContent(i)
  if (data > 0 and data > mc):
    error = g_eff[1].GetErrorYhigh(i-1)
  else:
    error = g_eff[1].GetErrorYlow(i-1)
  if mc>0:
    ratio.SetBinError(i,error/mc)
  else:
    ratio.SetBinError(i,0)


ratio.SetFillColor(0)
ratio.SetLineColor(kBlack)
ratio.SetMarkerColor(kBlack)
ratio.SetMarkerSize(0.9)
ratio.GetXaxis().SetRangeUser(minX_mass, maxX_mass)
ratio.GetYaxis().SetRangeUser(0.8, 1.2)
ratio.GetYaxis().SetNdivisions(405, kTRUE)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitle("data / MC")
ratio.GetXaxis().SetTitle("Dijet Mass [GeV]")
ratio.GetXaxis().SetTitleSize(0.2)
ratio.GetXaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetTitleSize(0.15)
ratio.GetYaxis().SetTitleOffset(0.5)
#ratio.GetXaxis().SetTitleOffset(0.8)
ratio.Draw("p")

#RedrawAxis
#pad1.cd()
#gPad.RedrawAxis()
pad2.cd()
gPad.RedrawAxis()

c2.SaveAs(outputDir+"/deltaeta_check_withRatio_"+tag+".png")
c2.SaveAs(outputDir+"/deltaeta_check_withRatio_"+tag+".pdf")
c2.Write()
ratio.Write()

file_out.Close()
