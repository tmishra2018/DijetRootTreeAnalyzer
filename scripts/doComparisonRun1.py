#!usr/bin/python

from setTDRStyle import setTDRStyle
import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi

import math
import optparse

CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0

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
parser.add_option("--input1",action="store",type="string",dest="input1",default='DijetMassSpectrum_CMS_13TeV_36pb-1.txt')
parser.add_option("--input2",action="store",type="string",dest="input2",default='DijetMassSpectrum_CMS_8TeV_19p7fb-1.txt')
parser.add_option("--name1",action="store",type="string",dest="name1",default='13 TeV 36 pb^{-1}')
parser.add_option("--name2",action="store",type="string",dest="name2",default='8 TeV 19.7 fb^{-1}')
parser.add_option("--outdir", action="store", type="string", dest="outdir", default="./")

(options, args) = parser.parse_args()

name1 = options.name1
name2 = options.name2
input1 = options.input1
input2 = options.input2
outdir = options.outdir

gROOT.Reset()
setTDRStyle()
gROOT.ForceStyle()
gROOT.SetStyle('tdrStyle')
#####################

minX_mass = 1118
maxX_mass = 5058 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
    
massBins = array("d",massBins_list)

ins1 = open(input1)
ins2 = open(input2)
ins_qcd1 = open("qcd_mass_cteq6ll_lhc13.out")
ins_qcd2 = open("qcd_mass_cteq6ll_lhc8.out")


#---- read the inputfiles -----------------

x1=[]
x1_el=[]
x1_eh=[]
y1=[]
y1_el=[]
y1_eh=[]

x2=[]
x2_el=[]
x2_eh=[]
y2=[]
y2_el=[]
y2_eh=[]

y1_qcd=[]
y1_qcd_el=[]
y1_qcd_eh=[]
y2_qcd=[]
y2_qcd_el=[]
y2_qcd_eh=[]
ratio_qcd=[]
ratio_qcd_el=[]
ratio_qcd_eh=[]


ratio=[]
err0=[]
err1_over_xsec1 = []
err2_over_xsec2 = []
err_ratio=[]

ii=0
for line in ins1:
  line.rsplit()
  if not line.startswith("#"):
    bin_low = float(line.split()[0])
    bin_high = float(line.split()[1])
    bin_content = float(line.split()[2])
    bin_err_low = float(line.split()[3])
    bin_err_high = float(line.split()[4])
    x1.append((bin_low+bin_high)/2)
    x1_el.append((bin_high-bin_low)/2)
    x1_eh.append((bin_high-bin_low)/2)
    y1.append(bin_content)
    y1_el.append(bin_err_low)
    y1_eh.append(bin_err_high)
    if not bin_content==0:
      err1_over_xsec1.append( (bin_err_low+bin_err_high)/(2*bin_content))
    else:
      err1_over_xsec1.append(0)
#print "%s  %s  %s  %s  %s" %(bin_low, bin_high, bin_content, bin_err_low, bin_err_high)
    ii+=1

ii=0
for line in ins2:
  line.rsplit()
  if not line.startswith("#"):
    bin_low = float(line.split()[0])
    bin_high = float(line.split()[1])
    bin_content = float(line.split()[2])
    bin_err_low = float(line.split()[3])
    bin_err_high = float(line.split()[4])
    x2.append((bin_low+bin_high)/2)
    x2_el.append((bin_high-bin_low)/2)
    x2_eh.append((bin_high-bin_low)/2)
    y2.append(bin_content)
    y2_el.append(bin_err_low)
    y2_eh.append(bin_err_high)
    if not bin_content==0:
      ratio.append( y1[ii]/y2[ii])
      err2_over_xsec2.append( (bin_err_low+bin_err_high)/(2*bin_content))
      err_ratio.append(ratio[ii]*math.sqrt(math.pow(err1_over_xsec1[ii],2)+math.pow(err2_over_xsec2[ii],2)))
    else: 
      ratio.append( 0)   
      err2_over_xsec2.append(0)
      err_ratio.append(0)
    print "%f  %f  %f  %f  %f" %(x1[ii],x2[ii],y1[ii],y2[ii],ratio[ii])
    err0.append(0) 
    ii+=1

##------ read files qcd cross sec
ii=0
for line in ins_qcd1:
  line.rsplit()
  if not line.startswith("#"):
    bin_content = float(line.split()[1])
    y1_qcd.append(bin_content)
    ii+=1

ii=0
for line in ins_qcd2:
  line.rsplit()
  if not line.startswith("#"):
    bin_content = float(line.split()[1])
    y2_qcd.append(bin_content)
    if ii==0:
      ratio_qcd_0=y1_qcd[0] / y2_qcd[0]  
      ratio_qcd.append(ratio_qcd_0)
      print "%f  %f" %  (ratio_qcd_0,ratio[0])
      #ratio_qcd.append(ratio[0])
    else: 
      #norm = ratio_qcd_0/ratio[0]
      norm = 1
      ratio_qcd.append(y1_qcd[ii] / y2_qcd[ii] / norm )
    ii+=1



v_x1 = array("d",x1)
v_x1_el = array("d",x1_el)
v_x1_eh = array("d",x1_eh)
v_y1 = array("d",y1)
v_y1_el = array("d",y1_el)
v_y1_eh =array("d",y1_eh)

v_x2 = array("d",x2)
v_x2_el = array("d",x2_el)
v_x2_eh = array("d",x2_eh)
v_y2 = array("d",y2)
v_y2_el = array("d",y2_el)
v_y2_eh =array("d",y2_eh)
##assuming that the number of bins correspond in g1 and g2
v_ratio = array("d",ratio)
v_err0= array("d",err0)
v_err_ratio= array("d",err_ratio)

v_ratio_qcd = array("d",ratio_qcd)


print "%d" % len(x1)
print v_x1
print v_y1
print v_x1_el
print v_x1_eh
print v_y1_el
print v_y1_eh
g1 = TGraphAsymmErrors(len(x1),v_x1,v_y1,v_x1_el,v_x1_eh,v_y1_el,v_y1_eh)
g2 = TGraphAsymmErrors(len(x2),v_x2,v_y2,v_x2_el,v_x2_eh,v_y2_el,v_y2_eh)
g_ratio = TGraphAsymmErrors(len(x1),v_x1,v_ratio,v_x1_el,v_x1_eh,v_err_ratio,v_err_ratio)
g_ratio_qcd = TGraphAsymmErrors(len(x1),v_x1,v_ratio_qcd,v_x1_el,v_x1_eh,v_err0,v_err0)

g_ratio_qcd.SetMarkerSize(0.9)
g_ratio_qcd.SetMarkerStyle(24)
g_ratio_qcd.SetMarkerColor(2)
g_ratio_qcd.SetLineColor(2)
g_ratio.SetMarkerSize(0.9)
g1.SetMarkerSize(0.9)
g2.SetMarkerSize(0.9)
g2.SetMarkerStyle(24)

g1.SetName("g1")
g2.SetName("g2")
g_ratio.SetName("ratio_1_2")
g_ratio_qcd.SetName("ratioqcd_1_2")
##----- Drawing  and save on file -----------------------

outFile = TFile(outdir+"comparison_"+name1+"_"+name2+".root", "recreate")
outFile.cd()
g1.Write()
g2.Write()
g_ratio.Write()
g_ratio_qcd.Write()

c1 = TCanvas('c1','',800,600) 
leg1 = TLegend(0.6, 0.7, 0.85, 0.85)
leg1.SetLineColor(0)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetBorderSize(0)
leg1.AddEntry(g_ratio, "data ratio 13/8", "p")
leg1.AddEntry(g2, "qcd ratio 13/8", "p")

c1.cd()
g_ratio.GetXaxis().SetRangeUser(1118,5058)
g_ratio.GetYaxis().SetRangeUser(0., 25.)
g_ratio.GetYaxis().SetNdivisions(8, kTRUE)
g_ratio.GetYaxis().SetTitle("xsec ratio 13TeV / 8TeV")
g_ratio.Draw("AP")
g_ratio_qcd.Draw("PL SAME")
leg1.Draw()
##draw the lumi text on the canvas
CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
c1.SaveAs(outdir+'/ratio_run2_run1_qcd.png')

c = TCanvas('c','',600,600)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(g1, name1, "p")
leg.AddEntry(g2, name2, "p")

c.cd()

##----- pad 1 -----------
pad1 = TPad("pad1", "pad1",0.01,0.13,0.99,0.99)  
pad1.SetRightMargin(0.1)

pad1.Draw()
pad1.cd()

pad1.SetLogy()

xtitle = "Dijet mass (GeV)"
g1.GetXaxis().SetTitle(xtitle)
title_y = "d#sigma / dm (pb / GeV)"
g1.GetYaxis().SetTitleFont(42)
g1.GetYaxis().SetTitleSize(0.04)
g1.GetYaxis().SetTitle(title_y)
g1.GetYaxis().SetRangeUser(0.0000004,10)
g1.GetXaxis().SetRangeUser(1118,5058)
g2.GetXaxis().SetRangeUser(1118,5058)

g1.Draw("APE0")
g2.Draw("PE0 SAME")
leg.Draw()
gPad.RedrawAxis()
##draw the lumi text on the canvas
CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

c.cd()
##-------pad 2------
pad2 = TPad("pad2", "pad2",0.01,0.001,0.99,0.25)
pad2.SetGrid()
	      
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.4)
pad2.SetRightMargin(0.1)
pad2.Draw()	       
pad2.cd()

g_ratio.GetXaxis().SetRangeUser(1118,5058)
g_ratio.Draw("AP")
g_ratio.GetYaxis().SetRangeUser(0., 15.)
g_ratio.GetYaxis().SetNdivisions(4, kTRUE)
g_ratio.GetYaxis().SetTitleFont(42)
g_ratio.GetYaxis().SetTitle("#frac{"+name1+"}{"+name2+"}")
g_ratio.GetXaxis().SetTitleSize(0.2)
g_ratio.GetXaxis().SetLabelSize(0.16)
g_ratio.GetYaxis().SetLabelSize(0.16)
g_ratio.GetYaxis().SetTitleSize(0.15)
g_ratio.GetXaxis().SetTitleOffset(0.8)

#RedrawAxis
pad1.cd()
gPad.RedrawAxis()
pad2.cd()
gPad.RedrawAxis()

c.Write()   
c.SaveAs(outdir+'/ratio_run2_run1_err.png')
