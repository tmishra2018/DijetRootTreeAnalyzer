#!usr/bin/python

import sys, os, subprocess, string, re
from ROOT import *
from array import array
import CMS_lumi
import optparse
from setTDRStyle import setTDRStyle


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
parser.add_option("--var",action="store",type="string",dest="var",default='ptHat')
parser.add_option("--xmin",action="store",type="float",dest="xmin",default=1)
parser.add_option("--xmax",action="store",type="float",dest="xmax",default=1)
parser.add_option("--xtitle",action="store",type="string",dest="xtitle",default='')
parser.add_option("--bins",action="store",type="int",dest="bins",default=11111111111)
parser.add_option("--rebin",action="store",type="int",dest="rebin",default=1)
parser.add_option("--units",action="store",type="string",dest="units",default='')
parser.add_option("--logy",action="store_true",default=False,dest="logy")
parser.add_option("--outputDir",action="store",type="string",default="./",dest="outputDir")
parser.add_option("--inputList_mc",action="store",type="string",default="list_mc.txt",dest="inputList_mc")
parser.add_option("--inputList_data",action="store",type="string",default="list_data.txt",dest="inputList_data")
parser.add_option("--lumi",action="store",type="float",default="1000.",dest="lumi")
parser.add_option("--runMin",action="store",type="float",default="-1",dest="runMin")
parser.add_option("--runMax",action="store",type="float",default="100000000000",dest="runMax")
parser.add_option("--golden",action="store_true",default=False,dest="golden")
parser.add_option("--plotSig",action="store_true",default=False,dest="plotSig")

(options, args) = parser.parse_args()

var = options.var
xmin = options.xmin
xmax = options.xmax
bins = options.bins
xtitle = options.xtitle
rebin = options.rebin
units = options.units
logy = options.logy
outputDir = options.outputDir
inputList_mc = options.inputList_mc
inputList_data = options.inputList_data
lumi = options.lumi
runMin = options.runMin
runMax = options.runMax
golden = options.golden
plotSig = options.plotSig
#############################

CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = str(int(options.lumi))+" pb^{-1} (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
iPeriod = 0
#######################
#minX_mass = 526.
minX_mass = 1181.
#maxX_mass = 5877. 
maxX_mass = 7866. 

massBins_list = [1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7060, 7320, 7589, 7866, 8152, 8447, 8752, 9067, 9391, 9726, 10072, 10430, 10798, 11179, 11571, 11977, 12395, 12827, 13272, 13732, 14000]
    
massBins = array("d",massBins_list)


#############################
#fileNames = ['QCD_Pt-301to470','QCD_Pt-470to600','QCD_Pt-600to800', 'QCD_Pt-800to1000', 'QCD_Pt-1000to1400', 'QCD_Pt-1400to1800', 'QCD_Pt-1800to2400', 'QCD_Pt-2400to3200', 'QCD_Pt-3200']
#xsections = [7475 ,587., 167, 28.25, 8.195, 0.7346, 0.102, 0.00644, 0.000163]
#colorF    = [ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9, ROOT.kBlue-9,ROOT.kBlue-9]
#colorL    = [ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack, ROOT.kBlack,ROOT.kBlack]
hist_allCuts      = []

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


#---- open the files --------------------
#var1 = ""
#if var=="pTWJ_j1" : var1 = "pT_j1"
#elif var=="pTWJ_j2" : var1 = "pT_j2"
#elif var=="etaWJ_j1" : var1 = "eta_j2"
#elif var=="etaWJ_j2" : var1 = "eta_j2"
#elif var=="phiWJ_j1" : var1 = "phi_j1"
#elif var=="phiWJ_j2" : var1 = "phi_j2"
#else : var1 = var


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
  tree = inf.Get('rootTupleTree/tree')

  #no deta
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<2.6 && mjj > '+str(minX_mass))

  #standard
  tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > '+str(minX_mass))
  # blinded 4 TeV
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > '+str(minX_mass)+' && mjj<4000')
  # "peak mjj 4 TeV"
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<3 && mjj > 3558 && mjj < 4509.')
  #MET 200-300 GeV
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > 1118 && MET > 200 && MET < 300')
  #EE
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > 1118 &&  TMath::Abs(etaWJ_j1)>1.4 &&  TMath::Abs(etaWJ_j2)>1.4')
  #BB
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > 1118 &&  TMath::Abs(etaWJ_j1)<1.4 &&  TMath::Abs(etaWJ_j2)<1.4')
  #EB or BE
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > 1118 && ( (TMath::Abs(etaWJ_j1)>1.4 &&  TMath::Abs(etaWJ_j2)<1.4) || (TMath::Abs(etaWJ_j1)<1.4 &&  TMath::Abs(etaWJ_j2)>1.4) )')
  ##eta <2
  #tree.Project(h_allCuts.GetName(), var,'deltaETAjj<1.3 && mjj > 1118 &&  TMath::Abs(etaWJ_j1)<2. &&  TMath::Abs(etaWJ_j2)<2.')
  Npassed = h_allCuts.GetEntries()
  eff = float(Npassed)/Nev
  print('eff : %f' % eff)
  print('(not using efficiency in the weight)')
  if not (i_f == 9):
    wt = options.lumi*float(xsecs[i_f])/Nev
  else:
    wt = options.lumi*float(xsecs[i_f])/(Nev*eff)  ###for the signal the eff is already in the xsec

  print('weight : %f' % wt)
  h_allCuts.Scale(wt)
  h_allCuts.SetDirectory(0)
  h_allCuts.SetFillColor(kBlue-9)
  h_allCuts.SetLineColor(kBlue-9)
  h_allCuts.SetMarkerColor(kBlue-9)
  hist_allCuts.append(h_allCuts)
  print "entries: %d" % h_allCuts.GetEntries()
  print "integral: %f" % h_allCuts.Integral()
   
  i_f += 1

#--------- data  --------------
chain = TChain("rootTupleTree/tree")
for i in range(0,len(fileNames_data)):
  chain.Add(fileNames_data[i])
  print fileNames_data[i]

h_dat = TH1F("h_dat", "", bins, xmin, xmax)
h_dat.Sumw2()
#passJSON
#chain.Project("h_dat",var,'deltaETAjj<1.3 && mjj > '+str(minX_mass)+' && PassJSON==1')
#standard
chain.Project("h_dat",var,'deltaETAjj<1.3 && mjj > '+str(minX_mass))
# "peak mjj 4 TeV"
#chain.Project(h_dat.GetName(), var,'deltaETAjj<3 && mjj>3558 && mjj<4509.')
#no deta
#chain.Project("h_dat",var,'deltaETAjj<2.6 && mjj > '+str(minX_mass))

#if golden:
#  chain.Project("h_dat", var,'deltaETAjj<1.3 && mjj > '+str(minX_mass)+' && run >= '+ str(runMin) + '&& run <= ' + str(runMax) +' && PassJSON==1')
#else:
#  chain.Project("h_dat", var,'deltaETAjj<1.3 && mjj > '+str(minX_mass)+' && run >= '+ str(runMin) + ' && run <= ' + str(runMax) )
  

#h_dat = hist_allCuts[9].Clone()
#h_dat.SetName("h_dat")
#h_dat.SetLineStyle(2)
#h_dat.SetLineWidth(2)
h_dat.SetMarkerColor(kBlack)
h_dat.SetLineColor(kBlack)

                   
#--------------------------		   
NQCD_allCuts = hist_allCuts[0].Integral()

#for i in range(0,len(fileNames)) :
for i in range(1,9) :
  NQCD_allCuts += hist_allCuts[i].Integral()
    
hist_allCutsQCD = hist_allCuts[0].Clone('hist_allCutsQCD')

#for i in range(1,len(fileNames)):
for i in range(1,9):
  hist_allCutsQCD.Add(hist_allCuts[i]) 
  print "hist_allCutsQCD entries: %f" % hist_allCutsQCD.Integral()

hsQCD_allCuts = THStack('QCD_allCuts','QCD_allCuts')

#for i in range(0,len(fileNames)) :
for i in range(0,9) :
  hsQCD_allCuts.Add(hist_allCuts[i])

#--- signal ---
h_sig = hist_allCuts[9]
h_sig.SetLineColor(2)
h_sig.SetLineWidth(2)
h_sig.SetFillStyle(0)

NDAT = h_dat.GetEntries()
NQCD = hist_allCutsQCD.Integral(0,hist_allCutsQCD.GetNbinsX()+1)
NSIG = h_sig.Integral(0,hist_allCutsQCD.GetNbinsX()+1)
## k factor calculated including overflow and underflow
kFactor = NDAT/NQCD
#kFactor = 1.
print ("kFactor set to %f" % kFactor)
print("NQCD = "+str(NQCD))
print("NDAT = "+str(NDAT))
print("NSIG = "+str(NSIG))
print ("NDAT / NQCD = = %f" % (NDAT/NQCD))

hist_allCutsQCD.Scale(kFactor)

print ("---- After scaling signal to bkg (if not plotting mjj) -----")
print ("bkg integral all cuts = %f" % NQCD_allCuts)


#----- Drawing  and save on file -----------------------

outFile = TFile(outputDir+"histo_data_"+var+"_fromTree.root", "recreate")
outFile.cd()
#scale to data
integral_data = h_dat.Integral()
integral_qcd = hist_allCutsQCD.Integral()
#if not integral_qcd==0:
#  hist_allCutsQCD.Scale(integral_data/integral_qcd)
#else:
#  hist_allCutsQCD.Scale(0)
#  print("QCD scaled to 0!")

#rebin only for plot
hist_allCutsQCD.Write()
h_dat.Write()
h_sig.SetName("h_sig")
h_sig.Write()

#Rebin only for the plots, add last bin overflow only for the plot
#hist_allCutsQCD.Rebin(rebin)
#h_dat.Rebin(rebin)
if (var=="mjj" or var=="Dijet_MassAK4") and rebin==-1:
  h_sig_rebin = h_sig.Rebin(len(massBins_list)-1,var+"_rebin",massBins)
  hist_allCutsQCD_rebin = hist_allCutsQCD.Rebin(len(massBins_list)-1,var+"_rebin",massBins)
  h_dat_rebin = h_dat.Rebin(len(massBins_list)-1,var+"_data_rebin",massBins)
  hist_allCutsQCD_rebin.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
  h_dat_rebin.GetXaxis().SetRangeUser(minX_mass,maxX_mass)  
  ## add last bin overflow
  lastbin = h_dat_rebin.FindBin(maxX_mass-1.)
  h_dat_rebin.SetBinContent(lastbin,h_dat_rebin.Integral(lastbin,len(massBins_list)))
  hist_allCutsQCD_rebin.SetBinContent(lastbin,hist_allCutsQCD_rebin.Integral(lastbin,len(massBins_list)))
  filetxt_data = open(outputDir+"/DijetMassSpectrum_data_events_CMS_13TeV_"+str(options.lumi)+"pb-1.txt","w+")  
  filetxt_mc = open(outputDir+"/DijetMassSpectrum_mc_events_CMS_13TeV_"+str(options.lumi)+"pb-1.txt","w+")  
  filetxt_data.write("########################################################################################\n")
  filetxt_data.write("#binLow(GeV)    binHigh(GeV)    events    errLow (evt)    errUp (evt)\n")
  filetxt_data.write("########################################################################################\n")
  filetxt_mc.write("########################################################################################\n")
  filetxt_mc.write("#binLow(GeV)    binHigh(GeV)    events    errLow (evt)    errUp (evt)\n")
  filetxt_mc.write("########################################################################################\n")
  for i in range(1,len(massBins_list)):
    filetxt_data.write(str(h_dat_rebin.GetXaxis().GetBinLowEdge(i))+"     "+str(h_dat_rebin.GetXaxis().GetBinUpEdge(i))+"      "+str(h_dat_rebin.GetBinContent(i))+"        "+str(h_dat_rebin.GetBinError(i))+"      "+str(h_dat_rebin.GetBinError(i))+"\n")
    filetxt_mc.write(str(hist_allCutsQCD_rebin.GetXaxis().GetBinLowEdge(i))+"     "+str(hist_allCutsQCD_rebin.GetXaxis().GetBinUpEdge(i))+"      "+str(hist_allCutsQCD_rebin.GetBinContent(i))+"        "+str(hist_allCutsQCD_rebin.GetBinError(i))+"      "+str(hist_allCutsQCD_rebin.GetBinError(i))+"\n")
  filetxt_data.close() 
  filetxt_mc.close() 

else :
  h_sig_rebin = h_sig.Rebin(rebin)
  hist_allCutsQCD.Rebin(rebin)
  h_dat.Rebin(rebin)
  hist_allCutsQCD_rebin = hist_allCutsQCD.Clone(var+"_rebin")
  h_dat_rebin = h_dat.Clone(var+"_data_rebin")
  ## add last bin overflow
  h_dat_rebin.SetBinContent(
      h_dat_rebin.GetNbinsX(),
      h_dat_rebin.GetBinContent(h_dat_rebin.GetNbinsX()) + h_dat_rebin.GetBinContent(h_dat_rebin.GetNbinsX()+1)
      )
  hist_allCutsQCD_rebin.SetBinContent(
      hist_allCutsQCD_rebin.GetNbinsX(), 
      hist_allCutsQCD_rebin.GetBinContent(hist_allCutsQCD_rebin.GetNbinsX()) + hist_allCutsQCD_rebin.GetBinContent(hist_allCutsQCD_rebin.GetNbinsX()+1)
      )


can_allCuts = TCanvas('can_allCuts_'+var,'can_allCuts_'+var,600,650)

leg = TLegend(0.6, 0.7, 0.85, 0.85)
leg.SetLineColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.AddEntry(hist_allCutsQCD_rebin, "QCD", "f")
if plotSig:
  leg.AddEntry(h_sig, "q* (4 TeV)", "l")
leg.AddEntry(h_dat_rebin, "data", "p")

can_allCuts.cd()

#----- pad 1 -----------
pad1 = TPad("pad1", "pad1",0,0.15,1,1)  
#pad1.SetRightMargin(0.1)

pad1.Draw()
pad1.Clear()
pad1.cd()

max = hist_allCutsQCD_rebin.GetBinContent(hist_allCutsQCD_rebin.GetMaximumBin())

if logy:
  gPad.SetLogy(1)
  hist_allCutsQCD_rebin.SetMaximum(100.*max)  
  h_dat_rebin.SetMaximum(100.*max)
  if not (var=="mjj" or var=="Dijet_MassAK4"):
    hist_allCutsQCD_rebin.SetMinimum(0.1)
    h_dat_rebin.SetMinimum(0.1)
else:
  gPad.SetLogy(0)
  hist_allCutsQCD_rebin.SetMaximum(max + 0.5*max)
  h_dat_rebin.SetMaximum(max + 0.5*max)
  h_sig_rebin.Scale(h_dat_rebin.Integral()/h_sig_rebin.Integral())
  if not (var=="mjj" or var=="Dijet_MassAK4"):
    hist_allCutsQCD_rebin.SetMinimum(0.)
    h_dat_rebin.SetMinimum(0.)


#hist_allCutsQCD.Reset()
hist_allCutsQCD_rebin.GetXaxis().SetTitle(xtitle)
h_dat_rebin.GetXaxis().SetTitle(xtitle)
binwidth = hist_allCutsQCD_rebin.GetBinWidth(1)
if not (var=="mjj" or var=="Dijet_MassAK4"):
  hist_allCutsQCD_rebin.GetYaxis().SetTitle("Events / %.2f %s" % (binwidth, units))
  h_dat_rebin.GetYaxis().SetTitle("Events / %2f %s" )
  h_sig_rebin.GetYaxis().SetTitle("Events / %2f %s")
else:
  hist_allCutsQCD_rebin.GetYaxis().SetTitle("Events / GeV")
  h_dat_rebin.GetYaxis().SetTitle("Events / GeV" )
  h_sig_rebin.GetYaxis().SetTitle("Events / GeV")

#maximumBin = array('f',  [hist_allCutsQCD_rebin.GetBinContent(hist_allCutsQCD_rebin.GetMaximumBin()), h_dat_rebin.GetBinContent(h_dat_rebin.GetMaximumBin())])
#max = TMath.MaxElement(2, maximumBin)
#hist_allCutsQCD_rebin.SetMaximum(1.2*max)
hist_allCutsQCD_rebin.SetMinimum(0.0002)
h_dat_rebin.Draw("p")
hist_allCutsQCD_rebin.Draw("hist")
if plotSig:
  h_sig_rebin.Draw("hist same")
h_dat_rebin.Draw("p same")
leg.Draw()

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

gPad.RedrawAxis()

#-------pad 2------
can_allCuts.cd()
pad2 = TPad("pad2", "pad2",0.,0.,1,0.26)
pad2.SetGrid()
	      
pad2.SetTopMargin(0)
pad2.SetBottomMargin(0.4)
#pad2.SetRightMargin(0.1)
pad2.Draw()	       
pad2.cd()


ratio = h_dat_rebin.Clone("ratio")
ratio.Divide(hist_allCutsQCD_rebin)
ratio.SetFillColor(0)
ratio.SetLineColor(kBlack)
ratio.SetMarkerColor(kBlack)
ratio.GetYaxis().SetRangeUser(0., 2.)
ratio.GetYaxis().SetNdivisions(405, kTRUE)
ratio.GetYaxis().SetTitleFont(42)
ratio.GetYaxis().SetTitle("data / MC")
ratio.GetXaxis().SetTitleSize(0.2)
ratio.GetXaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetLabelSize(0.16)
ratio.GetYaxis().SetTitleSize(0.15)
ratio.GetYaxis().SetTitleOffset(0.5)
#ratio.GetXaxis().SetTitleOffset(0.8)
ratio.Draw("p")

#RedrawAxis
pad2.cd()
gPad.RedrawAxis()

can_allCuts.Write()   
if(logy):
  can_allCuts.SaveAs(outputDir+var+'_allCuts_logy.C')
  can_allCuts.SaveAs(outputDir+var+'_allCuts_logy.png')
  can_allCuts.SaveAs(outputDir+var+'_allCuts_logy.pdf')
if(rebin == -1):
  can_allCuts.SaveAs(outputDir+var+'_allCuts_varBin.C')
  can_allCuts.SaveAs(outputDir+var+'_allCuts_varBin.png')
  can_allCuts.SaveAs(outputDir+var+'_allCuts_varBin.pdf')


else:
  can_allCuts.SaveAs(outputDir+var+'_allCuts.C')
  can_allCuts.SaveAs(outputDir+var+'_allCuts.png')
  can_allCuts.SaveAs(outputDir+var+'_allCuts.pdf')


can_allCuts.Write()   
can_allCuts.Close()

outFile.Close()

##----- keep the GUI alive ------------
#if __name__ == '__main__':
#  rep = ''
#  while not rep in ['q','Q']:
#    rep = raw_input('enter "q" to quit: ')
#    if 1 < len(rep):
#      rep = rep[0]
