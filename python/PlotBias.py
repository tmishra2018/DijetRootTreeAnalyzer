from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import random
import sys
import math
from scipy.integrate import quad
from itertools import *
from operator import *


def setStyle():
    rt.gStyle.SetOptStat(0)
    rt.gStyle.SetOptTitle(0)
    rt.gStyle.SetPaintTextFormat("1.2g")
    rt.gROOT.SetBatch()
    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)

def getBiasHistos(varName,toyTree):        
    
    h = rt.TH1D('h_bias','h_bias',70,-4,4)
    toyTree.Project('h_bias','%s'%(varName))

    return h

def setHist(h_data,xTitle,yTitle,color=rt.kBlack):
    h_data.SetMarkerColor(color)
    h_data.SetLineColor(color)
    h_data.SetMarkerStyle(20)
    h_data.SetLineColor(color)
    h_data.GetXaxis().SetTitle(xTitle)
    h_data.GetYaxis().SetTitle(yTitle)
    #h_data.GetXaxis().SetLabelOffset(0.16)
    h_data.GetXaxis().SetLabelSize(0.05)
    h_data.GetYaxis().SetLabelSize(0.05)
    h_data.GetXaxis().SetTitleSize(0.05)
    h_data.GetYaxis().SetTitleSize(0.05)
    #h_data.GetXaxis().SetTitleOffset(0.8)
    #h_data.GetYaxis().SetTitleOffset(0.7)
    h_data.SetMaximum(90)
    #h_data.SetMinimum(2e-1)
    return h_data

    
def print1DBias(c,rootFile,h,func,printName,xTitle,yTitle,lumiLabel="",boxLabel="",plotLabel="",tLeg=None):

    #c.SetLogy(1)
    setHist(h,xTitle,yTitle)
    
    #h.SetLineWidth(2)    
    h.Draw("pe")
    if func!=None:
        #npoints = 1000
        #listx = []
        #listy = []
        #xmin = h.GetXaxis().GetXmin()
        #xmax = h.GetXaxis().GetXmax()
        #for ix in range(0,npoints+1):
        #    x = xmin+ix*(xmax-xmin)/npoints
        #    listx.append(x)
        #    listy.append(h.GetEntries()*func.Eval(x)/func.Integral(xmin,xmax))     
        #graph = rt.TGraph(npoints, array('d',listx), array('d',listy))
        #graph.SetLineColor(rt.kRed)
        #graph.SetLineWidth(2)
        #graph.Draw("same")
        h.Fit(func,"mle")
        func.Draw("same")
    

        
    tLeg = rt.TLegend(0.55,0.6,0.89,0.89)
    tLeg.SetLineColor(rt.kWhite)
    tLeg.SetLineWidth(0)
    tLeg.SetFillStyle(0)
    tLeg.AddEntry(h,"toy data, %s, %s GeV, #mu=%.2f"%(options.model, options.mass,options.r),"lep")
    tLeg.AddEntry(func,"gaus, mean = %.2f, std. = %.2f"%(func.GetParameter(1),func.GetParameter(2)),"l")
            
    tLeg.Draw("same")

    l = rt.TLatex()
    l.SetTextAlign(11)
    l.SetTextSize(0.06)
    l.SetTextFont(62)
    l.SetNDC()
    l.DrawLatex(0.12,0.91,"CMS")
    l.SetTextSize(0.05)
    l.SetTextFont(52)
    l.DrawLatex(0.23,0.91,"Preliminary")
    l.SetTextFont(42)
    l.DrawLatex(0.62,0.91,"%s"%lumiLabel)
    l.SetTextFont(52)
    l.SetTextSize(0.045)

    pdf_dict = {'modexp':'mod. exp.',
                'fourparam':'4-param. dijet',
                'fiveparam':'5-param. dijet'
                }
    l.DrawLatex(0.15,0.82,'gen. pdf = %s'%pdf_dict[options.genPdf])
    l.DrawLatex(0.15,0.77,'fit. pdf = %s'%pdf_dict[options.fitPdf])

    c.cd()
    
    c.Print(printName)
    c.Print(os.path.splitext(printName)[0]+'.C')
    cWrite = c.Clone(os.path.splitext(printName)[0].split('/')[-1])
    rootFile.cd()
    c.Write(os.path.splitext(printName)[0].split('/')[-1])

    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store results")
    parser.add_option('-l','--lumi',dest="lumi", default=3000.,type="float",
                  help="integrated luminosity in pb^-1")
    parser.add_option('-b','--box',dest="box", default="MultiJet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="qq",type="string",
                  help="signal model name")
    parser.add_option('--mass',dest="mass", default='750',type="string",
                  help="mass of resonance")
    parser.add_option('-t','--input-toy-file',dest="inputToyFile", default=None,type="string",
                  help="input toy file")
    parser.add_option('-r',dest="r",default=1,type="float",
                  help="expect signal r value")
    parser.add_option('--gen-pdf',dest="genPdf", default="modexp", choices=['modexp','fourparam','fiveparam'],
                  help="pdf for generating")
    parser.add_option('--fit-pdf',dest="fitPdf", default="fourparam", choices=['modexp','fourparam','fiveparam'],
                  help="pdf for fitting")
    
    (options,args) = parser.parse_args()
     
    setStyle()
    
    box = options.box
    lumi = options.lumi
    cfg = Config.Config(options.config)

    
    toyTree = None
    if options.inputToyFile is not None:
        toyFiles = options.inputToyFile.split(',')
        toyTree = rt.TChain("tree_fit_sb")
        for toyFile in toyFiles:
            toyTree.Add(toyFile)
    
    c = rt.TCanvas('c','c',500,400)    
    c.SetLeftMargin(0.12) 
    c.SetBottomMargin(0.12)
    rootFile = rt.TFile.Open(options.outDir + '/' + 'Plots_%s'%box + '.root','recreate')
    tdirectory = rootFile.GetDirectory(options.outDir)
    if tdirectory==None:
        rootFile.mkdir(options.outDir)
        tdirectory = rootFile.GetDirectory(options.outDir)

    dataString = "Sim. Data"

    eventsLabel = "Toy Datasets"
    
    x = array('d', cfg.getBinning(box)[0]) # mjj binning
    nBins = (len(x)-1)

    toyTree.Print()
        
    h_bias = getBiasHistos('(mu-%.3f)/muErr'%options.r,toyTree)     
    
    btagLabel = ""

    lumiLabel = "%.1f fb^{-1} (13 TeV)" % (lumi/1000.)
    boxLabel = ''

    
    
    gaus_func = rt.TF1("gaus_func","gaus(0)",-4,4)
    gaus_func.SetParameter(0,30)
    gaus_func.SetParameter(1,0)
    gaus_func.SetParameter(2,1)
    print1DBias(c,tdirectory,h_bias,gaus_func,options.outDir+"/bias_%s_%s_lumi-%.3f_r-%.3f_%s_%s_%s.pdf"%(options.model,options.mass,(options.lumi/1000.),options.r,box,options.genPdf,options.fitPdf),"(#hat{#mu} - #mu)/#sigma_{#mu}",eventsLabel,lumiLabel,boxLabel,'',None)
    
 
