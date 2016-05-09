from optparse import OptionParser
import ROOT as rt
import sys
import rootTools
import glob
from math import *
from framework import Config
import os
from array import *

def getFileName(hybridLimit, massPoint, box, model, lumi,  directory, method,t):
    if hybridLimit == "higgsCombineToys":
        fileName = "%s/%s%s_%s_lumi-%.3f_%s_%i.%s.mH120.root"%(directory,hybridLimit,model,massPoint,lumi,box,t,method)
    else:
        fileName = "%s/%s%s_%s_lumi-%.3f_%s.%s.mH120.root"%(directory,hybridLimit,model,massPoint,lumi,box,method)
    return fileName


def writeXsecTree(box, model, directory, massPoint, xsecULObs, xsecULExpPlus2, xsecULExpPlus, xsecULExp, xsecULExpMinus, xsecULExpMinus2):
    outputFileName = "%s/xsecUL_%s_%s_%s.root" %(directory, model, massPoint, box)
    print "INFO: xsec UL values being written to %s"%outputFileName
    fileOut = rt.TFile.Open(outputFileName, "recreate")
    
    xsecTree = rt.TTree("xsecTree", "xsecTree")
    myStructCmd = "struct MyStruct{Double_t mass;"
    ixsecUL = 0
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+0)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+1)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+2)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+3)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+4)
    myStructCmd+= "Double_t xsecUL%i;"%(ixsecUL+5)
    ixsecUL+=6
    myStructCmd += "}"
    rt.gROOT.ProcessLine(myStructCmd)
    from ROOT import MyStruct

    s = MyStruct()
    xsecTree.Branch("mass", rt.AddressOf(s,"mass"),'mass/D')
    
    
    s.mass = float(massPoint)
    
    ixsecUL = 0
    xsecTree.Branch("xsecULObs_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+0)),'xsecUL%i/D'%(ixsecUL+0))
    xsecTree.Branch("xsecULExpPlus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+1)),'xsecUL%i/D'%(ixsecUL+1))
    xsecTree.Branch("xsecULExpPlus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+2)),'xsecUL%i/D'%(ixsecUL+2))
    xsecTree.Branch("xsecULExp_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+3)),'xsecUL%i/D'%(ixsecUL+3))
    xsecTree.Branch("xsecULExpMinus_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+4)),'xsecUL%i/D'%(ixsecUL+4))
    xsecTree.Branch("xsecULExpMinus2_%s"%box, rt.AddressOf(s,"xsecUL%i"%(ixsecUL+5)),'xsecUL%i/D'%(ixsecUL+5))
    exec 's.xsecUL%i = xsecULObs[ixsecUL]'%(ixsecUL+0)
    exec 's.xsecUL%i = xsecULExpPlus2[ixsecUL]'%(ixsecUL+1)
    exec 's.xsecUL%i = xsecULExpPlus[ixsecUL]'%(ixsecUL+2)
    exec 's.xsecUL%i = xsecULExp[ixsecUL]'%(ixsecUL+3)
    exec 's.xsecUL%i = xsecULExpMinus[ixsecUL]'%(ixsecUL+4)
    exec 's.xsecUL%i = xsecULExpMinus2[ixsecUL]'%(ixsecUL+5)
    ixsecUL += 4

    xsecTree.Fill()

    fileOut.cd()
    xsecTree.Write()
    
    fileOut.Close()
    
    return outputFileName

if __name__ == '__main__':

    
    parser = OptionParser()
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")
    parser.add_option('-m','--model',dest="model", default="gg",type="string",
                  help="signal model name")
    parser.add_option('-l','--lumi',dest="lumi", default=1.981,type="float",
                  help="lumi in fb^-1, e.g.: 0.210")
    parser.add_option('--mass',dest="mass", default='750',type="string",
                  help="mass of resonance")
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Input/Output directory to store output")
    parser.add_option('--signif',dest="doSignificance",default=False,action='store_true',
                  help="for significance instead of limit")
    parser.add_option('--toys',dest="doHybridNew",default=False,action='store_true',
                  help="for toys instead of asymptotic")
    parser.add_option('--xsec',dest="xsec",default=1,type="float",
                  help="reference xsec")

    (options,args) = parser.parse_args()

    boxInput = options.box
    model = options.model
    lumi = options.lumi
    directory = options.outDir
    doHybridNew = options.doHybridNew
    doSignificance = options.doSignificance
    refXsec = options.xsec
    

    boxes = boxInput.split('_')

    box = boxInput       

    haddOutputs = []

    
    if len(options.mass.split(','))==1:
        massIterable = [options.mass]
    else:
        massIterable = list(eval(options.mass))
    for massPoint in massIterable:
        
        if doSignificance and doHybridNew:
            if not glob.glob(getFileName("higgsCombineSignif",massPoint,boxInput,model,lumi,directory,"HybridNew",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombineSignif",massPoint,boxInput,model,lumi,directory,"HybridNew",0))
            tFile = rt.TFile.Open(getFileName("higgsCombineSignif",massPoint,boxInput,model,lumi,directory,"HybridNew",0))
        elif doHybridNew: 
            if not glob.glob(getFileName("higgsCombineToys",massPoint,boxInput,model,lumi,directory,"HybridNew",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombineToys",massPoint,boxInput,model,lumi,directory,"HybridNew",0))
            tFile = rt.TFile.Open(getFileName("higgsCombineToys",massPoint,boxInput,model,lumi,directory,"HybridNew",0))
        elif doSignificance: 
            if not glob.glob(getFileName("higgsCombine",massPoint,boxInput,model,lumi,directory,"ProfileLikelihood",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombine",massPoint,boxInput,model,lumi,directory,"ProfileLikelihood",0))
            tFile = rt.TFile.Open(getFileName("higgsCombine",massPoint,boxInput,model,lumi,directory,"ProfileLikelihood",0))
        else:
            if not glob.glob(getFileName("higgsCombine",massPoint,boxInput,model,lumi,directory,"Asymptotic",0)): continue
            print "INFO: opening %s"%(getFileName("higgsCombine",massPoint,boxInput,model,lumi,directory,"Asymptotic",0))
            tFile = rt.TFile.Open(getFileName("higgsCombine",massPoint,boxInput,model,lumi,directory,"Asymptotic",0))

        try:
            if tFile.InheritsFrom("TFile") is False:
                continue
        except:
            continue
            
        limit = tFile.Get("limit")
        try:
            if limit.InheritsFrom("TTree") is False: 
                tFile.cd()
                tFile.Close()
                continue
        except:
            tFile.cd()
            tFile.Close()
            continue
        if doSignificance and limit.GetEntries() < 1: 
            tFile.cd()
            tFile.Close()
            continue
        if (not doSignificance) and limit.GetEntries() < 5: 
            tFile.cd()
            tFile.Close()
            continue
        limit.Draw('>>elist','','entrylist')
        elist = rt.gDirectory.Get('elist')
        entry = elist.Next()
        limit.GetEntry(entry)
        limits = []
        while True:
            if entry == -1: break
            limit.GetEntry(entry)
                
            if doSignificance:
                limits.append(max(0.0,limit.limit))
            else:                    
                limits.append(refXsec*limit.limit)
            entry = elist.Next()
        tFile.cd()
        tFile.Close()

        # no observed limit
        if len(limits)==5:
            limits.append(0)
            
        limits.reverse()
        print massPoint
        print limits
        
        haddOutput = writeXsecTree(boxInput, model, directory, massPoint, [limits[0]],[limits[1]],[limits[2]],[limits[3]],[limits[4]],[limits[5]])
        haddOutputs.append(haddOutput)


    #if doSignificance:
    #    c = plotSignificance(boxInput,model,sigHist,doHybridNew)
    #else:
    if doHybridNew:
        os.system("hadd -f %s/xsecUL_HybridNew_%s_%s.root %s"%(directory,model,boxInput," ".join(haddOutputs)))
    else:
        os.system("hadd -f %s/xsecUL_Asymptotic_%s_%s.root %s"%(directory,model,boxInput," ".join(haddOutputs)))
