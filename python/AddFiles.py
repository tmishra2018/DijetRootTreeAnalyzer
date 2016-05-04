from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from array import *
import os
import sys

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-d','--dir',dest="outDir",default="./",type="string",
                  help="Output directory to store output")
    parser.add_option('-l','--list',dest="list", default="lists/testlist.txt",type="string",
                  help="test list")
    

    (options,args) = parser.parse_args()

    tchain = rt.TChain('rootTupleTree/tree')
    tchain.SetBranchStatus('*',0)
    tchain.SetBranchStatus('mjj',1)
    tchain.SetBranchStatus('passHLT_CaloScoutingHT250',1)
    tchain.SetBranchStatus('deltaETAjj',1)
    tchain.SetBranchStatus('deltaPHIjj',1)
    tchain.SetBranchStatus('phiWJ_j1',1)
    tchain.SetBranchStatus('phiWJ_j2',1)
    tchain.SetBranchStatus('etaWJ_j1',1)
    tchain.SetBranchStatus('etaWJ_j2',1)
    tchain.SetBranchStatus('pTWJ_j1',1)
    tchain.SetBranchStatus('pTWJ_j2',1)
    tchain.SetBranchStatus('deltaETAjjAK4',1)
    tchain.SetBranchStatus('phiAK4_j1',1)
    tchain.SetBranchStatus('phiAK4_j2',1)
    tchain.SetBranchStatus('etaAK4_j1',1)
    tchain.SetBranchStatus('etaAK4_j2',1)
    tchain.SetBranchStatus('pTAK4_j1',1)
    tchain.SetBranchStatus('pTAK4_j2',1)
    
    h_mjj_1GeVbin = rt.TH1D('h_mjj_HLTpass_HT250_1GeVbin','h_mjj_HLTpass_HT250_1GeVbin',14000,0,14000)
    h_deltaETAjj = rt.TH1D('h_deltaETAjj_HLTpass_HT250','h_deltaETAjj_HLTpass_HT250',1000,0,5)
    h_deltaPHIjj = rt.TH1D('h_deltaPHIjj_HLTpass_HT250','h_deltaPHIjj_HLTpass_HT250',1000,0,4)
    h_etaWJ_j1 = rt.TH1D('h_etaWJ_j1_HLTpass_HT250','h_etaWJ_j1_HLTpass_HT250',1000,-5,5)
    h_etaWJ_j2 = rt.TH1D('h_etaWJ_j2_HLTpass_HT250','h_etaWJ_j2_HLTpass_HT250',1000,-5,5)
    h_phiWJ_j1 = rt.TH1D('h_phiWJ_j1_HLTpass_HT250','h_phiWJ_j1_HLTpass_HT250',1000,-4,4)
    h_phiWJ_j2 = rt.TH1D('h_phiWJ_j2_HLTpass_HT250','h_phiWJ_j2_HLTpass_HT250',1000,-4,4)
    h_pTWJ_j1 = rt.TH1D('h_pTWJ_j1_HLTpass_HT250','h_pTWJ_j1_HLTpass_HT250',10000,0,10000)
    h_pTWJ_j2 = rt.TH1D('h_pTWJ_j2_HLTpass_HT250','h_pTWJ_j2_HLTpass_HT250',10000,0,10000)    
    h_deltaETAjjAK4 = rt.TH1D('h_deltaETAjjAK4_HLTpass_HT250','h_deltaETAjjAK4_HLTpass_HT250',1000,0,5)
    h_deltaPHIjjAK4 = rt.TH1D('h_deltaPHIjjAK4_HLTpass_HT250','h_deltaPHIjjAK4_HLTpass_HT250',1000,0,4)
    h_etaAK4_j1 = rt.TH1D('h_etaAK4_j1_HLTpass_HT250','h_etaAK4_j1_HLTpass_HT250',1000,-5,5)
    h_etaAK4_j2 = rt.TH1D('h_etaAK4_j2_HLTpass_HT250','h_etaAK4_j2_HLTpass_HT250',1000,-5,5)
    h_phiAK4_j1 = rt.TH1D('h_phiAK4_j1_HLTpass_HT250','h_phiAK4_j1_HLTpass_HT250',1000,-4,4)
    h_phiAK4_j2 = rt.TH1D('h_phiAK4_j2_HLTpass_HT250','h_phiAK4_j2_HLTpass_HT250',1000,-4,4)
    h_pTAK4_j1 = rt.TH1D('h_pTAK4_j1_HLTpass_HT250','h_pTAK4_j1_HLTpass_HT250',10000,0,10000)
    h_pTAK4_j2 = rt.TH1D('h_pTAK4_j2_HLTpass_HT250','h_pTAK4_j2_HLTpass_HT250',10000,0,10000)
    
    cut = 'passHLT_CaloScoutingHT250&&abs(deltaETAjj)<1.3&&abs(etaWJ_j1)<2.5&&abs(etaWJ_j2)<2.5&&pTWJ_j1>60&&pTWJ_j2>30'
    #cut = 'passHLT_CaloScoutingHT250'
    
    histoExpr = { h_mjj_1GeVbin: 'mjj', 
                  h_deltaETAjj: 'deltaETAjj',
                  h_deltaPHIjj: 'deltaPHIjj',
                  h_etaWJ_j1: 'etaWJ_j1',
                  h_etaWJ_j2: 'etaWJ_j2',
                  h_phiWJ_j1: 'phiWJ_j1',
                  h_phiWJ_j2: 'phiWJ_j2',
                  h_pTWJ_j1: 'pTWJ_j1',
                  h_pTWJ_j2: 'pTWJ_j2',
                  h_deltaETAjjAK4: 'deltaETAjjAK4',
                  h_deltaPHIjjAK4: 'deltaPHIjjAK4',
                  h_etaAK4_j1: 'etaAK4_j1',
                  h_etaAK4_j2: 'etaAK4_j2',
                  h_phiAK4_j1: 'phiAK4_j1',
                  h_phiAK4_j2: 'phiAK4_j2',
                  h_pTAK4_j1: 'pTAK4_j1',
                  h_pTAK4_j2: 'pTAK4_j2',
                }
    
    f = open(options.list)
    for line in f:
        print line.replace('\n','')
        tchain.Add(line.replace('\n',''))

    def project(tree, h, var, cut):
        print 'projecting var: %s, cut: %s from tree: %s into hist: %s'%(var, cut, tree.GetName(), h.GetName())
        tree.Project(h.GetName(),var,cut)

    for hist,var in histoExpr.iteritems():
        #if var!='mjj': continue
        project(tchain, hist, var, cut)

    #tchain.Draw('>>elist',cut,'entrylist')
    #elist = rt.gDirectory.Get('elist')
    #entry = -1
    #while True:
    #    entry = elist.Next()
    #    if entry == -1: break
    #    tchain.GetEntry(entry)
    #    if entry%10000==0: print 'processing entry %i'%entry
    #    for hist,var in histoExpr.iteritems():
    #        hist.Fill(eval('tchain.%s'%var))

    output = rt.TFile(options.outDir+"/output.root",'recreate')
    output.cd()
    for h in histoExpr:
        h.Write()
    output.Close()
