import os
import ROOT as rt
from setTDRStyle import *

class Drawer():
    """Class to draw overlayed histos for data and signals"""

    def __init__(self, hData, hSignal):
        print "Drawer::init"

        self._dataHistos = {}
        self._sigHistos = {}

        for sample,opts in hData.items():
            self.loopfile(opts[0],sample)
        for sample,opts in hSignal.items():
            self.loopfile(opts[0],sample)

    def loopfile(self, infile, sample):
        print "Drawer::loopfile"
        hlist = []
        rootFile = rt.TFile(infile)
        hnames = [k.GetName() for k in rootFile.GetListOfKeys()]
        for name in hnames:
            hlist.append(rootFile.Get(name))
        self._dataHistos[sample] = hlist

    
    def setStyle(self):
        print "Drawer::setStyle"
        setTDRStyle()
        

    def addRatioBox(self, histo):
        print "Drawer::addRatioBox"


    def printPlots(self, outPath):
        print "Drawer::printCanvas"
        self.setStyle()


        for dataplot in self._dataHistos["data"]:
            corrCanvas = rt.TCanvas('c','c',500,500)
            corrCanvas.SetRightMargin(0.15)            

            print type(dataplot)

            myTH1 = dataplot
            print myTH1.GetName()
            #myTH1.Draw()

            #for mcplot in _sigHistos:
            
            corrCanvas.Print(outPath+'/'+dataplot.GetName()+'.pdf')

        #self.addRatioBox()



    #def getHisto
    #def drawSignificance
