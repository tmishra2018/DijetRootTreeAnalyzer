from optparse import OptionParser
import ROOT as rt
import rootTools
from framework import Config
from framework import Drawer
from array import *
from itertools import *
from operator import *
import os
import random
import sys
import math



if __name__ == '__main__':
    ###################################################################
    usage = "usage: python python/bTag_performance.py -c config/bTag.cfg -b BTag2016"
    parser = OptionParser(usage)
    parser.add_option('-c','--config',dest="config",type="string",default="config/run2.config",
                  help="Name of the config file to use")
    parser.add_option('-b','--box',dest="box", default="CaloDijet",type="string",
                  help="box name")

    rt.RooMsgService.instance().setGlobalKillBelow(rt.RooFit.FATAL)
    rt.gStyle.SetPaintTextFormat('+.2f')

    (options,args) = parser.parse_args()
    
    cfg = Config.Config(options.config)
    box = options.box

    lumi        = cfg.getVariables(box,"lumi")
    drawOptions = cfg.getVariables(box,"drawOptions")

    inputDataHistos    = cfg.getVariables(box,"inputDataHistos")
    inputSignalsHistos = cfg.getVariables(box,"inputSignalsHistos")
    ###################################################################


    print inputDataHistos
    print inputSignalsHistos

    plots = Drawer.Drawer(inputDataHistos,inputSignalsHistos)
    plots.printPlots('./plots')


