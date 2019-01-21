#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime
import os
parser = argparse.ArgumentParser(description='Submit jot to batch')
parser.add_argument('--config_file',type=str,dest="config_file",default=1, help=" config_file of the era")
parser.add_argument('--outputname',type=str,dest="outputname",default=1, help="Run era")
parser.add_argument('--list',type=str,dest="list",default=1, help="list of file to run")
parser.add_argument('--dirlist',type=str,dest="dirlist",default=1, help="directory where are lists of file to run")
args = parser.parse_args()
config_file = args.config_file #path to the processed lumi JSON file
outputname = args.outputname # which run
listtorun = args.list
dirlist = args.dirlist
outputdir = "/eos/user/{}/{}/private/Output_lsfBatch/".format(os.environ["USER"][0], os.environ["USER"])
cmd4="./main "+dirlist+listtorun+" "+config_file+" dijets/events "+outputdir+outputname+" "+outputdir+outputname
cmd1="cd /afs/cern.ch/work/${USER:0:1}/$USER/private/HT_Condor_output/DijetRootTreeAnalyzer/"
cmd2="export SCRAM_ARCH=slc6_amd64_gcc530"
cmd3="eval `scramv1 runtime -sh`"
cmd5="export X509_USER_PROXY=/afs/cern.ch/user/${USER:0:1}/$USER/.globus/gridproxy.cert" 
call(cmd1+" && "+cmd5+" && "+cmd2+" && "+cmd3+" && "+cmd4, shell=True )
print(cmd4)

