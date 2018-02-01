#!/usr/bin/python
import argparse, os, tempfile, shutil, sys,math,pickle,itertools
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from subprocess import call, PIPE, STDOUT, Popen
import argparse
import datetime



parser = argparse.ArgumentParser(description='Submit jot to batch')
parser.add_argument('--config_file',type=str,dest="config_file",default=1, help=" config_file of the era")
parser.add_argument('--whichrun',type=str,dest="whichrun",default=1, help="Run era")
parser.add_argument('--Nlist',type=int,dest="Nlist",default=1, help="number of list")
parser.add_argument('--dirlist',type=str,dest="dirlist",default=1, help="directory where to finf the lists")
args = parser.parse_args()
config_file = args.config_file #path to the processed lumi JSON file
Runsuffix = args.whichrun # which run
Nlist = args.Nlist
dirlist = args.dirlist

today = datetime.date.today()
today.strftime('%d-%m-%Y')

if not os.path.exists(config_file):
  print("Error: %s not found" % config_file)
  exit(1)

cmd_1="voms-proxy-init --voms cms --valid 100:00 -out $HOME/.globus/gridproxy.cert"
cmd_2="export X509_USER_PROXY=${HOME}/.globus/gridproxy.cert"

#os.system(cmd_1)
#os.system(cmd_2)
List_splited = ['xaa','xab','xac','xad','xae','xaf','xag','xah','xai','xaj','xak','xal','xam','xan','xao','xap','xaq','xar','xas','xat','xau','xav','xaw','xax','xay','xaz','xba','xbb','xbc','xbd','xbe','xbf','xbg','xbh','xbi','xbj','xbk','xbl','xbm','xbn','xbo','xbp','xbq','xbr','xbs','xbt','xbu','xbv','xbw','xbx','xby','xbz','xca','xcb','xcc','xcd','xce','xcf','xcg','xch','xci','xcj','xck','xcl','xcm','xcn','xco','xcp','xcq','xcr','xcs','xct','xcu','xcv','xcw','xcx','xcy','xcz']

for i in range(0, Nlist): 
	list_to_submit = List_splited[i]
	outputN = str(i+1)
	argforscript = "--config_file "+config_file+" --outputname "+Runsuffix+"_"+outputN+" --list "+list_to_submit+" --dirlist "+dirlist
	cmd_3="bsub  -q 1nd -o /afs/cern.ch/work/h/hlattaud/txt_output/"+Runsuffix+"_"+outputN+"_out.txt  /afs/cern.ch/work/h/hlattaud/private/CMSSW_8_0_25/src/CMSDIJET/responsecomputing/DijetRootTreeAnalyzer/Run_Analyser_short.py "+argforscript
	print(cmd_3)
	os.system(cmd_1+" && "+cmd_2+" && "+cmd_3)




