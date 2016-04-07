#!/bin/csh 

set STARTDIR = "/cmshome/santanas/CMS/Releases/CMSSW_7_4_15_DijetScouting/src/CMSDIJET/DijetRootTreeAnalyzer/"
set LISTTREE = "/cmshome/santanas/CMS/Releases/CMSSW_7_4_15_DijetScouting/src/CMSDIJET/DijetRootTreeAnalyzer/test/santanas/DijetScouting/list/ScoutingPFCommissioning_reduced.txt"
set TESTTREE = "/t3/users/santanas/Dijet13TeVScouting/rootTrees_reduced/ScoutingPFCommissioning__09_03_2016_20160309_225640/merged/rootfile_ScoutingPFCommissioning__Run2015D-v1__RAW_ScoutingPFCommissioning__09_03_2016_20160309_225640_reduced_skim.root"
set TREENAME = "rootTupleTree/tree"
set CODENAME = "analysisClass_mainDijetScoutingTrigger.C"
set CUTFILE = "config/cutFile_mainDijetScoutingTrigger.txt"
set OUTPUTDIR = "/cmshome/santanas/CMS/Releases/CMSSW_7_4_15_DijetScouting/src/CMSDIJET/DijetRootTreeAnalyzer/test/santanas/DijetScouting/trigger/"
set OUTPUTNAME = "outputTrigger"

cd  $STARTDIR
ls

if ($#argv == 0) then
    echo No arguments ... starting from scratch
    ./scripts/make_rootNtupleClass.sh -f $TESTTREE -t $TREENAME
    ln -sf $CODENAME src/analysisClass.C
    make clean
endif

make
./main $LISTTREE $CUTFILE $TREENAME $OUTPUTDIR/$OUTPUTNAME $OUTPUTDIR/$OUTPUTNAME
cd ..
