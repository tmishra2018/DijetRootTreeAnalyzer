#!/bin/bash

outputdir="plots_comparison_Run2015BC_50ns_DCSreprocessing/"
#outputdir="plots_comparison_50ns_vs_Run2015C/"
#outputdir="plots_comparison_50ns_vs_Run2015D_goldenJson/"
#outputdir="plots_comparison_Run2015C_vs_Run2015D_goldenJson/"
#outputdir="plots_comparison_50ns_vs_Run2015D_DCSonly/"
#outputdir="plots_comparison_25s_74X_vs_75X_noJEC_fullrange/"
#outputdir="plots_comparison_50ns_25ns/"
#outputdir="plots_comparison_50ns_JEC50nsV4_fixedJEC/"
list1="list_for_trigger_Run2015B_plus_Run2015C_50ns_goldenJson_29Aug2015_JEC-Summer15_50nsV4.txt"
#list2="list_for_trigger_Run2015C_25ns_goldenJson_2Sept2015.txt"
#list1="list_for_data_comparison_50ns_JEC_Summer2015_50ns_V4.txt"
#list2="list_for_data_comparison_25ns_JEC_Summer2015_25ns_V3.txt"
#list1="list_for_data_comparison_50ns_JEC_Summer2015_50ns_V4.txt"
#list1="list_for_data_comparison_50ns_JEC_Summer2015_50ns_V4_fixedJEC.txt"
#list2="list_for_data_comparison_25ns_JEC_Summer2015_25ns_V3.txt"
#list1="list_for_data_comparison_50ns_noJEC.txt"
#list2="list_for_data_comparison_25ns_74X_16pb-1_noJEC.txt"
#list2="list_for_data_comparison_25ns_Run2015D_goldenJson_noJEC.txt"
#list2="list_for_plots_Run2015C_CMSSW_7_5_3-75X_dataRun2_HLT_frozen_v2_RelVal_20150930_190058.txt"
list2="list_for_plots_data4T_DCSreprocessed_RunBC_50ns_merged.txt"

name1="Run2015BC_68invpb"
#name2="Run2015C_25ns_74x"
#name2="data_Run2015D"
#list1="list_for_data_comparison_25ns_74X_16pb-1.txt"
#list2="list_for_data_comparison_RelVal_25ns_75X.txt"
#name1="data_25ns_74X"
#name2="Run2015C_25ns_753"
name2="Run2015BC_71invpb"


lumi=1. #not used

mkdir -p $outputdir

python DrawFromTree_comparison_data_data.py --var mjj --xmin 1 --xmax 14000 --xtitle "Dijet Mass [GeV]" --bins 13999 --rebin -1 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var mjj --xmin 1 --xmax 14000 --xtitle "Dijet Mass [GeV]" --bins 13999 --rebin -1 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi 
#python DrawFromTree_comparison_data_data.py --var pTWJ_j1 --xmin 30 --xmax 6000 --xtitle "p_{T}(j1) [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var pTWJ_j2 --xmin 30 --xmax 6000 --xtitle "p_{T}(j2) [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var etaWJ_j1 --xmin -3 --xmax 3 --xtitle "#eta(j1)" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi
#python DrawFromTree_comparison_data_data.py --var etaWJ_j2 --xmin -3 --xmax 3 --xtitle "#eta(j2)" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi
#python DrawFromTree_comparison_data_data.py --var  phiWJ_j1  --xmin -3.1415 --xmax 3.1415  --xtitle "#phi (j1)" --bins 200 --rebin 5  --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi 
#python DrawFromTree_comparison_data_data.py --var  phiWJ_j2  --xmin -3.1415 --xmax 3.1415  --xtitle "#phi (j2)" --bins 200 --rebin 5  --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi 
#python DrawFromTree_comparison_data_data.py --var  massWJ_j1  --xmin 0 --xmax 1000  --xtitle "#m (j1) [GeV]" --bins 200 --rebin 5  --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy 
#python DrawFromTree_comparison_data_data.py --var  massWJ_j2  --xmin 0 --xmax 1000  --xtitle "#m (j2) [GeV]" --bins 200 --rebin 5  --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy 
#python DrawFromTree_comparison_data_data.py --var deltaETAjj --xmin 0 --xmax 2 --xtitle "#Delta#eta(jj)" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi
#python DrawFromTree_comparison_data_data.py --var deltaPHIjj --xmin 0 --xmax 3.14 --xtitle "#Delta#phi(jj)" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi
#python DrawFromTree_comparison_data_data.py --var deltaPHIjj --xmin 0 --xmax 3.14 --xtitle "#Delta#phi(jj)" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#
#python DrawFromTree_comparison_data_data.py --var Dijet_MassAK4 --xmin 1 --xmax 14000 --xtitle "Dijet Mass AK4 [GeV]" --bins 13999 --rebin -1 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var Dijet_MassAK4 --xmin 1 --xmax 14000 --xtitle "Dijet Mass AK4 [GeV]" --bins 13999 --rebin -1 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi 
#python DrawFromTree_comparison_data_data.py --var pTAK4_j1 --xmin 30 --xmax 6000 --xtitle "p_{T}(j1) AK4 [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var pTAK4_j2 --xmin 30 --xmax 6000 --xtitle "p_{T}(j2) AK4 [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var etaAK4_j1 --xmin -3 --xmax 3 --xtitle "#eta(j1) AK4" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi
#python DrawFromTree_comparison_data_data.py --var etaAK4_j2 --xmin -3 --xmax 3 --xtitle "#eta(j2) AK4" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi


#python DrawFromTree_comparison_data_data.py --var MET  --xmin 0 --xmax 1000 --xtitle "MET [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var metSig  --xmin 0 --xmax 1 --xtitle "MET / #Sigma E_{T}" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var nVtx --xmin 0 --xmax 50 --xtitle "nvtx" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var chargedHadEnFrac_j1 --xmin 0 --xmax 1 --xtitle "Charged hadron En. fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var neutrHadEnFrac_j1 --xmin 0 --xmax 1 --xtitle "Neutral hadron En.fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var photonEnFrac_j1  --xmin 0 --xmax 1 --xtitle "Photon En. fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var eleEnFract_j1 --xmin 0 --xmax 1 --xtitle "Electron En. fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy        
#python DrawFromTree_comparison_data_data.py --var muEnFract_j1 --xmin 0 --xmax 1 --xtitle "Muon En. fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var neutrElectromFrac_j1  --xmin 0 --xmax 1 --xtitle "neutr. EM En. fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var  chargedElectromFrac_j1  --xmin 0 --xmax 1 --xtitle "ch. EM En. fraction j1" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var chargedMult_j1 --xmin 0 --xmax 50 --xtitle "charged mult. j1" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy       
#python DrawFromTree_comparison_data_data.py --var neutrMult_j1  --xmin 0 --xmax 50 --xtitle "neutral mult. j1" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy        
#python DrawFromTree_comparison_data_data.py --var photonMult_j1 --xmin 0 --xmax 50 --xtitle "photon mult. j1" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy        
#python DrawFromTree_comparison_data_data.py --var jetPtAK4matchCaloJet_j1  --xmin 30 --xmax 5000 --xtitle "pT(j1) calo match [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var chargedHadEnFrac_j2 --xmin 0 --xmax 1 --xtitle "Charged hadron En. fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var neutrHadEnFrac_j2 --xmin 0 --xmax 1 --xtitle "Neutral hadron En.fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var photonEnFrac_j2  --xmin 0 --xmax 1 --xtitle "Photon En. fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var eleEnFract_j2 --xmin 0 --xmax 1 --xtitle "Electron En. fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy        
#python DrawFromTree_comparison_data_data.py --var muEnFract_j2 --xmin 0 --xmax 1 --xtitle "Muon En. fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var neutrElectromFrac_j2  --xmin 0 --xmax 1 --xtitle "neutr. EM En. fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var  chargedElectromFrac_j2  --xmin 0 --xmax 1 --xtitle "ch. EM En. fraction j2" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var chargedMult_j2 --xmin 0 --xmax 50 --xtitle "charged mult. j2" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy       
#python DrawFromTree_comparison_data_data.py --var neutrMult_j2  --xmin 0 --xmax 50 --xtitle "neutral mult. j2" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy        
#python DrawFromTree_comparison_data_data.py --var photonMult_j2 --xmin 0 --xmax 50 --xtitle "photon mult. j2" --bins 50 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy        
#python DrawFromTree_comparison_data_data.py --var jetPtAK4matchCaloJet_j2  --xmin 30 --xmax 5000 --xtitle "pT(j2)calo match [GeV]" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var CosThetaStarWJ --xmin -1 --xmax 1  --xtitle "cos #theta *" --bins 200 --rebin 5 --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi --logy
#python DrawFromTree_comparison_data_data.py --var  Nak4 --xmin 0 --xmax 15  --xtitle "N jets" --bins 15  --outputDir $outputdir --inputList_1 $list1 --inputList_2 $list2 --name1 $name1 --name2 $name2 --lumi $lumi 
