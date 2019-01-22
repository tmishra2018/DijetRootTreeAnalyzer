cd $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer
mkdir -p /afs/cern.ch/work/${USER:0:1}/$USER/private/HT_Condor_output/DijetRootTreeAnalyzer/output_txtfile/
mkdir -p log
voms-proxy-init --voms cms --valid 100:00 -out /afs/cern.ch/user/${USER:0:1}/$USER/.globus/gridproxy.cert
condor_submit output_dir=/afs/cern.ch/work/${USER:0:1}/$USER/private/HT_Condor_output/DijetRootTreeAnalyzer/output_txtfile/ config_File=config/cutFile_mainGammaplusJetSelection_Newsmearing.txt output_name=TEST_condor name_list=test_prefix dir_list=dir_list_test/ Njob=1 Submit_to_HTcondor.sub
# name_list = prefixe de la liste, voir avec la commande split -l10 -d long_list.txt
# dir_list = repertoire contenant les listes
# Njobs = nb de listes
condor_q
