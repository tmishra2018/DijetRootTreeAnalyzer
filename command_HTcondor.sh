cd $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer

## Define inputs
### repertoire contenant les listes
dir_list=$1
### prefixe de la liste, voir avec la commande split -l10 -d long_list.txt
name_list=x
### Nombre de listes
Njob=$(ls -l $dir_list | wc -l)

config_File=config/cutFile_mainGammaplusJetSelection_Newsmearing.txt

## Define outputs

output_name=TEST_condor

output_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/output_txtfile/$dir_list
errors_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/errors_txtfile/$dir_list
logs_dir=log

#_________________________________________________________________#
## Creating directories

echo 'Submission to HTcondor:'
echo '  Inputs in ' $dir_list
echo ' ' $Njob 'jobs.'

echo 'Creating output directories...'
for dir in $output_dir  $errors_dir  $logs_dir  ; do
    mkdir -p $dir
done

## VOMS

echo ''
echo 'Setting VOMS:'
echo ''

voms-proxy-init --voms cms --valid 100:00 -out ~/.globus/gridproxy.cert

## Submission

condor_submit output_dir=$output_dir errors_dir=$errors_dir logs_dir=$logs_dir config_File=$config_File output_name=$output_name name_list=$name_list dir_list=$dir_list Njob=$Njob  Submit_to_HTcondor.sub

condor_q
