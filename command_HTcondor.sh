cd $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer

## VOMS

echo ''
echo 'Setting VOMS:'
echo ''

voms-proxy-init --voms cms --valid 168:00 -out ~/.globus/gridproxy.cert
#voms-proxy-init --voms cms --valid 168:00


config_File=$1
#example: config/cutFile_mainGammaplusJetSelection_Newsmearing.txt
dir_data=$2
dir_mc=$3

for arg in $dir_data $dir_mc; do
    ## Define inputs
    ### repertoire contenant les listes
    dir_list=$arg/
    ### prefixe de la liste, voir avec la commande split -l10 -d long_list.txt
    name_list=x
    ### Nombre de listes
    Njob=$(ls -l $dir_list | wc -l)
    Njob=$(($Njob - 1))

    ## Define outputs
    
    if [ $dir_list == $dir_data ]; then output_name=DATA; fi
    if [ $dir_list == $dir_mc ]; then output_name=MC; fi

    thedate=$(date +"%Y-%m-%d")
    output_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/$dir_data/$thedate/output_txtfile/
    errors_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/$dir_data/$thedate/errors_txtfile/
    logs_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/$dir_data/$thedate/logs

    #_________________________________________________________________#
    ## Creating directories

    echo 'Submission to HTcondor:'
    echo '  Inputs in ' $dir_list
    echo '  Config file is ' $config_File
    echo ' ' $Njob 'jobs.'

    echo 'Creating output directories...'
    for dir in $output_dir  $errors_dir  $logs_dir  ; do
        mkdir -p $dir
    done

    ## Submission

    condor_submit output_dir=$output_dir errors_dir=$errors_dir logs_dir=$logs_dir config_File=$config_File output_name=$output_name name_list=$name_list dir_list=$dir_list Njob=$Njob  Submit_to_HTcondor.sub
done

condor_q
