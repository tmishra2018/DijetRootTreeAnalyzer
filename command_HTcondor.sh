cd $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer

## VOMS

echo ''
echo 'Setting VOMS:'
echo ''

voms-proxy-init --voms cms --valid 168:00 -out ~/.globus/gridproxy.cert
#voms-proxy-init --voms cms --valid 168:00

thedate=$(date +"%Y-%m-%d_%H-%M-%S")
args=($@)
for ((i=0;i<${#@};i+=3)); do
    ## Define inputs
    config_File=${args[$i]}
    ### repertoire contenant les listes
    dir_list=${args[$i+1]}
    ### fine naming
    specialname=${args[$i+2]}
    ### prefixe de la liste, voir avec la commande split -l10 -d long_list.txt
    name_list=x
    ### Nombre de listes
    Njob=$(ls -l $dir_list/ | wc -l)
    Njob=$(($Njob - 1))

    ## Define outputs
    
    output_name=$(basename $dir_list)'_'$specialname
    output_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/$dir_list/$specialname/$thedate/output_txtfile/
    errors_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/$dir_list/$specialname/$thedate/errors_txtfile/
    logs_dir=/afs/cern.ch/work/${USER:0:1}/$USER/JEC-task/HT_Condor_output/DijetRootTreeAnalyzer/$dir_list/$specialname/$thedate/logs

    #_________________________________________________________________#
    ## Creating directories

    echo ''
    echo ''
    echo 'Submission to HTcondor: ' $output_name
    echo '  Config file is ' $config_File
    echo '  Inputs in ' $dir_list
    echo ' ' $Njob 'jobs.'

    echo 'Creating output directories...'
    for dir in $output_dir  $errors_dir  $logs_dir  ; do
        mkdir -p $dir
    done

    ## Submission

    condor_submit output_dir=$output_dir errors_dir=$errors_dir logs_dir=$logs_dir config_File=$config_File output_name=$output_name name_list=$name_list dir_list=$dir_list/ Njob=$Njob  Submit_to_HTcondor.sub
done

echo ''
echo 'Done.'
echo ''
echo ''
condor_q
