cd /afs/cern.ch/work/l/ltortero/JEC-task/CMSSW_8_0_31/src/CMSDIJET/DijetRootTreeAnalyzer/lists_2017UL

for sample in GJets_HT-100To200_RunIISummer19MiniAOD-106X GJets_HT-200To400_RunIISummer19MiniAOD-106X GJets_HT-400To600_RunIISummer19MiniAOD-106X GJets_HT-40To100_RunIISummer19MiniAOD-106X GJets_HT-600ToInf_RunIISummer19MiniAOD-106X Run2017B-09Aug2019_UL2017-v1 Run2017C-09Aug2019_UL2017-v1 Run2017D-09Aug2019_UL2017-v1 Run2017E-09Aug2019_UL2017-v1 Run2017F-09Aug2019_UL2017-v1
do
    echo $sample
    rm -rf $sample*
    for full_sample_name in $(xrdfs lyogrid06.in2p3.fr ls /dpm/in2p3.fr//home/cms/data/store/user/ltortero/${sample})
    do
        #xrdfs lyogrid06.in2p3.fr ls ${full_sample_name}/crab_${sample}
        most_recent=$(xrdfs lyogrid06.in2p3.fr ls ${full_sample_name}/crab_${sample} | tail -1)
        echo "Use $(basename ${most_recent})"
        for N in $(xrdfs lyogrid06.in2p3.fr ls ${most_recent})
        do
            for file in $(xrdfs lyogrid06.in2p3.fr ls ${N} | grep .root)
            do
                echo 'root://lyogrid06.in2p3.fr/'$file >> ${sample}.txt
            done
        done
    done
    Nfiles=$(cat ${sample}.txt | wc -l)
    echo "Found ${Nfiles} files, splitting..."
    mkdir -p ${sample}
    cd ${sample}
    split -l 10 -a $((${#Nfiles}-1)) -d ../${sample}.txt
    i=0
    while [ $((${#Nfiles}-2)) -gt $i ]
    do
        i=$(( $i + 1 ))
        rename x0 x x*
    done
    cd ..
done

# rm -f Run2018B-17Sep2018-v1.txt
# for ((num=0; num<3; num++)); do 
#     for file in $(gfal-ls srm://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data/store/user/ltortero/Run2018B-17Sep2018-v1/EGamma/crab_Run2018B-17Sep2018-v1/190304_141244/000$num/ | grep .root) ; do 
#         echo 'root://lyogrid06.in2p3.fr//dpm/in2p3.fr/home/cms/data/store/user/ltortero/Run2018B-17Sep2018-v1/EGamma/crab_Run2018B-17Sep2018-v1/190304_141244/000'$num'/'$file >> Run2018B-17Sep2018-v1.txt
#     done
# done
# echo 'You should have 2616 root files, you got:'
# cat Run2018B-17Sep2018-v1.txt | wc -l
