old_JEC=$1
new_JEC=$2
old_JER=$3
new_JER=$4

period=Autumn18
runs=(A B C D)

cd $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer/config

if [ "$old_JEC" == "$new_JEC" ] ; then
    echo 'Not changing JEC version'
else
    for run in $runs ; do
        for lvl in withoutL2Res onlyL2Res L2L3Res JER ; do
            file_to_modify='cutFile_Run2018'$run'_'$lvl'.txt'
            touch $file_to_modify
            for frun in $runs ; do
                sed -i "s|data\/$period'_Run'$frun'_V'$old_JEC'_DATA'\/$period'_Run'$frun'_V'$old_JEC'_DATA'|data\/$period'_Run'$frun'_V'$new_JEC'_DATA'\/$period'_Run'$frun'_V'$new_JEC'_DATA'|g" $file_to_modify
                sed -i "s|data\/$period'_V'$old_JEC'_MC'\/$period'_V'$old_JEC'_MC'|data\/$period'_V'$new_JEC'_MC'\/$period'_V'$new_JEC'_MC'|g" $file_to_modify
            done
        done
    done
fi

if [ "$old_JER" == "$new_JER" ] ; then
    echo 'Not changing JER version'
else
    for run in $runs ; do
        for lvl in withoutL2Res onlyL2Res L2L3Res JER ; do
            file_to_modify='cutFile_Run2018'$run'_'$lvl'.txt'
            touch $file_to_modify
            for frun in ABC D ; do
                for dataMC in DATA MC ; do
                    sed -i "s|data\/JER_textfile\/$period'_Run'$frun'_V'$old_JEC'_'$dataMC\/$period'_Run'$frun'_V'$old_JEC'_'$dataMC|data\/JER_textfile\/$period'_Run'$frun'_V'$new_JEC'_'$dataMC\/$period'_Run'$frun'_V'$new_JEC'_'$dataMC|g" $file_to_modify
                done
            done
        done
    done
fi
