old_JEC=$1
new_JEC=$2
old_JER=$3
new_JER=$4

period=Summer19UL17

cd $CMSSW_BASE/src/CMSDIJET/DijetRootTreeAnalyzer/config

if [ "$old_JEC" == "$new_JEC" ] ; then
    echo 'Not changing JEC version'
else
    for run in B C D E F ; do
        for lvl in withoutL2Res onlyL2Res L2L3Res JER ; do
            for L1_approach in ComplexL1 SimpleL1 ; do
                file_to_modify=cutFile_Run2017UL_${run}_${L1_approach}_${lvl}.txt
                touch $file_to_modify
                echo 'Modify '$file_to_modify
                for frun in B C D E F ; do
                    sed -i "s|data/${period}_Run${frun}_V${old_JEC}_${L1_approach}_DATA/${period}_Run${frun}_V${old_JEC}_${L1_approach}_DATA|data/${period}_Run${frun}_V${new_JEC}_${L1_approach}_DATA/${period}_Run${frun}_V${new_JEC}_${L1_approach}_DATA|g" $file_to_modify
                    sed -i "s|data/${period}_V${old_JEC}_${L1_approach}_MC/${period}_V${old_JEC}_${L1_approach}_MC|data/${period}_V${new_JEC}_${L1_approach}_MC/${period}_V${new_JEC}_${L1_approach}_MC|g" $file_to_modify
                done
            done
        done
    done
fi

period=Fall17

if [ "$old_JER" == "$new_JER" ] ; then
    echo 'Not changing JER version'
else
    for run in B C D E F ; do
        for lvl in withoutL2Res onlyL2Res L2L3Res JER ; do
            for L1_approach in ComplexL1 SimpleL1 ; do
                file_to_modify=cutFile_Run2017UL_${run}_${lvl}.txt
                touch $file_to_modify
                echo 'Modify '$file_to_modify
                for dataMC in DATA MC ; do
                    sed -i "s|data/JER_textfile/${period}_V${old_JER}_${dataMC}/${period}_V${old_JER}_${dataMC}|data/JER_textfile/${period}_V${new_JER}_${dataMC}/${period}_V${new_JER}_${dataMC}|g" $file_to_modify
                done
            done
        done
    done
fi
