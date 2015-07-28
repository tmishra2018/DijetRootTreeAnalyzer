#!/bin/bash

DIR='/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/mc/output_Spring15_v2_JESminus5percent/'
#DIR='/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/data/25July2015-CertJson-251586-251883_JEC-Summer15_50nsV2_501dfb2/'
ENDPOINT="dcap://cmsrm-se01.roma1.infn.it/"
for i in `ls $DIR/*_reduced_skim.root|sed 's/[0-9]*_reduced_skim.root//'|sort -u`
do
    DATASET=`basename $i`
    TARGET="${DATASET}reduced_skim.root"
    #TARGET="${DIR}/merged/${DATASET}reduced_skim.root"
    COMMAND="hadd -f /tmp/${TARGET}"
    #COMMAND="hadd -f ${ENDPOINT}${TARGET}"
    for i in `ls $DIR/$DATASET*_reduced_skim.root`
    do
	COMMAND="$COMMAND ${ENDPOINT}$i"
	COMMAND2="dccp /tmp/${TARGET} $DIR/merged/${TARGET}"
	#COMMAND="$COMMAND $i"
    done
    echo $COMMAND
    $COMMAND
    $COMMAND2
done
