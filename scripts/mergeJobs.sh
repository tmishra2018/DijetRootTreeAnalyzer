#!/bin/bash

#DIR='/pnfs/roma1.infn.it/data/cms/store/user/roma-group1/Dijet/reducedTrees/mc/output_Spring15_v2_JESplus5percent'
#ENDPOINT="dcap://cmsrm-se01.roma1.infn.it/"
DIR='output_Spring15_v2_JECplus5percent/'
ENDPOINT=''
for i in `ls $DIR/*_reduced_skim.root|sed 's/[0-9]*_reduced_skim.root//'|sort -u`
do
    DATASET=`basename $i`
    TARGET="$DIR/merged/${DATASET}reduced_skim.root"
    COMMAND="hadd -f ${ENDPOINT}${TARGET}"
    for i in `ls $DIR/$DATASET*_reduced_skim.root`
    do
	#COMMAND="$COMMAND ${ENDPOINT}/$i"
	COMMAND="$COMMAND $i"
    done
    echo $COMMAND
    $COMMAND
done
