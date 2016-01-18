#!/bin/bash

DIR='/eos/cms/store/group/phys_exotica/dijet/Dijet13TeVScouting/rootTrees_reduced/ScoutingPFHT__17_01_2016_20160117_182422/'
eos="/afs/cern.ch/project/eos/installation/0.3.84-aquamarine/bin/eos.select"
$eos mkdir -p $DIR/merged/
ENDPOINT="root://eoscms/"
for i in `$eos ls $DIR/ | grep _reduced_skim.root |sed 's/[0-9]*_reduced_skim.root//'|sort -u`
do
    echo $i
    DATASET=`basename $i`
    TARGET="${DATASET}reduced_skim.root"
    #TARGET="${DIR}/merged/${DATASET}reduced_skim.root"
    COMMAND="hadd -f /tmp/${TARGET}"
    #COMMAND="hadd -f ${ENDPOINT}${TARGET}"
    for i in `$eos ls $DIR/$DATASET*_reduced_skim.root`
    do
	COMMAND="$COMMAND ${ENDPOINT}${DIR}$i"
	#COMMAND="$COMMAND $i"
    done
    COMMAND2="xrdcp /tmp/${TARGET} ${ENDPOINT}${DIR}/merged/${TARGET}"
    echo $COMMAND
    echo $COMMAND2
    $COMMAND
    $COMMAND2
done
