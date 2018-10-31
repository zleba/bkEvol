#!/bin/bash
address=$PWD/../
#/nfs/dust/cms/user/zlebcr/powheg/armando/pythia/pythiaLocal/share/Pythia8/examples
cd $TMP
echo $PWD


source $address/../setup62.sh

cp  $address/build/one .
#cp one  $address/farm/histos/${1}_${3}.root

./one $1 $2 $3

cp *.root  $address/farm/histos/${1}_${3}.root
