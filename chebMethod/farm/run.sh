#!/bin/bash
address=$PWD/../
#/nfs/dust/cms/user/zlebcr/powheg/armando/pythia/pythiaLocal/share/Pythia8/examples
cd $TMP
echo $PWD


source $address/../setup62.sh

cp  $address/one .

./one $1 $2

cp *.root  $address/farm/histos/Sol${2}.root
