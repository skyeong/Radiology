#! /bin/csh 

set subjname=$argv[1]
set curdir=$PWD
cd ${curdir}/${subjname}

# skull strip using FSL's bet
bet metastasis.nii metastasis_bet.nii -S -B
gunzip *.gz
cd ${curdir}
