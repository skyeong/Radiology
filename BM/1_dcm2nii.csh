#! /bin/csh 

set subjname=$argv[1]
set curdir=$PWD
cd ${curdir}/${subjname}

# converting images and remove unnecessary images
dcm2nii -g n -4 *.dcm
rm -rf co*.nii o*.nii
mv *.nii metastasis.nii

# compressing existing original raw files
mkdir raw
tar -cvzf raw/${subjname}_dicom.tar.gz *.dcm
rm -rf *.dcm
cd ${curdir}
