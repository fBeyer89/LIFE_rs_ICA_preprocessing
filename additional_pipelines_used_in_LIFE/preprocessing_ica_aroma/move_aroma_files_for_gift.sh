#!/bin/bash
rm -rf /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/slicesdir.sh
printf "slicesdir -o " >> /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/slicesdir.sh
chmod a+x /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/slicesdir.sh

while read subject
#for subject in LI00037838
do


if [ -f "/data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/${subject}/aroma_with_wmcsf/transforms/withCSFWM/aroma_csfwm_denoised_MNI.nii" ];
then
echo "exists ${subject}"

if [ -f "/data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/${subject}/aroma_with_wmcsf/input_gift_is_aroma_csfwm_denoised_MNI.nii" ];
then
echo "already copied"
else
echo "copying"
cp /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/${subject}/aroma_with_wmcsf/transforms/withCSFWM/aroma_csfwm_denoised_MNI.nii /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/${subject}/aroma_with_wmcsf/input_gift_is_aroma_csfwm_denoised_MNI.nii
fi
printf "/data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/${subject}/aroma_with_wmcsf/transforms/withCSFWM/aroma_csfwm_denoised_MNI.nii /afs/cbs.mpg.de/software/fsl/currentversion.ubuntu-precise-amd64/ubuntu-trusty-amd64/data/fsl-mni152-templates/MNI152_T1_2mm_brain.nii.gz " >> /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/slicesdir.sh

else
echo "${subject} not done yet"
fi



done < /home/raid1/fbeyer/Documents/LIFE/Subjects/OlderSubjects/Third_InclusionScheme_APOE_SOP/521Subjects.txt


#bash /data/pt_life/data_fbeyer/BMI_RSN_analysis/Subjects_512_191/slicesdir.sh
