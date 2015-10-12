# LIFE_rs_ICA_preprocessing

Scripts used to preprocess LIFE resting state data (based on J. Huntenburgs preprocessing scripts available on github 
(https://github.com/juhuntenburg/pipelines/tree/master/src/lsd_lemon))

Structural preprocessing (run_structural.py using structural.py) performs

- freesurfer segmentation if not yet done
- brain extraction based on freesurfer brain segmentation
- registration of the skull-stripped brain to MNI (2mm isotropic) space
- 
Functional preprocessing (run_lemon_resting.py using lemon_resting.py) performs

- deletion of first 5 volumes
- coregistration of mean EPI to fieldmap magnitude scan and calculation of unwarping field
- motion correction
- applying unwarping and motion correction to time series
- slice time correction
- registration to 3mm isotropic MNI space
- smoothing
- visualizing and saving of .png images (functional and anatomical images overlaid)

