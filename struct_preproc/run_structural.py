# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 11:10:25 2015

@author: fbeyer
"""

from structural import create_structural
import sys
'''
Meta script to run structural preprocessing
------------------------------------------
Can run in two modes:
python run_structural.py s {subject_id}
python run_structural.py f {text file containing list of subjects}
'''
mode=sys.argv[1]

if mode == 's':
    subjects=[sys.argv[2]]
    print subjects
elif mode == 'f':
    with open(sys.argv[2], 'r') as f:
        subjects = [line.strip() for line in f]

for subject in subjects:
    print 'Running subject '+subject
    working_dir = '/scr/kennedy2/data_fbeyer/genetics/Lemon_mod/'+subject+'/'
    
    
    data_dir = '/scr/kennedy2/data_fbeyer/genetics/subjects/'+subject+'/'
    out_dir = '/scr/kennedy2/data_fbeyer/genetics/subjects/' +subject+'/'
    
    #change depending on where freesurfer files are located
    freesurfer_dir = '/scr/kennedy2/LIFE//freesurfer_all/'
   
    
    ##########you have to change this!!
    standard_brain = '/usr/share/fsl/5.0/data/standard/MNI152_T1_2mm_brain.nii.gz'
    create_structural(subject=subject, working_dir=working_dir, data_dir=data_dir,
    freesurfer_dir=freesurfer_dir, out_dir=out_dir,
    standard_brain=standard_brain)