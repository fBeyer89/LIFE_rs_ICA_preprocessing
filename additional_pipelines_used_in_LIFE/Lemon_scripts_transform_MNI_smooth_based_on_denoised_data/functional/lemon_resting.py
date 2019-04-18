# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:26:20 2015

@author: fbeyer
"""

from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
from strip_rois import strip_rois_func
from transform_timeseries import create_transform_pipeline
from ants_registration import create_ants_registration_pipeline
from smoothing import create_smoothing_pipeline
from visualize import create_visualize_pipeline
from normalize import create_normalize_pipeline
from slicetiming_correction import create_slice_timing_pipeline

'''
Main workflow for lemon resting state preprocessing.
====================================================
Uses file structure set up by conversion script.
'''
def create_lemon_resting(subject, working_dir, data_dir, freesurfer_dir, out_dir,
    vol_to_remove, TR, epi_resolution, highpass, lowpass,
    echo_space, te_diff, pe_dir, standard_brain, standard_brain_resampled, standard_brain_mask, 
    standard_brain_mask_resampled, fwhm_smoothing):
    # set fsl output type to nii.gz
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
    # main workflow
    func_preproc = Workflow(name='lemon_resting')
    func_preproc.base_dir = working_dir
    func_preproc.config['execution']['crashdump_dir'] = func_preproc.base_dir + "/crash_files"
    # select files
    templates={
    #'func': 'func/EPI_t2.nii',
    #'fmap_fullwarp' : 'unwarp/B0_ph.nii',
    #'mat_moco': '',
    #'transform_ts': 'lemon_resting/transform_timeseries/merge/rest2anat.nii.gz'
    'anat_head' : 'structural/T1.nii.gz', #either with mod or without
    'anat_brain' : 'structural/brain.nii.gz', #new version with brain_extraction from freesurfer  #T1_brain_brain.nii.gz',
    'brain_mask' : 'structural/T1_brain_mask.nii.gz', #T1_brain_brain_mask.nii.gz',
    'ants_affine': 'structural/transforms2mni/transform0GenericAffine.mat',
    'ants_warp':   'structural/transforms2mni/transform1Warp.nii.gz',
    'transform_ts': 'resting_state/coregister/rest_coregistered_nativespace.nii.gz'
    }
    
    
    selectfiles = Node(nio.SelectFiles(templates,
    base_directory=data_dir),
    name="selectfiles")
    
    #only do transformation, slice-timing and warping to MNI
    #for FRANZ preprocessed data 
    
 
    #workflow to convert signal into percent signal change
    normalize=create_normalize_pipeline()
    normalize.inputs.inputnode.tr = TR 
 
    #workflow to transform timeseries to MNI
    ants_registration=create_ants_registration_pipeline()
    ants_registration.inputs.inputnode.ref=standard_brain_resampled    
    
    #workflow to smooth
    smoothing = create_smoothing_pipeline() 
    smoothing.inputs.inputnode.fwhm=fwhm_smoothing
   
   
    #visualize registration results
    visualize = create_visualize_pipeline()
    visualize.inputs.inputnode.mni_template=standard_brain_resampled 
    
    
    #sink to store files
    sink = Node(nio.DataSink(parameterization=False,
    base_directory=out_dir,
    substitutions=[
    ('rest_coregistered_nativespace_trans_smooth.nii', 'rest_mni_smoothed.nii'),
    ('rest_coregistered_nativespace_trans.nii.gz', 'rest_ants_tr')]
    ),
    name='sink')
    
    
    
    # connections
    func_preproc.connect([
    
    #for FRANZ preprocessed data
    (selectfiles, ants_registration,    [('transform_ts', 'inputnode.denoised_ts')]),
    
    #added by me
    (selectfiles, ants_registration, [('ants_affine', 'inputnode.ants_affine')]),
    (selectfiles, ants_registration, [('ants_warp', 'inputnode.ants_warp')]),

    (ants_registration, smoothing, [('outputnode.ants_reg_ts', 'inputnode.ts_transformed')]),

    (smoothing, visualize,  [('outputnode.ts_smoothed', 'inputnode.ts_transformed')]),

    (ants_registration, sink, [('outputnode.ants_reg_ts', 'ants.@antsnormalized')]),
    (smoothing, sink, [('outputnode.ts_smoothed', '@smoothed.FWHM6')]),
    (visualize, sink, [('outputnode.output_image', 'check.notnormalized.@registration')])
    ])


    func_preproc.run(plugin='MultiProc')#
    

