# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:30:06 2015

@author: fbeyer
"""


from nipype.pipeline.engine import MapNode, Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
import nipype.interfaces.afni as afni
import nipype.interfaces.utility as util
'''
Workflow to apply all spatial transformations to each volume
of a time series in a single interpolation
'''
def create_transform_pipeline(name='transfrom_timeseries'):
    # set fsl output type
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
    # initiate workflow
    transform_ts = Workflow(name='transform_timeseries')
    # inputnode
    inputnode=Node(util.IdentityInterface(fields=['orig_ts',
    'anat_head',
    'mat_moco',
    'fullwarp',
    'resolution',
    'brain_mask'
    ]),
    name='inputnode')
    # outputnode
    outputnode=Node(util.IdentityInterface(fields=['trans_ts',
    'trans_ts_mean',
    'trans_ts_masked',
    'resamp_brain',
    ]),
    name='outputnode')
    #resample anatomy
    resample = Node(fsl.FLIRT(datatype='float',
    out_file='T1_resampled.nii.gz'),
    name = 'resample_anat')
    transform_ts.connect([(inputnode, resample, [('anat_head', 'in_file'),
    ('anat_head', 'reference'),
    ('resolution', 'apply_isoxfm')
    ]),
    (resample, outputnode, [('out_file', 'resamp_brain')])
    ])
    # split timeseries in single volumes
    split=Node(fsl.Split(dimension='t',
    out_base_name='timeseries'),
    name='split')
    transform_ts.connect([(inputnode, split, [('orig_ts','in_file')])])
    
    # applymoco premat and fullwarpfield
    applywarp = MapNode(fsl.ApplyWarp(interp='spline',
    relwarp=True,
    out_file='rest2anat.nii.gz',
    datatype='float'),
    iterfield=['in_file', 'premat'],
    name='applywarp')
    transform_ts.connect([(split, applywarp, [('out_files', 'in_file')]),
    (inputnode, applywarp, 
    [('mat_moco', 'premat'),
    ('fullwarp','field_file')]),
    (resample, applywarp, [('out_file', 'ref_file')])
    ])
    # re-concatenate volumes
    merge=Node(fsl.Merge(dimension='t',
    merged_file='rest2anat.nii.gz'),
    name='merge')
    transform_ts.connect([(applywarp,merge,[('out_file','in_files')]),
    (merge, outputnode, [('merged_file', 'trans_ts')])])
    # calculate new mean
    tmean = Node(fsl.maths.MeanImage(dimension='T',
    out_file='rest_mean2anat_lowres.nii.gz'),
    name='tmean')
    transform_ts.connect([(merge, tmean, [('merged_file', 'in_file')]),
    (tmean, outputnode, [('out_file', 'trans_ts_mean')])
    ])
    
    # resample brain mask
    resample_brain = Node(afni.Resample(resample_mode='NN',
    outputtype='NIFTI_GZ',
    out_file='T1_brain_mask_lowres.nii.gz'),
    name = 'resample_brain')
    transform_ts.connect([(inputnode, resample_brain, [('brain_mask', 'in_file')]),
                          (tmean, resample_brain,     [('out_file', 'master')]),
                          (resample_brain, outputnode, [('out_file', 'brain_mask_resamp')])
                          ])
    
    #mask the transformed file
    mask = Node(fsl.ApplyMask(), name="mask")
    transform_ts.connect([(resample_brain,mask, [('out_file', 'mask_file')]),
                          (merge, mask, [('merged_file', 'in_file')]),
                          (mask, outputnode, [('out_file', 'trans_ts_masked')])])
    



    
    return transform_ts