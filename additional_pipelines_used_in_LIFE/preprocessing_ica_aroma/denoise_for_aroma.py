# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:34:01 2015

@author: fbeyer
"""

from nipype.pipeline.engine import Node, Workflow, MapNode
import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
import nipype.algorithms.rapidart as ra
from compcor import extract_noise_components
from normalize_timeseries import time_normalizer
from fix_header_tr import fix_TR_fs

'''
Main workflow for denoising
Largely based on https://github.com/nipy/nipype/blob/master/examples/
rsfmri_vol_surface_preprocessing_nipy.py#L261
use denoising in functional space (project csf/wm mask to functional)
only regress wm+csf compcor-components (not motion components)
using highpass instead of bandpassfiltering and no normalization of timeseries
'''


def create_denoise_pipeline(name='denoise'):
    # workflow
    denoise = Workflow(name='denoise')
    # Define nodes
    inputnode = Node(interface=util.IdentityInterface(fields=['anat_brain',
                                                              'brain_mask',
                                                              'flirt_mat',
                                                              'unwarped_mean',
                                                              'epi_coreg',
                                                              'highpass_sigma',
                                                              'tr']),
                     name='inputnode')
    outputnode = Node(interface=util.IdentityInterface(fields=['wmcsf_mask',
                                                               'brain2epi',
                                                               'wmcsf_mask2epi',
                                                               'combined_motion',
                                                               'comp_regressor',
                                                               'comp_F',
                                                               'comp_pF',
                                                               'out_betas',
                                                               'ts_fullspectrum',
                                                               'normalized_file']),
                      name='outputnode')
    # run fast to get tissue probability classes
    fast = Node(fsl.FAST(), name='fast')
    denoise.connect([(inputnode, fast, [('anat_brain', 'in_files')])])
    # functions to select tissue classes
    def selectindex(files, idx):
        import numpy as np
        from nipype.utils.filemanip import filename_to_list, list_to_filename
        return list_to_filename(np.array(filename_to_list(files))[idx].tolist())

    def selectsingle(files, idx):
        return files[idx]


    # binarize tissue classes
    binarize_tissue = MapNode(fsl.ImageMaths(op_string='-nan -thr 0.99 -ero -bin'),
                              iterfield=['in_file'],
                              name='binarize_tissue')
    denoise.connect([(fast, binarize_tissue, [(('partial_volume_files', selectindex, [0, 2]), 'in_file')]),
                     ])
    # combine tissue classes to noise mask
    wmcsf_mask = Node(fsl.BinaryMaths(operation='add',
                                      out_file='wmcsf_mask.nii'),
                      name='wmcsf_mask')
    denoise.connect([(binarize_tissue, wmcsf_mask, [(('out_file', selectsingle, 0), 'in_file'),
                                                    (('out_file', selectsingle, 1), 'operand_file')]),
                     (wmcsf_mask, outputnode, [('out_file', 'wmcsf_mask')])])
        
    
  
    
    # project wm_csf mask from anatomical to original epi space using inverse FLIRT-matrix
    invmat = Node(fsl.ConvertXFM(),
                         name='invmat')
    invmat.inputs.invert_xfm = True
    
    apply_inv=Node(fsl.ApplyXfm(), name='apply_inv')
    apply_inv.inputs.apply_xfm=True               
    denoise.connect([(inputnode, invmat, [('flirt_mat', 'in_file')]),
                     (invmat, apply_inv, [('out_file', 'in_matrix_file')]),
                     (inputnode, apply_inv, [('unwarped_mean', 'reference')]),
                     (wmcsf_mask, apply_inv, [('out_file', 'in_file')]),
                     (apply_inv, outputnode, [('out_file', 'wmcsf_mask2epi')])
                     ])
     #project brain to epi space as a checkup
    apply_inv_brain = Node(fsl.ApplyXfm(), name='apply_inv_brain')
    apply_inv_brain.inputs.apply_xfm=True
    denoise.connect([
                     (invmat, apply_inv_brain, [('out_file', 'in_matrix_file')]),
                     (inputnode, apply_inv_brain, [('unwarped_mean', 'reference')]),
                     (inputnode, apply_inv_brain, [('anat_brain', 'in_file')]),
                     (apply_inv_brain, outputnode, [('out_file', 'brain2epi')])
                     ])
                     
                     
    #no artifact detection and motion regression done because of AROMA
    
    # create filter with compcor components
    createfilter2 = Node(util.Function(input_names=['realigned_file', 'mask_file',
                                                    'num_components',
                                                    'extra_regressors'],
                                       output_names=['out_files'],
                                       function=extract_noise_components),
                         name='makecompcorfilter')
    createfilter2.inputs.num_components = 6
    createfilter2.inputs.extra_regressors=None
    createfilter2.plugin_args = {'submit_specs': 'request_memory = 17000'}
    denoise.connect([
                     (inputnode, createfilter2, [('epi_coreg', 'realigned_file')]),
                     (apply_inv, createfilter2, [('out_file', 'mask_file')]),
                     (createfilter2, outputnode, [('out_files', 'comp_regressor')]),
                     ])
    # regress compcor and other noise components
    filter2 = Node(fsl.GLM(out_f_name='F_noise.nii.gz',
                           out_pf_name='pF_noise.nii.gz',
                           out_res_name='rest2anat_denoised.nii.gz',
                           output_type='NIFTI_GZ',
                           demean=True),
                   name='filternoise')
    filter2.plugin_args = {'submit_specs': 'request_memory = 17000'}
    denoise.connect([(inputnode, filter2, [('epi_coreg', 'in_file')]),
                     (createfilter2, filter2, [('out_files', 'design')]),
                     (inputnode, filter2, [('brain_mask', 'mask')]),
                     (filter2, outputnode, [('out_f', 'comp_F'),
                                            ('out_pf', 'comp_pF'),
                                            ('out_file', 'out_betas')
                                            ])
                     ])



    # write TR into header again (glms remove it)
    # do not use mri_convert interface as it has a bug (already fixed in niyppe master)
    fix_tr = Node(util.Function(input_names=['in_file', 'TR_sec'], output_names=['out_file'], function=fix_TR_fs),
                  name='fix_tr')
    denoise.connect(inputnode, 'tr', fix_tr, 'TR_sec')
    denoise.connect(filter2, 'out_res', fix_tr, 'in_file')



    #use only highpass filter (because high-frequency content (otherwise filtered by lowpass is already considered in AROMA)) 
    highpass_filter = Node(fsl.TemporalFilter(out_file='rest_denoised_highpassed.nii'),
                           name='highpass_filter')
    highpass_filter.plugin_args = {'submit_specs': 'request_memory = 17000'}
    denoise.connect([(inputnode, highpass_filter, [('highpass_sigma', 'highpass_sigma')
                                                   ]),
                     (fix_tr, highpass_filter, [('out_file', 'in_file')]),
                     (fix_tr, outputnode, [('out_file', 'ts_fullspectrum')])
                     ])
    
    # time-normalize scans (could be set to percent change etc.  but here NO normalization is used
    #                 http://nipy.org/nitime/api/generated/nitime.fmri.io.html)
    normalize_time = Node(util.Function(input_names=['in_file', 'tr'],
                                        output_names=['out_file'],
                                        function=time_normalizer),
                          name='normalize_time')
    normalize_time.plugin_args = {'submit_specs': 'request_memory = 17000'}
    denoise.connect([(inputnode, normalize_time, [('tr', 'tr')]),
                     (highpass_filter, normalize_time, [('out_file', 'in_file')]),
                     (normalize_time, outputnode, [('out_file', 'normalized_file')])
                     ])
    return denoise
