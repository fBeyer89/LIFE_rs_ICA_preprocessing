# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:35:16 2016

@author: fbeyer
"""
from nipype.pipeline.engine import Node, Workflow, MapNode
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util

'''
Workflow for correcting distortion using fieldmaps
using mean_epi from motion correction as input,
coregistration with magnitude image of fieldmap,
calculating of shiftmap and applying it to moco-mean and timeseries 
calculating new mean image after unwarping
see:https://lcni.uoregon.edu/kb-articles/kb-0003
'''

def create_fmap_pipeline(name='fmap_coreg'):
    
    fmap_coreg = Workflow(name='fmap_coreg')

    # inputnode
    inputnode = Node(util.IdentityInterface(fields=['epi_mean',
                                                    'epi_coreg',
                                                    'mag',
                                                    'phase',
                                                    'echo_space',
                                                    'te_diff',
                                                    'pe_dir'
                                                    ]),
                     name='inputnode')
    # outputnode
    outputnode = Node(util.IdentityInterface(fields=['epi2fmap',
                                                     'fmap',
                                                     'unwarpfield_epi2fmap',
                                                     'unwarped_mean_epi2fmap',
                                                     'unwarped_epi',                                                     
                                                     'mean_from_unwarped_epi'
                                                    
                                                     ]),
                      name='outputnode')
    
    #correct for distortions using fieldmap    
    #### prepare fieldmap ####
    # split first magnitude image from mag input
    splitmag = Node(fsl.ExtractROI(t_min=0,  t_size=1),    name='splitmag')
    fmap_coreg.connect(inputnode, 'mag', splitmag, 'in_file')
    
    # strip magnitude image and erode even further
    bet = Node(fsl.BET(frac=0.5,    mask=True),    name='bet')
    fmap_coreg.connect(splitmag,'roi_file', bet,'in_file')   
    erode = Node(fsl.maths.ErodeImage(kernel_shape='sphere',    kernel_size=3,    args=''),    name='erode')
    fmap_coreg.connect(bet,'out_file', erode, 'in_file')
    
    # prepare fieldmap
    prep_fmap = Node(fsl.epi.PrepareFieldmap(),    name='prep_fmap')
    fmap_coreg.connect([(erode, prep_fmap, [('out_file', 'in_magnitude')]),   
                                            (inputnode, prep_fmap, [('phase', 'in_phase'),('te_diff', 'delta_TE')]),
                                            (prep_fmap, outputnode, [('out_fieldmap','fieldmap')])
    ])
   
   
   #### unmask fieldmap ####
    fmap_mask = Node(fsl.maths.MathsCommand(args='-abs -bin'),
    name='fmap_mask')
    
    unmask = Node(fsl.FUGUE(save_unmasked_fmap=True),
    name='unmask')
    fmap_coreg.connect([(prep_fmap, fmap_mask, [('out_fieldmap', 'in_file')]),
                        (fmap_mask, unmask, [('out_file', 'mask_file')]),
                        (prep_fmap, unmask,[('out_fieldmap','fmap_in_file')]),
                        (inputnode, unmask, [('pe_dir', 'unwarp_direction')]) 
                        ])
                        
    #### register epi to fieldmap ####
    epi2fmap = Node(fsl.FLIRT(dof=6,
    out_file='rest_mean2fmap.nii.gz',
    interp='spline'),
    name='epi2fmap')
    fmap_coreg.connect(
    [(inputnode,epi2fmap,[('epi_mean', 'in_file')]),
    (splitmag, epi2fmap, [('roi_file', 'reference')]),
    (epi2fmap, outputnode, [('out_file', 'epi2fmap')])
    ])
    
    #### first create shiftmap from fieldmap ####
    unwarp = Node(fsl.FUGUE(save_shift=True),
                  name='unwarp')
    fmap_coreg.connect([(epi2fmap, unwarp, [('out_file', 'in_file')]),
                        (unmask, unwarp, [('fmap_out_file', 'fmap_in_file')]),
                        (fmap_mask, unwarp, [('out_file', 'mask_file')]),
                        (inputnode, unwarp, [('echo_space', 'dwell_time'),
                                             ('pe_dir', 'unwarp_direction')]),
                        (unwarp, outputnode, [('shift_out_file', 'shiftmap')])
                        ])
    


    #### make warpfield combining func2fieldmap-FLIRT and shiftmap and apply to the mean image####
    convertwarp0 = Node(fsl.utils.ConvertWarp(out_relwarp=True,
                                              out_file='rest_mean2fmap_unwarpfield.nii.gz'),
                        name='convertwarp0')
    applywarp0 = Node(fsl.ApplyWarp(interp='spline',
                                    relwarp=True,
                                    out_file='rest_mean2fmap_unwarped.nii.gz',
                                    datatype='float'),
                      name='applywarp0')

    fmap_coreg.connect([(splitmag, convertwarp0, [('roi_file', 'reference')]),
                        (epi2fmap, convertwarp0, [('out_matrix_file', 'premat')]),#here the registration from mean epi to fieldmap is included
                        (unwarp, convertwarp0, [('shift_out_file', 'shift_in_file')]),
                        (inputnode, convertwarp0, [('pe_dir', 'shift_direction')]),
                        (inputnode, applywarp0, [('epi_mean', 'in_file')]),
                        (splitmag, applywarp0, [('roi_file', 'ref_file')]),
                        (convertwarp0, applywarp0, [('out_file', 'field_file')]),
                        (convertwarp0, outputnode, [('out_file', 'unwarpfield_epi2fmap')]),
                        (applywarp0, outputnode, [('out_file', 'unwarped_mean_epi2fmap')])
                        ])

    # apply warpfield to motion corrected volumes
    applywarp = Node(fsl.ApplyWarp(interp='spline',
                                      relwarp=True,
                                      out_file='rest2anat.nii.gz',
                                      datatype='float'),
                        #iterfield=['in_file'],
                        name='applywarp')
    fmap_coreg.connect([
                        (inputnode, applywarp, [('epi_coreg', 'in_file')]),
                        (convertwarp0, applywarp,
                           [('out_file', 'field_file')]),
                        (splitmag, applywarp, [('roi_file', 'ref_file')]),
                        (applywarp, outputnode, [('out_file', 'unwarped_epi')])
                        ])

    
    # calculate new mean
    tmean = Node(fsl.maths.MeanImage(dimension='T',
                                     out_file='mean_from_unwarped_epi.nii.gz'),
                 name='tmean')
    fmap_coreg.connect([(applywarp, tmean, [('out_file', 'in_file')]),
                          (tmean, outputnode, [('out_file', 'mean_from_unwarped_epi')])
                          ])                    
                        
                        
    return fmap_coreg