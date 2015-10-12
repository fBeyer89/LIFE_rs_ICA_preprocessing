# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 14:00:12 2015

@author: fbeyer
"""
from nipype.pipeline.engine import MapNode, Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants


def create_ants_registration_pipeline(name='ants_registration'):
    # set fsl output type
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
    # initiate workflow
    ants_registration = Workflow(name='ants_registration')
    # inputnode
    inputnode=Node(util.IdentityInterface(fields=['denoised_ts',
    'ants_affine',
    'ants_warp',
    'ref'
    ]),
    name='inputnode')
    # outputnode
    outputnode=Node(util.IdentityInterface(fields=['ants_reg_ts',
    ]),
    name='outputnode')

    #also transform to mni space
    collect_transforms = Node(interface = util.Merge(2),name='collect_transforms')    
    
    ants_reg = Node(ants.ApplyTransforms(input_image_type = 3, dimension = 3, interpolation = 'Linear'), name='ants_reg')
    
    
    
    
    ants_registration.connect([
                          (inputnode, ants_reg, [('denoised_ts', 'input_image')]),
                          (inputnode, ants_reg, [('ref', 'reference_image')]),
                          (inputnode, collect_transforms, [('ants_affine', 'in1')]),
                          (inputnode, collect_transforms, [('ants_warp', 'in2')]),
                          (collect_transforms, ants_reg,  [('out', 'transforms')]),
                          (ants_reg, outputnode, [('output_image', 'ants_reg_ts')])
                          ])
                          
    return ants_registration
    
    
