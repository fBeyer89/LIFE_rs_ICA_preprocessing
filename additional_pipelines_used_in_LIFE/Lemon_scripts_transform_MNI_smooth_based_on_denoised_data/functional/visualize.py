# -*- coding: utf-8 -*-
"""
Created on Wed Feb 11 08:58:48 2015

@author: fbeyer
"""

from nipype.pipeline.engine import MapNode, Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.ants as ants
import nipype.interfaces.utility as util
'''
Workflow to smooth functional image
'''
def create_visualize_pipeline(name='visualize'):

    # initiate workflow
    visualize = Workflow(name='visualize')
    # inputnode
    inputnode=Node(util.IdentityInterface(fields=['ts_transformed',
    'mni_template'
    ]),
    name='inputnode')
    # outputnode
    outputnode=Node(util.IdentityInterface(fields=['output_image'
    ]),
    name='outputnode')
    
    
    #apply smoothing
    slicer = Node(fsl.Slicer(sample_axial=6, image_width=750),name = 'smooth')
   
    
    visualize.connect([
    (inputnode, slicer, [('ts_transformed', 'in_file'),('mni_template', 'image_edges')]),     
    (slicer, outputnode,[('out_file', 'output_image')])
    ])
    
 
   
    return visualize