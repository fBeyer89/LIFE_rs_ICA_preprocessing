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
def create_smoothing_pipeline(name='smoothing'):
    # set fsl output type
    fsl.FSLCommand.set_default_output_type('NIFTI')
    # initiate workflow
    smoothing = Workflow(name='smoothing')
    # inputnode
    inputnode=Node(util.IdentityInterface(fields=['ts_transformed',
    'fwhm'
    ]),
    name='inputnode')
    # outputnode
    outputnode=Node(util.IdentityInterface(fields=['ts_smoothed'
    ]),
    name='outputnode')
    
    
    #apply smoothing
    smooth = Node(fsl.Smooth(),name = 'smooth')
   
    
    smoothing.connect([
    (inputnode, smooth, [
    ('ts_transformed', 'in_file'),
    ('fwhm', 'fwhm')]
    ), 
    (smooth, outputnode, [('smoothed_file', 'ts_smoothed')]
    )
    ])
    
 



    
    return smoothing