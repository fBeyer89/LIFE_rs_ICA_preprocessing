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
def create_slice_timing_pipeline(name='slicetiming'):
    # set fsl output type
    fsl.FSLCommand.set_default_output_type('NIFTI')
    # initiate workflow
    slicetiming = Workflow(name='slicetiming')
    # inputnode
    inputnode=Node(util.IdentityInterface(fields=['ts'
    ]),
    name='inputnode')
    # outputnode
    outputnode=Node(util.IdentityInterface(fields=['ts_slicetcorrected'
    ]),
    name='outputnode')
    
    
    #use FSL slicetiming (default ascending bottom to top)
    timer = Node(fsl.SliceTimer(),name = 'timer')
    timer.inputs.time_repetition=2.0
   
    
    slicetiming.connect([
    (inputnode, timer, [
    ('ts', 'in_file')]
    ), 
    (timer, outputnode, [('slice_time_corrected_file', 'ts_slicetcorrected')]
    )
    ])
    
 



    
    return slicetiming