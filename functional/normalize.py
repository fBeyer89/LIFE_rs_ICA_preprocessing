# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:34:01 2015

@author: fbeyer
"""

from nipype.pipeline.engine import Node, Workflow, MapNode
import nipype.interfaces.fsl as fsl
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.utility as util
import nipype.algorithms.rapidart as ra
import nipype.interfaces.afni as afni
from normalize_timeseries import time_normalizer
from nipype.utils.filemanip import list_to_filename
'''
Main workflow for denoising
Largely based on https://github.com/nipy/nipype/blob/master/examples/
rsfmri_vol_surface_preprocessing_nipy.py#L261
but denoising in anatomical space
'''
def create_normalize_pipeline(name='normalize'):
    # workflow
    normalize = Workflow(name='normalize')
    # Define nodes
    inputnode = Node(interface=util.IdentityInterface(fields=['epi_coreg',
    'tr']),
    name='inputnode')
    outputnode = Node(interface=util.IdentityInterface(fields=[
    'normalized_file']),
    name='outputnode')

    # time-normalize scans
    normalize_time=Node(util.Function(input_names=['in_file','tr'],
    output_names=['out_file'],
    function=time_normalizer),
    name='normalize_time')
    normalize_time.plugin_args={'submit_specs': 'request_memory = 17000'}
    normalize.connect([(inputnode, normalize_time, [('tr', 'tr')]),
    (inputnode, normalize_time, [('epi_coreg', 'in_file')]),
    (normalize_time, outputnode, [('out_file', 'normalized_file')])
    ])
    
    # time-normalize scans    
    
    return normalize