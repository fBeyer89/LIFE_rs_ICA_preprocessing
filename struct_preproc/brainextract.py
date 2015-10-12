# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 11:19:21 2015

@author: fbeyer
"""

from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.freesurfer as fs
import nipype.interfaces.fsl as fsl
'''
Workflow to extract relevant output from freesurfer directory
'''
def create_brainextract_pipeline(name='brainextract'):
    # workflow
    brainextract = Workflow(name='brainextract')
    #inputnode
    inputnode=Node(util.IdentityInterface(fields=['anat', 'fraction']),
                   name='inputnode')
    #outputnode
    outputnode=Node(util.IdentityInterface(fields=['anat_brain', 'anat_brain_mask']),
                    name='outputnode')
    #use bet brain extraction
    bet = Node(interface=fsl.BET(mask=True),
               name = 'bet')
  
    # connections
    brainextract.connect([(inputnode, bet, [('anat','in_file'),
    ('fraction', 'frac')]),
    (bet, outputnode, [('out_file', 'anat_brain')]),
    (bet, outputnode, [('mask_file', 'anat_brain_mask')])
    ])
    
    return brainextract
