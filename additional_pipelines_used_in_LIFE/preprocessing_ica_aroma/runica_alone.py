# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 14:30:53 2016

@author: fbeyer
"""
import os
from nipype.pipeline.engine import Node, Workflow
from nipype.interfaces.utility import Function
import nipype.interfaces.utility as util
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
from strip_rois import strip_rois_func
from moco import create_moco_pipeline
#from ICA_AROMA2 import create_ica_aroma
from smoothing import create_smoothing_pipeline
import ICA_AROMA_functions as aromafunc
    
    
runICA = Function(input_names=["fslDir", "outDir", "inFile", "melDirIn", "mask", "dim", "TR"],
                  output_names=["mdir"],
                  function=aromafunc.runICA)  
runICA.inputs.fslDir= '/usr/share/fsl/5.0/bin/'
#os.path.join(os.environ["FSLDIR"],'bin','')
runICA.inputs.dim=0
runICA.inputs.TR=2
runICA.inputs.melDirIn=""
runICA.inputs.outDir="/scr/lessing2/data_fbeyer/FTO_YFAS/WDR/"
runICA.inputs.inFile="/scr/lessing2/data_fbeyer/FTO_YFAS/Subjects/LI00037838/aroma_inputs/func_preproc_smoothed.nii"
runICA.inputs.mask="/scr/lessing2/data_fbeyer/FTO_YFAS/Subjects/LI00037838/aroma_inputs/func_brain_mask.nii"

runICA.run()