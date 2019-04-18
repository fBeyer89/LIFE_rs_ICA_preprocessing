# -*- coding: utf-8 -*-
"""
Created on Mon Feb  9 12:26:20 2015

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
from denoise_for_aroma import create_denoise_pipeline
from smoothing import create_smoothing_pipeline
from fmap_coreg import create_fmap_pipeline
import ICA_AROMA_functions as aromafunc


'''
Workflow to run 
-recommended preprocessing for AROMA,
(removal of five volumes, motion correction, fieldmap correction, 
intensity normalization, smoothing and calculation of registration parameters to MNI space)
-ICA-AROMA denoising and 
-postprocessing 
(regression of wm/csf noise using compcor, highpass-filtering
and registration to MNI-space using precomputed parameters)

for pre/postprocessing steps and general idea see http://www.sciencedirect.com/science/article/pii/S1053811915001822.
for ICA-AROMA code see: https://github.com/rhr-pruim/ICA-AROMA
====================================================
'''
def create_aroma_prep(subject, working_dir, data_dir, freesurfer_dir, out_dir,
    vol_to_remove, TR, epi_resolution, highpass,
    echo_space, te_diff, pe_dir, standard_brain, standard_brain_resampled, standard_brain_mask, 
    standard_brain_mask_resampled, fwhm_smoothing):
    # set fsl output type to nii.gz
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')
    # main workflow
    aroma_prep = Workflow(name='aroma_prep')
    aroma_prep.base_dir = working_dir
    aroma_prep.config['execution']['crashdump_dir'] = aroma_prep.base_dir + "/crash_files"
    

    #helper function to save array output of AROMA functions
    def small_save(filename, in_data):
        import numpy as np
        np.save(filename,in_data)
        return filename    
    
    # select files    
    templates={
    #'anat_head' : subject+'/preprocessed/mod/anat/T1.nii.gz', 
    'anat_brain' : subject+'/preprocessed/mod/anat/brain.nii.gz', 
    'brain_mask' : subject+'/preprocessed/mod/anat/T1_brain_mask.nii.gz', 
    'func': subject+'/func/EPI_t2.nii',
    'mag': subject+'/unwarp/B0_mag.nii',
    'phase': subject+'/unwarp/B0_ph.nii'
    }

    
    selectfiles = Node(nio.SelectFiles(templates,
    base_directory=data_dir),
    name="selectfiles")
    
    ##################preprocessing############################################
    # node to remove first volumes
    remove_vol = Node(util.Function(input_names=['in_file','t_min'],
    output_names=["out_file"],
    function=strip_rois_func),
    name='remove_vol')
    remove_vol.inputs.t_min = vol_to_remove    
    
    aroma_prep.connect([
    
    (selectfiles, remove_vol,    [('func', 'in_file')])])  
       
    #motion correction
    moco = create_moco_pipeline()
   
    #creating and applying the fieldmap
    fmap_coreg=create_fmap_pipeline()
    fmap_coreg.inputs.inputnode.echo_space = echo_space
    fmap_coreg.inputs.inputnode.te_diff = te_diff
    fmap_coreg.inputs.inputnode.pe_dir = pe_dir
    
    aroma_prep.connect([
    (selectfiles, fmap_coreg, [('mag', 'inputnode.mag')]),
    (selectfiles, fmap_coreg, [('phase', 'inputnode.phase')]),
    (moco, fmap_coreg, [('outputnode.epi_moco', 'inputnode.epi_coreg')]),
    (moco, fmap_coreg, [('outputnode.epi_mean','inputnode.epi_mean' )])
    
    ])
    
    #reorient to std
    reorient2std=Node(fsl.Reorient2Std(), name="reorient2std")
    #mean intensity normalization
    meanintensnorm = Node(fsl.ImageMaths(op_string= '-ing 10000'), name='meanintensnorm')
    
    aroma_prep.connect([
    #(unwarp, moco,    [('warped_file', 'inputnode.epi')]),
    (remove_vol, moco, [('out_file', 'inputnode.epi')]),
    (selectfiles, reorient2std, [('anat_brain', 'in_file')]),#reorient to standard to avoid registration issues seen previously
    (fmap_coreg, meanintensnorm, [('outputnode.unwarped_epi', 'in_file')])
    ])        
    
    
    #mask functional image
    betfunctional=Node(fsl.BET(frac=0.3), name='betfunctional')
    binmask = Node(fsl.ImageMaths(op_string= '-bin'), name='binmask')
    
    #smoothing (@6mm FWHM)
    smoothing = create_smoothing_pipeline() 
    smoothing.inputs.inputnode.fwhm=fwhm_smoothing
    
    
    aroma_prep.connect([
    (meanintensnorm, smoothing,  [('out_file', 'inputnode.ts_transformed')]),
    (fmap_coreg, betfunctional, [('outputnode.unwarped_mean_epi2fmap','in_file' )]),
    (betfunctional, binmask, [('out_file', 'in_file')])])
    
    # Func > Anat
    # register example func to high-resolution (use linear registration with 7 degrees of freedom and output
    #matrix example_func2highres.mat
    flirt = Node(fsl.FLIRT(cost_func='mutualinfo', interp='trilinear'), name='flirt')
    flirt.inputs.dof=7
#
    # Anat > Standard
    # register high-resolution to standard template - ###flirt## (as preparation for fnirt)
    flirt_prep = Node(fsl.FLIRT(cost_func='mutualinfo', interp='trilinear'), name='flirt_prep')
    flirt_prep.inputs.reference=standard_brain
    flirt_prep.inputs.interp='trilinear'
    flirt_prep.inputs.dof=12   
   
    fnirt=Node(fsl.FNIRT(), name='fnirt')
    fnirt.inputs.ref_file=standard_brain
    fnirt.inputs.field_file=True
    fnirt.inputs.fieldcoeff_file=True
    
    
    aroma_prep.connect([   
    (reorient2std, flirt,  [('out_file', 'reference')]), 
    (betfunctional, flirt, [('out_file', 'in_file')]),  
    (reorient2std, flirt_prep,  [('out_file', 'in_file')]),    
    (flirt_prep, fnirt,    [('out_matrix_file', 'affine_file')]),          
    (reorient2std,fnirt, [('out_file', 'in_file')]),
    ])

    ##################ICA-AROMA###############################################
    #Step 1) MELODIC
    runICA = Node(name="runICA", interface=Function(input_names=["fslDir", "outDir", "inFile", "melDirIn", "mask", "dim", "TR"],
                  output_names=["mdir", 'melICmix', 'melodic_FTmix'],
                  function=aromafunc.runICA))  
    runICA.inputs.fslDir= os.path.join(os.environ["FSLDIR"],'bin','')
    runICA.inputs.dim=0 #automatically estimates network number using MDL
    runICA.inputs.TR=2.0
    runICA.inputs.melDirIn=""
    runICA.inputs.outDir=out_dir 

    #Step 2) Automatic classification of the components
    #  - registering the spatial maps to MNI
    regMNI= Node(name="regMNI", interface=Function(input_names=["fslDir", "inFile", "outFile", "affmat", "warp"],
                              output_names=['melodic_IC_MNI2mm'],
            function=aromafunc.register2MNI))
            
    regMNI.inputs.fslDir= os.path.join(os.environ["FSLDIR"],'bin','')
    regMNI.inputs.outFile=out_dir+'melodic.ica/melodic_IC_thr_MNI2mm.nii.gz'
    
    
    #connect inputs to Melodic-ICA
    aroma_prep.connect([(binmask, runICA, [('out_file', 'mask')]),
    (smoothing, runICA, [('outputnode.ts_smoothed', 'inFile')]),
    #connect inputs to registration node
    (runICA, regMNI, [('mdir', 'inFile')]),
    (flirt,regMNI, [('out_matrix_file', 'affmat')]),
    (fnirt, regMNI, [('fieldcoeff_file', 'warp')])
    ])
     
    #extracting the Maximum RP correlation feature
    feature_time_series=Node(name="feature_time_series", interface=Function(input_names=["melmix", "mc"],
                             output_names=["maxRPcorr"],
                            function=aromafunc.feature_time_series))
    
    save_featts=Node(name="save_featts", interface=Function(input_names=["filename", "in_data"],
                     output_names=["filename"],
                     function=small_save))
    
    save_featts.inputs.filename=working_dir+'/aroma_prep/save_featts/maxRPcorr.npy' 

    aroma_prep.connect([
    (runICA,feature_time_series, [('melICmix', 'melmix')])
    ])

     #connect moco to time_series features
    aroma_prep.connect([(moco, feature_time_series, [('outputnode.par_moco', 'mc')]),
                        (feature_time_series, save_featts, [("maxRPcorr", "in_data")])
                        ])

    #extracting the High-frequency content feature                  
    feature_freq=Node(name="feature_freq", interface=Function(input_names=["melFTmix", "TR"],
                       output_names=["HFC"],
                       function=aromafunc.feature_frequency))
    feature_freq.inputs.TR=2.0
    
    
    save_featfreq=Node(name="save_featfreq", interface=Function(input_names=["filename", "in_data"],
                     output_names=["filename"],
                     function=small_save))
    save_featfreq.inputs.filename=working_dir+'/aroma_prep/feature_freq/HFC.npy' 
    
    aroma_prep.connect([
    (runICA, feature_freq, [('melodic_FTmix', 'melFTmix')]),
    (feature_freq, save_featfreq, [('HFC', 'in_data')])
    ])

    #extracting the CSF & Edge fraction features    
    feature_spatial=Node(name="feature_spatial", interface=Function(input_names=["fslDir", "tempDir", "aromaDir", "melIC"],
                                                                    output_names=["edgeFract", "csfFract"],
    function=aromafunc.feature_spatial))
    
    feature_spatial.inputs.fslDir= os.path.join(os.environ["FSLDIR"],'bin','')
    feature_spatial.inputs.tempDir=working_dir+"/aroma_prep/feature_spatial/"
    feature_spatial.inputs.aromaDir="/home/raid1/fbeyer/Documents/Scripts/ICA-AROMA/"
    
        
    save_featsp_edge=Node(name="save_featsp_edge", interface=Function(input_names=["filename", "in_data"],
                     output_names=["filename"],
                     function=small_save))
    save_featsp_edge.inputs.filename=working_dir+'/aroma_prep/feature_spatial/edge.npy'
 
    save_featsp_csf=Node(name="save_featsp_csf", interface=Function(input_names=["filename", "in_data"],
                     output_names=["filename"],
                     function=small_save))
    save_featsp_csf.inputs.filename=working_dir+'/aroma_prep/feature_spatial/csf.npy'
    
    aroma_prep.connect([
    (regMNI, feature_spatial, [('melodic_IC_MNI2mm','melIC')]),
    (feature_spatial, save_featsp_edge, [('edgeFract', 'in_data')]),
    (feature_spatial, save_featsp_csf, [('csfFract', 'in_data')])
    ])   

    #classification of features using predefined feature space
    classification=Node(name="classification", interface=Function(input_names=["outDir", "maxRPcorr", "edgeFract", "HFC", "csfFract"],
                        output_names=["motionICs"],
    function=aromafunc.classification)) 
    classification.inputs.outDir=out_dir+'/melodic.ica/'
    
    save_class_mic=Node(name="save_class_mic", interface=Function(input_names=["filename", "in_data"],
                     output_names=["filename"],
                     function=small_save))
    save_class_mic.inputs.filename=working_dir+'/aroma_prep/classification/motionICs.npy'
    
    aroma_prep.connect([
    (classification, save_class_mic, [('motionICs', 'in_data')])
    ])

    #connections of classification with all features
    aroma_prep.connect([    
    (feature_time_series, classification, [('maxRPcorr', 'maxRPcorr')]),
    (feature_spatial, classification, [('edgeFract', 'edgeFract')]),
    (feature_spatial, classification, [('csfFract', 'csfFract')]),
    (feature_freq, classification, [('HFC', 'HFC')])
    ])

    #Step 3) Data denoising' (using input data (=intensity normalized, motion corrected, smoothed fMRI in subject space))
    denoising_ICA=Node(name="denoising_ICA", interface=Function(input_names=["fslDir", "inFile", "outDir", "melmix", "denType", "denIdx"],
                    output_names=["denoised_func"],
                    function=aromafunc.denoising))
    denoising_ICA.inputs.fslDir= os.path.join(os.environ["FSLDIR"],'bin','')
    denoising_ICA.inputs.outDir=out_dir
    denoising_ICA.inputs.melmix=out_dir+'/melodic.ica/melodic_mix'
    denoising_ICA.inputs.denType=2 #1=aggr, 2=non-aggr, 3=both
     
    aroma_prep.connect([
    (classification, denoising_ICA, [('motionICs', 'denIdx')]),
    (smoothing, denoising_ICA, [('outputnode.ts_smoothed', 'inFile')])
    ])


    ##################post-processing###########################################
    #WM/CSF/linear trend removal using regression of compcor-components
    #highpass-filtering here!!
    #registration of denoised image to MNI space
    postprocess=create_denoise_pipeline()
    postprocess.inputs.inputnode.highpass_sigma = 1. / (2 * TR * highpass)
    # https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=ind1205&L=FSL&P=R57592&1=FSL&9=A&I=-3&J=on&d=No+Match%3BMatch%3BMatches&z=4
    postprocess.inputs.inputnode.tr = TR
    
    aroma_prep.connect([
    (binmask, postprocess,            [('out_file', 'inputnode.brain_mask')]),                              
    (denoising_ICA, postprocess,         [('denoised_func', 'inputnode.epi_coreg')]),
    (moco, postprocess,             [('outputnode.epi_mean', 'inputnode.unwarped_mean')]),
    (flirt, postprocess,            [('out_matrix_file', 'inputnode.flirt_mat')]),
    (reorient2std, postprocess,      [('out_file', 'inputnode.anat_brain')])
    ])   

    #register only-AROMA denoised file	into MNI space (without CSF/WM removal + highpass-filtering)
    apply_TF= Node(fsl.ApplyWarp(), name="apply_TF")
    apply_TF.inputs.ref_file = standard_brain

    #register AROMA and CSF/WM denoised file	into MNI space
    apply_TF_denoised= Node(fsl.ApplyWarp(), name="apply_TF_denoised")
    apply_TF_denoised.inputs.ref_file = standard_brain

    aroma_prep.connect([
    (postprocess, apply_TF_denoised, [('outputnode.normalized_file', 'in_file')]), #with bandpass-filtering
    (flirt, apply_TF_denoised, [('out_matrix_file', 'premat')]),
    (fnirt, apply_TF_denoised, [('fieldcoeff_file', 'field_file')]),
    (denoising_ICA, apply_TF, [('denoised_func', 'in_file')]),
    (flirt, apply_TF, [('out_matrix_file', 'premat')]),
    (fnirt, apply_TF, [('fieldcoeff_file', 'field_file')]),
    ]) 

    #sink to store files
    sink = Node(nio.DataSink(parameterization=False,
    base_directory=out_dir,
    substitutions=[
    ('brain_reoriented_warped.nii', 'brain2mni_warp_fnirt.nii'),
    ('brain_reoriented_flirt.nii', 'brain2mni_aff_flirt.nii'),
    ('brain_reoriented_fieldwarp.nii', 'brain2mni_warpcoeff_cout.nii'),
    ('brain_reoriented_field.nii', 'brain2mni_warpfield_fout.nii'),
    ('rest2anat_maths_smooth', 'func_preproc_smoothed'),
    ('rest_mean2fmap_unwarped_brain_maths', 'func_brain_mask'),   
    ('denoised_func_data_nonaggr_warp', 'aroma_denoised_MNI'),
    ('rest_denoised_highpassed_norm_warp.nii', 'aroma_csfwm_denoised_MNI.nii'),
    ('rest2anat_denoised', 'rest_denoised_fullspectrum'),
    ('rest_denoised_highpassed_norm', 'rest_denoised_highpass'),
    ('wmcsf_mask_lowres_flirt', 'wmcsf_mask2epi'),
    ('brain2mni_aff_flirt', 'brain2epi')]
    ),
    name='sink')
    
    
    ##all the output
    aroma_prep.connect([
    (moco, sink,[
    ('outputnode.epi_moco', 'realign.@realigned_ts'),
    ('outputnode.par_moco', 'realign.@par'),
    ('outputnode.rms_moco', 'realign.@rms'),
    ('outputnode.mat_moco', 'realign.MAT.@mat'),
    ('outputnode.epi_mean', 'realign.@mean'),
    ('outputnode.rotplot', 'realign.plots.@rotplot'),
    ('outputnode.transplot', 'realign.plots.@transplot'),
    ('outputnode.dispplots', 'realign.plots.@dispplots'),
    ('outputnode.tsnr_file', 'realign.@tsnr')]),
    (fmap_coreg, sink, [('outputnode.fmap', 'fmap.transforms2anat.@fmap'),
                        ('outputnode.unwarped_mean_epi2fmap', 'fmap.@unwarped_mean_epi2fmap'),
                        ('outputnode.epi2fmap', 'fmap.@epi2fmap'),
                        ('outputnode.unwarpfield_epi2fmap', 'fmap.@fmap_fullwarp'),
                        ('outputnode.mean_from_unwarped_epi', 'fmap.@mean_from_unwarped')
                        ]),                            
    (binmask, sink, [('out_file', 'aroma_inputs.@brainmask')]),
    (smoothing, sink, [('outputnode.ts_smoothed', 'aroma_inputs.@FWHM6')]),
    (flirt, sink,    [('out_file', 'aroma_inputs.reg.affine.@result')]),
    (flirt, sink,    [('out_matrix_file', 'aroma_inputs.reg.affine.@matrix')]),
    (fnirt, sink,    [('warped_file', 'aroma_inputs.reg.fnirt.@image' )]),
    (fnirt, sink,    [('fieldcoeff_file', 'aroma_inputs.reg.fnirt.@fieldcoeff' )]),
    (fnirt, sink,    [('field_file', 'aroma_inputs.reg.fnirt.@field' )]),
    (flirt_prep, sink,    [('out_matrix_file','aroma_inputs.reg.fnirt.@mat')]),
    (flirt_prep, sink,    [('out_file', 'aroma_inputs.reg.fnirt.@resultflirt')]),
    (save_featts, sink, [('filename', 'aroma_res.features.@rescorr')]),
    (save_featsp_edge, sink, [('filename', 'aroma_res.features.@edge')]),
    (save_featsp_csf, sink, [('filename', 'aroma_res.features.@csf')]),
    (save_featfreq, sink,  [('filename', 'aroma_res.features.@HFC')]),
    (denoising_ICA, sink,      [('denoised_func', 'aroma_res.@denoised_func')]),
    (postprocess, sink, [
            ('outputnode.wmcsf_mask', 'denoise.mask.@wmcsf_masks'),
            ('outputnode.combined_motion', 'denoise.artefact.@combined_motion'),
            ('outputnode.comp_regressor', 'denoise.regress.@comp_regressor'),
            ('outputnode.comp_F', 'denoise.regress.@comp_F'),
            ('outputnode.comp_pF', 'denoise.regress.@comp_pF'),
            ('outputnode.brain2epi', 'denoise.mask.@brain2epi'),
            ('outputnode.wmcsf_mask2epi', 'denoise.mask.@wmcsf_mask2epi'),
            ('outputnode.normalized_file', 'denoise.@normalized'),
            # FL added fullspectrum
            ('outputnode.ts_fullspectrum', 'denoise.@ts_fullspectrum')]),
    (apply_TF, sink,              [('out_file', 'transforms.withoutCSFWM.@denoised_result_MNI')]), #with only motion
    (apply_TF_denoised, sink,              [('out_file', 'transforms.withCSFWM.@denoised_result_MNI')]) #with motion, compcor CSF/WM + highpass
    ])


    aroma_prep.write_graph(dotfilename='aroma_prep.dot', graph2use='colored', format='pdf', simple_form=True)
    aroma_prep.run(plugin='MultiProc' )#plugin='CondorDAGMan', plugin_args={'initial_specs': 'request_memory = 1500'}
    

