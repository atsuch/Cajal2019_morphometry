#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:07:19 2019

A small workflow to perform surface-based morphometry analyses using Freesurfer
results and mri_glmfit

@author: tsuchida
"""

import os
import os.path as op
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces.freesurfer.utils import SurfaceSmooth, MRIsCalc
from nipype.interfaces.freesurfer.model import MRISPreproc, GLMFit
from ginnipi_tools.interfaces.custom import FDR
from ginnipi_tools.toolbox.flow import (getValFromDict, createTuple2, prependString,
                                        getElementFromList, getElementsFromTxtList)
from ginnipi_tools.toolbox.plotting_tools import plot_surf_stat, plot_surf_map



def genFreesurferSBMglmWF(name='fsSBM',
                          base_dir=op.abspath('.'),
                          group_sublist=[],
                          model_dir=None,
                          model_info={'Model1': 'dods'},
                          design_input='fsgd',
                          fs_subjects_dir='/data/analyses/work_in_progress/freesurfer/fsmrishare-flair6.0/',
                          fwhm=[0.0, 10.0],
                          measure_list=['thickness', 'area'],
                          target_atlas='fsaverage',
                          target_atlas_surfreg='sphere.reg',
                          correction_method='FDR'):
    
    wf = Workflow(name)
    wf.base_dir = base_dir

    
    # Node: model List
    modelList = Node(IdentityInterface(fields=['model_name'], mandatory_inputs=True),
                    name='modelList')
    modelList.iterables = ('model_name', list(model_info.keys()))

    
    # Grab fsgd or design mat and contrast mtx files from model_dir
    fileList_temp_args = {'contrast_files': [['model_name', '*.mtx']],
                          'contrast_sign_files':  [['model_name', '*.mdtx']]}
    if design_input == 'fsgd':
        fileList_temp_args['fsgd_file'] = [['model_name', '*.fsgd']]
    elif design_input == 'design_mat':
        fileList_temp_args['design_mat'] = [['model_name', 'X.mat']]
        
    fileList = Node(DataGrabber(infields=['model_name'],
                                outfields=list(fileList_temp_args.keys())), 
                    name="fileList")
    fileList.inputs.base_directory = model_dir
    fileList.inputs.ignore_exception = False
    fileList.inputs.raise_on_empty = True
    fileList.inputs.sort_filelist = True
    fileList.inputs.template = '%s/%s'
    fileList.inputs.template_args =  fileList_temp_args
    wf.connect(modelList, "model_name", fileList, "model_name")


    # preproc for each hemisphere to produce concatenated file for glmfit and 
    # also a mean map
    
    # Define a few other iterables
    measList = Node(IdentityInterface(fields=['measure'],
                                      mandatory_inputs=True),
                    name='measList')
    measList.iterables = ('measure', measure_list)
    
    smoothList = Node(IdentityInterface(fields=['fwhm'],
                                        mandatory_inputs=True),
                      name='smoothList')
    smoothList.iterables = ('fwhm', fwhm)
    
    surfaces = ['inflated', 'pial']
    plotSurfList = Node(IdentityInterface(fields=['surf']),
                      name='plotSurfList')
    plotSurfList.iterables = ('surf', surfaces)
    
    
    # MRI_preproc
    lhSBMpreproc = MapNode(MRISPreproc(),
                           name='lhSBMpreproc',
                           iterfield=['args', 'out_file'])
    lhSBMpreproc.inputs.subjects_dir = fs_subjects_dir
    lhSBMpreproc.inputs.target = target_atlas
    lhSBMpreproc.inputs.hemi = 'lh'
    lhSBMpreproc.inputs.args = ['', '--mean']
    lhSBMpreproc.inputs.out_file = ['{}.lh.{}.mgh'.format(out_name, target_atlas) for out_name in ['stacked', 'mean']]
    lhSBMpreproc.inputs.subjects = group_sublist
    wf.connect(measList, "measure", lhSBMpreproc, "surf_measure")
    
    rhSBMpreproc = MapNode(MRISPreproc(),
                           name='rhSBMpreproc',
                           iterfield=['args', 'out_file'])
    rhSBMpreproc.inputs.subjects_dir = fs_subjects_dir
    rhSBMpreproc.inputs.target = target_atlas
    rhSBMpreproc.inputs.hemi = 'rh'
    rhSBMpreproc.inputs.args = ['', '--mean']
    rhSBMpreproc.inputs.out_file = ['{}.rh.{}.mgh'.format(out_name, target_atlas) for out_name in ['stacked', 'mean']]
    rhSBMpreproc.inputs.subjects = group_sublist
    wf.connect(measList, "measure", rhSBMpreproc, "surf_measure")
    
    
    # Create smoothed mean maps for each non-zero fwhm
    non_zero_fwhm = [val for val in fwhm if val != 0.0]
    lhSmoothMean = MapNode(SurfaceSmooth(),
                           name='lhSmoothMean',
                           iterfield=['fwhm', 'out_file'])
    lhSmoothMean.inputs.subject_id = target_atlas
    lhSmoothMean.inputs.hemi = 'lh'
    lhSmoothMean.inputs.subjects_dir = fs_subjects_dir
    lhSmoothMean.inputs.fwhm = non_zero_fwhm
    lhSmoothMean.inputs.cortex = True
    lhSmoothMean.inputs.out_file = ['mean.lh.fwhm{}.{}.mgh'.format(str(int(val)), target_atlas) for val in non_zero_fwhm]
    wf.connect(lhSBMpreproc, ('out_file', getElementFromList, 1), lhSmoothMean, 'in_file')
    
    rhSmoothMean = MapNode(SurfaceSmooth(),
                           name='rhSmoothMean',
                           iterfield=['fwhm', 'out_file'])
    rhSmoothMean.inputs.subject_id = target_atlas
    rhSmoothMean.inputs.hemi = 'rh'
    rhSmoothMean.inputs.subjects_dir = fs_subjects_dir
    rhSmoothMean.inputs.fwhm = non_zero_fwhm
    rhSmoothMean.inputs.cortex = True
    rhSmoothMean.inputs.out_file = ['mean.rh.fwhm{}.{}.mgh'.format(str(int(val)), target_atlas) for val in non_zero_fwhm]
    wf.connect(rhSBMpreproc, ('out_file', getElementFromList, 1), rhSmoothMean, 'in_file')

    
    # For each concatenated surfaces produced by the SBMpreproc, run glmfit
    
    if correction_method == 'FDR':
        save_res = False
    elif correction_method == 'perm':
        save_res = True
    
    if design_input == 'fsgd': 
        fsgdInput = Node(Function(input_names=['item1', 'item2'],
                                  output_names=['out_tuple'],
                                  function=createTuple2),
                         name='fsgdInput')
        wf.connect(fileList, 'fsgd_file', fsgdInput, 'item1')
        wf.connect(modelList, ('model_name', getValFromDict, model_info),
                   fsgdInput, 'item2')
    
    lhSBMglmfit = Node(GLMFit(),
                       name='lhSBMglmfit')
    lhSBMglmfit.inputs.subjects_dir = fs_subjects_dir
    lhSBMglmfit.inputs.surf = True
    lhSBMglmfit.inputs.subject_id = target_atlas
    lhSBMglmfit.inputs.hemi = 'lh'
    lhSBMglmfit.inputs.cortex = True
    lhSBMglmfit.inputs.save_residual = save_res
    wf.connect(smoothList, 'fwhm', lhSBMglmfit, 'fwhm')
    wf.connect(lhSBMpreproc, ('out_file', getElementFromList, 0), lhSBMglmfit, 'in_file')
    if design_input == 'fsgd':
        wf.connect(fsgdInput, 'out_tuple', lhSBMglmfit, 'fsgd')
    elif design_input == 'design_mat':
        wf.connect(fileList, 'design_mat', lhSBMglmfit, 'design')
    wf.connect(fileList, 'contrast_files', lhSBMglmfit, 'contrast')
    
    rhSBMglmfit = Node(GLMFit(),
                       name='rhSBMglmfit')
    rhSBMglmfit.inputs.subjects_dir = fs_subjects_dir
    rhSBMglmfit.inputs.surf = True
    rhSBMglmfit.inputs.subject_id = target_atlas
    rhSBMglmfit.inputs.hemi = 'rh'
    rhSBMglmfit.inputs.cortex = True
    rhSBMglmfit.inputs.save_residual = save_res
    wf.connect(smoothList, 'fwhm', rhSBMglmfit, 'fwhm')
    wf.connect(rhSBMpreproc, ('out_file', getElementFromList, 0), rhSBMglmfit, 'in_file')
    if design_input == 'fsgd':
        wf.connect(fsgdInput, 'out_tuple', rhSBMglmfit, 'fsgd')
    elif design_input == 'design_mat':
        wf.connect(fileList, 'design_mat', rhSBMglmfit, 'design')
    wf.connect(fileList, 'contrast_files', rhSBMglmfit, 'contrast')


    # perfrom FDR correction if 'FDR' is chosen
    if correction_method == 'FDR':
        
        mriFDR = MapNode(FDR(),
                         iterfield=['in_file1', 'in_file2', 'fdr_sign'],
                         name='mriFDR')
        mriFDR.inputs.fdr = 0.05
        mriFDR.inputs.out_thr_file = 'fdr_threshold.txt'
        mriFDR.inputs.out_file1 = 'lh.sig_corr.mgh'
        mriFDR.inputs.out_file2 = 'rh.sig_corr.mgh'
        wf.connect(lhSBMglmfit, 'sig_file', mriFDR, 'in_file1')
        wf.connect(lhSBMglmfit, 'mask_file', mriFDR, 'in_mask1')
        wf.connect(rhSBMglmfit, 'sig_file', mriFDR, 'in_file2')
        wf.connect(rhSBMglmfit, 'mask_file', mriFDR, 'in_mask2')
        wf.connect(fileList, ('contrast_sign_files', getElementsFromTxtList),
                   mriFDR, 'fdr_sign')
        

    # perform Permutation if 'perm' is chosen
    elif correction_method == 'perm':
        
#        glmSim = MapNode(GLMFitSim(),
#                         iterfield=['glm_dir', 'permutation'],
#                         name='glmSim')
#        glmSim.inputs.spaces = '2spaces'
        
         raise NotImplementedError
     
        
    ### Plotting ###
    lh_bg_map = op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.sulc')
    rh_bg_map = op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.sulc')
    
    # Plot the mean map
    plotMeanMaps = MapNode(Function(input_names=['lh_surf', 'lh_surf_map', 'lh_bg_map',
                                                 'rh_surf', 'rh_surf_map', 'rh_bg_map',
                                                 'out_fname'],
                                    output_name=['out_file'],
                                    function=plot_surf_map),
                           iterfield=['lh_surf_map', 'rh_surf_map', 'out_fname'],
                           name='plotMeanMaps')
    plotMeanMaps.inputs.lh_bg_map = lh_bg_map
    plotMeanMaps.inputs.rh_bg_map = rh_bg_map
    plotMeanMaps.inputs.out_fname = ['mean_fwhm{}.png'.format(s) for s in non_zero_fwhm]
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.')),
               plotMeanMaps, 'lh_surf')
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.')),
               plotMeanMaps, 'rh_surf')
    wf.connect(lhSmoothMean, 'out_file', plotMeanMaps, 'lh_surf_map')
    wf.connect(rhSmoothMean, 'out_file', plotMeanMaps, 'rh_surf_map')
    
      
    # Plot uncorrected maps
    plot_stat_inputs = ['lh_surf', 'lh_stat_map', 'lh_bg_map',
                        'rh_surf', 'rh_stat_map', 'rh_bg_map',
                        'out_fname', 'cmap', 'upper_lim', 'threshold']
    
    plotUncorrectedG = MapNode(Function(input_names=plot_stat_inputs,
                                        output_name=['out_file'],
                                        function=plot_surf_stat),
                               iterfield=['lh_stat_map', 'rh_stat_map', 'out_fname'],
                               name='plotUncorrectedG')
    plotUncorrectedG.inputs.lh_bg_map = lh_bg_map
    plotUncorrectedG.inputs.rh_bg_map = rh_bg_map
    plotUncorrectedG.inputs.cmap = 'jet'
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.')),
               plotUncorrectedG, 'lh_surf')
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.')),
               plotUncorrectedG, 'rh_surf')
    wf.connect(fileList, ('contrast_files',  makeFStringElementFromFnameList, '.mtx', '_uncorrected_gamma_map.png', True),
               plotUncorrectedG, 'out_fname')
    wf.connect(lhSBMglmfit, 'gamma_file', plotUncorrectedG, 'lh_stat_map')
    wf.connect(rhSBMglmfit, 'gamma_file', plotUncorrectedG, 'rh_stat_map')
    
    plotUncorrectedP = MapNode(Function(input_names=plot_stat_inputs,
                                        output_name=['out_file'],
                                        function=plot_surf_stat),
                               iterfield=['lh_stat_map', 'rh_stat_map', 'out_fname'],
                               name='plotUncorrectedP')
    plotUncorrectedP.inputs.lh_bg_map = lh_bg_map
    plotUncorrectedP.inputs.rh_bg_map = rh_bg_map
    plotUncorrectedP.inputs.upper_lim = 10.0
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.')),
               plotUncorrectedP, 'lh_surf')
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.')),
               plotUncorrectedP, 'rh_surf')
    wf.connect(fileList, ('contrast_files',  makeFStringElementFromFnameList, '.mtx', '_uncorrected_p_map.png', True),
               plotUncorrectedP, 'out_fname')
    wf.connect(lhSBMglmfit, 'sig_file', plotUncorrectedP, 'lh_stat_map')
    wf.connect(rhSBMglmfit, 'sig_file', plotUncorrectedP, 'rh_stat_map')
    
    # Plot the corrected map
    
    # For gamma first create gamma masked by corrected p
    lhMaskGamma = MapNode(MRIsCalc(),
                          iterfield = ['in_file1', 'in_file2'],
                          name='lhMaskGamma')
    lhMaskGamma.inputs.action = 'masked'
    lhMaskGamma.inputs.out_file = 'lh.masked_gamma.mgh'
    wf.connect(lhSBMglmfit, 'gamma_file', lhMaskGamma, 'in_file1')
    if correction_method == 'FDR':
        wf.connect(mriFDR, 'out_file1', lhMaskGamma, 'in_file2')
        
    rhMaskGamma = MapNode(MRIsCalc(),
                          iterfield = ['in_file1', 'in_file2'],
                          name='rhMaskGamma')
    rhMaskGamma.inputs.action = 'masked'
    rhMaskGamma.inputs.out_file = 'rh.masked_gamma.mgh'
    wf.connect(rhSBMglmfit, 'gamma_file', rhMaskGamma, 'in_file1')
    if correction_method == 'FDR':
        wf.connect(mriFDR, 'out_file2', rhMaskGamma, 'in_file2')
        
    # Plot masked gamma 
    plotCorrectedG = MapNode(Function(input_names=plot_stat_inputs,
                                      output_name=['out_file'],
                                      function=plot_surf_stat),
                             iterfield=['lh_stat_map', 'rh_stat_map', 'out_fname'],
                             name='plotCorrectedG')
    plotCorrectedG.inputs.lh_bg_map = lh_bg_map
    plotCorrectedG.inputs.rh_bg_map = rh_bg_map
    plotCorrectedG.inputs.cmap = 'jet'
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.')),
               plotCorrectedG, 'lh_surf')
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.')),
               plotCorrectedG, 'rh_surf')
    wf.connect(fileList, ('contrast_files',  makeFStringElementFromFnameList, '.mtx', '_masked_gamma_map.png', True),
               plotCorrectedG, 'out_fname')
    wf.connect(lhMaskGamma, 'out_file', plotCorrectedG, 'lh_stat_map')
    wf.connect(rhMaskGamma, 'out_file', plotCorrectedG, 'rh_stat_map')
    
    # Plot thresholded P
    plotCorrectedP = MapNode(Function(input_names=plot_stat_inputs,
                                      output_name=['out_file'],
                                      function=plot_surf_stat),
                             iterfield=['lh_stat_map', 'rh_stat_map',
                                        'threshold', 'out_fname'],
                             name='plotCorrectedP')
    plotCorrectedP.inputs.lh_bg_map = op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.sulc')
    plotCorrectedP.inputs.rh_bg_map = op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.sulc')
    plotCorrectedP.inputs.upper_lim = 10.0
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'lh.')),
               plotCorrectedP, 'lh_surf')
    wf.connect(plotSurfList, ('surf', prependString, op.join(fs_subjects_dir, target_atlas, 'surf', 'rh.')),
               plotCorrectedP, 'rh_surf')
    wf.connect(lhSBMglmfit, 'sig_file', plotCorrectedP, 'lh_stat_map')
    wf.connect(rhSBMglmfit, 'sig_file', plotCorrectedP, 'rh_stat_map')
    wf.connect(fileList, ('contrast_files', makeFStringElementFromFnameList, '.mtx', '_corrected_p_map.png', True),
               plotCorrectedP, 'out_fname')
    if correction_method == 'FDR':
        wf.connect(mriFDR, 'out_thr_file', plotCorrectedP, 'threshold')
    
    
#    # Datasink
#    datasink = Node(DataSink(base_directory=base_dir,
#                             container='%sSink' % name),
#                    name='Datasink')
#    
#    glm_outputs = ['gamma_file', 'gamma_var_file', 'sig_file', 'ftest_file']
#    for out in glm_outputs:
#        wf.connect(lhSBMglmfit, out, datasink, 'lhSBMglm_{}'.format(out))
#        wf.connect(rhSBMglmfit, out, datasink, 'rhSBMglm_{}'.format(out))
#    
#    if correction_method == 'FDR':
#        wf.connect(mriFDR, 'out_file1', datasink, 'lhSBM_fdr_corrected_sig')
#        wf.connect(mriFDR, 'out_file2', datasink, 'rhSBM_fdr_corrected_sig')
#        
#    wf.connect(lhMaskGamma, 'out_file', datasink, 'lhSBM_masked_gamma')
#    wf.connect(rhMaskGamma, 'out_file', datasink, 'rhSBM_masked_gamma')
#    
#    wf.connect(plotMeanMaps, 'out_file', datasink, 'mean_map_png')  
#    wf.connect(plotUncorrectedG, 'out_file', datasink, 'uncorrected_gamma_png')
#    wf.connect(plotUncorrectedP, 'out_file', datasink, 'uncorrected_p_png')
#    wf.connect(plotCorrectedG, 'out_file', datasink, 'masked_gamma_png')
#    wf.connect(plotCorrectedP, 'out_file', datasink, 'corrected_p_png')
#    
    return wf
    
  