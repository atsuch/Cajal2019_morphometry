#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue 23 July 2019

This is a script to run FreesurferSBMglmWF in Freesurfer_SBM.py to perform
a simple set of GLM analyses on surface using fsgd and contrast files set up in 
/homes_unix/tsuchida/MRiSHARE/morphometry/SBM_main/Simple_Models.

It scans for the 3 groups in Simple_Models dir and runs the WF for each group
separately.

Tested with py3_ls (nipype 1.1.9 with py3)

@author: tsuchida
"""
import os
import os.path as op
import pandas as pd
import glob
import subprocess as sp
from nipype import config, logging
from ginnipi.workflows.Freesurfer_SBM import genFreesurferSBMglmWF    

MODEL_DIR = "/homes_unix/tsuchida/MRiSHARE/morphometry/SBM_main/Simple_Models"

def main(model_dir=MODEL_DIR, subid_label='mrishare_id', design_input='fsgd'):
    '''
    Runs the WF for all the groups found in model dir. 
    '''

    wd = "/beegfs_data/scratch/tsuchida-SBM"
    fs_subdir = '/data/analyses/work_in_progress/freesurfer/fsmrishare-flair6.0/'

    # find the groups in group dirs
    group_csv_glob = glob.glob(op.join(model_dir, '*', 'group_info.csv'))
    print ('Found {} group_info.csv...'.format(len(group_csv_glob)))
    
    if not group_csv_glob:
        raise Exception('No group info found in the model dir')
    
    for group_info_path in group_csv_glob:
        group_info = pd.read_csv(group_info_path)
        group_name = group_info_path.split('/')[-2]
        group_dir = op.dirname(group_info_path)
        
        # Copy the models to wd
        wd_indir = op.join(wd, 'input_dir')
        os.makedirs(wd_indir, exist_ok=True)
        group_indir = op.join(wd_indir, group_name)
        
        sp.call(['rsync', '-avh', '{}/'.format(group_dir), group_indir])
            
        # get the model name list
        model_dirs = glob.glob(op.join(group_dir, 'Model*/'))
        
        if not model_dirs:
            print ('No Model dirs found for the group {}'.format(group_name))
            break
        
        model_names = [m.split('/')[-2]for m in sorted(model_dirs)]
 
        print ('Found {} following models fround for the group {}'.format(len(model_names), group_name))
        print ('Gathering contrast info for each model...')
        
        model_info = {}
        for model_name, model_dir in zip(model_names, model_dirs):
            contrast_files = glob.glob(op.join(model_dir, '*.mtx'))
            if contrast_files:
                contrast_names = [op.basename(f).replace('.mtx', '') for f in contrast_files]
                print ('{}: {}'.format(model_name, contrast_names))

                for f, name in zip(contrast_files, contrast_names):
                    sign_file = f.replace('.mtx', '.mdtx')
                    if not op.exists(sign_file):
                        raise Exception ('Could not find corresponding sign file for contrast {}'.format(name))
                    else:
                        model_info[model_name] = 'dods'

            
        # log dir 
        log_dir = op.join(os.getcwd(), 'log_dir', group_name)
        os.makedirs(log_dir, exist_ok=True)
        
            
        # WF for the group
        group_sublist = group_info[subid_label].values.tolist()
        wf = genFreesurferSBMglmWF(name='SBM_{}'.format(group_name),
                                   base_dir=wd,
                                   group_sublist=group_sublist,
                                   model_dir=group_indir,
                                   model_info=model_info,
                                   design_input=design_input,
                                   fs_subjects_dir=fs_subdir,
                                   fwhm=[0.0, 10.0],
                                   measure_list=['thickness', 'area'],
                                   target_atlas='fsaverage',
                                   target_atlas_surfreg='sphere.reg',
                                   correction_method='FDR')
        
        config.update_config({'logging': {'log_directory': log_dir, 'log_to_file': True}})
        logging.update_logging(config)
        config.set('execution','job_finished_timeout','20.0')
        config.set('execution','keep_inputs', 'true')

        wf.run(plugin='SLURM',
               plugin_args = {#'sbatch_args': '--partition=gindev',
                               'dont_resubmit_completed_jobs': True,
                               'max_jobs': 50})


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description='generates New Freesurfer SBM wf for models specified')
    parser.add_argument('-d', '--model_dir', action='store',
                        default=MODEL_DIR,
                        help='directory containing groups, and models for each group')
    parser.add_argument('-s', '--subid_label', action='store',
                        default='mrishare_id',
                        help='column name for the subject id')
    parser.add_argument('--design_input', action='store',
                        default='fsgd', const='fsgd', nargs='?', choices=['fsgd', 'design_mat'],
                        help='mri_glmfit design input')

    args = vars(parser.parse_args())

    main(**args)
