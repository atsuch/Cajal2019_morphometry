#!/usr/bin/env python
"""
A simplified pipeline of ABACI Anatomical processing pipeline for demo.

@Author: Ami Tsuchida
@Author: Alexandre Laurent
"""

import os
import os.path as op
import sys

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import DataGrabber
from nipype.interfaces.base import Undefined
from nipype.algorithms.misc import Gunzip
from nipype.interfaces.freesurfer.preprocess import ReconAll
from nipype.interfaces.freesurfer.utils import MakeAverageSubject
from ginnipi.interfaces.custom import AbcReconAll, AbcFreeSurferSource
from ginnipi.toolbox.flow import (getElementFromList, getElementFromSplit,
                                  createList8)
from ginnipi.toolbox.utilities import checkfsVersion, get_scanlist_BIDS, get_scanlist
from ginnipi.workflows.subworkflows.basic_structural import genBasicStructuralPipeline
from ginnipi.workflows.subworkflows.anat_supplQC import genVbmFsCrossQCWorkflow


def genIshareAnat(name='iShareAnat',
                  base_dir=os.path.abspath('morphometry_tutorial/workdir'),
                  use_FLAIR=True,
                  do_recon_all=True,
                  do_fast_first=False,
                  fs_subjects_dir=os.path.abspath('morphometry_tutorial/subjects_dir'),
                  input_dir=os.path.abspath('morphometry_tutorial/data'),
                  subjects=None,
                  sinks=False,
                  spm=op.join(os.getenv('HOME'),'matlab', 'spm12'),
                  spm_standalone=None,
                  mcr=None):
    '''
    Generates a nipype workflow for the structural parts : Freesurfer, SPM
    NewSegment, optionally FSL FAST/FIRST.
    
    Either just process T1 or use FLAIR with use_FLAIR.
        
    Inputs:
        - name: name of the workflow
        - base_dir: the working directory
        - use_FLAIR:
        - do_recon_all:
        - do_fast_first:
        - fs_subjects_dir: the freesurfer SUBJECTS_DIR to use
        - input_dir: the directory were the input files downloaded from XNAT are written
        - subjects: a subject list (names should correspond to the folders inside the input_dir
        - sinks: create output sink nodes
        - spm_standalone: True if using the SPM standalone version and matlab MCR.
        - mcr: path to matlab compiled runtime folder
        - mcr2015:
        
    '''
    ## Construction of the workflow
    ishanat = Workflow(name)
    ishanat.base_dir = op.abspath(base_dir)

    ## List of subject
    if subjects:
        subject_list = subjects
    else:
        subject_list = [subject_dir.split("/")[-1] for subject_dir in os.listdir(input_dir) if op.isdir(op.join(input_dir,subject_dir))]
    
    # List --> Nipype iterable
    subjectList = Node(IdentityInterface(fields=['subject_id'],
                                         mandatory_inputs=True),
                       name="subjectList")    
    subjectList.iterables = ('subject_id', subject_list)
    

    ## List of scans
    outfields = ['T1']
    if use_FLAIR:
        outfields.append('FLAIR')
        
    scanList = Node(DataGrabber(infields=['subject_id'],
                                outfields=outfields),
                    name="scanList")
    scanList.inputs.base_directory = input_dir
    scanList.inputs.ignore_exception = False
    scanList.inputs.raise_on_empty = True
    scanList.inputs.sort_filelist = True
    scanList.inputs.template = '%s/%s/%s'
        ishanat.connect(subjectList, ("subject_id", get_scanlist_BIDS, input_dir, use_FLAIR, False, False),
                        scanList, "template_args")
    else:     
        scanList.inputs.template = '%s/%s/*.dcm'
        scanList.inputs.template_args = get_scanlist(use_FLAIR, False, False) 
    ishanat.connect(subjectList, "subject_id", scanList, "subject_id")

    ## Recon-all has already been run --> need to set the subject directory and the subject ID before structural processing
    if not do_recon_all:

        # we substitute an fsSource node (the recon-all has already been run)
        fsReconAll = Node(AbcFreeSurferSource(), name="fsReconAll")
        
        if op.isdir(subjects_dir):
            # Set subject directory
            fsReconAll.inputs.subjects_dir = subjects_dir 
            
            # Set subject ID
            if op.isdir(subjects_dir + '/'+ subject_list[0]):
                ishanat.connect(subjectList, "subject_id", fsReconAll, "subject_id")
            elif (not BIDS) and op.isdir(subjects_dir + '/' + subject_list[0].split('_')[0]):
                ishanat.connect(subjectList, ("subject_id",getElementFromSplit,'_',0), fsReconAll, "subject_id")
            elif BIDS and op.isdir(subjects_dir + '/' + subject_list[0].split('sub-')[1]):
                ishanat.connect(subjectList, ("subject_id",getElementFromSplit,'sub-',1), fsReconAll, "subject_id")
            else:
                if BIDS:
                    print('There is no subject directory named ' + subject_list[0] + ' or ' + subject_list[0].split('sub-')[1] + ' in recon-all directory ' + subjects_dir)
                else:
                    print('There is no subject directory named ' + subject_list[0] + ' or ' + subject_list[0].split('_')[0] + ' in recon-all directory ' + subjects_dir)
                sys.exit(2)
                
        else:
            print('Recon-all directory ' + subjects_dir + ' does not exist. Please make sure your configuration tab is properly settled')
            sys.exit(2)


    ## Structural processing
    # Main structural workflow
    struct = genBasicStructuralPipeline('struct',
                                        spm_standalone=spm_standalone,
                                        mcr=mcr,
                                        spmpath=spm,
                                        TPM_file=tpm,
                                        brainmask=brainmask,
                                        multichannel_seg=use_FLAIR,
                                        seg_channel_info=seg_channel_info,
                                        seg_sampling_distance=seg_sampling_distance,
                                        subjects_dir=subjects_dir,
                                        do_recon_all=do_recon_all,
                                        do_fast_first=do_fast_first,
                                        BIDS=BIDS)
    if do_recon_all:
        ishanat.connect(subjectList, "subject_id", struct, "inputNode.subject_id")
    else:
        ishanat.connect(fsReconAll, "subject_id", struct, "inputNode.subject_id")
        ishanat.connect(fsReconAll, "nu", struct, "inputNode.nu")
        if use_FLAIR:
            ishanat.connect(fsReconAll, "FLAIR", struct, "inputNode.FLAIR")
    ishanat.connect(scanList, "main", struct,'inputNode.main')
    if use_FLAIR:
        ishanat.connect(scanList, "acc", struct,'inputNode.acc')
    
    if do_recon_all:
        # two version of recon-all : fs5 and fs6 commands
        if fsversion == "freesurfer5":
            
            if not checkfsVersion("v5"):
                raise Exception("If I may .. you are using freesurfer 6.0 and requested a freesurfer5 command of recon-all")
            
            # we actually do the recon-all
            fsReconAll = Node(AbcReconAll(), name="fsReconAll5")
            fsReconAll.inputs.args = ' -3T -contrasurfreg -qcache -no-isrunning'
            fsReconAll.inputs.directive = 'all'
            fsReconAll.inputs.environ = {}
            fsReconAll.ignore_exception = False
            fsReconAll.inputs.subjects_dir = subjects_dir
            fsReconAll.terminal_output = 'stream'
            fsReconAll.inputs.use_FLAIR = use_FLAIR
            # no need to connect the T1 and T2flair since the process was already initialized in basic_structural
            ishanat.connect(struct, "outputNode.subject_id", fsReconAll, "subject_id")
            
        else:
            
            if not checkfsVersion("v6"):
                raise Exception("If I may .. you are using some freesurfer 5 version and requested a freesurfer6 command of recon-all")
            
            # Node: fs.fsReconAll
            fsReconAll = Node(AbcReconAll(), name="fsReconAll6")
            fsReconAll.inputs.args = ' -brainstem-structures -3T -contrasurfreg -qcache -no-isrunning'
            fsReconAll.inputs.directive = 'all'
            fsReconAll.inputs.environ = {}
            fsReconAll.ignore_exception = True # may throw an exception that cannot find result files ==> dummy
            fsReconAll.inputs.subjects_dir = subjects_dir
            fsReconAll.inputs.use_FLAIR = use_FLAIR
            # no need to connect the T1 and T2flair since the process was already initialized in basic_structural
            ishanat.connect(struct, "outputNode.subject_id", fsReconAll, "subject_id")   
            
            # Node: fs.fsReconAllHippoT1 pour calculer l'hippocampe sans le FLAIR
            fsReconAllHippoT1 = Node(AbcReconAll(), name="fsHipT1")
            fsReconAllHippoT1.plugin_args={'sbatch_args': '-t 20:00:00'}
            fsReconAllHippoT1.inputs.args = ' -hippocampal-subfields-T1 '
            fsReconAllHippoT1.inputs.subjects_dir = subjects_dir
            ishanat.connect(fsReconAll, "subject_id", fsReconAllHippoT1, "subject_id")
            
            if use_FLAIR:
                # Node: fs.fsReconAllHippoT1T2 pour calculer l'hippocampe AVEC le FLAIR
                fsReconAllHippoT1T2 = Node(AbcReconAll(), name="fsHipT1T2")
                fsReconAllHippoT1T2.plugin_args={'sbatch_args': '--mem 7000'}
                fsReconAllHippoT1T2.inputs.subjects_dir = subjects_dir
                ishanat.connect(fsReconAllHippoT1, "subject_id", fsReconAllHippoT1T2, "subject_id")
                ishanat.connect(struct, "outputNode.nifti_acc", fsReconAllHippoT1T2, "hippo_file")

#    # Conversion of FS outputs
#    fsconv = genPostFsPipeline(name='fsconv', fsversion=fsversion, pipeline='connectomics',no_myelin=not use_FLAIR)
#    if fsversion == "freesurfer6" and do_recon_all:
#        if use_FLAIR:
#            ishanat.connect(fsReconAllHippoT1T2, "subject_id", fsconv, "inputNode.subject_id")
#            ishanat.connect(fsReconAllHippoT1T2, "subjects_dir", fsconv, "inputNode.subjects_dir")
#        else:
#            ishanat.connect(fsReconAllHippoT1, "subject_id", fsconv, "inputNode.subject_id")
#            ishanat.connect(fsReconAllHippoT1, "subjects_dir", fsconv, "inputNode.subjects_dir")
#    else:
#        ishanat.connect(fsReconAll, "subject_id", fsconv, "inputNode.subject_id")
#        ishanat.connect(fsReconAll, "subjects_dir", fsconv, "inputNode.subjects_dir")
    
    
    ### Compute myelin maps if use_FLAIR
    if use_FLAIR:
        myelin = genMyelinMapPipeline(name='myelin')
        ishanat.connect(fsReconAll, "subject_id", myelin, "inputNode.subject_id")
        ishanat.connect(fsReconAll, "subjects_dir", myelin, "inputNode.subjects_dir")
        ishanat.connect(struct, ("outputNode.bias_corrected_images", getElementFromList, 0), myelin, "inputNode.T1w")
        ishanat.connect(struct, ("outputNode.bias_corrected_images", getElementFromList, 0), myelin, "inputNode.coregT2w")
        
    ### Perform cross check with FS segmentation
        
    SpmFsCrossCheck = genVbmFsCrossQCWorkflow(name='SpmFsCrossCheck', fsHipp=False)
    ishanat.connect(fsReconAll, "subject_id", SpmFsCrossCheck, "inputNode.subject_id")
    ishanat.connect(struct, ("outputNode.bias_corrected_images", getElementFromList, 0), SpmFsCrossCheck, "inputNode.ref_main")
    ishanat.connect(struct, ("outputNode.native_class_images", getElementFromList, 0), SpmFsCrossCheck, "inputNode.vbm_native_gm")
    ishanat.connect(struct, ("outputNode.native_class_images", getElementFromList, 1), SpmFsCrossCheck, "inputNode.vbm_native_wm")
    ishanat.connect(struct, ("outputNode.native_class_images", getElementFromList, 2), SpmFsCrossCheck, "inputNode.vbm_native_csf")
 
    if sinks:
        from nipype.interfaces.io import DataSink
        
        # XNAT assessor sinks
        assessorSink = Node(DataSink(), name='assessorSink')  
        assessorSink.inputs.container = 'data'            
        assessorSink.inputs.parameterization = False 
        ishanat.connect([
                     (struct, assessorSink,[('outputNode.xnat_assessor', 'vbmAssessorData')]),
                     #(fsconv, assessorSink,[('outputNode.xnat_assessor', 'fsAssessorData')]),
                    ]) 

        # Structural results sink: VBM
        anatSink = Node(DataSink(), name='anatSink')
        anatSink.inputs.container = 'data'
        anatSink.inputs.parameterization = False
        ishanat.connect([(struct, anatSink, [("outputNode.native_class_images", "native_class_images"),
                                             ("outputNode.bias_corrected_images", "bias_corrected_images"),
                                             ("outputNode.dartel_input_images", "dartel_input_images"),
                                             ("outputNode.forward_deformation_field", "forward_deformation_field"),
                                             ("outputNode.inverse_deformation_field", "inverse_deformation_field"),
                                             ("outputNode.newsegment_stereotaxic_registration", "newsegment_stereotaxic_registration"),
                                             ("outputNode.new_modulated_class_images", "SPM12_modulated_class_images"),
                                             ("outputNode.structural_MNI_111", "structural_MNI_111"),
                                             ("outputNode.structural_MNI_222", "structural_MNI_222"),
                                             ("outputNode.tissue_classes_MNI_111", "normalized_class_images_MNI_111"),
                                             ("outputNode.tissue_classes_MNI_222", "normalized_class_images_MNI_222"),
                                             ("outputNode.jacobian_map", "jacobian_map_MNI_111"),
                                             ("outputNode.jacobian_modulated_gm_image", "jacobian_modulated_gm_image_MNI_111"),
                                             ("outputNode.jacobian_modulated_wm_image", "jacobian_modulated_wm_image_MNI_111"),
                                             ("outputNode.jacobian_modulated_csf_image", "jacobian_modulated_csf_image_MNI_111")]),
                   
                         ])

        # Structural results sink: coreg
        coregSink = Node(DataSink(), name='coregSink')
        coregSink.inputs.container = 'data'
        coregSink.inputs.parameterization = False
        ishanat.connect([(struct, coregSink,[("outputNode.cropped_BFcorrected_main", "cropped_main"),
                                             ("outputNode.cropped_skullstripped_BFcorrected_main", "skullstripped_nuc_cropped_main")])
                        ])
        if use_FLAIR:
            ishanat.connect([(struct, coregSink,[("outputNode.coregistered_acc_FS", "coregistered_acc_FS"),
                                                 ("outputNode.cropped_coregistered_BFcorrected_acc", "coregistered_cropped_acc")]),
                             (myelin, coregSink, [("outputNode.raw_myelin_ratio_map_vol", "raw_myelin_ratio_map_vol"),
                                                  ("outputNode.smoothed_myelin_ratio_map_surf", "smoothed_myelin_ratio_map_surf")])
                                             
                            ])

        # Structural results sink: freesurfer 2 nifti conversions ; possibly add some gii?
#        fsSink = Node(DataSink(), name='fsSink')
#        fsSink.inputs.container = 'data'
#        fsSink.inputs.parameterization = False
#        ishanat.connect([(fsconv, fsSink, [("outputNode.ventricle_mask", "fs_ventricle_mask"),
#                                       ("outputNode.aseg_image", "fs_aseg_image"),
#                                       ("outputNode.parcellation2005_image", "fs_parcellation2005_image"),
#                                       ("outputNode.parcellation2009_image", "fs_parcellation2009_image"),
#                                       ("outputNode.brain_mask", "fs_brain_mask"),
#                                       ("outputNode.binary_brain_mask", "fs_binary_brain_mask"),
#                                       ("outputNode.fs_preprocessed_T1", "fs_preprocessed_T1"),
#                                       ("outputNode.white_matter_parcellation", "fs_white_matter_parcellation"),
#                                       ("outputNode.cortical_ribbon", "fs_cortical_ribbon")
#                                      ])
#                    ])
#        if use_FLAIR:
#            ishanat.connect([(fsconv, fsSink, [("outputNode.smoothed_T1byT2_ratio_map", "fs_smoothed_T1byT2_ratio_map"),
#                                           ("outputNode.smoothed_T1byT2_ratio_map_corrected", "fs_smoothed_T1byT2_ratio_map_corrected")
#                                          ])
#                         ])
        
        # Quality Control sink
        qcSink = Node(DataSink(), name='qcSink')
        qcSink.inputs.container = 'data'
        qcSink.inputs.parameterization = False
        ishanat.connect([(struct, qcSink,[("outputNode.spm_seg_QC_images_native", "struct_spm_seg_QCplot_native"),
                                          ("outputNode.spm_seg_QC_images_mod", "struct_spm_seg_QCplot_mod"),
                                          #("outputNode.registration_isocontours", "struct_T1NormContoursPlot"),
                                          ("outputNode.brainmask_QC", "struct_brainmask_QCplot")])
                        ])
        if use_FLAIR:
            ishanat.connect([(struct, qcSink,[("outputNode.coregistration_isocontours", "struct_T2CoregContoursPlot"),
                                              ("outputNode.coregistration_CostFunction", "struct_T2CoregCostFunction")])
                            ])
          
        # Additional structural results sink: FAST and FIRST segmentations
        if do_fast_first:
            fastfirstSink = Node(DataSink(), name='fastfirstSink')
            fastfirstSink.inputs.container = 'data'
            fastfirstSink.inputs.parameterization = False
            ishanat.connect([(struct, fastfirstSink, [("outputNode.fast_class_images", "T1fast_class_images"),
                                                      ("outputNode.fast_class_map", "T1fast_class_map"),
                                                      ("outputNode.fast_pv_images", "T1fast_pv_images"),
                                                      ("outputNode.fast_pv_map", "T1fast_pv_map"),
                                                      ("outputNode.first_bvars", "first_bvars"),
                                                     ("outputNode.first_segmentations_3D", "first_segmentations_3D"),
                                                     ("outputNode.first_segmentation_4D", "first_segmentation_4D"),
                                                     ("outputNode.first_vtk_surfaces", "first_vtk_surfaces")])
                             ])
            ishanat.connect([(struct, qcSink, [("outputNode.fsl_seg_QC_images", "struct_fsl_seg_QCplot")])])



    ishanat.write_graph(dotfilename='ishare_anat', graph2use='flat', format='png', simple_form=True)       
    
    return ishanat

    
def sinkDict():
    '''
    A dictionary of sinks with annotations regarding their upload in XNAT
    '''
    sinks=  {
                'anatSink':('resources',
                            'ABC_Anatomical_Data'),
                'coregSink':('resources',
                             'ABC_Anatomical_Coregistration_Data'),
                'fastfirstSink':('resources',
                                 'ABC_Fast_First_Data'),
                'fsSink':('resources',
                          'ABC_Freesurfer_Data_v6.0'),
                'qcSink':  ('resources', 
                            'ABC_Quality_Control_Data'),
                'assessorSink': ('assessor_data', 'ABC_Assessors')

            }
    return sinks

