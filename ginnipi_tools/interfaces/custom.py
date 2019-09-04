# -*- coding: utf-8 -*-
'''
ABC pipeline interfaces and customizations of the original nipype interfaces
@author: Pierre-Yves Hervé 2015-16
'''
import warnings
import os
from glob import glob
import numpy as np
from scipy.io import savemat
import subprocess
from nipype.interfaces.base import (Directory, TraitedSpec, traits,
                                    isdefined, File, InputMultiPath)
from nipype.interfaces.base import (CommandLineInputSpec, CommandLine,
                                    OutputMultiPath)
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec
from nipype.utils.filemanip import (split_filename, list_to_filename)
from nipype.interfaces.freesurfer.base import FSTraitedSpec, FSCommand
from nipype.interfaces.freesurfer.preprocess import (ReconAllOutputSpec,
                                                     FreeSurferSource)
from nipype.interfaces.spm.base import (SPMCommand, SPMCommandInputSpec,
                                        scans_for_fname)

warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)


#################################################
##### Custom interface for Freesurfer tools
#################################################

class AbcReconAllInputSpec(CommandLineInputSpec):
    subject_id = traits.Str("recon_all", argstr='-subjid %s',
                            desc='subject name', usedefault=True)
    directive = traits.Enum('all', 'autorecon1', 'autorecon2', 'autorecon2-cp',
                            'autorecon2-wm', 'autorecon2-inflate1',
                            'autorecon2-perhemi', 'autorecon3', 'localGI',
                            'qcache', argstr='-%s', desc='process directive',
                            usedefault=False, position=0)
    hemi = traits.Enum('lh', 'rh', desc='hemisphere to process',
                       argstr="-hemi %s")
    T1_files = InputMultiPath(File(exists=True), argstr='-i %s...',
                              desc='name of T1 file to process')
    T2_file = File(exists=True, argstr="-T2 %s", min_ver='5.3.0',
                   desc='Specify a T2 image to refine the cortical surface')
    FLAIR_file = File(exists=True, argstr="-FLAIR %s", min_ver='5.3.0',
                   desc='Specify a FLAIR file to refine the cortical surface',
                   xor='T2_file')
    hippo_file = File(exists=True, argstr="-hippocampal-subfields-T1T2 %s FLAIR", min_ver='6.0',
                   desc='Specify a T2/DP/FLAIR file to generate an automated segmentation of the hippocampal subfields')
    use_FLAIR = traits.Bool(argstr="-FLAIRpial", min_ver='5.3.0',
                   desc='Use a FLAIR image to refine the cortical surface')
    openmp = traits.Int(argstr="-openmp %d",
                        desc="Number of processors to use in parallel")
    subjects_dir = Directory(exists=True, argstr='-sd %s', hash_files=False,
                             desc='path to subjects directory', genfile=True)
    flags = traits.Str(argstr='%s', desc='additional parameters')

class AbcReconAllOutputSpec(ReconAllOutputSpec):
    FLAIR = File(
        exists=True,
        desc='Raw T2/FLAIR image conformed to Freesurfer space',
        loc='mri',
        altkey='orig/FLAIRraw')


class AbcFreeSurferSource(FreeSurferSource):
    output_spec = AbcReconAllOutputSpec

    def _list_outputs(self):
        subjects_dir = self.inputs.subjects_dir
        subject_path = os.path.join(subjects_dir, self.inputs.subject_id)
        output_traits = self._outputs()
        outputs = output_traits.get()
        for k in list(outputs.keys()):
            if not k in ['subjects_dir','subject_id']:
                val = self._get_files(subject_path, k,
                                      output_traits.traits()[k].loc,
                                      output_traits.traits()[k].altkey)
                if val:
                    outputs[k] = list_to_filename(val)
        outputs['subject_id'] = self.inputs.subject_id
        outputs['subjects_dir'] = subjects_dir
        return outputs

class AbcReconAll(CommandLine):
    """
    Uses recon-all to generate surfaces and parcellations of structural data
    from anatomical images of a subject. This customization of the standard interface
    allows using FLAIR or T2 extra files.

    Examples
    --------

    >>> from nipype.interfaces.freesurfer import ReconAll
    >>> reconall = ReconAll()
    >>> reconall.inputs.subject_id = 'foo'
    >>> reconall.inputs.directive = 'all'
    >>> reconall.inputs.subjects_dir = '.'
    >>> reconall.inputs.T1_files = 'structural.nii'
    >>> reconall.cmdline
    'recon-all -all -i structural.nii -subjid foo -sd .'
    """

    _cmd = 'recon-all'
    _additional_metadata = ['loc', 'altkey']
    input_spec = AbcReconAllInputSpec
    output_spec = AbcReconAllOutputSpec
    _can_resume = True


    _steps = [
        #autorecon1
        ('motioncor', ['mri/rawavg.mgz', 'mri/orig.mgz']),
        ('talairach', ['mri/transforms/talairach.auto.xfm',
                       'mri/transforms/talairach.xfm']),
        ('nuintensitycor', ['mri/nu.mgz']),
        ('normalization', ['mri/T1.mgz']),
        ('skullstrip',
         ['mri/brainmask.auto.mgz',
          'mri/brainmask.mgz']),
        #autorecon2
        ('gcareg', ['mri/transforms/talairach.lta']),
        ('canorm', ['mri/norm.mgz']),
        ('careg', ['mri/transforms/talairach.m3z']),
        ('careginv', ['mri/transforms/talairach.m3z.inv.x.mgz',
                      'mri/transforms/talairach.m3z.inv.y.mgz',
                      'mri/transforms/talairach.m3z.inv.z.mgz']),
        ('rmneck', ['mri/nu_noneck.mgz']),
        ('skull-lta', ['mri/transforms/talairach_with_skull_2.lta']),
        ('calabel',
         ['mri/aseg.auto_noCCseg.mgz', 'mri/aseg.auto.mgz', 'mri/aseg.mgz']),
        ('normalization2', ['mri/brain.mgz']),
        ('maskbfs', ['mri/brain.finalsurfs.mgz']),
        ('segmentation', ['mri/wm.asegedit.mgz', 'mri/wm.mgz']),
        ('fill', ['mri/filled.mgz']),
        ('tessellate', ['surf/lh.orig.nofix', 'surf/rh.orig.nofix']),
        ('smooth1', ['surf/lh.smoothwm.nofix', 'surf/rh.smoothwm.nofix']),
        ('inflate1', ['surf/lh.inflated.nofix', 'surf/rh.inflated.nofix']),
        ('qsphere', ['surf/lh.qsphere.nofix', 'surf/rh.qsphere.nofix']),
        ('fix', ['surf/lh.orig', 'surf/rh.orig']),
        ('white',
         ['surf/lh.white',
          'surf/rh.white',
          'surf/lh.curv',
          'surf/rh.curv',
          'surf/lh.area',
          'surf/rh.area',
          'label/lh.cortex.label',
          'label/rh.cortex.label']),
        ('smooth2', ['surf/lh.smoothwm', 'surf/rh.smoothwm']),
        ('inflate2',
         ['surf/lh.inflated',
          'surf/rh.inflated',
          'surf/lh.sulc',
          'surf/rh.sulc',
          'surf/lh.inflated.H',
          'surf/rh.inflated.H',
          'surf/lh.inflated.K',
          'surf/rh.inflated.K']),
        #autorecon3
        ('sphere', ['surf/lh.sphere', 'surf/rh.sphere']),
        ('surfreg', ['surf/lh.sphere.reg', 'surf/rh.sphere.reg']),
        ('jacobian_white', ['surf/lh.jacobian_white',
                            'surf/rh.jacobian_white']),
        ('avgcurv', ['surf/lh.avg_curv', 'surf/rh.avg_curv']),
        ('cortparc', ['label/lh.aparc.annot', 'label/rh.aparc.annot']),
        ('pial',
         ['surf/lh.pial',
          'surf/rh.pial',
          'surf/lh.curv.pial',
          'surf/rh.curv.pial',
          'surf/lh.area.pial',
          'surf/rh.area.pial',
          'surf/lh.thickness',
          'surf/rh.thickness']),
        ('cortparc2', ['label/lh.aparc.a2009s.annot',
                       'label/rh.aparc.a2009s.annot']),
        ('parcstats2',
         ['stats/lh.aparc.a2009s.stats',
          'stats/rh.aparc.a2009s.stats',
          'stats/aparc.annot.a2009s.ctab']),
        ('cortribbon', ['mri/lh.ribbon.mgz', 'mri/rh.ribbon.mgz',
                        'mri/ribbon.mgz']),
        ('segstats', ['stats/aseg.stats']),
        ('aparc2aseg', ['mri/aparc+aseg.mgz', 'mri/aparc.a2009s+aseg.mgz']),
        ('wmparc', ['mri/wmparc.mgz', 'stats/wmparc.stats']),
        ('balabels', ['BA.ctab', 'BA.thresh.ctab'])]
        #('label-exvivo-ec', ['label/lh.entorhinal_exvivo.label',
        #                     'label/rh.entorhinal_exvivo.label'])]

    def _gen_subjects_dir(self):
        return os.getcwd()

    def _gen_filename(self, name):
        if name == 'subjects_dir':
            return self._gen_subjects_dir()
        return None

    def _list_outputs(self):
        """
        See io.FreeSurferSource.outputs for the list of outputs returned
        """
        if isdefined(self.inputs.subjects_dir):
            subjects_dir = self.inputs.subjects_dir
        else:
            subjects_dir = self._gen_subjects_dir()

        if isdefined(self.inputs.hemi):
            hemi = self.inputs.hemi
        else:
            hemi = 'both'
        outputs = self._outputs().get()
        outputs.update(FreeSurferSource(subject_id=self.inputs.subject_id,
                                        subjects_dir=subjects_dir,
                                        hemi=hemi)._list_outputs())
        outputs['subject_id'] = self.inputs.subject_id
        outputs['subjects_dir'] = subjects_dir
        return outputs

    def _is_resuming(self):
        subjects_dir = self.inputs.subjects_dir
        if not isdefined(subjects_dir):
            subjects_dir = self._gen_subjects_dir()
        if os.path.isdir(os.path.join(subjects_dir, self.inputs.subject_id,
                                      'mri')):
            return True
        return False

    def _format_arg(self, name, trait_spec, value):
        if name == 'T1_files':
            if self._is_resuming():
                return ''
        return super(AbcReconAll, self)._format_arg(name, trait_spec, value)

    @property
    def cmdline(self):
        cmd = super(AbcReconAll, self).cmdline
        if not self._is_resuming():
            return cmd
        subjects_dir = self.inputs.subjects_dir
        if not isdefined(subjects_dir):
            subjects_dir = self._gen_subjects_dir()
        #cmd = cmd.replace(' -all ', ' -make all ')
        # iflogger.info('Overriding recon-all directive')
        # this part used to be commented because nipype added unwanted (?) flags in the fs command
        flags = []
        whole_procedure=True
        if isdefined(self.inputs.args):
            if 'hippocampal' in self.inputs.args:
                whole_procedure=False 
        if isdefined(self.inputs.hippo_file):
            whole_procedure=False
        if whole_procedure:
            directive = 'all'
            for idx, step in enumerate(self._steps):
                step, outfiles = step
                if all([os.path.exists(os.path.join(subjects_dir,
                                                    self.inputs.subject_id,f)) for
                        f in outfiles]):
                    flags.append('-no%s'%step)
                    if idx > 4:
                        directive = 'autorecon2'
                    elif idx > 23:
                        directive = 'autorecon3'
                else:
                    flags.append('-%s'%step)
            cmd = cmd.replace(' -%s ' % self.inputs.directive, ' -%s ' % directive)
            cmd += ' ' + ' '.join(flags)
        # iflogger.info('resume recon-all : %s' % cmd)
        return cmd


 
class FDRInputSpec(FSTraitedSpec):
    # required
    in_file1 = File(argstr="--i %s %s %s", exists=True, mandatory=True, copyfile=True,
                   desc="source volume or surface overlay")
    in_mask1 = File(exists=True,
                    requires=['in_file1'],
                    desc="optional mask for in_file1")
    out_file1 = File(genfile=True,
                     requires=['in_file1'],
                     desc="optional thresholded output file for in_file1")
    
    # optional second input
    in_file2 = File(argstr="--i %s %s %s", exists=True, mandatory=False, copyfile=True,
                   desc="source volume or surface overlay")
    in_mask2 = File(exists=True,
                    requires=['in_file2'],
                    desc="optional mask for in_file2")
    out_file2 = File(genfile=True,
                     requires=['in_file2'],
                     desc="optional thresholded output file for in_file2")
    
    # other optional args
    fdr = traits.Float(argstr="--fdr %f", default=0.05, usedefault=True,
                       desc="value between 0 and 1, typically .05")
    fdr_sign = traits.Enum('pos', 'neg', 'abs',
                           argstr="--%s", 
                          desc="only consider positve (pos, negative (neg), or all (abs) voxels")
    out_thr_file = File(argstr="--thfile %s", 
                        desc="write threshold to text file")

class FDROutputSpec(TraitedSpec):
    out_file1 = File(desc="Output file for in_file1")
    out_file2 = File(desc="Output file for in_file2")
    out_thr_file = File(desc="Text file with a threshold")

class FDR(FSCommand):
    """ This program runs mri_fdr.
    Examples
    ========
    >>>from ginnipi.interfaces.custom import FDR
    >>>mriFDR = FDR()
    >>>mriFDR.inputs.fdr=0.05
    >>>mriFDR.inputs.fdr_sign = 'abs'
    >>>mriFDR.inputs.out_thr_file = 'fdr_thr.txt'
    >>>mriFDR.inputs.in_file1 = 'lh.sig.mgh'
    >>>mriFDR.inputs.mask_file1 = 'lh.mask.mgh'
    >>>mriFDR.inputs.out_file1 = 'lh.sig_corr_abs.mgh'
    
    >>>mriFDR.cmdline
    'mri_fdr --fdr 0.05 --thfile fdr_thr.txt -i lh.sig.mgh lh.mask.mgh lh.sig_cor_abs.mgh'
    
    """

    _cmd = 'mri_fdr'
    input_spec = FDRInputSpec
    output_spec = FDROutputSpec

    def _format_arg(self, opt, spec, val):
        if opt == 'in_file1':
            _si = self.inputs
            mask_file = _si.in_mask1 if isdefined(_si.in_mask1) else 'nomask'
            out_file = _si.out_file1 if isdefined(_si.out_file1) else 'nooutput'
            return spec.argstr % (_si.in_file1, mask_file, out_file)
        if opt == 'in_file2':
            _si = self.inputs
            mask_file = _si.in_mask2 if isdefined(_si.in_mask2) else 'nomask'
            out_file = _si.out_file2 if isdefined(_si.out_file2) else 'nooutput'
            return spec.argstr % (_si.in_file2, mask_file, out_file)

        return super(FDR, self)._format_arg(opt, spec, val)


    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file1) and self.inputs.out_file1:
            outputs['out_file1'] = os.path.abspath(self.inputs.out_file1)
        if isdefined(self.inputs.out_file2) and self.inputs.out_file2:
            outputs['out_file2'] = os.path.abspath(self.inputs.out_file2)
        if isdefined(self.inputs.out_thr_file) and self.inputs.out_thr_file:
            outputs['out_thr_file'] = os.path.abspath(self.inputs.out_thr_file)
        return outputs


# Custom interface for mri_glmfit-sim. Note that it wraps args required for permutation
# testing for now.
class GLMFitSimInputSpec(FSTraitedSpec):
    glm_dir = Directory(argstr='--glmdir %s', exists=True,
                        desc="glm directory from glm_fit")
    permutation = traits.Tuple(
                    traits.Int(default=1000),
                    traits.Float(default=4.0),
                    traits.Enum('pos', 'neg', 'abs'),
                    argstr='--perm %d %f %s',
                    desc="number of permutation to run, vertex-wise, cluster-forming thr, and sign")
    cwp = traits.Float(argstr="--cwp %f", default=0.05, usedefault=True,
                       desc="Keep clusters that have the specified cluster-wise p-values")
    spaces = traits.Enum('2spaces', '3spaces', argstr='--%s',
                         desc="correct for number of spaces bing tested")
    bg = traits.Int(default=1, desc="number of threads to use")
    #copy_inputs = traits.Bool(desc="If running as a node, set this to True " +
    #                          "this will copy some implicit inputs to the " + "node directory.")
    
class GLMFitSimOutputSpec(TraitedSpec):

    glm_dir = Directory(exists=True, desc="output directory")
    pdf_dat = File(exists=True, desc="probability distribution function of clusterwise correction")
    sig_corrected_map = File(exists=True, desc="cluster-wise corrected map (overlay)")
    sig_cluster_file = File(exists=True, desc="summary of clusters (text)")
    sig_masked_map = File(exists=True, desc="corrected sig values masked by the clusters that survive correction")
    ocn_annot = File(desc="output cluster number (annotation of clusters)")
    ocn_map = File(desc="output cluster number (segmentation showing where each numbered cluster is)")
    ocn_dat = File(desc="the average value of each subject in each cluster")

class GLMFitSim(FSCommand):
    """Use FreeSurfer v6.0 mri_glmfit-sim to run permutation.
    Examples
    --------
    >>> glmsim = GLMFitSim()
    >>> glmsim.inputs.glm_dir = 'lh.gender_age.glmdir'
    >>> glmsim.inputs.permutation = (1000, 4.0, 'abs')
    >>> glmsim.inputs.cwp = 0.05
    >>> glmsim.inputs.spaces = '2spaces'
    >>> glmsim.cmdline
    'mri_glmfit-sim --glmdir lh.gender_age.glmdir --perm 1000 4.0 abs --2spaces'
 
    """

    _cmd = 'mri_glmfit-sim'
    input_spec = GLMFitSimInputSpec
    
#    def run(self, **inputs):
#         
#        if self.inputs.copy_inputs:
#            node_glm_dir = os.getcwd()
#            files_to_copy = ['mri_glmfit.log', ]
#            glmfitlog = os.path.join(self.inputs.glm_dir, 'mri_glmfit.log')
#            self.inputs.glm_dir = os.getcwd()

    output_spec = GLMFitSimOutputSpec

    def _format_arg(self, name, spec, value):

        return super(GLMFitSim, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        # Get the top-level output directory
        if not isdefined(self.inputs.glm_dir):
            glmdir = os.getcwd()
        else:
            glmdir = os.path.abspath(self.inputs.glm_dir)
        outputs["glm_dir"] = glmdir

        # Assign the output files that always get created
        outputs["pdf_dat"] = os.path.join(glmdir, "perm.*.pdf.dat")
        outputs["sig_corrected_map"] = os.path.join(glmdir, "perm.*.sig.cluster.mgh")
        outputs["sig_cluster_file"] = os.path.join(glmdir, "perm.*.sig.cluster.summary")
        outputs["sig_masked_map"] = os.path.join(glmdir, "perm.*.sig.masked.mgh")
        outputs["ocn_annot"] = os.path.join(glmdir, "perm.*.sig.ocn.annot")
        outputs["ocn_map"] = os.path.join(glmdir, "perm.*.sig.ocn.mgh")
        outputs["ocn_dat"] = os.path.join(glmdir, "perm.*.sig.ocn.dat")


        return outputs
    
    
#################################################
##### Custom interface for SPM tools
#################################################
class AbcNewSegmentInputSpec(SPMCommandInputSpec):
    channel_files = InputMultiPath(File(exists=True),
                              desc="A list of files to be segmented",
                              field='channel', copyfile=False, mandatory=True)
    channel_info = traits.List(traits.Tuple(traits.Float(), traits.Float(),
                                traits.Tuple(traits.Bool, traits.Bool)),
                                desc="""A list of tuples with the following fields:
            - bias regularisation (0-10)
            - FWHM of Gaussian smoothness of bias
            - which maps to save (Corrected, Field) - a tuple of two boolean values""",
            field='channel')

    tissues = traits.List(traits.Tuple(traits.Tuple(File(exists=True), traits.Int()), traits.Int(),
                                       traits.Tuple(traits.Bool, traits.Bool), traits.Tuple(traits.Bool, traits.Bool)),
                         desc="""A list of tuples (one per tissue) with the following fields:
            - tissue probability map (4D), 1-based index to frame
            - number of gaussians
            - which maps to save [Native, DARTEL] - a tuple of two boolean values
            - which maps to save [Unmodulated, Modulated] - a tuple of two boolean values""",
            field='tissue')

    affine_regularization = traits.Enum('mni', 'eastern', 'subj', 'none', field='warp.affreg',
                      desc='mni, eastern, subj, none ')

    warping_regularization = traits.List(traits.Float(),
                                         field='warp.reg',
                                          desc='Aproximate distance between sampling points.',
                                          minlen=5,maxlen=5,default=[0, 0.001, 0.5000, 0.0500, 0.2000])
    mrf = traits.Int(default=1,
                     field='warp.mrf',
                     desc='Number of iterations for Markov Random Field cleanup. Set to zero to disable')

    cleanup =traits.Int(default=1,
                        field='warp.cleanup',
                        desc='Brain extraction procedure (0:no clean-up,1: light clean-up, 2:thorough clean-up')

    smoothness = traits.Float(default=0.0,
                              field='warp.fwhm',
                              desc='Smoothness of the images: 0.0 for MRI, 5.0 for PET or SPECT')

    sampling_distance = traits.Float(default= 1.0, field='warp.samp',
                                     desc='Sampling distance on data for parameter estimation')

    write_deformation_fields = traits.List(traits.Bool(), minlen=2, maxlen=2, field='warp.write',
                                           desc="Which deformation fields to write:[Inverse, Forward]")

    write_bounding_box = traits.List(traits.List(traits.Float(),
                                                 minlen=3, maxlen=3),
                                     field='warp.bb', minlen=2, maxlen=2,
                                     desc=('3x2-element list of lists representing '
                                           'the bounding box (in mm) to be written'))

    write_voxel_size = traits.Float(field='warp.vox',
                                    default=1.0,
                                    desc=('Isotropic voxel size (in mm) of the written normalised images'))

    niter = traits.Int(default = 1,
                       field = 'iterations',
                       desc = 'Number of iterations')


class AbcNewSegmentSPMCommand(SPMCommand):
    """
    SPMCommand is a Pseudo prototype class, meant to be subclassed.
    So this is what we do here.
    This subclass overloads the _make_matlab_command method
    to bypass spm_jobman_run at the end and replace it directly
    with the spatial job launcher:
    job=jobs{1}.spm.spatial.preproc;
    spm_preproc_run(job)
    This solves the problem with the NewSegment matlabbatch interface
    where you can't select a voxel size for writing images,
    and can't get modulated images through Normalise('write') because
    preserve volume option is absent from the matlabbatch interface.
    WARNING: have to be careful with neuro/radio convention of outputs
    """

    def _make_matlab_command(self, contents, postscript=None):
        """Generates a mfile to build job structure
        Parameters
        ----------

        contents : list
            a list of dicts generated by _parse_inputs
            in each subclass

        cwd : string
            default os.getcwd()

        Returns
        -------
        mscript : string
            contents of a script called by matlab

        """
        cwd = os.getcwd()
        mscript = """
        %% Generated by nipype.interfaces.spm
        if isempty(which('spm')),
             throw(MException('SPMCheck:NotFound', 'SPM not in matlab path'));
        end
        [name, version] = spm('ver');
        fprintf('SPM version: %s Release: %s\\n',name, version);
        fprintf('SPM path: %s\\n', which('spm'));
        spm('Defaults','fMRI');

        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
           spm_jobman('initcfg');
           spm_get_defaults('cmdline', 1);
        end\n
        """
        if self.mlab.inputs.mfile:
            if isdefined(self.inputs.use_v8struct) and self.inputs.use_v8struct:
                mscript += self._generate_job('jobs{1}.spm.%s.%s' %
                                              (self.jobtype, self.jobname),
                                              contents[0])
            else:
                if self.jobname in ['st', 'smooth', 'preproc', 'preproc8',
                                    'fmri_spec', 'fmri_est', 'factorial_design',
                                    'defs']:
                    # parentheses
                    mscript += self._generate_job('jobs{1}.%s{1}.%s(1)' %
                                                  (self.jobtype, self.jobname),
                                                  contents[0])
                else:
                    # curly brackets
                    mscript += self._generate_job('jobs{1}.%s{1}.%s{1}' %
                                                  (self.jobtype, self.jobname),
                                                  contents[0])
        else:
            jobdef = {'jobs': [{self.jobtype:
                                [{self.jobname:
                                  self.reformat_dict_for_savemat(contents[0])}]
                                }]}
            savemat(os.path.join(cwd, 'pyjobs_%s.mat' % self.jobname), jobdef)
            mscript += "load pyjobs_%s;\n\n" % self.jobname
#          mscript += """
#         job=jobs{1}.spm.spatial.preproc;\n
#         spm_preproc_run(job);\n
#         % ABACI custom interface,bypassing this: spm_jobman(\'run\', jobs);\n
#         """
        mscript += """
        job=jobs{1}.spm.%s.%s;\n
        """ %(self.jobtype, self.jobname)
        mscript += """
        spm_preproc_run(job);\n
        % ABACI custom interface,bypassing this: spm_jobman(\'run\', jobs);\n
        """

        if self.inputs.use_mcr:
            mscript += """
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close(\'all\', \'force\');
        end;
            """
        if postscript is not None:
            mscript += postscript
        return mscript


class AbcNewSegmentOutputSpec(TraitedSpec):
    native_class_images = traits.List(File(exists=True), desc='native space probability maps')
    dartel_input_images = traits.List(File(exists=True), desc='dartel imported class images')
    normalized_class_images = traits.List(File(exists=True), desc='normalized class images')
    modulated_class_images = traits.List(File(exists=True), desc='modulated+normalized class images')
    transformation_mat = OutputMultiPath(File(exists=True), desc='Normalization transformation')
    bias_corrected_images = OutputMultiPath(File(exists=True), desc='bias corrected images')
    bias_field_images = OutputMultiPath(File(exists=True), desc='bias field images')
    forward_deformation_field = OutputMultiPath(File(exists=True))
    inverse_deformation_field = OutputMultiPath(File(exists=True))

class AbcNewSegment(AbcNewSegmentSPMCommand):
    """
    Use spm_preproc8 (New Segment) to separate structural images into different
    tissue classes. Supports multiple modalities.

    NOTE: This interface handles multichannel segmentation, unlike the default nipype one.

    http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=43

    Examples
    --------
    >>> import nipype.interfaces.spm as spm
    >>> seg = spm.NewSegment()
    >>> seg.inputs.channel_files = 'structural.nii'
    >>> seg.inputs.channel_info = (0.0001, 60, (True, True))
    >>> seg.run() # doctest: +SKIP

    For VBM pre-processing [http://www.fil.ion.ucl.ac.uk/~john/misc/VBMclass10.pdf],
    TPM.nii should be replaced by /path/to/spm8/toolbox/Seg/TPM.nii

    >>> seg = NewSegment()
    >>> seg.inputs.channel_files = 'structural.nii'
    >>> tissue1 = (('TPM.nii', 1), 2, (True,True), (False, False))
    >>> tissue2 = (('TPM.nii', 2), 2, (True,True), (False, False))
    >>> tissue3 = (('TPM.nii', 3), 2, (True,False), (False, False))
    >>> tissue4 = (('TPM.nii', 4), 2, (False,False), (False, False))
    >>> tissue5 = (('TPM.nii', 5), 2, (False,False), (False, False))
    >>> seg.inputs.tissues = [tissue1, tissue2, tissue3, tissue4, tissue5]
    >>> seg.run() # doctest: +SKIP
    """

    input_spec = AbcNewSegmentInputSpec
    output_spec = AbcNewSegmentOutputSpec

    def __init__(self, **inputs):
        _local_version = SPMCommand().version
        if _local_version and '12.' in _local_version:
            self._jobtype = 'spatial'
            self._jobname = 'preproc'
        else:
            self._jobtype = 'tools'
            self._jobname = 'preproc8'

        SPMCommand.__init__(self, **inputs)

    def _format_arg(self, opt, spec, val):
        """
        Convert input to appropriate format for spm
        """

        if opt in ['channel_files', 'channel_info']:
            # structure have to be recreated, because of some weird traits error
            channel_list=[]
            if not isdefined(self.inputs.channel_info):
                for channel in self.inputs.channel_files:
                    new_channel = {}
                    new_channel['vols'] = scans_for_fname(channel)
                    channel_list.append(new_channel)
            else:
                c=0
                for channel in self.inputs.channel_files:
                    new_channel = {}
                    new_channel['vols'] = scans_for_fname(channel)
                    info = self.inputs.channel_info[c]
                    new_channel['biasreg'] = info[0]
                    new_channel['biasfwhm'] = info[1]
                    new_channel['write'] = [int(info[2][0]), int(info[2][1])]
                    channel_list.append(new_channel)
                    c=c+1
            return channel_list
        elif opt == 'tissues':
            new_tissues = []
            for tissue in val:
                new_tissue = {}
                new_tissue['tpm'] = np.array([','.join([tissue[0][0], str(tissue[0][1])])], dtype=object)
                new_tissue['ngaus'] = tissue[1]
                new_tissue['native'] = [int(tissue[2][0]), int(tissue[2][1])]
                new_tissue['warped'] = [int(tissue[3][0]), int(tissue[3][1])]
                new_tissues.append(new_tissue)
            return new_tissues
        elif opt == 'write_deformation_fields':
            return super(AbcNewSegment, self)._format_arg(opt, spec, [int(val[0]), int(val[1])])
        else:
            return super(AbcNewSegment, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['native_class_images'] = []
        outputs['dartel_input_images'] = []
        outputs['normalized_class_images'] = []
        outputs['modulated_class_images'] = []
        outputs['transformation_mat'] = []
        outputs['bias_corrected_images'] = []
        outputs['bias_field_images'] = []
        outputs['inverse_deformation_field'] = []
        outputs['forward_deformation_field'] = []

        n_classes = 5
        if isdefined(self.inputs.tissues):
            n_classes = len(self.inputs.tissues)
        pth, base, ext = split_filename(self.inputs.channel_files[0])
        if isdefined(self.inputs.tissues):
            for i, tissue in enumerate(self.inputs.tissues):
                if tissue[2][0]:
                    outputs['native_class_images'].append(os.path.join(pth, "c%d%s.nii" % (i+1, base)))
                if tissue[2][1]:
                    outputs['dartel_input_images'].append(os.path.join(pth, "rc%d%s.nii" % (i+1, base)))
                if tissue[3][0]:
                    outputs['normalized_class_images'].append(os.path.join(pth, "wc%d%s.nii" % (i+1, base)))
                if tissue[3][1]:
                    outputs['modulated_class_images'].append(os.path.join(pth, "mwc%d%s.nii" % (i+1, base)))
        else:
            for i in range(n_classes):
                outputs['native_class_images'][i].append(os.path.join(pth, "c%d%s.nii" % (i+1, base)))

        outputs['transformation_mat'].append(os.path.join(pth, "%s_seg8.mat" % base))

        if isdefined(self.inputs.write_deformation_fields):
            if self.inputs.write_deformation_fields[0]:
                outputs['inverse_deformation_field'].append(os.path.join(pth, "iy_%s.nii" % base))
            if self.inputs.write_deformation_fields[1]:
                outputs['forward_deformation_field'].append(os.path.join(pth, "y_%s.nii" % base))

        if isdefined(self.inputs.channel_info):
            if len(set([ os.path.split(fpth)[1] for fpth in self.inputs.channel_files]))==1 & len(self.inputs.channel_files)>1:
                channel_info=self.inputs.channel_info[0]
                if channel_info[2][0]:
                    outputs['bias_corrected_images'].append(os.path.join(pth, "m%s.nii" % (base)))
                if channel_info[2][1]:
                    outputs['bias_field_images'].append(os.path.join(pth, "BiasField_%s.nii" % (base)))
                for c in range(1,len(self.inputs.channel_info)):
                    channel_info=self.inputs.channel_info[c]
                    if channel_info[2][0]:
                        outputs['bias_corrected_images'].append(os.path.join(pth, "m%s_c%s.nii" % (base,str(c).zfill(4))))
                    if channel_info[2][1]:
                        outputs['bias_field_images'].append(os.path.join(pth, "BiasField_%s_c%s.nii" % (base,str(c).zfill(4))))
            else:
                c=0
                for channel_info in self.inputs.channel_info:
                    pth, base, ext = split_filename(self.inputs.channel_files[c])
                    if channel_info[2][0]:
                        outputs['bias_corrected_images'].append(os.path.join(pth, "m%s.nii" % (base)))
                    if channel_info[2][1]:
                        outputs['bias_field_images'].append(os.path.join(pth, "BiasField_%s.nii" % (base)))
                    c=c+1
        return outputs


class ABCApplyDeformationFieldInputSpec(SPMCommandInputSpec):
    in_files = InputMultiPath(exists=True, mandatory=True, field='out{1}.pull.fnames')
    deformation_field = File(exists=True, mandatory=True, field='comp{1}.def')
    reference_volume = File(exists=True, mandatory=True, field='comp{2}.id.space')
    interp = traits.Range(
        low=0,
        high=7,
        field='out{1}.pull.interp',
        desc='degree of b-spline used for interpolation')


class ABCApplyDeformationFieldOutputSpec(TraitedSpec):
    out_files = OutputMultiPath(File(exists=True))


class ABCApplyDeformations(SPMCommand):
    input_spec = ABCApplyDeformationFieldInputSpec
    output_spec = ABCApplyDeformationFieldOutputSpec

    _jobtype = 'util'
    _jobname = 'defs'

    def _make_matlab_command(self, contents, postscript=None):
        """Generates a mfile to build job structure
        """
        cwd = os.getcwd()
        mscript = """
        %% Generated by nipype.interfaces.spm
        if isempty(which('spm')),
             throw(MException('SPMCheck:NotFound', 'SPM not in matlab path'));
        end
        [name, version] = spm('ver');
        fprintf('SPM version: %s Release: %s\\n',name, version);
        fprintf('SPM path: %s\\n', which('spm'));
        spm('Defaults','fMRI');

        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
           spm_jobman('initcfg');
           spm_get_defaults('cmdline', 1);
        end\n
        """

        mscript += 'jobs{1}.spm.util.defs.comp{1}.def = { '+ "'{0}'".format(self.inputs.deformation_field) + ' }'
        mscript += '\n'
        mscript += "jobs{1}.spm.util.defs.comp{2}.id.space = {" + "'{0}'".format(self.inputs.reference_volume) + ' }'
        mscript += '\n'
        mscript += "jobs{1}.spm.util.defs.out{1}.pull.fnames = {" + "'{0}'".format("' \n '".join(self.inputs.in_files)) + ' }'
        mscript += '\n'
        mscript += "jobs{1}.spm.util.defs.out{1}.pull.interp = " + "{0}".format(self.inputs.interp)
        mscript += '\n'

        mscript += """
        spm_jobman(\'run\', jobs);\n
        """
        if self.inputs.use_mcr:
            mscript += """
        if strcmp(name, 'SPM8') || strcmp(name(1:5), 'SPM12'),
            close(\'all\', \'force\');
        end;
            """
        if postscript is not None:
            mscript += postscript
        return mscript

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_files'] = []
        for filename in self.inputs.in_files:
            _, fname = os.path.split(filename)
            outputs['out_files'].append(os.path.realpath('w%s' % fname))
        return outputs


class JacobianInputSpec(SPMCommandInputSpec):
    deformation_file = traits.File(exists=True,
                               field='comp{1}.def',
                               desc='The deformation field to use',mandatory=False)

    reference = traits.File(exists=True,
                            field='comp{2}.id.space',
                            desc='the image specifying write direction, bounding box, voxel sizes (i.e. template)',
                            mandatory=False)

    out_fname = traits.File(default='jacobian.nii',exists=False,usedefault=True,
                                  field='out{1}.savejac.ofname',
                                  desc="write to working directory",
                                  mandatory=False)

    write_to_workdir = traits.Int(default=1,
                                  field='out{1}.savejac.savedir.savepwd',
                                  desc="write to working directory",
                                  usedefault=True,
                                  mandatory=True)
    
    ## usage of Dartel flow
    dartel_flow = traits.File(exists=True,
                               field='comp{1}.dartel.flowfield',
                               desc='The deformation field to use. The same field can be used for both forward and backward deformations.',mandatory=False)
    
    dartel_dir = traits.List(traits.Int(), 
                             field='comp{1}.dartel.times',
                             minlen=2, maxlen=2,
                             desc='The direct of the Dartel flow',mandatory=False)
    
    dartel_timesteps = traits.Int(default=6,
                               field='comp{1}.dartel.K',
                               desc='The number of time points used for solving the partial differential equations.',mandatory=False)
    


class JacobianOutputSpec(TraitedSpec):
    jacobian_map = File(desc = "The deformation's jacobian determinant map image", exists = True)

class Jacobian(SPMCommand):

    input_spec = JacobianInputSpec
    output_spec = JacobianOutputSpec

    def __init__(self, **inputs):
        _local_version = SPMCommand().version
        if _local_version and '12.' in _local_version:
            self._jobtype = 'util'
            self._jobname = 'defs'
        SPMCommand.__init__(self, **inputs)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['jacobian_map']=os.path.abspath('j_'+self.inputs.out_fname)
        return outputs

    def _format_arg(self, opt, spec, val):
        """
        Convert input to appropriate format for spm
        Need cell strings for def and id fields
        """
        if opt in ['deformation_file', 'reference','dartel_flow']:
            new_input = np.array([val], dtype=object)
            return new_input
        else:
            return val


########################################################
##### Custom interface for FSL tools
#######################################################
class fslFDRInputSpec(FSLCommandInputSpec):
    in_file = File(
        exists=True,
        argstr='-i %s',
        mandatory=True,
        position=1,
        desc='p-value file name (3D/4D image file)')
    out_file = File(
        name_template="%s_fdr_adjusted",
        argstr='-a %s',
        desc=('FDR adjusted output file name'),
        name_source="in_file",
        keep_extension=True)
    mask = File(
        exists=True,
        argstr='-m %s',
        desc=('mask image file name'))
    q_thresh = traits.Float(
        argstr='-q %.2f', desc='q-value threshold')
    out_thresh_file = File(
        name_template="%s_fdr_thresholded",
        argstr='--othresh=%s',
        desc='output a thresholded p val image')
    oneminusp = traits.Bool(
        argstr='--oneminusp',
        desc=('treat input as 1-p'))
    
class fslFDROutputSpec(TraitedSpec):
    out_file = File(
        exists=True, desc=('FDR adjusted output file'))
    out_thresh_file = File(
        exists=True, desc=('FDR-thresholded p-val file'))
    
class fslFDR(FSLCommand):
    """ This program runs fdr (FSL)."""
    
    _cmd = 'fdr'
    input_spec = fslFDRInputSpec
    output_spec = fslFDROutputSpec
    
    def _list_outputs(self):
        outputs = self.output_spec().get()
        if isdefined(self.inputs.out_file) and self.inputs.out_file:
            outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        else:
            outputs['out_file'] = glob(self._gen_fname('*_fdr_adjusted.nii'))[0]
        if isdefined(self.inputs.out_thresh_file) and self.inputs.out_thresh_file:
            outputs['out_thresh_file'] = os.path.abspath(self.inputs.out_thresh_file)
       
        return outputs

# Randomise is wrapped in nipype, but output are not ordered.
# So fix that in this custom version.    
class RandomiseInputSpec(FSLCommandInputSpec):
    in_file = File(
        exists=True,
        desc='4D input file',
        argstr='-i %s',
        position=0,
        mandatory=True)
    base_name = traits.Str(
        'randomise',
        desc='the rootname that all generated files will have',
        argstr='-o "%s"',
        position=1,
        usedefault=True)
    design_mat = File(
        exists=True, desc='design matrix file', argstr='-d %s', position=2)
    tcon = File(
        exists=True, desc='t contrasts file', argstr='-t %s', position=3)
    fcon = File(exists=True, desc='f contrasts file', argstr='-f %s')
    mask = File(exists=True, desc='mask image', argstr='-m %s')
    x_block_labels = File(
        exists=True, desc='exchangeability block labels file', argstr='-e %s')
    demean = traits.Bool(
        desc='demean data temporally before model fitting', argstr='-D')
    one_sample_group_mean = traits.Bool(
        desc=('perform 1-sample group-mean test instead of generic '
              'permutation test'),
        argstr='-1')
    show_total_perms = traits.Bool(
        desc=('print out how many unique permutations would be generated '
              'and exit'),
        argstr='-q')
    show_info_parallel_mode = traits.Bool(
        desc='print out information required for parallel mode and exit',
        argstr='-Q')
    vox_p_values = traits.Bool(
        desc='output voxelwise (corrected and uncorrected) p-value images',
        argstr='-x')
    tfce = traits.Bool(
        desc='carry out Threshold-Free Cluster Enhancement', argstr='-T')
    tfce2D = traits.Bool(
        desc=('carry out Threshold-Free Cluster Enhancement with 2D '
              'optimisation'),
        argstr='--T2')
    f_only = traits.Bool(desc='calculate f-statistics only', argstr='--f_only')
    raw_stats_imgs = traits.Bool(
        desc='output raw ( unpermuted ) statistic images', argstr='-R')
    p_vec_n_dist_files = traits.Bool(
        desc='output permutation vector and null distribution text files',
        argstr='-P')
    num_perm = traits.Int(
        argstr='-n %d',
        desc='number of permutations (default 5000, set to 0 for exhaustive)')
    seed = traits.Int(
        argstr='--seed=%d',
        desc='specific integer seed for random number generator')
    var_smooth = traits.Int(
        argstr='-v %d', desc='use variance smoothing (std is in mm)')
    c_thresh = traits.Float(
        argstr='-c %.1f', desc='carry out cluster-based thresholding')
    cm_thresh = traits.Float(
        argstr='-C %.1f', desc='carry out cluster-mass-based thresholding')
    f_c_thresh = traits.Float(
        argstr='-F %.2f', desc='carry out f cluster thresholding')
    f_cm_thresh = traits.Float(
        argstr='-S %.2f', desc='carry out f cluster-mass thresholding')
    tfce_H = traits.Float(
        argstr='--tfce_H=%.2f', desc='TFCE height parameter (default=2)')
    tfce_E = traits.Float(
        argstr='--tfce_E=%.2f', desc='TFCE extent parameter (default=0.5)')
    tfce_C = traits.Float(
        argstr='--tfce_C=%.2f', desc='TFCE connectivity (6 or 26; default=6)')


class RandomiseOutputSpec(TraitedSpec):
    tstat_files = traits.List(
        File(exists=True), desc='t contrast raw statistic')
    fstat_files = traits.List(
        File(exists=True), desc='f contrast raw statistic')
    t_p_files = traits.List(
        File(exists=True), desc='f contrast uncorrected p values files')
    f_p_files = traits.List(
        File(exists=True), desc='f contrast uncorrected p values files')
    t_corrected_p_files = traits.List(
        File(exists=True),
        desc='t contrast FWE (Family-wise error) corrected p values files')
    f_corrected_p_files = traits.List(
        File(exists=True),
        desc='f contrast FWE (Family-wise error) corrected p values files')


class Randomise(FSLCommand):
    """FSL Randomise: feeds the 4D projected FA data into GLM
    modelling and thresholding
    in order to find voxels which correlate with your model
    Example
    -------
    >>> import nipype.interfaces.fsl as fsl
    >>> rand = fsl.Randomise(in_file='allFA.nii', mask = 'mask.nii', tcon='design.con', design_mat='design.mat')
    >>> rand.cmdline
    'randomise -i allFA.nii -o "randomise" -d design.mat -t design.con -m mask.nii'
    """

    _cmd = 'randomise'
    input_spec = RandomiseInputSpec
    output_spec = RandomiseOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['tstat_files'] = glob(
            self._gen_fname('%s_tstat*.nii' % self.inputs.base_name))
        outputs['fstat_files'] = glob(
            self._gen_fname('%s_fstat*.nii' % self.inputs.base_name))
        prefix = False
        if self.inputs.tfce or self.inputs.tfce2D:
            prefix = 'tfce'
        elif self.inputs.vox_p_values:
            prefix = 'vox'
        elif self.inputs.c_thresh or self.inputs.f_c_thresh:
            prefix = 'clustere'
        elif self.inputs.cm_thresh or self.inputs.f_cm_thresh:
            prefix = 'clusterm'
        if prefix:
            outputs['t_p_files'] = sorted(glob(
                self._gen_fname('%s_%s_p_tstat*' % (self.inputs.base_name,
                                                    prefix))))
            outputs['t_corrected_p_files'] = sorted(glob(
                self._gen_fname('%s_%s_corrp_tstat*.nii' %
                                (self.inputs.base_name, prefix))))

            outputs['f_p_files'] = sorted(glob(
                self._gen_fname('%s_%s_p_fstat*.nii' % (self.inputs.base_name,
                                                        prefix))))
            outputs['f_corrected_p_files'] = sorted(glob(
                self._gen_fname('%s_%s_corrp_fstat*.nii' %
                                (self.inputs.base_name, prefix))))
            
        return outputs

########################################################
##### Custom interface for custom-made commandline tools
#######################################################
class CoregQCInputSpec(CommandLineInputSpec):
    coreg_im_file = traits.Str(mandatory=True,
                               desc='Coregistered image to be plotted',
                               argstr='%s',
                               position=1)
    ref_file = traits.Str(mandatory=True,
                            desc ='Ref image file (typically T1 brain) whose edge will be shown',
                            argstr='%s',
                            position=2)
    mask_file = traits.Str(mandatory=True,
                           desc='Brain mask',
                           argstr='%s',
                           position=3)
    out_basename = traits.Str(mandatory=True,
                              desc='Output basename for png and txt',
                              argstr='%s',
                              position=4)

class CoregQCOutputSpec(TraitedSpec):
    out_plot = File(exists=True)
    out_txt = File(exists=True)

class CoregQC(CommandLine):
    """
    Creates multi-slice axial plot with Slicer showing coregistation IsoContours
    and also computes FLIRT costfunction.
    """
    _cmd = 'coreg_QC.sh'
    input_spec = CoregQCInputSpec
    output_spec = CoregQCOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_plot'] = "{}.png".format(os.path.abspath(self.inputs.out_basename))
        outputs['out_txt'] = "{}.txt".format(os.path.abspath(self.inputs.out_basename))
        return outputs


class MaskOverlayQCplotInputSpec(CommandLineInputSpec):
    bg_im_file = traits.Str(mandatory=True,
                            desc ='Background ref image file (typically T1 brain)',
                            argstr='%s',
                            position=1)
    mask_file = traits.Str(mandatory=True,
                           desc='Brain mask',
                           argstr='%s',
                           position=2)
    transparency = traits.Enum(0, 1,
                              argstr='%d',
                              desc='Set transparency (0: solid, 1:transparent)',
                              mandatory=True,
                              position=3)
    out_file = traits.Str('mask_overlay.png',
                           mandatory=False,
                           desc='Output png filename',
                           argstr='%s',
                           position=4)
    bg_max = traits.Float(argstr="%.3f",
                          mandatory=False,
                          desc='Optionally specifies the bg img intensity range as a percentile',
                          position=5)

class MaskOverlayQCplotOutputSpec(TraitedSpec):
    out_file = File(exists=True)

class MaskOverlayQCplot(CommandLine):
    """
    Creates multi-slice axial plot with Slicer showing mask overlaied on
    background image.
    """
    _cmd = 'mask_overlay_QC_images.sh'
    input_spec = MaskOverlayQCplotInputSpec
    output_spec = MaskOverlayQCplotOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_file'] = os.path.abspath(self.inputs.out_file)
        return outputs

 
    
class VbmQCplotInputSpec(CommandLineInputSpec):
    bg_file = traits.Str(mandatory=True,
                         desc ='background image file to plot',
                         argstr='%s',
                         position=1)
    
    gm_file = traits.Str(mandatory=True,
                         desc ='GM image file to overlay',
                         argstr='%s',
                         position=2)
    
    wm_file = traits.Str(mandatory=True,
                         desc ='WM image file to overlay',
                         argstr='%s',
                         position=3)
    
    use_jet = traits.Enum(0, 1,
                          mandatory=False,
                          usedefault=True,
                          desc='Use jet colorscheme',
                          argstr='%d',
                          position=4)
    
    out_name = traits.Str('VBM',
                          mandatory=False,
                          usedefault=True,
                          desc='Prefix to the output png filename',
                          argstr='%s',
                          position=5)
    
    bg_max = traits.Float(99.9,
                          argstr="%.3f",
                          mandatory=False,
                          usedefault=True,
                          desc='Optionally specifies the bg img intensity range as a percentile',
                          position=6)
    
    gm_max = traits.Float(1.0,
                          argstr="%.3f",
                          mandatory=False,
                          usedefault=True,
                          desc='Optionally specifies the max val of overlaying GM',
                          position=7)
    
    wm_max = traits.Float(1.0,
                          argstr="%.3f",
                          mandatory=False,
                          usedefault=True,
                          desc='Optionally specifies the max val of overlaying WM',
                          position=8)


class VbmQCplotOutputSpec(TraitedSpec):
    output_images = OutputMultiPath(File(exists=True), desc='QC images in png format')

class VbmQCplot(CommandLine):
    """
    Plots mid-saggital slices per frame for viewing slice-to-volume
    movement artefacts in DWI data.
    """
    _cmd = 'create_vbm_QC_images.sh'
    input_spec = VbmQCplotInputSpec
    output_spec = VbmQCplotOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['output_images']= []
        outputs['output_images'].append(os.path.abspath('{}_gm_overlay_sag.png'.format(self.inputs.out_name)))
        outputs['output_images'].append(os.path.abspath('{}_gm_overlay_cor.png'.format(self.inputs.out_name)))
        outputs['output_images'].append(os.path.abspath('{}_gm_overlay_ax.png'.format(self.inputs.out_name)))
        outputs['output_images'].append(os.path.abspath('{}_wm_overlay_sag.png'.format(self.inputs.out_name)))
        outputs['output_images'].append(os.path.abspath('{}_wm_overlay_cor.png'.format(self.inputs.out_name)))
        outputs['output_images'].append(os.path.abspath('{}_wm_overlay_ax.png'.format(self.inputs.out_name)))
        
        return outputs

