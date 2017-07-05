function epiqa = tbx_scfg_mriq_epi_qa
% (Sub)configuration file for the MRIquality toolbox, partim EPI/QA.
% Quality assurance includes SNR and stability estimates.
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%==========================================================================
% Written by Evelyne Balteau, Cyclotron Research Centre, June 2017
%==========================================================================

%==========================================================================
% Unprocessed EPI series: can be either in DICOM or nii format. Nifti data
% are axpected to be associated with JSON metadata
%==========================================================================
EPIseries           = cfg_files;
EPIseries.tag       = 'EPIseries';
EPIseries.name      = 'EPI time series';
EPIseries.help      = {['Select EPI images (either DICOM or nifti format), ' ...
    'discarding the first few volumes to avoid T1 saturation effect.']};
EPIseries.filter    = '^*\.(IMA|dcm|nii)$';
EPIseries.ufilter   = '.*';
EPIseries.num       = [0 Inf];
EPIseries.val       = {''};

%==========================================================================
% Noise images if available (EPI volumes acquired without RF)
%==========================================================================
NOISEseries           = cfg_files;
NOISEseries.tag       = 'NOISEseries';
NOISEseries.name      = 'Noise images';
NOISEseries.help      = {'EPI images acquired without RF.', ...
    ['The most reliable way to estimate SNR is to estimate noise from pure ' ...
    'noise images, i.e. acquired without RF pulses. If not available, the ' ...
    'noise level will be estimated either on a noise slice (i.e. a slice ' ...
    'located out of the phantom) if any, or gathering out-of-phantom voxels ' ...
    'based on volume masking. The former is only appropriate for single-slice ' ...
    'acquisitons, to avoid any contamination of the noise slice by signal coming ' ...
    'from other slices. The latter can be strongly biased by artefacts ' ...
    'superimposed on the background noise.']};
NOISEseries.filter    = '^*\.(IMA|dcm|nii)$';
NOISEseries.ufilter   = '.*';
NOISEseries.num       = [0 Inf];
NOISEseries.val       = {''};

%==========================================================================
% Images
%==========================================================================
series              = cfg_branch;
series.tag       = 'series';
series.name      = 'Input images';
series.help      = {['Input EPI images (time series and noise images if ' ...
    'available) either in DICOM or NIFTI format.']};
series.val       = {EPIseries NOISEseries};

%==========================================================================
% Output directory for summary results
%==========================================================================
indir         = cfg_entry;
indir.tag     = 'indir';
indir.name    = 'Input directory';
indir.help    = {'Output files will be written to the same folder ',...
    'as each corresponding input file.'};
indir.strtype = 's';
indir.num     = [1 Inf];
indir.val     = {'yes'};

%==========================================================================
% outdir Output directory
%==========================================================================
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'User-defined output directory';
outdir.help    = {'Select a directory where output files will be written to.'};
outdir.filter = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

%==========================================================================
% output Output choice
%==========================================================================
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output directory';
output.help    = {'Output directory for summary results.', ...
    ['The output directory can be either the input directory or ' ...
    'any user-refined directory. In the former case, results will be saved ' ...
    'in a ''stab'' directory located within the input directory.'], ...
    ['Note that all other results will be saved in the ' ...
    '<input directory>/stab directory.']};
output.values  = {indir outdir };
output.val = {indir};


%==========================================================================
% Signal plane
%==========================================================================
sigplane         = cfg_entry;
sigplane.tag     = 'sigplane';
sigplane.name    = 'Signal plane';
sigplane.help    = {'Enter the number of the slice used to estimate the signal level.', ...
    ['Typically, the slice at the centre of the phantom is used. If the number ' ...
    'is set to 0, the middle slice of the volume will be used.']};
sigplane.strtype = 'i';
sigplane.num     = [1 1];
sigplane.val     = {0};

%==========================================================================
% Noise plane
%==========================================================================
noiplane         = cfg_entry;
noiplane.tag     = 'noiplane';
noiplane.name    = 'Noise plane';
noiplane.help    = {'Enter the number of the slice used to estimate the noise level.', ...
    ['Typically, if a noise slice has been acquired, it is the outmost (last) ' ...
    'one in the volume. If the noise plane number is set to 0, it is assumed ' ...
    'that no noise plane has been acquired. The noise level ' ...
    'will be estimated based on noise series (data acquired without ' ...
    'RF, if available) or ''background'' voxels selected using a ' ...
    'mask. The latter method is the least robust since it can be strongly biased ' ...
    'by artefacts contaminating the background voxels.']};
noiplane.strtype = 'i';
noiplane.num     = [1 1];
noiplane.val     = {0};

%==========================================================================
% Noise plane
%==========================================================================
roisize         = cfg_entry;
roisize.tag     = 'roisize';
roisize.name    = 'ROI size';
roisize.help    = {['Enter the length (number of voxels) of rectangular ROI ' ...
    '(central ROI within signal plane) for quantitative ROI analysis.']};
roisize.strtype = 'i';
roisize.num     = [1 1];
roisize.val     = {21};

%==========================================================================
% Comment
%==========================================================================
comment         = cfg_entry;
comment.tag     = 'comment';
comment.name    = 'Comment';
comment.help    = {'Enter a comment or short description of the data analysed.'};
comment.strtype = 's';
comment.num     = [0 Inf];
comment.val     = {''};

%==========================================================================
% Processing parameters
% - noise slice number
% - signal slice number
% - comment
% (- tag: will be defined based on series number)
% - N_max: max length of rectangular ROI (central ROI within signal plabe)
% for quantitative ROI analysis 
%==========================================================================
procpar              = cfg_branch;
procpar.tag       = 'procpar';
procpar.name      = 'Processing parameters';
procpar.help      = {['List of input parameters that can be filled in by ' ...
    'the user. Note that not all the parameters need to be specified. Read ' ...
    'each parameter''s help for details.']};
procpar.val       = {sigplane noiplane roisize comment};


%==========================================================================
% Tools for EPI quality assurance
%==========================================================================
epiqa         = cfg_exbranch;
epiqa.tag     = 'epiqa';
epiqa.name    = 'Quality assurance tools for EPI';
epiqa.val     = {series output procpar};
epiqa.help    = {'Tools for EPI quality assurance, including SNR and stability estimates.'};
epiqa.prog    = @mriq_run_epi_qa;
% epiqa.vout    = @vout_create;

end
%----------------------------------------------------------------------


% ========================================================================
%% VOUT & OTHER SUBFUNCTIONS
% ========================================================================
% The RUN function:
% - out = hmri_run_create(job)
% is defined separately.
%_______________________________________________________________________

function dep = vout_create(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

if ~isfield(job, 'subj') % Many subjects
    dep(1) = cfg_dep;
    dep(1).sname = 'R1 Maps';
    dep(1).src_output = substruct('.','R1','()',{':'});
    dep(1).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    
    dep(2) = cfg_dep;
    dep(2).sname = 'R2s Maps';
    dep(2).src_output = substruct('.','R2s','()',{':'});
    dep(2).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    
    dep(3) = cfg_dep;
    dep(3).sname = 'MT Maps';
    dep(3).src_output = substruct('.','MT','()',{':'});
    dep(3).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    
    dep(4) = cfg_dep;
    dep(4).sname = 'A Maps';
    dep(4).src_output = substruct('.','A','()',{':'});
    dep(4).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    
    dep(5) = cfg_dep;
    dep(5).sname = 'T1w Maps';
    dep(5).src_output = substruct('.','T1w','()',{':'});
    dep(5).tgt_spec = cfg_findspec({{'filter','image','strtype','e'}});
    
else
    k=1;
    cdep(5*numel(job.subj),1) = cfg_dep;
    for i=1:numel(job.subj)
        
        cdep(k)            = cfg_dep;
        cdep(k).sname      = sprintf('R1_subj%d',i);
        cdep(k).src_output = substruct('.','subj','()',{i},'.','R1','()',{':'});
        cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        
        k=k+1;
        
        cdep(k)            = cfg_dep;
        cdep(k).sname      = sprintf('R2s_subj%d',i);
        cdep(k).src_output = substruct('.','subj','()',{i},'.','R2s','()',{':'});
        cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        
        k=k+1;
        
        cdep(k)            = cfg_dep;
        cdep(k).sname      = sprintf('MT_subj%d',i);
        cdep(k).src_output = substruct('.','subj','()',{i},'.','MT','()',{':'});
        cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        
        k=k+1;
        
        cdep(k)            = cfg_dep;
        cdep(k).sname      = sprintf('A_subj%d',i);
        cdep(k).src_output = substruct('.','subj','()',{i},'.','A','()',{':'});
        cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        
        k=k+1;
        
        cdep(k)            = cfg_dep;
        cdep(k).sname      = sprintf('T1w_subj%d',i);
        cdep(k).src_output = substruct('.','subj','()',{i},'.','T1w','()',{':'});
        cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        
        k=k+1;
    end
    dep = cdep;
    
end
end
%_______________________________________________________________________

function c = unlimit(c)
try
    if isa(c, 'cfg_files')
        c.num = [0 Inf];
    end
catch e %#ok<*NASGU>
end
try
    for i=1:numel(c.val)
        c.val{i} = unlimit(c.val{i});
    end
catch e
end
end
%_______________________________________________________________________
