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

% defaultoutdir = mriq_get_defaults('epi_qa.paths.stab');
% defaultarchiv = mriq_get_defaults('epi_qa.paths.arch');


%==========================================================================
% Unprocessed EPI series: can be either in DICOM or nii format. Nifti data
% are axpected to be associated with JSON metadata
%==========================================================================
EPIseries           = cfg_files;
EPIseries.tag       = 'EPIseries';
EPIseries.name      = 'EPI time series';
EPIseries.help      = {['Select EPI images (either DICOM or nifti format), ' ...
    'discarding the first few volumes to avoid T1 saturation effect.']};
EPIseries.dir       = mriq_get_defaults('epi_qa.paths.input');
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
NOISEseries.dir       = mriq_get_defaults('epi_qa.paths.input');
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
% outdir Output directory for summary results
%==========================================================================
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
% outdir.val{1}  = {''};
outdir.help    = {'Output directory for summary results.'
    ['If no directory is given, summary results will be saved ' ...
    'in a ''stab'' directory located within the input directory.']
    ['Note that all other (non summary) results will be saved in the ' ...
    '<input directory>/stab directory.']};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];
outdir.def     = @(val)mriq_get_defaults('epi_qa.paths.stab',val{:}); % NOTE: mriq_get_defaults must return a cellstr!

%==========================================================================
% Archiving options: OFF
%==========================================================================
archOFF         = cfg_entry;
archOFF.tag     = 'archOFF';
archOFF.name    = 'No archiving/compression';
archOFF.help    = {'Original files, intermediate files and results are kept untouched.'};
archOFF.strtype = 's';
archOFF.num     = [1 Inf];
archOFF.val     = {'Archiving disabled'};

%==========================================================================
% Archiving options: ON
%==========================================================================
archON          = cfg_files;
archON.tag      = 'archON';
archON.name     = 'Archiving enabled';
archON.help     = {'Select a directory where original files will be compressed/archived.'};
archON.filter   = 'dir';
archON.ufilter  = '.*';
archON.num      = [1 1];
archON.def      = @(val)mriq_get_defaults('epi_qa.paths.arch', val{:}); % NOTE: mriq_get_defaults must return a cellstr!


%==========================================================================
% Archiving options
%==========================================================================
archive         = cfg_choice;
archive.tag     = 'archive';
archive.name    = 'Archiving options';
archive.help    = {['If enabled, the input images (NIFTI or DICOM format) ' ...
    'are tar/gz compressed and saved into the selected archiving directory. ' ...
    'All intermediate files are deleted.'], ...
    ['If archiving is disabled, original images, intermediate files and ' ...
    'results will be kept untouched.']};
archive.values  = {archON archOFF};
archive.val = {archON};

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
epiqa.name    = 'Quality assurance';
epiqa.val     = {series outdir archive procpar};
epiqa.help    = {'Tools for EPI quality assurance, including SNR and stability estimates.'};
epiqa.prog    = @mriq_run_epi_qa;
% epiqa.vout    = @vout_epi_qa; % not implemented

end