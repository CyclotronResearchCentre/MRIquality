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
EPIseries.help      = {'Select EPI images (either DICOM or nifti format)'};
EPIseries.dir       = mriq_get_defaults('path_input');
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
NOISEseries.dir       = mriq_get_defaults('path_input');
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
% Archiving options
%==========================================================================
archive        = cfg_menu;
archive.tag    = 'archive';
archive.name   = 'Archiving';
archive.help   = {'1. Enabled:', ...
    ['Original input data are compressed and archived ' ...
    '(see archiving directory in the [Configure] module of the toolbox). ' ...
    'Intermediate files are cleaned up and results are saved in the output ' ...
    'directory (see output directory in the [Configure] module of the toolbox).'], ...
    '2. Disabled:', ...
    'Original files, intermediate files and results are kept untouched.'};
archive.labels = {
               'Enabled'
               'Disabled'}';
archive.values = {1 0};
archive.val    = {1};



%==========================================================================
% Number of volumes to discard
%==========================================================================
dummies         = cfg_entry;
dummies.tag     = 'dummies';
dummies.name    = 'Dummies';
dummies.help    = {'Enter the number of dummy volumes. ', ...
    ['To discard the first few volumes of the time series, ' ...
    'and avoid T1 saturation effects impacting the QA results.']};
dummies.strtype = 'i';
dummies.num     = [1 1];
dummies.val     = {0};

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
procpar.val       = {dummies sigplane noiplane roisize comment};


%==========================================================================
% Tools for EPI quality assurance
%==========================================================================
epiqa         = cfg_exbranch;
epiqa.tag     = 'epiqa';
epiqa.name    = 'Quality assurance';
epiqa.val     = {series archive procpar};
epiqa.help    = {'Tools for EPI quality assurance, including SNR and stability estimates.'};
epiqa.prog    = @mriq_run_epi_qa;
epiqa.vout    = @vout_epi_qa; 

end


%% VOUT & OTHER SUBFUNCTIONS
% ========================================================================
% The RUN function:
% - out = hmri_run_create(job)
% is defined separately.
%_______________________________________________________________________

function dep = vout_epi_qa(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

% All output from a given dataset have a common filename formed as follows:
% <Scanner>_<DateYYYYMMDD>_stud0000_ser0000 (see the asterisk * below). The
% file name is appended with suffix and the file extension varies according
% to the output. 
% For each run of the QA batch, the following output are produced:  
% - *.json file with the results
% - *.txt file with same information but a bit easier to read
% - *.png screenshot of the general results window
% - *_batch.mat job saved for the current run
% - *_MEAN.nii mean volume
% - *_NOISE_DISTRIB.png screenshot of the noise distribution figure
% - *_SD.nii standard deviation volume
% - *_SNR_noRF.nii SNR volume (calculated with noise image)
% - *_tSNR.nii tSNR volume

dep = [];
% k=1;
% cdep(1,5*numel(job.subj)) = cfg_dep;
% for i=1:numel(job.subj)
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('R1_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','R1','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('R2s_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','R2s','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('MT_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','MT','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('A_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','A','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('T1w_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','T1w','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('MTw_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','MTw','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('PDw_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','PDw','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
% end
% dep = cdep;
    
end