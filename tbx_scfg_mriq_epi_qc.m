function epiqc = tbx_scfg_mriq_epi_qc
%==========================================================================
% (Sub)configuration file for the MRIquality toolbox, partim EPI/QC.
% Quality control for EPI data currently includes:
% - sequential check: view a series of data to visually detect any obvious
%   artefact or big subject motion.
% - spike check (and correct) - could be extended for outlier detection 
%   (in DWI data from Allegra scanner, some DWI volumes came out oddly
%   signal-less)...  
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%==========================================================================
% Written by Evelyne Balteau - 2013-2018
% Cyclotron Research Centre, University of Liege, Belgium
% e.balteau@uliege.be
%==========================================================================

% ---------------------------------------------------------------------
% EPI images
% ---------------------------------------------------------------------
EPIimages           = cfg_files;
EPIimages.tag       = 'EPIimages';
EPIimages.name      = 'EPI images';
EPIimages.help      = {'Select EPI images'
    'Input the whole series of EPI images (from fMRI or DWI acquistions).'};
EPIimages.filter    = 'image';
EPIimages.ufilter   = '.*';
EPIimages.num       = [1 Inf];
EPIimages.val       = {''};

% ---------------------------------------------------------------------
% outdir Output directory for summary results
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.val{1}  = {''};
outdir.help    = {'Output directory for results.'
    ['If no directory is given, summary results will be saved ' ...
    'in the EPI images input directory.']};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [0 1];


% % % % % ---------------------------------------------------------------------
% % % % % parameters for sequential check - centre of windowing
% % % % % ---------------------------------------------------------------------
% % % % seqcentre         = cfg_entry;
% % % % seqcentre.tag     = 'seqcentre';
% % % % seqcentre.name    = 'Centre for windowing';
% % % % seqcentre.help    = {
% % % %                      'The images will be displayed sequentially with a constant windowing defined by a centre and width. You can adjust these values to improve the efficiency of the visual inspection of your data, and re-run the quality control as many times as necessary.'
% % % %                      'Typical values:'
% % % %                      '- for background visualisation (recommended): centre = 200 / width = 200'
% % % %                      '- for normal visualisation: centre = 700 / width = 700.'
% % % %                     }';
% % % % seqcentre.strtype = 'e';
% % % % seqcentre.num     = [1 1];
% % % % seqcentre.val     = {50}; 
% % % % 
% % % % % ---------------------------------------------------------------------
% % % % % parameters for sequential check - centre of windowing
% % % % % ---------------------------------------------------------------------
% % % % seqwidth         = cfg_entry;
% % % % seqwidth.tag     = 'seqwidth';
% % % % seqwidth.name    = 'Width for windowing';
% % % % seqwidth.help    = {
% % % %                      'The images will be displayed sequentially with a constant windowing defined by a centre and width. You can adjust these values to improve the efficiency of the visual inspection of your data, and re-run the quality control as many times as necessary.'
% % % %                      'Typical values:'
% % % %                      '- for background visualisation (recommended): centre = 200 / width = 200'
% % % %                      '- for normal visualisation: centre = 700 / width = 700.'
% % % %                     }';
% % % % seqwidth.strtype = 'e';
% % % % seqwidth.num     = [1 1];
% % % % seqwidth.val     = {50}; 
% % % % 
% % % % % ---------------------------------------------------------------------
% % % % % parameters for sequential check - inter-frame delay and rate
% % % % % ---------------------------------------------------------------------
% % % % seqrate         = cfg_entry;
% % % % seqrate.tag     = 'seqrate';
% % % % seqrate.name    = 'Number of frames per second';
% % % % seqrate.help    = {'To determine the rate of scrolling through the images. Enter 0 if you want to scroll manually.'};
% % % % seqrate.strtype = 'e';
% % % % seqrate.num     = [1 1];
% % % % seqrate.val     = {5}; 
% % % % 
% % % % % ---------------------------------------------------------------------
% % % % % parameters for sequential check
% % % % % ---------------------------------------------------------------------
% % % % params_seq         = cfg_branch;
% % % % params_seq.tag     = 'params_seq';
% % % % params_seq.name    = 'Processing parameters';
% % % % params_seq.val     = {seqcentre seqwidth seqrate};
% % % % params_seq.help    = {'Processing parameters for sequential check'};

% ---------------------------------------------------------------------
% parameters for spike detect and correct - do try and correct or not?
% ---------------------------------------------------------------------
do_correct         = cfg_menu;
do_correct.tag     = 'do_correct';
do_correct.name    = 'Spike correction';
do_correct.help    = {'Whether spike correction is performed (if possible) or not...'};
do_correct.labels  = {'No'
                      'Yes'
                      };
do_correct.values  = {0 1};
do_correct.val     = {0}; 

% ---------------------------------------------------------------------
% parameters for spike detect and correct - threshold for correction
% ---------------------------------------------------------------------
correction_threshold         = cfg_entry;
correction_threshold.tag     = 'correction_threshold';
correction_threshold.name    = 'Correction threshold';
correction_threshold.help    = {['Threshold determining whether a data set ' ...
    'affected by spikes can or cannot be corrected.']
    ['Is considered as "spike correctable" a data set where the noise ' ...
    'variation is not greater than corrThresh*(average noise level). ' ...
    'Default value = 0.10. Increase this value if you want to force ' ...
    'the spike correction.']};
correction_threshold.strtype = 'e';
correction_threshold.num     = [1 1];
correction_threshold.val     = {0.1}; 

% ---------------------------------------------------------------------
% parameters for spike detect and correct - threshold for spike detection
% ---------------------------------------------------------------------
spike_threshold         = cfg_entry;
spike_threshold.tag     = 'spike_threshold';
spike_threshold.name    = 'Spike threshold';
spike_threshold.help    = {'Threshold determining whether a slice is spiky or not. '
    ['Is considered as spiky a slice where the noise level is higher ' ...
    'than spikeThresh*(standard deviation of noise). Default value = 4.']};
spike_threshold.strtype = 'e';
spike_threshold.num     = [1 1];
spike_threshold.val     = {4}; 

% ---------------------------------------------------------------------
% parameters for spike detect and correct
% ---------------------------------------------------------------------
params_spike         = cfg_branch;
params_spike.tag     = 'params_spike';
params_spike.name    = 'Spike detect&correct parameters';
params_spike.val     = {do_correct correction_threshold spike_threshold};
params_spike.help    = {'Parameters for spike detection and correction'};

% ---------------------------------------------------------------------
% sequential_check Check sequentially the EPI series by displaying the 
% images as mosaic with interframe delay deltat and windowing and
% centring centwin = [centre width]
% ---------------------------------------------------------------------
sequential_check        = cfg_exbranch;
sequential_check.tag    = 'qc_check_sequential';
sequential_check.name   = 'Sequential check';
sequential_check.val    = {EPIimages}; % {params_seq};
sequential_check.help   = {'To visually check for spikes, artefacts and movements.'
    ['Load a set of EPI images and run the batch. ' ...
    'The images are displayed as mosaic. ' ...
    'The interactive display window allows you to scroll through ' ...
    'the EPI series and adjust the windowing to properly detect spikes ' ...
    'head motion and artefacts (just as you should have done at the ' ...
    'scanner during or straight after data acquisition!!).']};
sequential_check.prog   = @qc_check_sequential;

% ---------------------------------------------------------------------
% spike_check Spike check & correct (optional)
% ---------------------------------------------------------------------
spike_check          = cfg_exbranch;
spike_check.tag      = 'qc_spike_check';
spike_check.name     = 'Spike check (and correct)';
spike_check.val      = {EPIimages outdir params_spike};
spike_check.help     = {
    'To automatically detect spikes and correct if possible (correction is optional).'
    ['WARNING: only apply correction to data that are affected by spikes. ' ...
    'Use this tool to check whether there are spikes indeed and how bad ' ...
    'they are affecting your data. If your data are affected by spikes, ' ...
    'check carefully the way the correction has succeeded (or not) before ' ...
    'you decide to work with corrected data. Using badly corrected data ' ...
    '(e.g. because too many spikes and resulting poor correction) will ' ...
    'affect your results, too. ']
    'Please report any bugs or problems to e.balteau@uliege.be.'
    };
spike_check.prog     = @qc_spike_check;
spike_check.vout     = @vout_spike_check;

% ---------------------------------------------------------------------
% epi QC tools for EPI acquisitions
% ---------------------------------------------------------------------
epiqc         = cfg_choice;
epiqc.tag     = 'check_epi';
epiqc.name    = 'Quality control';
epiqc.help    = {
                'Quality control tools for EPI images (including fMRI and diffusion imaging).'
              }';
epiqc.values  = {sequential_check spike_check};

end


function dep = vout_spike_check(job)
% To manage dependencies for next steps in processing pipeline.
% If spike correction is enabled by user (do_correct==1), spike corrected
% images can be used at the next step. 
% WARNINGS:
% 1- this shouldn't save you to check the spike corrected data before you
% decide to include them into your analysis!!!
% 2- if the correction is not possible, the original images are passed as
% outputs to the next step!

dep = [];
if job.params_spike.do_correct
    % create new dependency:
    dep            = cfg_dep;
    % name it:
    dep.sname      = 'Spike Corrected Images';
    % substructure of the output containing the list of generated
    % files (it'll be in out.spcfiles):
    dep.src_output = substruct('.','spcfiles');
    % filter the list of generated files to keep images only:
    dep.tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end

end