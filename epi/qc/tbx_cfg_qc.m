%=========================================================================%
% This file is part of the Quality Control Toolbox (TCQ)
% Copyright (C) 2013 - Cyclotron Research Centre
% University of Liege, Belgium
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%=========================================================================%

function qc = tbx_cfg_qc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: configuration file for the Quality Control Toolbox (QCT).
% This toolbox is a collection of quality control tools for MRI
% acquisitions. 
%--------------------------------------------------------------------------
% Warning and disclaimer: This software is for research use only. 
% Do not use it for clinical or diagnostic purposes.
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - 2013
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENT CONTENT
% QCT - EPI - sequential check
% QCT - EPI - spike check (and correct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','QCT')); end

% ---------------------------------------------------------------------
% raws Raw Images
% ---------------------------------------------------------------------
raws           = cfg_files;
raws.tag       = 'raws';
raws.name      = 'EPI images';
raws.help      = {'Input the whole series of EPI images (from fMRI or DWI acquistions).'}; 
raws.filter    = 'image';
raws.ufilter   = '.*';
raws.num       = [1 Inf];

% ---------------------------------------------------------------------
% subj Subject
% --------------------------------------------------------------------
subj           = cfg_branch;
subj.tag       = 'subj';
subj.name      = 'Subject';
subj.val       = {raws};
subj.help      = {'Specify a subject.'};

% ---------------------------------------------------------------------
% data Data
% ---------------------------------------------------------------------
sdata          = cfg_repeat;
sdata.tag      = 'data';
sdata.name     = 'Data';
sdata.val      = {subj};
sdata.help     = {'Specify the number of subjects.'};
sdata.values   = {subj};
sdata.num      = [1 Inf];

% ---------------------------------------------------------------------
% indir Input directory as output directory
% ---------------------------------------------------------------------
indir          = cfg_menu;
indir.tag      = 'indir';
indir.name     = 'Input directory';
indir.help     = {'Output files will be written to the same folder as each corresponding input file.'};
indir.labels   = {'Yes'};
indir.values   = {1};
indir.val      = {1};

% ---------------------------------------------------------------------
% outdir Output directory
% ---------------------------------------------------------------------
outdir         = cfg_files;
outdir.tag     = 'outdir';
outdir.name    = 'Output directory';
outdir.help    = {'Select a directory where output files will be written to.'};
outdir.filter  = 'dir';
outdir.ufilter = '.*';
outdir.num     = [1 1];

% ---------------------------------------------------------------------
% output Output choice
% ---------------------------------------------------------------------
output         = cfg_choice;
output.tag     = 'output';
output.name    = 'Output choice';
output.help    = {'Output directory can be the same as the input directory for each input file or user selected.'};
output.values  = {indir outdir};
output.val  = {indir};

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
correction_threshold.help    = {'Threshold determining whether a data set affected by spikes can or cannot be corrected. Is considered as "spike correctable" a data set where the noise variation is not greater than corrThresh*(average noise level). Default value = 0.10. Increase this value if you want to force the spike correction.'};
correction_threshold.strtype = 'e';
correction_threshold.num     = [1 1];
correction_threshold.val     = {0.1}; 

% ---------------------------------------------------------------------
% parameters for spike detect and correct - threshold for spike detection
% ---------------------------------------------------------------------
spike_threshold         = cfg_entry;
spike_threshold.tag     = 'spike_threshold';
spike_threshold.name    = 'Spike threshold';
spike_threshold.help    = {'Threshold determining whether a slice is spiky or not. Is considered as spiky a slice where the noise level is higher than spikeThresh*(standard deviation of noise). Default value = 4.'};
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
subj.val                = {raws};
sdata.val               = {subj};
sdata.values            = {subj};
% sequential_check.val    = {sdata params_seq};
sequential_check.val    = {sdata};
sequential_check.help   = {'To visually check for spikes, artefacts and movements. Load a set of EPI images and run the batch. The images are displayed as mosaic. The interactive display window allows you to scroll through the EPI series and adjust the windowing to properly detect spikes and artefacts (just as you should have done at the scanner straight after acquiring the data!!).'};
sequential_check.prog   = @qc_check_sequential;

% ---------------------------------------------------------------------
% spike_check Spike check & correct (optional)
% ---------------------------------------------------------------------
spike_check          = cfg_exbranch;
spike_check.tag      = 'qc_spike_check';
spike_check.name     = 'Spike check (and correct)';
subj.val             = {output raws};
sdata.val            = {subj};
sdata.values         = {subj};
spike_check.val      = {sdata params_spike};
spike_check.help     = {
    'To automatically detect spikes and correct if possible (correction is optional).'
    ''
    'WARNING: only apply correction to data that are affected by spikes. Use this tool to check whether there are spikes indeed and how bad they are affecting your data. If your data are affected by spikes, check carefully the way the correction has succeeded (or not) before you decide to work with corrected data. Using badly corrected data (e.g. because too many spikes and resulting poor correction) will affect your results, too. Please report any bugs or problems to e.balteau@ulg.ac.be.'
    };
spike_check.prog     = @qc_spike_check;
spike_check.vout     = @vout_spike_check;

% ---------------------------------------------------------------------
% epi QC tools for EPI acquisitions
% ---------------------------------------------------------------------
epi         = cfg_choice;
epi.tag     = 'check_epi';
epi.name    = 'EPI series';
epi.help    = {
                'Quality control tools for EPI images (including fMRI and diffusion imaging).'
              }';
epi.values  = {sequential_check spike_check};

% ---------------------------------------------------------------------
% QC Toolbox
% ---------------------------------------------------------------------
qc         = cfg_choice;
qc.tag     = 'QCT';
qc.name    = 'Quality Control Toolbox';
qc.help    = {
                  'This toolbox allows you to easily control the quality of your data.'
                  'This toolbox should be considered as a beta (trial) version, and will include a number of (as yet unspecified) extensions in future updates. Please report any bugs or problems to e.balteau@ulg.ac.be.'
              }';
qc.values  = {epi};
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
    for k=1:numel(job.subj)
        % create new dependency:
        cdep(1)            = cfg_dep;
        % name it:
        cdep(1).sname      = sprintf('Spike Corrected Images (Subj %d)', k); 
        % substructure of the output containing the list of generated
        % files (it'll be in out.subj(k).spcfiles): 
        cdep(1).src_output = substruct('.','subj','()',{k},'.','spcfiles');
        % filter the list of generated files to keep images only:
        cdep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
        % append to the list of dependencies
        dep = [dep cdep];
    end
end

end