function mriq = tbx_cfg_mriq
% Configuration file for the MRIquality toolbox.
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%==========================================================================
% Written by Evelyne Balteau, Cyclotron Research Centre, June 2017
%==========================================================================

if ~isdeployed, addpath(genpath(fileparts(mfilename('fullpath')))); end

% The following subsections are available:
% - tbx_scfg_mriq_epi: analysis of EPI series, including stability
%   and SNR parameters (SNR, tSNR, RDC, ...), detection of DWI outliers,
%   and more...
% - tbx_scfg_mriq_coil: coil utility scripts for coil quality assessment,
%   including g-maps, noise correlation, ... 
% - tbx_scfg_mriq_struct: nothing implemented yet here...
% - tbx_scfg_mriq_visual: for systematic visual check -> refer to the hMRI
%   implementation...

% ---------------------------------------------------------------------
% mriq Tools
% ---------------------------------------------------------------------
mriq         = cfg_choice;
mriq.tag     = 'mriq';
mriq.name    = 'MRI Quality';
mriq.help    = {
    ['This toolbox has been initially implemented for general QA and QC ', ...
    'procedures at the Cyclotron Research Centre, Liege. Feel free to ', ... 
    'adapt and use it for your own research.']
    ['This toolbox should be considered as only a beta (trial) version, ',...
    'and will include a number of (as yet unspecified) extensions in ',...
    'future updates. Please report any bugs or problems to e.balteau@uliege.be.']
    }';
mriq.values  = {tbx_scfg_mriq_epi};% tbx_scfg_mriq_coil};
end

