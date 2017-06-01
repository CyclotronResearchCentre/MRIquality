function mriq = tbx_cfg_mriq
% Configuration file for the MRIquality toolbox including QA procedures and
% processing for EPI stability and coil specifications.
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
%   and SNR parameters (SNR, tSNR, RDC, ...)
% - tbx_scfg_mriq_coil: coil utility scripts and procesing including
%   g-maps, noise correlation, ... 

% ---------------------------------------------------------------------
% hmri hMRI Tools
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
    'future updates.  Please report any bugs or problems to e.balteau@ulg.ac.be.']
    }';
mriq.values  = {tbx_scfg_mriq_epi};% tbx_scfg_mriq_coil};
end

