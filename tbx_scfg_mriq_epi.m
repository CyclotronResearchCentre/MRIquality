function epi = tbx_scfg_mriq_epi
% (Sub)configuration file for the MRIquality toolbox, partim EPI data.
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%==========================================================================
% Written by Evelyne Balteau, Cyclotron Research Centre, June 2017
%==========================================================================

% The following subsections are available:
% - tbx_scfg_mriq_epi_qa: quality assurance of EPI data, i.e. analysis of
%   EPI series acquired (mainly) on phantom, checking for SNR and stability
%   parameters to detect early any deviation from good quality.  
% - tbx_scfg_mriq_epi_qc: quality control of  utility scripts and procesing including
%   g-maps, noise correlation, ... 

% ---------------------------------------------------------------------
% hmri hMRI Tools
% ---------------------------------------------------------------------
epi         = cfg_choice;
epi.tag     = 'epi';
epi.name    = 'EPI Tools';
epi.help    = {
    ['Coolection of tools for quality assurance and quality control ', ...
    'of EPI data.']
    }';
epi.values  = {tbx_scfg_mriq_epi_qa};% tbx_scfg_mriq_epi_qc};
end

