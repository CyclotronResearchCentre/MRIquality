function mriq = tbx_cfg_mriq_redirect
% Diversion config file for the MRI quality toolbox
%
% PURPOSE
% The present config file redirects towards the full implementation of
% the toolbox (tbx_cfg_mriq) only if present in the Matlab path. The
% toolbox can therefore be stored in a directory independent from the SPM
% implementation and synch'd with the main MRIquality repository whenever
% needed. If tbx_cfg_mriq is not found in the Matlab path, the MRI quality 
% tools are listed in the SPM Batch GUI but not available. A brief help
% section provides the user with instructions for MRIq Tools installation.
%
% USAGE
% Copy this file into a directory in the SPM toolbox directory (e.g.
% <path-to-SPM>/toolbox/MRIquality). Add the MRIquality directory
% (containing the full implementation of the toolbox) to the Matlab path.
% Restart SPM and the Batch GUI. The MRI quality tools will be available in
% the SPM>Tools menu.
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%__________________________________________________________________________
% Cyclotron Research Centre - University of Liège
% Evelyne Balteau - June 2017
%==========================================================================

if ~isdeployed, addpath(fileparts(mfilename('fullpath'))); end

try
    % MRIquality is available
    mriq = tbx_cfg_mriq;
catch %#ok<CTCH>
    % No MRIquality toolbox found
    mriq         = cfg_exbranch;
    mriq.tag     = 'MRIq';
    mriq.name    = 'MRI Quality - not available!';
    mriq.help    = {
        ['The MRI quality toolbox does not seem to be available on this computer. ',...
        'In order to use the toolbox in SPM, the directory containing the ' ...
        'toolbox implementation does not need to be in the SPM/toolbox ' ...
        'directory, but it should be added to the Matlab path.']
        'Contact Evelyne Balteau if you need help with the toolbox (e.balteau@ulg.ac.be).'
        }';
    mriq.val  = {};
end

end

