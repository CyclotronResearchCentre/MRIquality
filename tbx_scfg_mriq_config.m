function mriq_config = tbx_scfg_mriq_config
% (Sub)configuration file for the MRIquality toolbox
% -> Dealing with local defaults ("Configure Toolbox" branch)
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%==========================================================================
% Written by Evelyne Balteau, Cyclotron Research Centre, June 2019
%==========================================================================

%==========================================================================
% mriq_setdef Select MRIquality default parameter file
%==========================================================================
mriq_setdef         = cfg_files;
mriq_setdef.tag     = 'mriq_setdef';
mriq_setdef.name    = 'Defaults parameters';
mriq_setdef.help    = {['You can either stick with standard default parameters ' ...
    'from [mriq_defaults.m] or select your own customised default file.' ...
    'Some parameters may also be manually defined within this Batch GUI.']};
mriq_setdef.filter  = 'm';
mriq_setdef.dir     = fullfile(fileparts(mfilename('fullpath')),'config');
mriq_setdef.ufilter = '^mriq_.*\.m$';
mriq_setdef.num     = [1 1];
mriq_setdef.def     = @(val)mriq_get_defaults('def_file', val{:});

%==========================================================================
% path_arch Directory for archiving original data 
%==========================================================================
path_arch         = cfg_files;
path_arch.tag     = 'path_arch';
path_arch.name    = 'Archiving directory';
path_arch.val{1}  = {''};
path_arch.help    = {['Directory for archiving original data. ', ...
    '[When applicable and] when archiving is enabled, original input data ' ...
    'are comressed and archived in this directory.']};
path_arch.filter  = 'dir';
path_arch.ufilter = '.*';
path_arch.num     = [0 1];

%==========================================================================
% path_output Directory for output results
%==========================================================================
path_output         = cfg_files;
path_output.tag     = 'path_output';
path_output.name    = 'Output directory';
path_output.val{1}  = {''};
path_output.help    = {'Directory for output results.'};
path_output.filter  = 'dir';
path_output.ufilter = '.*';
path_output.num     = [0 1];

%==========================================================================
% path_tmp Directory where intermediate files are created and stored
% during processing. These files are deleted after processing if archiving
% is enabled.
% Directory for output results
% Directory for archiving original data 
%==========================================================================
path_tmp         = cfg_files;
path_tmp.tag     = 'path_tmp';
path_tmp.name    = 'Temporary directory';
path_tmp.val{1}  = {''};
path_tmp.help    = {['Directory where intermediate files are created and ' ...
    'stored during processing. [When applicable] these files are deleted ' ...
    'after processing if archiving and/or cleanup is enabled.']};
path_tmp.filter  = 'dir';
path_tmp.ufilter = '.*';
path_tmp.num     = [0 1];

%==========================================================================
% path_input Directory where original input data are stored
%==========================================================================
path_input         = cfg_files;
path_input.tag     = 'path_input';
path_input.name    = 'Input directory';
path_input.val{1}  = {''};
path_input.help    = {'Directory where original input data are stored.'};
path_input.filter  = 'dir';
path_input.ufilter = '.*';
path_input.num     = [0 1];

%==========================================================================
% Manually-defined defaults (within the GUI). These default parameters will
% overwrite all the other previously set defaults (from either standard or
% customised default files).
% WARNING: the job.mriq_manualdef fields must match the mriq_def ones!
%==========================================================================
mriq_manualdef           = cfg_branch;
mriq_manualdef.tag       = 'mriq_manualdef';
mriq_manualdef.name      = 'Customise defaults';
mriq_manualdef.help      = {'To further customise the default parameters.', ...
    ['The parameters defined here will overwrite the values defined in ' ...
    'either the standard or customised default file above. Parameters ' ...
    'below may be left blank - default values defined in the standard ' ...
    'or customised default file will be used.']};
mriq_manualdef.val       = {path_input path_tmp path_output path_arch};

%==========================================================================
% Configure the MRIquality toolbox - load local, user-defined defaults 
% file or manually selected defaults, and overwrite standard defaults 
%==========================================================================
mriq_config         = cfg_exbranch;
mriq_config.tag     = 'mriq_config';
mriq_config.name    = 'Configure';
mriq_config.val     = { mriq_setdef mriq_manualdef };
mriq_config.help    = {'Configure default parameters', ...
    ['Customised default parameters can be set here by selecting ' ...
    'a customised default file to replace ' ...
    '[mriq_defaults.m] and/or by defining manually '...
    'a series of parameters within this Configure module.']};
mriq_config.prog    = @mriq_run_config;

end
%----------------------------------------------------------------------

% =========================================================================
% (VOUT &) RUN SUBFUNCTION(S)
% =========================================================================
function out = mriq_run_config(job)
%==========================================================================
% PURPOSE
% Load standard defaults and overwrite them by customised values. 
%==========================================================================

% 1. reinitialise defaults to standard ones
mriq_defaults;

% 2. load and run selected default file 
% (may be identical to the standard one)
deffnam = job.mriq_setdef;
spm('Run',deffnam);

% 3. overwrite any customised defaults
% WARNING: the job.mriq_manualdef fields must match the mriq_def ones!
f = fieldnames(job.mriq_manualdef);
for cf = 1:length(f)
    if ~isempty(cell2mat(job.mriq_manualdef.(f{cf})))
        mriq_get_defaults(f{cf},job.mriq_manualdef.(f{cf}));
    end
end
out = mriq_get_defaults;

end
%_______________________________________________________________________
