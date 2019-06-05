function mriq_defaults
% Sets the defaults which are used by the MRI quality toolbox.
%
% FORMAT mriq_defaults
%_______________________________________________________________________
%
% This file is given as a template to set default parameters according to
% your own needs. Default parameters include parameters such as output
% directory location for QA results, temporary directory for intermediate
% steps that can be deleted at the end of the procedure, etc.  
% Individual users can make copies which can be stored on their own
% matlab path, and loaded using the Configure toolbox module.
%
% It is recommended not to modify this file, but to save a copy with a
% meaningful name and then modify it according to your own needs. 
%
% The structure and content of this file are largely inspired by the
% equivalent file in SPM.
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Written by ebalteau, 2019.
% Cyclotron Research Centre, University of Liege, Belgium

%%
global mriq_def

% This very default file
mriq_def.def_file = {[mfilename('fullpath') '.m']};

% Directory where original data are initially stored (input images)
mriq_def.path_input = {fullfile('C:','MRIquality','data')}; 

% Directory where intermediate files are created and stored during
% processing. These files are deleted after processing if archiving is
% enabled.
mriq_def.path_tmp = {fullfile('C:','MRIquality','tmp')}; 

% Directory for results (stability or else)
mriq_def.path_output = {fullfile('C:','MRIquality','results')}; 

% Directory for archiving original images 
mriq_def.path_arch = {fullfile('C:','MRIquality','archives')}; 

end