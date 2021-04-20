function mriq_defaults_prisma_autoqa
% Sets the defaults which are used by the MRI quality toolbox for automated
% QA data processing.
% FORMAT mriq_defaults_prisma_autoqa
%_______________________________________________________________________
%
% This file has been adapted to process QA data from the Prisma scanner at
% CRC, assuming the data to be processed are initially in
% \\fmriserver2\mriQA which has been mounted as "W:" on the computer where
% the automated QA is run...
%
% The structure and content of this file are largely inspired by the
% equivalent file in SPM.
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Written by ebalteau, 2020.
% Cyclotron Research Centre, University of Liege, Belgium

%%
global mriq_def

% This very default file
mriq_def.def_file = {[mfilename('fullpath') '.m']};

% Directory where original data are initially stored (input images)
mriq_def.path_input = {fullfile('W:')}; 

% Directory where intermediate files are created and stored during
% processing. These files are deleted after processing if archiving is
% enabled.
mriq_def.path_tmp = {fullfile('C:','MRIquality','tmp')}; 

% Directory for results (stability or else)
mriq_def.path_output = {fullfile('W:','MRIquality','results')}; 

% Directory for archiving original images 
mriq_def.path_arch = {fullfile('W:','MRIquality','archives')}; 

end