function mriq_defaults
% Sets the defaults which are used by the MRI quality toolbox.
%
% FORMAT mriq_defaults
%_______________________________________________________________________
%
% This file can be customised to any the site/person own setup.
% Individual users can make copies which can be stored on their own
% matlab path. Make sure then that your 'mriq_defaults' is the first one
% found in the path. See matlab documentation for details on setting path.
%
% Care must be taken when modifying this file!
%
% The structure and content of this file are largely inspired by the
% equivalent file in SPM.
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Written by C. Phillips, 2013.
% Cyclotron Research Centre, University of Liege, Belgium

%%
global mriq_def

%========================= EPI STABILITY QA ===============================
% Directory where original data are intitially stored (input images)
mriq_def.epi_qa.paths.input = {'D:\home\logistic\mri\qa\data\input'}; 
% Directory where intermediate files are created and stored during
% processing. These files are deleted after processing if archiving is
% enabled.
mriq_def.epi_qa.paths.temp = {'D:\home\logistic\mri\qa\data\temp'}; 
% Directory for stability results (tSNR & SNR maps, drift, RDC, etc...)
mriq_def.epi_qa.paths.stab = {'D:\home\logistic\mri\qa\results\stability'}; 
% Directory for archiving original images 
mriq_def.epi_qa.paths.arch = {'D:\home\logistic\mri\qa\data\archives'}; 


% call local defaults to check for available protocols
if exist('mriq_defaults_local.m','file')
    mriq_defaults_local;
end

end