%==========================================================================
% (VOUT &) RUN SUBFUNCTION(S)
%==========================================================================
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
