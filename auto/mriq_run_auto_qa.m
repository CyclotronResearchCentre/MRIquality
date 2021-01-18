function out = mriq_run_auto_qa(dicomdir)

% Main script for running automated QA. The argument is the directory
% containing the DICOM (*.IMA) images as transferred after QA acquisition.
% If not provided, it'll retrieve the default path from the default
% parameter file (see below).

% the main directory where data arrive
if (nargin==0)
    dicomdir = 'W:\inpQA'; % for current development ant testing... eb 10/12/2020
    % dicomdir = cell2mat(mriq_get_defaults('path_input'));
end

% default for DICOM to NIfTI conversion
json = struct('extended',false,'separate',true);
        
% specify temporary directory for processing the data
% WARNING: the default tmp directory is defined as
% fullfile(pwd,'MRIquality','tmp'). Might need to re-run mriq_defaults
% manually if current working directory has been wrongly set and needs to
% be changed. 
PARAMS.paths.process = cell2mat(mriq_get_defaults('path_tmp'));
if ~exist(PARAMS.paths.process,'dir')
    try
        mkdir(PARAMS.paths.process);
        fprintf(1,['\nWARNING: %s' ...
            '\nThe specified temporary directory for intermediate files produced ' ...
            '\nduring data processing does not exist and has been automatically created.' ...
            '\n'],PARAMS.paths.process);
    catch
        error(['\nERROR: %s' ...
            '\nThe specified temporary directory for intermediate files produced ' ...
            '\nduring data processing does not exist and could not be created.' ...
            '\nPlease use the Configure module to define the temporary directory.' ...
            '\n'],PARAMS.paths.process);
    end
end

% check any new files in dicomdir
dicomfilelist = spm_select('FPList',dicomdir,'^.*\.IMA$');

% if not empty...
if ~isempty(dicomfilelist)
    % 0) copy all files to tmp directory
    for cf = 1:size(dicomfilelist,1)
        copyfile(dicomfilelist(cf,:),PARAMS.paths.process); % must be changed to MOVE in final version
    end
    dicomfilelist = spm_select('FPList',PARAMS.paths.process,'^.*\.IMA$');

    % 1) read headers 
    hdr = spm_dicom_headers(dicomfilelist);
    
    % 2) since EPI and NOISE data to be processed together must have been
    % acquired during the same session (i.e. same study date and ID), we
    % move each file in the corresponding <acquisition date>_<StudyID>
    % subdirectory. 
    % Format from header: - StudyDate: 737512
    %                     - StudyID: '50 '
    for chdr = 1:length(hdr)
        cdate = datestr(hdr{chdr}.StudyDate,'yyyymmdd');
        cid = sprintf('%0.4d',str2double(hdr{chdr}.StudyID));
        cser = sprintf('%0.4d',hdr{chdr}.SeriesNumber);
        cdir = fullfile(PARAMS.paths.process, sprintf('%s_%s', cdate, cid), sprintf('ser%s',cser));
        if ~exist(cdir,'dir')
            mkdir(cdir);
        end
        movefile(hdr{chdr}.Filename, cdir);
    end
    
    % 4) list the created directories
    dirlist = spm_select('FPList', PARAMS.paths.process,'dir');
    
    % 5) Proceed one directory (i.e. one session) at a time
    for cdir = 1:size(dirlist,1)
        % 5.1) list the series of the current session
        serieslist = spm_select('FPList', dirlist(cdir,:),'dir');
        % 5.2) DICOM2NII convert each series and retrieve acq params
        for cser = 1:size(serieslist,1)
            dicomfilelist = spm_select('FPList',serieslist(cser,:),'^.*\.IMA$');
            % convert into NIfTI
            hdr = spm_dicom_headers(dicomfilelist);
            niifilelist{cser} = spm_dicom_convert(hdr,'all','flat','nii',dirlist(cdir,:),json); %#ok<*AGROW>
            % gather header information for the current series
            niifilelist{cser}.nnii = length(niifilelist{cser}.files);
            niifilelist{cser}.sernum = get_metadata_val(niifilelist{cser}.files{1},'SeriesNumber');
            niifilelist{cser}.filelist = [];
            niifilelist{cser}.refampl = get_metadata_val(niifilelist{cser}.files{1},'flReferenceAmplitude');
            if isempty(niifilelist{cser}.refampl)
                niifilelist{cser}.refampl = 0.0;
            end
            tmp = get_metadata_val(niifilelist{cser}.files{1},'ImageType'); % 'ORIGINAL\PRIMARY\M\MB\ND\MOSAIC '
            niifilelist{cser}.acqparams.ImageType = deblank(tmp{1});
            niifilelist{cser}.acqparams.ScanningSequence = deblank(get_metadata_val(niifilelist{cser}.files{1},'ScanningSequence')); % 'EP'
            niifilelist{cser}.acqparams.SequenceVariant = deblank(get_metadata_val(niifilelist{cser}.files{1},'SequenceVariant')); % 'SK\SS '
            niifilelist{cser}.acqparams.ScanOptions = deblank(get_metadata_val(niifilelist{cser}.files{1},'ScanOptions')); % 'FS'
            niifilelist{cser}.acqparams.MRAcquisitionType = deblank(get_metadata_val(niifilelist{cser}.files{1},'MRAcquisitionType')); % '2D'
            niifilelist{cser}.acqparams.SequenceName = deblank(get_metadata_val(niifilelist{cser}.files{1},'SequenceName')); % 'epfid2d1_72 '
            niifilelist{cser}.acqparams.SliceThickness = get_metadata_val(niifilelist{cser}.files{1},'SliceThickness'); % 3
            niifilelist{cser}.acqparams.RepetitionTime = get_metadata_val(niifilelist{cser}.files{1},'RepetitionTime'); % 1170
            niifilelist{cser}.acqparams.EchoTime = get_metadata_val(niifilelist{cser}.files{1},'EchoTime'); % 30
            niifilelist{cser}.acqparams.ImagingFrequency = 0.001*round(1000*get_metadata_val(niifilelist{cser}.files{1},'ImagingFrequency')); % 123.2542
            niifilelist{cser}.acqparams.SpacingBetweenSlices = get_metadata_val(niifilelist{cser}.files{1},'SpacingBetweenSlices'); % 3.7500
            niifilelist{cser}.acqparams.NumberOfPhaseEncodingSteps = get_metadata_val(niifilelist{cser}.files{1},'NumberOfPhaseEncodingSteps'); % 72
            tmp = get_metadata_val(niifilelist{cser}.files{1},'EchoTrainLength'); % 72
            niifilelist{cser}.acqparams.EchoTrainLength = tmp{1};
            niifilelist{cser}.acqparams.PercentSampling = get_metadata_val(niifilelist{cser}.files{1},'PercentSampling'); % 100
            niifilelist{cser}.acqparams.PercentPhaseFieldOfView = get_metadata_val(niifilelist{cser}.files{1},'PercentPhaseFieldOfView'); % 100
            niifilelist{cser}.acqparams.PixelBandwidth = get_metadata_val(niifilelist{cser}.files{1},'PixelBandwidth'); % 2570
            niifilelist{cser}.acqparams.DeviceSerialNumber = deblank(get_metadata_val(niifilelist{cser}.files{1},'DeviceSerialNumber')); % '66021 '
            tmp = get_metadata_val(niifilelist{cser}.files{1},'AcquisitionMatrix'); % [72 0 0 72]
            niifilelist{cser}.acqparams.AcquisitionMatrix = tmp{1};
            niifilelist{cser}.acqparams.InPlanePhaseEncodingDirection = deblank(get_metadata_val(niifilelist{cser}.files{1},'InPlanePhaseEncodingDirection')); % 'COL '
            niifilelist{cser}.acqparams.FlipAngle = get_metadata_val(niifilelist{cser}.files{1},'FlipAngle'); % 65
            niifilelist{cser}.acqparams.ImagePositionPatient = get_metadata_val(niifilelist{cser}.files{1},'ImagePositionPatient'); % [3x1 double]
            niifilelist{cser}.acqparams.ImageOrientationPatient = get_metadata_val(niifilelist{cser}.files{1},'ImageOrientationPatient'); % [6x1 double]
            niifilelist{cser}.acqparams.ImageComments = deblank(get_metadata_val(niifilelist{cser}.files{1},'ImageComments')); % 'Unaliased MB2/PE2 '
            tmp = get_metadata_val(niifilelist{cser}.files{1},'CoilString');
            niifilelist{cser}.acqparams.CoilString = deblank(tmp{1});
            niifilelist{cser}.acqparams.AccelFactorPE = get_metadata_val(niifilelist{cser}.files{1},'AccelFactorPE');
            niifilelist{cser}.acqparams.MultiBandFactor = get_metadata_val(niifilelist{cser}.files{1},'MultiBandFactor');
        end
        
        % 5.3 Search for EPIseries-NOISEseries pairs, excludes non-EPI
        % series from the analysis (THIS MIGHT CHANGE IN THE FUTURE FOR
        % OTHER TYPE OF QA PROCEDURES, E.G. STRUCTURAL IMAGES, B0 or B1
        % IMAGES, etc...)
        % For each EPIseries, search for the corresponding NOISEseries. If
        % none, use empty NOISEseries. To determine whether a NOISEseries
        % is to be processed together with a given EPIseries, we compare
        % all the acqparams. They must match. 
        for cEPIser = 1:size(serieslist,1)
            if (niifilelist{cEPIser}.refampl ~= 0.0) && strcmp(niifilelist{cEPIser}.acqparams.ScanningSequence, 'EP')
                EPIseries = niifilelist{cEPIser}.files;
                fieldlist = fieldnames(niifilelist{cEPIser}.acqparams);
                NOISEseries = [];
                for cNOISEser = 1:size(serieslist,1)
                    if isempty(NOISEseries) && (niifilelist{cNOISEser}.refampl == 0.0) && strcmp(niifilelist{cNOISEser}.acqparams.ScanningSequence, 'EP')
                        NOISEok = true;
                        cfield = 1;
                        while NOISEok && cfield<length(fieldlist)+1
                            NOISEok = all(niifilelist{cNOISEser}.acqparams.(fieldlist{cfield}) == niifilelist{cEPIser}.acqparams.(fieldlist{cfield}));
                            if ~NOISEok
                                fprintf(1,'\nEPIseries = #%d - NOISEseries = #%d', niifilelist{cEPIser}.sernum, niifilelist{cNOISEser}.sernum);
                                fprintf(1,'\n DO NOT MATCH because of parameter %s:',fieldlist{cfield});
                                disp([niifilelist{cEPIser}.acqparams.(fieldlist{cfield}) niifilelist{cNOISEser}.acqparams.(fieldlist{cfield})]);
                            end
                            cfield = cfield + 1;
                        end
                        if NOISEok
                            NOISEseries = niifilelist{cNOISEser}.files;
                        end
                    end
                end
                
                fprintf(1,'\nEPIseries (%s etc...)\n',EPIseries{1});
                disp(niifilelist{cEPIser}.acqparams);
                if ~isempty(NOISEseries)
                    fprintf(1,'\nNOISEseries (%s etc...)\n',NOISEseries{1});
                    disp(niifilelist{cNOISEser}.acqparams);
                else
                    fprintf(1,'\nEmpty NOISEseries...\n');
                end
                out = mriq_run_epi_qa_wrapper(EPIseries, NOISEseries);
            end
        end
    end
end

end

%=========================================================================%
% WRAPPER TO BE CALLED TO RUN THE QA OVER PRESELECTED EPI & NOISE DATA
%=========================================================================%
function out = mriq_run_epi_qa_wrapper(EPIseries, NOISEseries)

% NOTE: we assume MRIQ_CONFIG has been run beforehand or that standard
% defaults are all right. This includes paths to various directories (data,
% tmp, archives, results).

% matlabbatch{1}.spm.tools.mriq.mriq_config.mriq_setdef = {'D:\home\git\MRIquality\config\mriq_defaults.m'};
% matlabbatch{1}.spm.tools.mriq.mriq_config.mriq_manualdef.path_input = {fileparts(EPIseries{1})}; % actually not used
% matlabbatch{1}.spm.tools.mriq.mriq_config.mriq_manualdef.path_tmp = {''}; % keep defaults
% matlabbatch{1}.spm.tools.mriq.mriq_config.mriq_manualdef.path_output = {flags.LogFile.LogDir};
% matlabbatch{1}.spm.tools.mriq.mriq_config.mriq_manualdef.path_arch = {'W:\archQA'};

clear matlabbatch;
cb = 1;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.series.EPIseries = EPIseries;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.series.NOISEseries = NOISEseries;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.archive = 1;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.procpar.dummies = 0;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.procpar.sigplane = 0;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.procpar.noiplane = 0;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.procpar.roisize = 21;
matlabbatch{cb}.spm.tools.mriq.epi.epiqa.procpar.comment = 'Automated QA';

spm ('defaults', 'fmri');      
% spm_jobman('initcfg'); NOTE 20210118 - init jobman must have been run
% beforehand or it'll be run here anyway. The consequence of running it
% here is that the defaults are set to the spm/config directory :/... To be
% tested: run spm_jobman('initcfg') at the beginning of the
% mriq_run_auto_qa script if "~deployed"... or in mriq_run_config... To be
% continued... The problem does not occur if mriq_run_auto_qa run from
% Batch GUI (not currently implemented) since in that case SPM is
% deployed...
out = spm_jobman('run', matlabbatch);

end