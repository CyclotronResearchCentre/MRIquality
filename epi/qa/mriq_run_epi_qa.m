function out = mriq_run_epi_qa(job)
%==========================================================================
% USAGE: mriq_run_epi_qa(job)
% QA tools based on and modified from the FBIRN QA and others...
% Friedman and Glover, JMRI 23:827-839 (2006)
% Implemented as part of the MRI quality toolbox (partim EPI)
%==========================================================================
% Written by Evelyne Balteau 
% Cyclotron Research Centre - March 2016
%==========================================================================

% Set paths
PARAMS.paths.input = fileparts(job.series.EPIseries{1});
PARAMS.paths.matlab = fileparts(mfilename('fullpath')); % directory containing the present script
PARAMS.paths.output = cell2mat(mriq_get_defaults('path_output'));
if ~exist(PARAMS.paths.output,'dir')
    try
        mkdir(PARAMS.paths.output);
        fprintf(1,['\nWARNING: %s' ...
            '\nThe specified output directory does not exist and has been automatically created.' ...
            '\n'],PARAMS.paths.output);
    catch
        error(['\nERROR: %s' ...
            '\nThe specified output directory does not exist and could not be created.' ...
            '\nPlease use the Configure module to define the output directory.' ...
            '\n'],PARAMS.paths.output);
    end
end

PARAMS.archive = job.archive;
if PARAMS.archive
    PARAMS.paths.archive = cell2mat(mriq_get_defaults('path_arch'));
    fprintf(1,['\nINFO: Data will be cleaned up at the end of the processing, including ' ...
        '\noriginal images compression and archiving to the following directory:' ...
        '\n\t%s\n'],PARAMS.paths.archive);
    if ~exist(PARAMS.paths.archive,'dir')
        try
            mkdir(PARAMS.paths.archive);
            fprintf(1,['\nWARNING: %s' ...
                '\nThe specified archive directory does not exist and has been automatically created.' ...
                '\n'],PARAMS.paths.archive);
        catch
            error(['\nERROR: %s' ...
                '\nThe specified archive directory does not exist and could not be created.' ...
                '\nPlease use the Configure module to define the archive directory.' ...
                '\n'],PARAMS.paths.archive);
        end
    end
else
    fprintf(1,'\nData won''t be cleaned up, compressed nor archived.\n');
end

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
            '\nThe specified temporary directory intermediate files produced ' ...
            '\nduring data processing does not exist and could not be created.' ...
            '\nPlease use the Configure module to define the temporary directory.' ...
            '\n'],PARAMS.paths.process);
    end
end


% The SPM Graphics window must be available (or maybe not, but this is
% doing no harm... will be used to display the realign results)...
fg = spm_figure('FindWin','Graphics');
if isempty(fg)
    spm_figure('Create','Graphics','Graphics','on')
end


% detect whether we have DICOM or NIFTI input files
Ninim = char(job.series.EPIseries);
Ninno = char(job.series.NOISEseries);
if strcmp('nii',spm_str_manip(job.series.EPIseries{1},'e'))
    % NIFTI CASE: nothing to do
    PARAMS.intype = 'NIFTI';
    % check whether metadata are available!! 
    metadatatmp = get_metadata(Ninim(1,:));
    if isempty(metadatatmp{1})
        error(['No metadata associated with the input NIFTI files. Cannot proceed. ' ...
            'Please use the hMRI toolbox to convert your DICOM images into NIFTI format.']);
    end
else
    % DICOM CASE: conversion into NIFTI first (with metadata)
    PARAMS.intype = 'DICOM';
    json = struct('extended',false,'separate',true); 
    hdr = spm_dicom_headers(Ninim);
    Ninim = spm_dicom_convert(hdr,'all','flat','nii',PARAMS.paths.input,json); 
    Ninim = char(Ninim.files);
    if ~isempty(Ninno)
        hdr = spm_dicom_headers(char(job.series.NOISEseries));
        Ninno = spm_dicom_convert(hdr,'all','flat','nii',PARAMS.paths.input,json);
        Ninno = char(Ninno.files);
    end
end    

% comments and tags
PARAMS.comment = job.procpar.comment;

% retrieve values from headers
hdrim = get_metadata(Ninim(1,:));
% if ~isempty(Ninno);hdrno = get_metadata(Ninno(1,:));end

% define and create temporary working directory
PARAMS.date = datestr(get_metadata_val(hdrim{1}, 'StudyDate'),'yyyymmdd');
PARAMS.series = get_metadata_val(hdrim{1}, 'SeriesNumber');
PARAMS.studyID = str2double(get_metadata_val(hdrim{1}, 'StudyID'));
PARAMS.scanner = deblank(get_metadata_val(hdrim{1},'ManufacturerModelName'));
PARAMS.resfnam = sprintf('%s_%s_stud%0.4d_ser%0.4d%', PARAMS.scanner, PARAMS.date, PARAMS.studyID, PARAMS.series);

% Dummy volumes to be discarded (T1 saturation effect)
PARAMS.dummies = job.procpar.dummies;
if PARAMS.dummies>size(Ninim,1)
    fprintf(1,['\nWARNING: The number of dummy volumes exceeds the number ' ...
            '\nof volumes in the time series. No dummy volume will be discarded.' ...
            '\n']);
    PARAMS.dummies = 0;
end
Ninim = Ninim(1+PARAMS.dummies:end,:);

% copy NIFTI files to processing directory (always keep
% input images untouched, i.e. original NIFTI files must be copied to
% processing directory, while NIFTI files converted from DICOM are moved to
% processing directory. NB: copy/move json files as well!
NinimTmp = [];
for cf = 1:size(Ninim,1)
    [~,NAME,EXT] = fileparts(Ninim(cf,:));
    NinimTmp = [NinimTmp; fullfile(PARAMS.paths.process, [NAME EXT])]; %#ok<AGROW>
    switch PARAMS.intype
        case 'DICOM'
            try movefile(Ninim(cf,:),NinimTmp(cf,:)); end
            try movefile([spm_str_manip(Ninim(cf,:),'s') '.json'],[spm_str_manip(NinimTmp(cf,:),'s') '.json']); end %#ok<*TRYNC>
        case 'NIFTI'
            try copyfile(Ninim(cf,:),NinimTmp(cf,:)); end
            try copyfile([spm_str_manip(Ninim(cf,:),'s') '.json'],[spm_str_manip(NinimTmp(cf,:),'s') '.json']); end
    end
end
NinnoTmp = [];
if ~isempty(Ninno)
    for cf = 1:size(Ninno,1)
        [~,NAME,EXT] = fileparts(Ninno(cf,:));
        NinnoTmp = [NinnoTmp; fullfile(PARAMS.paths.process, [NAME EXT])]; %#ok<AGROW>
        switch PARAMS.intype
            case 'DICOM'
                try movefile(Ninno(cf,:),NinnoTmp(cf,:)); end
                try movefile([spm_str_manip(Ninno(cf,:),'s') '.json'],[spm_str_manip(NinnoTmp(cf,:),'s') '.json']); end
            case 'NIFTI'
                try copyfile(Ninno(cf,:),NinnoTmp(cf,:)); end 
                try copyfile([spm_str_manip(Ninno(cf,:),'s') '.json'],[spm_str_manip(NinnoTmp(cf,:),'s') '.json']); end
        end
    end
end
Ninno = NinnoTmp;
Ninim = NinimTmp;

% clear a few variables
clear NinimTmp NinnoTmp;

% gather a few more acquisition parameters for reccord and processing
PARAMS.TR = get_metadata_val(hdrim{1}, 'RepetitionTime');
tmp = get_metadata_val(hdrim{1},'ReferenceAmplitude');
PARAMS.ref_ampl = tmp{1};
PARAMS.rf_freq = get_metadata_val(hdrim{1},'Frequency');
PARAMS.field_strength = get_metadata_val(hdrim{1},'FieldStrength');
tmp = get_metadata_val(hdrim{1},'SAR');
PARAMS.SAR = tmp{1};
V = spm_vol(Ninim(1,:));
PARAMS.dim = V.dim;

PARAMS.MB = get_metadata_val(hdrim{1},'MultiBandFactor');
% PARAMS.MB = input(sprintf('\nINPUT REQUIRED - Enter MB factor (%d): ', PARAMS.MB));
    
PARAMS.PAT = 1;
try 
    PARAMS.PAT = get_metadata_val(hdrim{1},'AccelFactorPE')*get_metadata_val(hdrim{1},'AccelFactor3D');
end
% PARAMS.PAT = input(sprintf('\nINPUT REQUIRED - Enter PAT factor (PEx3D = %d): ', PARAMS.PAT));

PARAMS.coils = get_metadata_val(hdrim{1},'ImaCoilString');
if isempty(PARAMS.coils)
%     PARAMS.coils = input('\nINPUT REQUIRED - Enter name of Rx coil: ','s');
    PARAMS.coils = 'Unknown';
end

PARAMS.ncha = 0;
try
    tmp = get_metadata_val(hdrim{1},'aRxCoilSelectData');
    tmp = tmp{1};
    for ccha = 1:length(tmp(1).asList);
        if isfield(tmp(1).asList(ccha),'lElementSelected')
            if (tmp(1).asList(ccha).lElementSelected == 1)
                PARAMS.ncha = PARAMS.ncha+1;
            end
        end
    end
end
% PARAMS.ncha = input(sprintf('\nINPUT REQUIRED - Enter number of coil elements (%d): ', PARAMS.ncha));

PARAMS.nvols = size(Ninim,1);
PARAMS.scfac = 1;
% if ~isempty(Ninno) % the ratio between noise and image scale factors
%     tmpscim = get_metadata_val(hdrim{1},'sCoilSelectUI');
%     tmpscno = get_metadata_val(hdrno{1},'sCoilSelectUI');
%     PARAMS.scfac = tmpscno.dOverallImageScaleFactor/tmpscim.dOverallImageScaleFactor; 
% end
% PARAMS.scfac = input(sprintf(['Identical scaling factor is assumed between EPI time series and noise images.' ...
%     '\nIf so, enter 1. otherwise, enter the ratio noise scfac / EPI scfac : ']));

% Signal plane and noise plane
PARAMS.signalplane = job.procpar.sigplane;
if PARAMS.signalplane == 0 % automatically defined as being the mid-volume slice
    PARAMS.signalplane = round(hdrim{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sSliceArray.lSize/2);
end
PARAMS.noiseplane = job.procpar.noiplane;
% NB: if PARAMS.noiseplane==0, it is assumed that no noise plane was
% acquired. Automated masking is applied to select noise voxels.
    
% write general information about the acquisition
fid = fopen(fullfile(PARAMS.paths.output, [PARAMS.resfnam '.txt']),'a');
fprintf(fid,'%s - %s\n',PARAMS.comment, PARAMS.date);
fprintf(fid,'\nACQUISITION PARAMETERS\n');
fprintf(fid,'    Scanner: %s\n', PARAMS.scanner);
fprintf(fid,'    Coils: %s (%d channels)\n', PARAMS.coils, PARAMS.ncha);
fprintf(fid,'    Acceleration: MB%d + PAT%d\n', PARAMS.MB, PARAMS.PAT);
fprintf(fid,'    TR = %5.0f ms\n', PARAMS.TR);
fprintf(fid,'    B0 = %5.4f T\n', PARAMS.field_strength);
fprintf(fid,'    Frequency = %9.6f MHz\n', PARAMS.rf_freq/1000000);
fprintf(fid,'    RefAmpl = %6.3f V\n', PARAMS.ref_ampl);
fprintf(fid,'    Signal plane = %d\n', PARAMS.signalplane);
fprintf(fid,'    Noise plane = %d\n', PARAMS.noiseplane);
fclose(fid);

% realign and reslice in the PE(y) direction to account for B0 drifts if
% number of image files is >1 
RES.SPAT = [];
rNinim = Ninim;
if size(Ninim,1)>1
    % REALIGN
    flags1 = struct('quality',1,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,'PW','','lkp',2);
    [rNinim, Params] = spm_realign_fbirn(Ninim,flags1);
    
    % RESLICE
    flags2 = struct('interp',4,'mask',1,'mean',1,'which',2,'wrap',[0 0 0]');
    spm_reslice(rNinim,flags2);
    RES.SPAT.y_motion = Params(:,2); % spatial drift in mm
    RES.SPAT.max_drift = max(RES.SPAT.y_motion)-min(RES.SPAT.y_motion); % maximal drift in mm
    [p,n,e] = fileparts(rNinim(1).fname);
    SOURCE = fullfile(p, ['mean' n,e]);
    RES.SPAT.meanim = [PARAMS.resfnam '_MEAN.nii'];
    copyfile(SOURCE,fullfile(PARAMS.paths.output, RES.SPAT.meanim),'f')
    disp(RES.SPAT);
    
    % RETRIEVE FILE NAMES
    clear rNinim
    for i=1:size(Ninim,1)
        [p,n,e] = fileparts(Ninim(i,:));
        rNinim(i,:) = [p, '/r' n,e];
    end
    fid = fopen(fullfile(PARAMS.paths.output, [PARAMS.resfnam '.txt']),'a');
    fprintf(fid,'\nSPATIAL PROCESSING\n');
    fprintf(fid,'    Maximum spatial drift in mm: %5.4f mm\n', RES.SPAT.max_drift);
    fclose(fid);
end

% define central ROI for quantitative ROI analysis
N_max = job.procpar.roisize; % maximal length of rectangular ROI edge
xg = ceil(PARAMS.dim(2)/2)+1;
yg = ceil(PARAMS.dim(1)/2)+1;
PARAMS.x_roi = (-floor(N_max/2):floor(N_max/2)) + xg;
PARAMS.y_roi = (-floor(N_max/2):floor(N_max/2)) + yg;

% % Select center of ROI (x,y AND z):
% figure('color',[1 1 1],'units','centimeters','position',[1 1 20 20]);colormap(gray);
% montage(reshape(Y1,[V1.dim(1),V1.dim(2),1,V1.dim(3)]));
% caxis([min(Y1(:)), max(Y1(:))]);
% title('Select center of the ROI','fontname',def_fontname,'fontsize',title_fontsize);
% [xg yg] = ginput(1);
% QA.signalplane = ceil(xg/V1.dim(1)) + floor(yg/V1.dim(2))*7;
% xg = xg - floor(xg/V1.dim(1))*V1.dim(1);
% yg = yg - floor(yg/V1.dim(2))*V1.dim(2);
% close gcf;

% % Select center of ROI (x,y from the default 'signalplane'):
% figure('color',[1 1 1],'units','centimeters','position',[1 1 20 20]);
% colormap(gray);
% subplot(1,1,1);
% imshow(Y1(:,:,QA.signalplane),[]);
% set(gca,'units','centimeters','position',[1 1 18 18]);
% set(gcf,'position',[1 1 20 20]);
% caxis([min(Y1(:)), max(Y1(:))]);
% title('Select center of the ROI','fontname',def_fontname,'fontsize',title_fontsize*1.8);
% [xg yg] = ginput(1);
% close gcf;

% estimate noise distribution, calculate SNR map and central ROI average SNR
RES.SNR = qa_snr_mb(rNinim, Ninno, PARAMS);
disp(RES.SNR);

% estimate other FBIRN parameters
RES.FBIRN = qa_fbirn_mb(rNinim, PARAMS);
disp(RES.FBIRN);

% display summary figure
qa_stab_mb_epi_display_results(RES, PARAMS);

% save all results and parameters
spm_jsonwrite(fullfile(PARAMS.paths.output, [PARAMS.resfnam '.json']),struct('PARAMS',PARAMS,'RES',RES),struct('indent','\t'));
% save the job:
clear matlabbatch;
matlabbatch{1}.spm.tools.mriq.epi.epiqa = job;
save(fullfile(PARAMS.paths.output, [PARAMS.resfnam '_batch.mat']), 'matlabbatch');

% define outputs - TO BE IMPLEMENTED
out = [];

if PARAMS.archive
    % tidy up a bit...
    % - tar.gzip input images only and archive in defined directory
    files2tar = job.series.EPIseries;
    if ~isempty(char(job.series.NOISEseries))
        files2tar = cat(1,files2tar,job.series.NOISEseries);
    end
    % if NIFTI files as inputs, must add JSON metadata
    if strcmp(PARAMS.intype,'NIFTI')
        files2tarJSON = cell(length(files2tar),1);
        for cf = 1:size(files2tar,1)
            files2tarJSON{cf} = [spm_str_manip(files2tar{cf},'s') '.json']; 
        end
        files2tar = cat(1, files2tar, files2tarJSON);
    end
    tar(fullfile(PARAMS.paths.archive, [PARAMS.resfnam '.tar.gz']), files2tar);

    % - delete all other intermediate files from temp directory
    files2del = spm_select('FPList',PARAMS.paths.process);
    for cf = 1: size(files2del, 1)
        delete(strtrim(files2del(cf,:)));
    end
end

end