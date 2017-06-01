function mriq_run_epi_qa(job)
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
if isfield(job.output,'indir')
    PARAMS.paths.output = fullfile(PARAMS.paths.input,'stab');
    if ~exist(PARAMS.paths.output,'dir');mkdir(PARAMS.paths.output);end
else
    PARAMS.paths.output = job.output.ourdir{1};
end

% detect whether we have DICOM or NIFTI input files
if strcmp('nii',spm_str_manip(job.series.EPIseries{1},'e'))
    % NIFTI CASE: nothing to do
    PARAMS.intype = 'NIFTI';
    Ninim = char(job.series.EPIseries);
    Ninno = char(job.series.NOISEseries);
    %% check whether metadata are available!! error otherwise or...?
else
    % DICOM CASE: conversion into NIFTI first (wiht metadata)
    PARAMS.intype = 'DICOM';
    %% print DICOM conversion process...
    json = struct('extended',false,'separate',true); 
    hdr = spm_dicom_headers(char(job.series.EPIseries));
    Ninim = spm_dicom_convert(hdr,'all','flat','nii',PARAMS.paths.input,json); 
    Ninim = char(Ninim.files);
    hdr = spm_dicom_headers(char(job.series.NOISEseries));
    Ninno = spm_dicom_convert(hdr,'all','flat','nii',PARAMS.paths.input,json); 
    Ninno = char(Ninno.files);
end    

% comments and tags
PARAMS.comment = job.procpar.comment;

% retrieve values from headers
hdrim = get_metadata(Ninim(1,:));
if ~isempty(Ninno);hdrno = get_metadata(Ninno(1,:));end %#ok<NASGU>

% define and create temporary working directory
PARAMS.date = datestr(get_metadata_val(hdrim{1}, 'StudyDate'),'yyyymmdd');
PARAMS.series = get_metadata_val(hdrim{1}, 'SeriesNumber');
PARAMS.studyID = str2double(get_metadata_val(hdrim{1}, 'StudyID'));
PARAMS.resfnam = sprintf('%s_stud%0.4d_ser%0.4d%', PARAMS.date, PARAMS.studyID, PARAMS.series);
PARAMS.paths.process = fullfile(PARAMS.paths.input,PARAMS.resfnam);
[SUCCESS,MESSAGE,~] = mkdir(PARAMS.paths.process);
if ~SUCCESS; error(MESSAGE); end

% copy NIFTI files to processing directory (always keep input images
% untouched, i.e. original NIFTI files must be copied to processing
% directory, while NIFTI files converted from DICOM are moved to processing
% directory. NB: copy/move json files as well!
NinimTmp = [];
for cf = 1:size(Ninim,1)
    [~,NAME,EXT] = fileparts(Ninim(cf,:));
    NinimTmp = [NinimTmp; fullfile(PARAMS.paths.process, [NAME EXT])]; %#ok<AGROW>
    switch PARAMS.intype
        case 'DICOM'
            movefile(Ninim(cf,:),NinimTmp(cf,:));
            try movefile([spm_str_manip(Ninim(cf,:),'s') '.json'],[spm_str_manip(NinimTmp(cf,:),'s') '.json']); end %#ok<*TRYNC>
        case 'NIFTI'
            copyfile(Ninim(cf,:),NinimTmp(cf,:));
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
                movefile(Ninno(cf,:),NinnoTmp(cf,:));
                try movefile([spm_str_manip(Ninno(cf,:),'s') '.json'],[spm_str_manip(NinnoTmp(cf,:),'s') '.json']); end
            case 'NIFTI'
                copyfile(Ninno(cf,:),NinnoTmp(cf,:));
                try copyfile([spm_str_manip(Ninno(cf,:),'s') '.json'],[spm_str_manip(NinnoTmp(cf,:),'s') '.json']); end
        end
    end
end
Ninno = NinnoTmp;
Ninim = NinimTmp;
cd(PARAMS.paths.process);

% clear a few variables
clear NinimTmp NinnoTmp;

% gather a few more acquisition parameters for reccord and processing
PARAMS.nslices = get_metadata_val(get_metadata_val(hdrim{1}, 'sSliceArray'),'lSize');
PARAMS.TR = get_metadata_val(hdrim{1}, 'RepetitionTime');
tmp = get_metadata_val(hdrim{1},'AcquisitionMatrix');
PARAMS.PE_lin = tmp{1}(1);
PARAMS.RO_col = tmp{1}(4);
tmp = get_metadata_val(hdrim{1},'ReferenceAmplitude');
PARAMS.ref_ampl = tmp{1};
PARAMS.rf_freq = get_metadata_val(hdrim{1},'Frequency');
PARAMS.field_strength = get_metadata_val(hdrim{1},'FieldStrength');
tmp = get_metadata_val(hdrim{1},'SAR');
PARAMS.SAR = tmp{1};
PARAMS.MB = get_metadata_val(hdrim{1},'lMultiBandFactor');
if isempty(PARAMS.MB)
    PARAMS.MB = 1;
end
PARAMS.PAT = get_metadata_val(hdrim{1},'lAccelFactPE')*get_metadata_val(hdrim{1},'lAccelFact3D');
PARAMS.coils = get_metadata_val(hdrim{1},'ImaCoilString');
PARAMS.ncha = 0;
tmp = get_metadata_val(hdrim{1},'aRxCoilSelectData');
tmp = tmp{1};
for ccha = 1:length(tmp(1).asList);
    if isfield(tmp(1).asList(ccha),'lElementSelected')
        if (tmp(1).asList(ccha).lElementSelected == 1)
            PARAMS.ncha = PARAMS.ncha+1;
        end
    end
end
PARAMS.nvols = size(Ninim,1);

% Signal plane and noise plane are hardcoded
PARAMS.signalplane = job.procpar.sigplane;
if PARAMS.signalplane==0 % automatically defined as being the mid-volume slice
    PARAMS.signalplane = round(hdrim{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sSliceArray.lSize/2);
end
PARAMS.noiseplane = job.procpar.noiplane;
if PARAMS.signalplane==0 % automatically defined as being the last slice
    PARAMS.noiseplane = hdrim{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sSliceArray.lSize;
end
% NB: if PARAMS.signalplane<0, automated masking is applied to select noise
% voxels.
    
% write general information about the acquisition
fid = fopen(fullfile(PARAMS.paths.output, [PARAMS.resfnam '.txt']),'a');
fprintf(fid,'%s - %s\n',PARAMS.comment, PARAMS.date);
fprintf(fid,'\nACQUISITION PARAMETERS\n');
fprintf(fid,'    Coils: %s (%d channels)\n', PARAMS.coils, PARAMS.ncha);
fprintf(fid,'    Acceleration: MB%d + PAT%d\n', PARAMS.MB, PARAMS.PAT);
fprintf(fid,'    TR = %5.0f ms\n', PARAMS.TR);
fprintf(fid,'    Matrix = %dx%d\n', PARAMS.PE_lin, PARAMS.RO_col);
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
    flags1 = struct('quality',1,'fwhm',5,'sep',4,'interp',2,'wrap',[0 0 0],'rtm',0,'PW','','graphics',1,'lkp',2);
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

% PARAMS.signalplane = 30;
% PARAMS.noiseplane = 60;

% define central ROI for quantitative ROI analysis
N_max = 21; % maximal length of rectangular ROI edge
xg = ceil(PARAMS.RO_col/2)+1;
yg = ceil(PARAMS.PE_lin/2)+1;
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
save(fullfile(PARAMS.paths.output, [PARAMS.resfnam '.mat']),'RES','PARAMS');

% % % % % % % tidy up a bit...
% % % % % % % - delete all nifti files
% % % % % % % - tar.gzip the current data directory with dicom data only
% % % % % % % - delete the directory
% % % % % % delete('*.nii');
% % % % % % cd ..;
% % % % % % tar([PARAMS.paths.process(1:end-1) '.tar.gz'],PARAMS.paths.process);
% % % % % % rmdir(PARAMS.paths.process,'s');
% % % % % % 
% % % % % % cd(PARAMS.paths.matlab);
% % % % % % rmpath(PARAMS.paths.matlab);
% % % % % % clear all;