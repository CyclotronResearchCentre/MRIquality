function qa_stab_mb_epi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: qa_stab_mb_epi
% Written by Evelyne Balteau 
% Cyclotron Research Centre - March 2016
% 
% QA tool based on and modified from the FBIRN QA and others...
% Friedman and Glover, JMRI 23:827-839 (2006)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% loads default paths
% DIRS.dicom = ['W:' filesep]; % directory where to find the input dicom files
% DIRS.dicom = 'D:\home\logistic\mri\qa';
% DIRS.dicom = 'D:\home\data\20170124_terra_erlangen';
DIRS.dicom = 'D:\home\data\20170419_phantom_mb';
DIRS.matlab = fileparts(mfilename('fullpath')); % directory containing the present script
% DIRS.output = fullfile(fileparts(DIRS.matlab),'qa_stability', filesep);
% DIRS.output = fullfile('D:\home','logistic','mri','qa','qa_stability', filesep);
DIRS.output = 'D:\home\data\20170419_phantom_mb\stab';
% DIRS.data = fullfile(fileparts(DIRS.matlab),'qa_data', filesep);
% DIRS.data = fullfile('D:\home','logistic','mri','qa','qa_data', filesep);
DIRS.data = DIRS.output;
PARAMS.outdir = DIRS.output;
addpath(DIRS.matlab);

% select nii files
Ninim = spm_select(Inf, '^*\.nii$','select QA images (time series discarding first few images)',{},DIRS.dicom);
Ninno = spm_select(Inf, '^*\.nii$','select QA NOISE images (0, 1 or more)',{},DIRS.dicom);

% comments and tags
PARAMS.comment = input('Enter comment if required: ','s');
PARAMS.tag = input('Enter tag if several runs today: ','s');
% PARAMS.comment = 'MB2/PAT1 data acquired with 20-channel head coil (HE1-4)';

% retrieve values from headers
hdrim = get_metadata(Ninim(1,:));
if ~isempty(Ninno);hdrno = get_metadata(Ninno(1,:));end

% define and create temporary working directory
PARAMS.date = datestr(get_metadata_val(hdrim{1}, 'StudyDate'),'yyyymmdd');
PARAMS.series = get_metadata_val(hdrim{1}, 'SeriesNumber');
PARAMS.studyID = str2double(get_metadata_val(hdrim{1}, 'StudyID'));
PARAMS.resfnam = sprintf('%s_stud%0.4d_ser%0.4d%', PARAMS.date, PARAMS.studyID, PARAMS.series);
DIRS.current = fullfile(DIRS.data,PARAMS.resfnam, filesep);
[SUCCESS,MESSAGE,~] = mkdir(DIRS.current);
if ~SUCCESS; error(MESSAGE); end

% copy DICOM files to temp directory
NinimTmp = [];
for cf = 1:size(Ninim,1)
    [~,NAME,EXT] = fileparts(Ninim(cf,:));
    NinimTmp = [NinimTmp; fullfile(DIRS.current, [NAME EXT])];
    copyfile(Ninim(cf,:),NinimTmp(cf,:));
end
NinnoTmp = [];
if ~isempty(Ninno)
    for cf = 1:size(Ninno,1)
        [~,NAME,EXT] = fileparts(Ninno(cf,:));
        NinnoTmp = [NinnoTmp; fullfile(DIRS.current, [NAME EXT])];
        copyfile(Ninno(cf,:),NinnoTmp(cf,:));
    end
end
Ninno = NinnoTmp;
Ninim = NinimTmp;
cd(DIRS.current);

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
PARAMS.signalplane = round(hdrim{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sSliceArray.lSize/2);
PARAMS.noiseplane = hdrim{1}.acqpar.CSASeriesHeaderInfo.MrPhoenixProtocol.sSliceArray.lSize;

% write general information about the acquisition
fid = fopen(fullfile(DIRS.output, [PARAMS.resfnam '.txt']),'a');
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
    copyfile(SOURCE,fullfile(DIRS.output, RES.SPAT.meanim),'f')
    disp(RES.SPAT);
    % RETRIEVE FILE NAMES
    clear rNinim
    for i=1:size(Ninim,1)
        [p,n,e] = fileparts(Ninim(i,:));
        rNinim(i,:) = [p, '/r' n,e];
    end
    fid = fopen(fullfile(DIRS.output, [PARAMS.resfnam '.txt']),'a');
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
save(fullfile(DIRS.output, [PARAMS.resfnam '.mat']),'RES','PARAMS');

% % % % % % % tidy up a bit...
% % % % % % % - delete all nifti files
% % % % % % % - tar.gzip the current data directory with dicom data only
% % % % % % % - delete the directory
% % % % % % delete('*.nii');
% % % % % % cd ..;
% % % % % % tar([DIRS.current(1:end-1) '.tar.gz'],DIRS.current);
% % % % % % rmdir(DIRS.current,'s');
% % % % % % 
% % % % % % cd(DIRS.matlab);
% % % % % % rmpath(DIRS.matlab);
% % % % % % clear all;