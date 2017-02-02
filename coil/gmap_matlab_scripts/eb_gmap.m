function eb_gmap(varargin)
%---------------------------------------------------------------------
% function to generate sensitivity maps, noise correlation matrices
% and g-maps from complex raw data.
% !! quite intensive computing, using a lot of memory !!
% compilation making use of several routines by Chloe Hutton, 
% Last modif: July 2008 - Evelyne Balteau
% Last modif: October 2013 - Evelyne Balteau - adapted for Prisma data
%---------------------------------------------------------------------

% addpath(genpath('/home/ebalteau/Documents/implem/spm2'));

% filenameraw = spm_select(1,'^*\.dat$','Select image raw data','',0)
% filenamenoise = spm_select(1,'^*\.dat$','Select noise raw data','',0)

filenameraw = 'D:\home\data\20170124_terra_erlangen\raw\meas_MID00454_FID00710_gre.dat';
filenamenoise = 'D:\home\data\20170124_terra_erlangen\raw\meas_MID00456_FID00712_gre.dat';

% filenameraw = 'D:\home\data\20170124_terra_erlangen\raw\meas_MID00385_FID00664_gre.dat';
% filenamenoise = 'D:\home\data\20170124_terra_erlangen\raw\meas_MID00387_FID00666_gre.dat';

% filenameraw = '/media/Elements/lausanne/test_liege_20131014_raw/meas_MID00081_FID01480_gre_sensitivity.dat';
% filenamenoise = '/media/Elements/lausanne/test_liege_20131014_raw/meas_MID00082_FID01481_gre_sensitivity.dat';

% filenameraw =   '/home/vivelyne/Documents/brols/lausanne/toworkonthetrain/meas_MID00081_FID01480_gre_sensitivity.dat';
% filenamenoise = '/home/vivelyne/Documents/brols/lausanne/toworkonthetrain/meas_MID00082_FID01481_gre_sensitivity.dat';

% filenameraw =   'E:\data_archive\new3Ttests\prisma\test_liege_20131014_raw\meas_MID00081_FID01480_gre_sensitivity.dat';
% filenamenoise = 'E:\data_archive\new3Ttests\prisma\test_liege_20131014_raw\meas_MID00082_FID01481_gre_sensitivity.dat';

%---------------------------------------------------------------------
% read the headers of the raw data (*.asc file):
%---------------------------------------------------------------------
[maindir,fname,e] = fileparts(filenameraw);
header = eb_read_protocol([maindir filesep fname '.dat']);

% some hacking to get it work with the Prisma data
% MID00089 & MID00091 - BC (2ch receivers)
% MDH_END = [120 200];
header.asc{4} = header.asc{2};
header.Meas{2} = header.Dicom{1};

% MID00080 & MID00081 - 20ch head-neck (HC1-4 + NE1,2
% MDH_END = [112 104];

% MID00054 & MID00055 - body18 + SP8 
% MDH_END = [128 120];

% MID00036 & MID00037 - 64ch head-neck (HC1-7 + NE1,2
% MDH_END = [168 152];

% MID00075 & MID00076 - 64ch head-neck (HC1-7 + NE1,2
% MDH_END = [200 176];

% MID00385 & MID00387 - 1Tx/32Rx 7T Erlangen
% MDH_END = [72 64];

% MID00454 & MID00456 - 8Tx/32Rx 7T Erlangen
MDH_END = [120 120];

ncoil = 0;
for i=1:length(header.asc{end}.sCoilSelectMeas.aRxCoilSelectData(1).asList)
    ncoil = ncoil + header.asc{end}.sCoilSelectMeas.aRxCoilSelectData(1).asList(i).lElementSelected;
end
ADClen = header.asc{end}.sKSpace.lBaseResolution;
ADCdwell = header.asc{end}.sRXSPEC.alDwellTime(1);
TR = header.asc{end}.alTR(1);
TE = header.asc{end}.alTE(1);
slthick = header.asc{end}.sSliceArray.asSlice(1).dThickness;
roFOV = header.asc{end}.sSliceArray.asSlice(1).dReadoutFOV;
peFOV = header.asc{end}.sSliceArray.asSlice(1).dPhaseFOV;
%nslices = header.asc{end}.sGroupArray.asGroup(1).nSize;
nslices = header.asc{end}.sSliceArray.lSize;
if isfield(header.asc{end}.sGroupArray.asGroup(1),'dDistFact')
    distfac = header.asc{end}.sGroupArray.asGroup(1).dDistFact;
else
    distfac = 0;
end
%resol = header.Meas{2}.DICOM.lBaseResolution;
resol = header.asc{end}.sKSpace.lBaseResolution;
%pelines = header.Meas{2}.DICOM.lPhaseEncodingLines;
pelines = header.asc{end}.sKSpace.lPhaseEncodingLines;
overs = header.Meas{2}.DICOM.flReadoutOSFactor; % RO oversampling factor
ser_interl = 1; % series interleaved - see sSliceArray.ucMode = 0x1 for ascending, 0x2 for descending, and 0x4 for interleaved
sli_interl = 1; % multislice mode interleaved

% TRICK TO DETERMINE THE LENGTH OF THE LAST HEADER AND RETRIEVE THE
% STARTING POINT OF THE DATA IN THE RAW.DAT FILE:
% Zoom into the displayed figure and count how many float points after the
% last data point - this will be the value to assign to MDH_END (first
% value for raw data with RF and second value for raw data without RF).
nread = ADClen*overs;
necho = nslices*pelines;
precis = sprintf('%d*float32=>double',nread*2);
skipn = 0;
ftp = fopen(filenameraw,'r');
% go to 1024 floats before the end of file
fseek(ftp,-4096,1);
% read the last 1024 floats
rawend = fread(ftp,nread*2*necho,precis,skipn);
% display to visually determine how many floats after the last echo
figure;plot(abs(rawend),'marker','.');
fclose(ftp);

ftp = fopen(filenamenoise,'r');
% go to 1024 floats before the end of file
fseek(ftp,-4096,1);
% read the last 1024 floats
rawend = fread(ftp,nread*2*necho,precis,skipn);
% display to visually determine how many floats after the last echo
figure;plot(abs(rawend),'marker','.');
fclose(ftp);

% header.asc{4}.sSliceArray.alSliceAcqOrder
% header.asc{4}.sSliceArray.lSize
% header.asc{4}.sKSpace.lBaseResolution
% header.asc{4}.sKSpace.lPhaseEncodingLines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the raw data (*.out file):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('reading the raw data... ');
 
%---------------------------------------------------------------------
% reorganize the raw data in a 3D array (pelines x nread x nslices):
% (c'est la qu'on rigole quand les series et les coupes sont en 
% interleaved !!!)
% NB interleaved mode:
% -	first the even slices, then the odd slices, if even number of slices
% -	first the odd slices, then the even slices, if odd number of slices !!!
%---------------------------------------------------------------------

% acquisition order of the slices:
if sli_interl
    if mod(nslices,2)==1
        sli_order = [1:2:nslices 2:2:nslices];
    else
        sli_order = [2:2:nslices 1:2:nslices];
    end
else
    sli_order = [1:nslices];
end
% sli_order = header.asc{end}.sSliceArray.alSliceAcqOrder + 1;
  
im = zeros(pelines,resol,ncoil,nslices);
figure('position',[50 50 800 800],'color',[1 1 1],'name','Magnitude images');
for ccha=1:ncoil
    
    disp(['   coil #' num2str(ccha)]);
    disp('      - lecture des donnees');
    raw = eb_readmeas(filenameraw(1,:),ADClen*overs,nslices*pelines,ncoil,ccha,MDH_END(1));
    raw = raw(1:2:end) + sqrt(-1)*raw(2:2:end);

    disp('      - reorganisation des donnees');
    data = zeros(resol*overs,pelines,nslices);
    if ser_interl % first loop over slices (order!) then over pelines
        for k=1:pelines
            for j=1:nslices
                start_dp = (k-1)*resol*overs*nslices + (j-1)*resol*overs + 1;
                stop_dp  = (k-1)*resol*overs*nslices + j*resol*overs;
                data(:,k,sli_order(j)) = raw(start_dp:stop_dp);
            end
        end
    else % first loop over pelines then over slices (order!) 
        for j=1:nslices
            for k=1:pelines
                start_dp = (k-1)*resol*overs + (j-1)*resol*overs*pelines + 1;
                stop_dp  = k*resol*overs + (j-1)*resol*overs*pelines;
                data(:,k,sli_order(j)) = raw(start_dp:stop_dp);
            end
        end
    end

    disp('      - reconstruction des images (complexes)');
    for j=1:nslices
        tmp = fftshift(ifft2(fftshift(data(:,:,j))));
        % reduce FOV if oversampling:
        if overs==2; im(:,:,ccha,j) = tmp(resol/2+1:resol*3/2,:); 
        else im(:,:,ccha,j) = tmp; 
        end
    end

    subplot(ceil(sqrt(ncoil)),ceil(sqrt(ncoil)),ccha);
    imshow(abs(squeeze(im(:,:,ccha,round(nslices/2)))),[]);
    text(resol*0.05,resol*0.9,num2str(ccha), 'color',[1 1 1])
    
    disp(['      - save data from coil']);
    eb_save_NIFTI(abs(squeeze(im(:,:,ccha,:))), header, [maindir filesep fname  num2str(ccha, '_%0.3d_abs') '.nii']);
    eb_save_NIFTI(real(squeeze(im(:,:,ccha,:))), header, [maindir filesep fname  num2str(ccha, '_%0.3d_real') '.nii']);
    eb_save_NIFTI(imag(squeeze(im(:,:,ccha,:))), header, [maindir filesep fname  num2str(ccha, '_%0.3d_imag') '.nii']);
    
    clear raw;
    clear data;
end
set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng',[maindir filesep fname '_abs.png']);
save([maindir filesep fname '_orig_cplx.mat'],'im');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the sensitivity map. The input should be reconstructed complex
% data with dimensions (resol, pelines, ncoil, nslices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensmap = make_sense_pad(datain, 1);
disp('calculate the sensitivity map...');
sensmap = make_sense_nosmooth(im, 1);
clear im;

figure('position',[50 50 800 800],'color',[1 1 1],'name','Sensitivity maps');
for ccha=1:ncoil
    eb_save_NIFTI(abs(squeeze(sensmap(:,:,ccha,:))), header, [maindir filesep fname  num2str(ccha, '_%0.3d_sensmap') '.nii']);
    subplot(ceil(sqrt(ncoil)),ceil(sqrt(ncoil)),ccha);
    imshow(abs(squeeze(sensmap(:,:,ccha,round(nslices/2)))),[0 1]);
    %title(['coil #' num2str(ccha)]);
    text(resol*0.05,resol*0.9,num2str(ccha), 'color',[1 1 1])
end
set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng',[maindir filesep fname '_sensmap.png']);
save([maindir filesep fname '_sensmap_cplx.mat'],'sensmap');
clear sensmap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate receiver noise matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('generate receiver noise matrix...');

%---------------------------------------------------------------------
% read the raw data (*.out file):
%---------------------------------------------------------------------
disp('read the noise raw data...');

noise = zeros(pelines,resol,ncoil,nslices);
for ccha=1:ncoil
    disp(['   coil #' num2str(ccha)]);
    disp(['      - lecture des donnees']);
    raw = eb_readmeas(filenamenoise(1,:),ADClen*overs,nslices*pelines,ncoil,ccha,MDH_END(2));
    raw = raw(1:2:end) + sqrt(-1)*raw(2:2:end);

    disp(['      - reorganisation des donnees']);
    data = zeros(resol*overs,pelines,nslices);
    if ser_interl % first loop over slices (order!) then over pelines
        for k=1:pelines
            for j=1:nslices
                start_dp = (k-1)*resol*overs*nslices + (j-1)*resol*overs + 1;
                stop_dp  = (k-1)*resol*overs*nslices + j*resol*overs;
                data(:,k,sli_order(j)) = raw(start_dp:stop_dp);
            end
        end
    else % first loop over pelines then over slices (order!) 
        for j=1:nslices
            for k=1:pelines
                start_dp = (k-1)*resol*overs + (j-1)*resol*overs*pelines + 1;
                stop_dp  = k*resol*overs + (j-1)*resol*overs*pelines;
                data(:,k,sli_order(j)) = raw(start_dp:stop_dp);
            end
        end
    end

    disp(['      - reconstruction des images (complexes)']);
    for j=1:nslices
        tmp = fftshift(ifft2(fftshift(data(:,:,j))));
        % reduce FOV if oversampling:
        if overs==2; noise(:,:,ccha,j) = tmp(resol/2+1:resol*3/2,:); 
        else noise(:,:,ccha,j) = tmp; 
        end
    end
 
    disp(['      - save data from coil']);
    eb_save_NIFTI(abs(squeeze(noise(:,:,ccha,:))), header, [maindir filesep fname  num2str(ccha, '_%0.3d_noise') '.nii']);
    
    clear raw;
    clear data;
end
save([maindir filesep fname '_noise_cplx.mat'],'noise');

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate psi: the psi matrix should be close to eye(ncoil)...
% the function eb_calc_psi returns a matrix normalized to ncoil...
% the results obtained with psi = eye(ncoil) are similar to the one
% obtained with the calculated psi...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('calculate the noise matrix psi...');

%---------------------------------------------------------------------
% calculate psi (NB: possible conflict of variable "psi"
% with toolbox/matlab/specfun/psi.mexglx :o(
%---------------------------------------------------------------------
psi_calc = eb_calc_psi(noise);
use_psi=1;
if use_psi==0
    psi_calc=eye(ncoil);
end
clear noise;

save([maindir filesep 'psi.mat'],'psi_calc');

figure('color',[1 1 1],'Position',[50 50 400 400],'name','Noise matrix PSI');
subplot(1,1,1);
set(gca,'position',[0.1 0.1 0.8 0.8])
imshow(abs(psi_calc),[]);colorbar;
set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng',[maindir filesep fname '_psi.png']);

figure('color',[1 1 1],'Position',[50 50 400 400],'name','Noise matrix PSI');
subplot(1,1,1);
set(gca,'position',[0.1 0.1 0.8 0.8])
imshow(abs(psi_calc),[0 1]);colorbar;
set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng',[maindir filesep fname '_psi1.png']);

%---------------------------------------------------------------------
% smooth the sensitivity maps.
%---------------------------------------------------------------------
disp('smooth the sensitivity maps...');
load([maindir filesep fname '_sensmap_cplx.mat']);
nvox_smooth = 2;
rsmoothvol = zeros(size(squeeze(sensmap(:,:,1,:))));
ismoothvol = zeros(size(squeeze(sensmap(:,:,1,:))));
figure('position',[50 50 800 800],'color',[1 1 1],'name','Smoothed sensitivity maps');
for ccha=1:ncoil
    sensevol=squeeze(sensmap(:,:,ccha,:));    
    spm_smooth(real(sensevol),rsmoothvol,nvox_smooth);
    spm_smooth(imag(sensevol),ismoothvol,nvox_smooth);
    sensmap(:,:,ccha,:) = rsmoothvol + sqrt(-1)*ismoothvol;
    subplot(ceil(sqrt(ncoil)),ceil(sqrt(ncoil)),ccha);
    imshow(abs(squeeze(sensmap(:,:,ccha,round(nslices/2)))),[0 1]);
    text(resol*0.05,resol*0.9,num2str(ccha), 'color',[1 1 1])
end
smoosensmap = sensmap;
save([maindir filesep fname '_smoosensmap_cplx.mat'],'smoosensmap');
set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng',[maindir filesep fname '_smoosensmap.png']);
clear sensmap;

%---------------------------------------------------------------------
% calculate the gmap.
%---------------------------------------------------------------------
% pour gagner du temps en cours de calcul
psi_inv = inv(psi_calc);

figure('position',[50 50 800 800],'color',[1 1 1],'name','G-MAPS');
for R=2:8 % considers acceleration in A>>P direction only
    disp(['R = ' num2str(R)]);
    fullfov = resol;
    redfov = floor(resol/R);
    for sl = 1:nslices
        disp(['slice number ' num2str(sl)]);
        for x=1:redfov
            for y=1:fullfov
                px = [x y];
                for j=2:R
                    npx = [px(j-1,1)+(redfov) px(1,2)];
                    px  = [px;npx];
                end

                smat = zeros(ncoil,R);
                for ccha=1:ncoil
                    for j=1:R
                        smat(ccha,j)=smoosensmap(px(j,2),px(j,1),ccha,sl);
                    end
                end

                U1=(smat'*psi_inv*smat);
                pU1=pinv(U1);
                for j=1:R
                    gmap(px(j,2),px(j,1),sl)=sqrt(pU1(j,j)*U1(j,j));
                end
            end
        end
    end

    subplot(3,3,R);
    imshow(abs(gmap(:,:,round(nslices/2))),[]);
    colorbar;
    text(resol*0.05,resol*0.9,['R' num2str(R)], 'color',[1 1 1])

    disp(['      - save gmap (R=' num2str(R) ')']);
    eb_save_NIFTI(abs(gmap), header, [maindir filesep fname '_gmap_R' num2str(R) '.nii']);

end

set(gcf,'PaperPositionMode','auto')
print(gcf,'-dpng',[maindir filesep fname '_gmaps.png']);
