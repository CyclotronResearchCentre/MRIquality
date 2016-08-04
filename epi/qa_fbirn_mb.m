function QA = qa_fbirn_mb(P, PARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QA tool based on and modified from the FBIRN QA
% Friedman and Glover, JMRI 23:827-839 (2006)
%
% The metrics estimated here are
%    - SFNR = signal to fluctuation noise ratio = temporal SNR
%    - Percent Drift = signal change over time (max-min) fitted by second
%      order polynomial (divided by mean signal intensity)
%    - Percent fluctuation = SD of residuals after removing drift (divided
%      by mean signal intensity)
%    - RDC = radius of decorrelation, deviation from theoretical value in
%      Weiskoff plot 
%
% N. Weiskopf, FIL, London 09/10/06
% updated by C. Hutton and N. Weiskopf, 03/01/07, for online application
% updated by E. Balteau 2014 for compatibility with Matlab R2012b.
% updated by E. Balteau 2016 for multi-band EPI acquisitions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N_max = length(PARAMS.x_roi);

% load images
VIM = spm_vol(P);
YIM = spm_read_vols(VIM);
dim = VIM(1).dim;
num_vols = size(YIM,4);

% Calculate SD volume over the series
accsum  = sum(YIM,4);
accsqr  = sum(YIM.^2,4);
varvol  = accsqr/num_vols - (accsum/num_vols).^2;
stdvol  = sqrt(varvol);
% Save SD volume
QA.sdmap   = [PARAMS.resfnam '_SD.nii'];
dm         = VIM(1).dim;
dt         = [spm_type('float32'),spm_platform('bigend')];
Ni         = nifti;
Ni.mat     = VIM(1).mat;
Ni.mat0    = VIM(1).mat;
Ni.descrip = 'voxel-wise standard deviation over the time series';
Ni.dat     = file_array(fullfile(PARAMS.outdir,QA.sdmap), dm, dt, 0, 1, 0);
create(Ni);
Ni.dat(:,:,:) = stdvol;


% Extract voxel time series from ROI
% Data, detrended data, and 2nd order polynomial drift fit
signalslicearray = squeeze(YIM(:,:,PARAMS.signalplane,:)); % signal plane
signalroiarray = signalslicearray(PARAMS.x_roi,PARAMS.y_roi,:);% ROI in the signal plane
detrslicearray     = zeros([N_max N_max num_vols]);
driftslicearray    = zeros([N_max N_max num_vols]);
for x_nr=1:N_max
    for y_nr=1:N_max
        p = polyfit(1:num_vols,squeeze(signalroiarray(x_nr,y_nr,:))',2);
        yfit = polyval(p, 1:num_vols);
        detrslicearray(x_nr,y_nr,:) = squeeze(signalroiarray(x_nr,y_nr,:))'-yfit;
        driftslicearray(x_nr,y_nr,:) = yfit;
    end
end

% tSNR CALCULATION -> requires mean and std after detrending
Ydetrend = YIM*0.0;
mean_im = mean(YIM,4);
mask = mean_im;
threshold = (prctile(mask(:),98)-prctile(mask(:),2))*0.2+prctile(mask(:),2);
mask = (mask>threshold);
mean_im = mean_im.*mask;
for i=1:dim(1)
    for j=1:dim(2)
        for k=1:dim(3)
            if (mask(i,j,k))
                p = polyfit(1:num_vols,squeeze(YIM(i,j,k,:))',2);
                yfit = polyval(p, 1:num_vols);
                Ydetrend(i,j,k,:) = squeeze(YIM(i,j,k,:))'-yfit;
            end
        end
    end
end
sigt_im = std(Ydetrend,0,4);
tsnr_im = (mean_im./sigt_im).*mask;
QA.tsnrmap   = [PARAMS.resfnam '_tsnrmap.nii'];
dm         = VIM(1).dim;
dt         = [spm_type('float32'),spm_platform('bigend')];
Ni         = nifti;
Ni.mat     = VIM(1).mat;
Ni.mat0    = VIM(1).mat;
Ni.descrip = 'tSNR volume';
Ni.dat     = file_array(fullfile(PARAMS.outdir,[PARAMS.resfnam '_tSNR.nii']), dm, dt, 0, 1, 0);
create(Ni);
Ni.dat(:,:,:) = tsnr_im;

% Intensity drift analysis
QA.data_drift_plot = squeeze(mean(mean(signalroiarray)));
QA.fit_drift_plot = squeeze(mean(mean(driftslicearray)));
QA.perc_drift = (max(QA.fit_drift_plot)-min(QA.fit_drift_plot))/mean(signalroiarray(:))*100;

% Fluctuation plot
QA.perc_fluct_plot = squeeze(mean(mean(detrslicearray)))/mean(signalroiarray(:))*100;
% according to FBIRN the fluctuation of the ROI time series
QA.perc_fluct = std(mean(mean(detrslicearray,1),2))/mean(signalroiarray(:))*100; 

% FFT of residuals
FT_resid = fftshift(fft(detrslicearray/mean(signalroiarray(:))*100,[],3));
QA.FT_resid = squeeze(mean(mean(abs(FT_resid),1),2));
QA.fcoord = ((1:size(FT_resid,3))-(size(FT_resid,3)+1)/2)/(0.001*PARAMS.TR*size(FT_resid,3)');

% Temporal SFNR, estimated as the mean SFNR across the voxels in the ROI
SFNR_voxel = mean(signalroiarray,3)./std(detrslicearray,[],3);
QA.SFNR_voxel = squeeze(mean(mean(SFNR_voxel))); % average SFNR value across voxels as defined by FBIRN

% Weiskoff test and RDC value
for i_perm = 1:10
    r_perm1 = randperm(N_max);
    r_perm2 = randperm(N_max);
    x_roi_p = r_perm1;
    y_roi_p = r_perm2;
    for i = 1:N_max
        mean_roi_val = mean(mean(mean(signalroiarray(x_roi_p(1:i),y_roi_p(1:i),:))));
        froi_val = std(mean(mean(detrslicearray(x_roi_p(1:i),y_roi_p(1:i),:),1),2));
        SFNR(i,i_perm) = mean_roi_val/froi_val;
    end
end
QA.SFNR = mean(SFNR,2);
QA.rdc = QA.SFNR(N_max)/QA.SFNR(1);

% Estimate x and y diameter of object from slice
meanvol = accsum/num_vols;

interp_fact=10;

% figure;
maxi=max(meanvol(:));
dia=0;
for i=1:size(meanvol,3);
    for j=1:size(meanvol,2);
        if ~isempty(find(meanvol(:,j,i)>maxi/2));
            meani=mean(meanvol(PARAMS.x_roi,j,i));
            interpx=interp(meanvol(:,j,i),interp_fact);
            n=find(interpx>meani/3);
            %plot(interpx);hold on;
            if ((n(end)-n(1))>dia);
                dia=(n(end)-n(1));
            end;
        end
    end
end
QA.x_diam = dia/interp_fact*sqrt(sum((VIM(1).mat*[1 0 0 0]').^2));

dia=0;
for i=1:size(meanvol,3);
    for j=1:size(meanvol,1);
        if ~isempty(find(meanvol(j,:,i)>maxi/2));
            meani=mean(meanvol(j,PARAMS.y_roi,i));
            interpy=interp(meanvol(j,:,i),interp_fact);
            n=find(interpy>meani/2);
            %plot(interpy);hold on;
            if ((n(end)-n(1))>dia);
                dia=(n(end)-n(1));
            end;
        end
    end
end
QA.y_diam = dia/interp_fact*sqrt(sum((VIM(1).mat*[1 0 0 0]').^2));

% Write results in file
fid = fopen(fullfile(PARAMS.outdir, [PARAMS.resfnam '.txt']),'a');
fprintf(fid,'\nSTABILITY RESULTS\n');
fprintf(fid,'    SD map: %s\n', QA.sdmap);
fprintf(fid,'    tSNR map: %s\n', QA.tsnrmap);
fprintf(fid,'    Intensity: drift = %6.5f | fluctuation = %6.5f (percent)\n', QA.perc_drift, QA.perc_fluct);
fprintf(fid,'    tSNR in central ROI: %5.2f\n', QA.SFNR_voxel);
fprintf(fid,'    RDC = %5.2f\n', QA.rdc);
fprintf(fid,'    Estimated phantom diameter: %4.1f x %4.1f mm2\n', QA.x_diam, QA.y_diam);
fclose(fid);


return;


