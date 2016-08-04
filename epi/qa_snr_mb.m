function SNR = qa_snr_mb(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USAGE: eb_snr_mb(imfiles, noisefiles, PARAMS)
%    imfiles     list of file names of the images. The first few volume
%                have already been discarded to avoid T1 saturation
%                effects. The images in the time series are already  
%                realigned and resliced.
%    noisefiles  zero (empty), one or more noise images acquired without RF
%                to estimate the noise distribution.
%
% Written by Evelyne Balteau - Cyclotron Research Centre - January 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global N fitn fitx % global variables used for histogram fit only

%% GET STARTED
if nargin == 0
    imfiles = spm_select(Inf,'image','Select images (acquired with RF)')
    noisefiles = spm_select(Inf,'image','Select noise images (acquired without RF)')
    PARAMS.MB = input('Enter MB acceleration factor: ');
    PARAMS.PAT = input('Enter PAT acceleration factor: ');
    PARAMS.comment = input('Enter comment if required: ','s');
    PARAMS.tag = input('Enter tag if several runs today: ','s');
    PARAMS.resfnam = 'results.txt';
    PARAMS.TR = input('Enter TR (ms): ');
    PARAMS.ncha = input('Enter number of channels: ');
    PARAMS.MB = input('Enter MultiBand factor: ');
    PARAMS.PAT = input('Enter PAT factor (PEx3D): ');
    PARAMS.outdir = '.';
    PARAMS.signalplane = input('Enter slice of refence (signal): ');
    PARAMS.noiseplane = input('Enter slice of refence (noise): ');
else
    imfiles = varargin{1};
    noisefiles = varargin{2};
    PARAMS = varargin{3};
end
PARAMS.datadir = fileparts(imfiles(1,:));


%% LOAD DATA
% load images
VIM = spm_vol(imfiles);
YIM = spm_read_vols(VIM);
dim = VIM(1).dim;
num_vols = size(YIM,4);
% load noise images
if ~isempty(noisefiles)
    VNO = spm_vol(noisefiles);
    YNO = spm_read_vols(VNO);
end
% data subsets
noiseslicearray = squeeze(YIM(:,:,PARAMS.noiseplane,:)); % noise plane
noiseroiarray = noiseslicearray(PARAMS.x_roi,PARAMS.y_roi,:); % ROI in the noise plane
signalslicearray = squeeze(YIM(:,:,PARAMS.signalplane,:)); % signal plane
signalroiarray = signalslicearray(PARAMS.x_roi,PARAMS.y_roi,:);% ROI in the signal plane

%% AVERAGE SIGNAL IN ROI
signal = mean(signalroiarray(:));

%% ESTIMATE NOISE AND SNR FROM noRF IMAGES
if ~isempty(noisefiles)
    noisy_vox = YNO(:);
    noisy_vox = noisy_vox(:);
    noisy_vox = noisy_vox(noisy_vox~=0); % to get rid of the bunch of zero values from edges of volume in Siemens images
    
    % plot histogram
    maxi = max(noisy_vox(:));
    x = 0:round(maxi); % siemens reconstruction has integer values only
    n = hist(noisy_vox(:),x);
    figure('color',[1 1 1],'position',[50 50 1000 800]);
    bar(x,n);hold on;
    
    % histogram fitting with central Chi distribution
    fitn = n; fitx = x;
    % first estimate of the parameters (see notebook 28/01/2016):
    N = PARAMS.ncha/PARAMS.MB/PARAMS.PAT;
    sig = x(n==max(n(:)))/sqrt(2*N-1);
    A = max(n(:))*gamma(N)*sig*2^(N-1)*exp(N-0.5)*(2*N-1)^(0.5-N);
    % adjust all parameters (3-params fit)
    [par3, fval3] = fminsearch(@minfun_CC3,[sig A N]);
    
    % plot the result of the fit
    xplot = 0:0.1:maxi;
    fit3 = central_Chi_PDF3(par3,xplot);
    plot(xplot,fit3,'color',[0 1 1]);
    
    % add title, axes labels, adjust axes...
    set(gca,'Xlim',[0 max(xplot)]);
    maxY = max(fit3);
    maxX = max(xplot);
    title(sprintf('Central Chi distribution (%d channels, MB%d, PAT%d)', PARAMS.ncha, PARAMS.MB, PARAMS.PAT));
    text2plot = [sprintf('Parameters estimated: \n    sigma = %5.2f \n    A = %6.0f \n    Neff = %5.2f \n    Residuals = %4.3e', par3(1), par3(2), par3(3), fval3) ...
        sprintf('\n\nAverage signal in central ROI = %5.2f', signal) ...
        sprintf('\n\nSNR = %5.2f', signal/par3(1))];
    text(maxX*0.04,0.95*maxY,text2plot,'VerticalAlignment','top');
    legend('data','fit');
    
    % save figure
    set(gcf,'PaperPositionMode','auto')
    print(gcf,'-dpng',fullfile(PARAMS.outdir,[PARAMS.resfnam '_NOISE_DISTRIB.png']));
    
    % save results in structure SNR
    SNR.sigma_noRF = par3(1);
    SNR.nchaeff_noRF = par3(3);
    SNR.snr_noRF = signal/par3(1);
    SNR.snrmap = [PARAMS.resfnam '_SNR_noRF.nii'];
    
    % save SNR volume
    dm         = VIM(1).dim;
    dt         = [spm_type('float32'),spm_platform('bigend')];
    Ni         = nifti;
    Ni.mat     = VIM(1).mat;
    Ni.mat0    = VIM(1).mat;
    Ni.descrip = 'SNR volume based on noRF noise estimation';
    Ni.dat     = file_array(fullfile(PARAMS.outdir,SNR.snrmap_noRF), dm, dt, 0, 1, 0);
    create(Ni);
    Ni.dat(:,:,:) = mean(YIM,4)/SNR.sigma_noRF;
end

%% VARIOUS OTHER SNR CALCULATIONS
n = PARAMS.ncha;
if ~isempty(noisefiles); neff = SNR.nchaeff_noRF; end
std_noise = zeros(1,num_vols);

% DIETRICH1: O. Dietrich et al. MRI 2008;26:754-762.
% correction factor for noncentral chi-distribution for n channels (0.6551
% for n = 1, i.e. Rayleigh distribution) 
for kk = 1:num_vols 
    tmp = noiseslicearray(PARAMS.x_roi,PARAMS.y_roi,kk);
    std_noise(kk) = std(tmp(:));
end
corr_noise_distrib = sqrt(2*n-(prod(1:2:2*n-1)/(2^(n-1)*prod(1:n-1)))^2*pi/2);
SNR.snr_DIETRICH1 = mean(signalroiarray(:))/sqrt(mean(std_noise.^2))*corr_noise_distrib;
% using the effective number of coils caluclated above
if ~isempty(noisefiles)
    corr_noise_distrib = sqrt(2*neff-(prod(1:2:2*neff-1)/(2^(neff-1)*prod(1:neff-1)))^2*pi/2);
    SNR.snr_DIETRICH1_neff = mean(signalroiarray(:))/sqrt(mean(std_noise.^2))*corr_noise_distrib;
end

% CONSTANTINIDES et al. 1997...
% sigma = sqrt(mean(noise_pixels.^2)/2n) where n = number of coil elements.
corr_noise_distrib = 1;
% SNR Constantinides calculated from slice out of phantom
SNR.snr_CONSTANTINI = mean(signalroiarray(:))/sqrt(mean(noiseroiarray(:).^2)/2/n)*corr_noise_distrib;
if ~isempty(noisefiles)
    SNR.snr_CONSTANTINI_neff = mean(signalroiarray(:))/sqrt(mean(noiseroiarray(:).^2)/2/neff)*corr_noise_distrib;
end

% and the FBIRN way (Friedman and Glover 2006): 
% DIFF = sumODD - sumEVEN (calculated over the signal slice only)
% NB: requires even number of images in order to have the sumODD and
% sumEVEN operating on identical number of images!
DIFFimage = sum(signalroiarray(:,:,1:2:end-1),3) - sum(signalroiarray(:,:,2:2:end),3);
SNR.snr_FRIEDMAN = mean(signalroiarray(:))/sqrt(var(DIFFimage(:))/num_vols);

% DIETRICH2: Dietrich et al. 2007
tmp = zeros(1,num_vols-1);
for kk = 1:num_vols-1
    S = signalroiarray(:,:,kk) + signalroiarray(:,:,kk+1);
    D = signalroiarray(:,:,kk) - signalroiarray(:,:,kk+1);  
    tmp(kk) = mean(S(:))/std(D(:))/sqrt(2);
end
SNR.snr_DIETRICH2 = mean(tmp);

% Write results in file
fid = fopen(fullfile(PARAMS.outdir, [PARAMS.resfnam '.txt']),'a');
if ~isempty(noisefiles);
    fprintf(fid,'\nNOISE RESULTS IN ROI\n');
    fprintf(fid,'Central Chi distribution (%d channels, MB%d, PAT%d)\n', PARAMS.ncha, PARAMS.MB, PARAMS.PAT);
    fprintf(fid,'Parameters estimated: \n    sigma = %5.2f \n    A = %6.0f \n    Neff = %5.2f \n    Residuals = %4.3e\n', par3(1), par3(2), par3(3), fval3);
end
fprintf(fid,'\nSNR RESULTS\n');
if ~isempty(noisefiles);fprintf(fid,'    noRF estimate:         %5.2f\n', SNR.snr_noRF);end
fprintf(fid,'    DIETRICH1:             %5.2f\n', SNR.snr_DIETRICH1);
if ~isempty(noisefiles);fprintf(fid,'    DIETRICH1 (Neff):      %5.2f\n', SNR.snr_DIETRICH1_neff);end
fprintf(fid,'    DIETRICH2:             %5.2f\n', SNR.snr_DIETRICH2);
fprintf(fid,'    FRIEDMAN:              %5.2f\n', SNR.snr_FRIEDMAN);
fprintf(fid,'    CONSTANTINIDES:        %5.2f\n', SNR.snr_CONSTANTINI);
if ~isempty(noisefiles);fprintf(fid,'    CONSTANTINIDES (Neff): %5.2f\n', SNR.snr_CONSTANTINI_neff);end
fclose(fid);


%% Additional functions for histogram fit
function f = central_Chi_PDF(par,m)
% non-central Chi distribution
global N
f = par(2)*m.^(2*N-1).*exp(-m.^2./(2*par(1)^2))./(2^(N-1)*par(1)^(2*N)*gamma(N));

function f = central_Chi_PDF3(par,m) % 3 parameters to be estimated incl. N
% non-central Chi distribution
f = par(2)*m.^(2*par(3)-1).*exp(-m.^2./(2*par(1)^2))./(2^(par(3)-1)*par(1)^(2*par(3))*gamma(par(3)));

function f = central_Chi_PDF4(par,m) % 4 parameters to be estimated incl. translation along x-axis
% non-central Chi distribution
f = par(2)*(m-par(4)).^(2*par(3)-1).*exp(-(m-par(4)).^2./(2*par(1)^2))./(2^(par(3)-1)*par(1)^(2*par(3))*gamma(par(3)));
f((m-par(4))<0) = 0;

function f = central_Chi_PDF3b(par,m) % sigma, ampl and translation (N=1)
% non-central Chi distribution
global N
f = par(2)*(m-par(3)).^(2*N-1).*exp(-(m-par(3)).^2./(2*par(1)^2))./(2^(N-1)*par(1)^(2*N)*gamma(N));
f((m-par(3))<0) = 0;

function f = minfun_CC(par)
global fitn fitx
f = sum((central_Chi_PDF(par,fitx) - fitn).^2);

function f = minfun_CC3(par)
global fitn fitx
f = sum((central_Chi_PDF3(par,fitx) - fitn).^2);

function f = minfun_CC4(par)
global fitn fitx
f = sum((central_Chi_PDF4(par,fitx) - fitn).^2);

function f = minfun_CC3b(par)
global fitn fitx
f = sum((central_Chi_PDF3b(par,fitx) - fitn).^2);

function f = non_central_Chi_PDF(par,m) 
% estimation includes sigma, A (constant amplitude) and eta (some
% background true signal level!?)
% non-central Chi distribution
global N
f = par(2)*m.^N.*par(3)^(1-N).*exp(-(m.^2+par(3)^2)./(2*par(1)^2)).*besseli(N-1,(m*par(3)/par(1)^2));

function f = non_central_Chi_PDF3(par,m)  % see notebook 28/01/2016
% same as above, estimating N (number of channel) as well
% estimation includes sigma, A (constant amplitude) and eta (some
% background true signal level!?)
% non-central Chi distribution
f = par(2)*m.^par(4).*par(3)^(1-par(4)).*exp(-(m.^2+par(3)^2)./(2*par(1)^2)).*besseli(par(4)-1,(m*par(3)/par(1)^2));

function f = minfun_nCC(par)
global fitn fitx
f = sum((non_central_Chi_PDF(par,fitx) - fitn).^2);

function f = minfun_nCC3(par)
global fitn fitx
f = sum((non_central_Chi_PDF3(par,fitx) - fitn).^2);

