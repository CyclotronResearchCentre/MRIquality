%=========================================================================%
% This file is part of the Quality Control Toolbox (TCQ)
% Copyright (C) 2013 - Cyclotron Research Centre
% University of Liege, Belgium
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%=========================================================================%

function out = qc_spike_check(job)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: automated procedure to detect spikes and correct for them
% when possible, applying a simple correction strategy to reduce the
% spurious variance introduced by spikes in an fMRI time series. 
% Always check the output data before using them for further analysis.
% Never assume it worked properly ;)!! 
%--------------------------------------------------------------------------
% Warning and disclaimer: This software is for research use only. 
% Do not use it for clinical or diagnostic purposes.
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - 2013
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Methods: the noise level is calculated for each slice and for each volume
% using a percentile-based mask eroded with a disk of radius 1.
% The time course of the noise level for each slice is considered.
% If the std of the time course is >10% of the mean noise level then the
% data are too strongly affected by spikes to be corrected!
% Otherwise, the spike_correc detects slices that have a noise level above
% the average noise level + 4 SD and replace them by a linear
% interpolation between the previous and the next non affected slices.
%
% Note that the reason(s) why an error "too strongly affected by spikes" is
% returned might be different than actually stated, i.e. not due to spikes
% but to e.g. subject movement or other system instabilities. Examine the
% output display with care when using this function.
%
% In order to check for miscorrection of spikes, please look carefully at
% the output screen shots:
% - the "spike_detected_*.png" images show the various spikes detected.
% Check whether these were spikes indeed!
% - the noise time course shows the background noise time course with the
% average +/- SD outlined in red and black, respectively. Spikes are
% detected based on these data whenever the background noise exceeds the
% average value + 4 SDs. Check consistency between these plots and the
% detected spikes. If any suspicious drift or jump in the background noise
% is observed, please notify the physics group for support.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for in=1:numel(job.subj)
    local_job.subj = job.subj(in);
    local_job.params_spike  = job.params_spike;
    out_loc = qc_spike_check_local(local_job);
    out.subj(in).spcfiles = out_loc.spcfiles;
end
end

function out_loc = qc_spike_check_local(job)

disp('----- Automated spike detection and correction tool -----');
fprintf(['\nWARNING: Always check the output of this procedure before going ahead\n'...
         '           with your analysis. Never assume it worked properly ;)!!\n\n']);

P = char(job.subj.raws);

% to retrieve output directory and a few other useful parameters
if isfield(job.subj(1).output,'indir') && job.subj(1).output.indir == 1
    outdir = fileparts(P(1,:));
else
    outdir = job.subj(1).output.outdir{1};
end

% whether correction should be applied or not (if possible)
bCorrect = job.params_spike.do_correct;
% is considered as "spike correctable" a data set where the noise variation
% is not greater than corrThresh*(average noise level)
corrThresh = job.params_spike.correction_threshold;
% is considered as spiky a slice where the noise level is higher then
% spikeThresh*(standard deviation of noise)
spikeThresh = job.params_spike.spike_threshold;

disp('----- Loading EPI series -----');

ndisc = 2; % discard the first two volumes to avoid any bias
P = P(ndisc+1:end,:);
V = spm_vol(P);
Y = spm_read_vols(V);

% Y is a 4D-matrix containing the time series of functional images
[ni,nj,nk,nl] = size(Y);
[p,n] = fileparts(V(1).fname);
mosaic_size = ceil(sqrt(nk));

% calculate mean volume:
Ymean = squeeze(mean(Y,4));

% use histograms to define threshold and erode the resulting mask to
% slightly extend the amount of excluded voxels, multiply mask over slices,
% more restrictive but safer!
threshold = (prctile(Ymean(:),98)-prctile(Ymean(:),2))*0.02+prctile(Ymean(:),2);
Ymask = Ymean<threshold;
se = strel('disk',1);
Ymask = imerode(Ymask,se);
maskprod = prod(double(Ymask),3);
Ymask = repmat(maskprod,[1 1 nk]);

disp('----- Analysis of the general noise level -----');

% let's calculate the noise level of each slice:
noise_slices = zeros(nk,nl);
for cslice = 1:nk
    current_mask = squeeze(Ymask(:,:,cslice));
    for nvol = 1:nl
        current_slice = squeeze(Y(:,:,cslice,nvol));
        noise_slices(cslice,nvol) = std(current_slice(find(current_mask==1)));
    end
    % %     % remove quadratic trend
    % %     pfit = polyfit(1:nl,squeeze(noise_slices(cslice,:)),2);
    % %     trend = pfit(1)*[1:nl].^2 + pfit(2)*[1:nl];
    % %     noise_slices(cslice,:) = squeeze(noise_slices(cslice,:))-trend;
end

average_noise = squeeze(mean(noise_slices,2));
variation_of_noise = squeeze(std(noise_slices,0,2));

% display figure with noise level time course for each slice
figure('color',[1 1 1],'position',[50 50 800 800]);
for cslice = 1:nk
    subplot(mosaic_size,mosaic_size,cslice);
    plot(noise_slices(cslice,:));
    line([0 nl],[1 1]*average_noise(cslice),'color',[1 0 0])
    line([0 nl],[1 1]*(average_noise(cslice)+variation_of_noise(cslice)),'color',[0 0 0]);
    line([0 nl],[1 1]*(average_noise(cslice)-variation_of_noise(cslice)),'color',[0 0 0]);
    set(gca,'XLim',[0 nl],'FontSize',5);
end
print(gcf,'-dpng',[outdir '\' n '_noise_time_course.png']);

% % display figure with average and std of noise for each slice
% figure('color',[1 1 1],'position',[50 50 800 800]);
% subplot(1,1,1);
% errorbar(1:nk,average_noise,average_noise);
% xlabel('slice number');
% title('Average noise amplitude +/- standard deviation for each slice');
% print(gcf,'-dpng',[outdir '\' n '_noise_level.png']);

% Spike correction if required (and possible)
if (bCorrect); textdisp = '(required by user)';
else textdisp = '(not required by user - just checking)';end
disp(['----- Spike correction ' textdisp ' -----']);

if (find(variation_of_noise./average_noise>corrThresh))
    bPossible = 0;
    fprintf(['\nSORRY...\n'...
        '\tThe data set is too strongly affected by spikes and can''t be corrected,\n'...
        '\tor spikes are not the only problem in this data set (subject moving, coil\n'...
        '\tinstabilities, other artefacts...). You may try and change the threshold \n'...
        '\tfor spike correction.\n\n']);
else bPossible = 1;
end

spike_count = 0;
spike_list = [];
spcY = Y;

for cslice = 1:nk
    % NB: the default threshold of 4*SD might sometimes not be appropriate...
    spiky = find(noise_slices(cslice,:)>average_noise(cslice)+spikeThresh*variation_of_noise(cslice));
    spike_count = spike_count + length(spiky);
    for m=1:length(spiky)
        spike_list = [spike_list;cslice spiky(m)+ndisc];

        if (bCorrect && bPossible)
            % determine the slices used to correct the affected slice
            ok_slice_before = spiky(m)-1;
            ok_slice_after = spiky(m)+1;
            while (~isempty(find(spiky==ok_slice_before,1)) && (ok_slice_before>0))
                ok_slice_before = ok_slice_before-1;
            end
            while (~isempty(find(spiky==ok_slice_after,1)) && (ok_slice_after<nl+1))
                ok_slice_after = ok_slice_after+1;
            end
            if (ok_slice_before==0); ok_slice_before=ok_slice_after; end
            if (ok_slice_after==nl+1); ok_slice_after=ok_slice_before; end

            % linear interpolation between ok_slice_before and ok_slice_after:
            slice_before = squeeze(Y(:,:,cslice,ok_slice_before));
            slice_after = squeeze(Y(:,:,cslice,ok_slice_after));
            if (ok_slice_before == ok_slice_after)
                spcY(:,:,cslice,spiky(m)) = slice_before;
            else
                slope = (slice_after-slice_before)/(ok_slice_after-ok_slice_before);
                spcY(:,:,cslice,spiky(m)) = slope*spiky(m) + slice_before - slope*ok_slice_before;
            end
        end
    end
end

if (bCorrect && bPossible)
    % save the volumes with "spc" for spike corrected
    if (ndisc)
        fprintf(['\nWARNING: The first %d images were discarded to avoid the interference of the \n'...
         '\tT1 saturation effect with the spike check procedure. These images haven''t been spike\n'...
         '\tcorrected and won''t be considered in the next steps of the analysis.\n\n'], ndisc);
    end
    for nvol = 1:nl
        [cpath,cnam,cext] = fileparts(V(nvol).fname);
        V(nvol).fname = [outdir '/spc' cnam cext];
        spm_write_vol(V(nvol),spcY(:,:,:,nvol));
    end
end

% save log file with processing information and pictures for visual
% inspection of the spiky volumes
[cpath,cnam] = fileparts(P(1,:));
fid = fopen([outdir '\' cnam '_spike_check_log.txt'],'w');
fprintf(fid,'LOG FILE FOR SPIKE CHECK\n');
fprintf(fid,'------------------------\n\n');
fprintf(fid,'Correction required: %d\n',bCorrect);
fprintf(fid,'Correction possible: %d\n',bPossible);
fprintf(fid,'Correction threshold: %4.1f\n',corrThresh);
fprintf(fid,'Spike detection threshold: %4.1d\n\n',spikeThresh);
if (isempty(spike_list))
    fprintf(fid,'No spike detected!\n');
else
    figure('color',[1 1 1],'position',[50 50 600 600]);
    fprintf(fid,'%d spikes detected:\n',size(spike_list,1));
    for cvol=1:nl
        % detect whether there are spikes in the current volume
        spiky = find(spike_list(:,2)-ndisc == cvol);
        if (~isempty(spiky))
            imshow(get_mosaic(squeeze(Y(:,:,:,cvol))),[0 50]);
            list4title = '';
            ns = length(spiky);
            [cpath,cnam] = fileparts(P(cvol,:));
            for cs=1:ns;
                cslice = spike_list(spiky(cs),1);
                list4title = [list4title num2str(cslice) ','];
                spikint = (noise_slices(cslice,cvol)-average_noise(cslice))/variation_of_noise(cslice);
                fprintf(fid,'%s - slice %d \t- intensity %5.2f\n',cnam,spike_list(spiky(cs),1),spikint);
            end
            title([num2str(ns) ' spike(s) detected in ' cnam ' (' list4title(1:end-1) ')']);
            print(gcf,'-dpng',[outdir '\' cnam '_' num2str(ns) '_spike(s)_detected.png']);
        end
    end
end
fclose(fid);

% generate the list of output files (dependencies)
% if correction not required by user, return empty list
if (bCorrect)
    out_loc.spcfiles = cell(length(job.subj.raws)-ndisc, 1);
    % if the correction was required but couldn't be performed, the non-
    % corrected files are listed as output files (to avoid the batch to
    % crash if this happens). In that case, there won't be any spc prefix
    % in the file names at the subsequent processing steps.
    if (bPossible); prefx='spc';else prefx='';end
    for i=ndisc+1:numel(job.subj.raws),
        [pth,nam,ext,num] = spm_fileparts(job.subj.raws{i});
        out_loc.spcfiles{i-ndisc} = fullfile(pth,[prefx, nam, ext, num]);
    end
else
    out_loc.spcfiles = {'Spike correction not enabled.'};
end
end