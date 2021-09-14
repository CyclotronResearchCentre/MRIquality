function out = mriq_run_epi_qareport(job)
%==========================================================================
% USAGE: mriq_run_epi_qareport(job)
% QA tools for longitudinal visualization of the results generated over a
% given time interval.
% Implemented as part of the MRI quality toolbox (partim EPI)
%==========================================================================
% Written by Evelyne Balteau 
% Cyclotron Research Centre - May 2021
%==========================================================================

% Set paths
mriq_defaults; % TOBEDELETED!!! Must be handled by the Config module...
PATHS.input = fileparts(job.RESlist{1});
PATHS.output = cell2mat(mriq_get_defaults('path_output'));
if ~exist(PATHS.output,'dir')
    try
        mkdir(PATHS.output);
        fprintf(1,['\nWARNING: %s' ...
            '\nThe specified output directory does not exist and has been automatically created.' ...
            '\n'],PATHS.output);
    catch
        error(['\nERROR: %s' ...
            '\nThe specified output directory does not exist and could not be created.' ...
            '\nPlease use the Configure module to define the output directory.' ...
            '\n'],PATHS.output);
    end
end

% The SPM Graphics window must be available (or maybe not, but this is
% doing no harm... will be used to display the realign results)...
fg = spm_figure('FindWin','Graphics');
if isempty(fg)
    spm_figure('Create','Graphics','Graphics','on')
end

DATA.logList = char(job.RESlist);
% Filters and tags
COILfilter = char(job.RESfilter_coil);
SCANNERfilter = char(job.RESfilter_scanner);
REStag = char(job.REStag);

%==========================================================================
% Initialize structure to collect the data from each time point
%==========================================================================
nlog = size(DATA.logList,1);
initcellfield = cell(nlog, 1);
initcellfield(:) = {''}; 
DATA.toPlot = struct( ... 
    'ACQ',struct( ...
        'date', zeros(nlog, 1), ...              % log.PARAMS.date: "20190131",
        'scanner', zeros(nlog, 1), ...           % log.PARAMS.scanner: "Prisma",
        'ref_ampl', zeros(nlog, 1), ...         % log.PARAMS.ref_ampl: 222.432723999,
        'rf_freq', zeros(nlog, 1), ...          % log.PARAMS.rf_freq: 123254542,
        'field_strength', zeros(nlog, 1), ...   % log.PARAMS.field_strength: 2.89362001419,
        'SAR', zeros(nlog, 1), ...              % log.PARAMS.SAR: 0.02550401926618,
        'coils', zeros(nlog, 1), ...             % log.PARAMS.coils: "HC1-6",
        'ncha', zeros(nlog, 1)), ...            % log.PARAMS.ncha: 44
    'RES', struct( ...
    	'sigma_noRF', zeros(nlog, 1), ...           % log.SNR.sigma_noRF: 35.81144175649234,
		'nchaeff_noRF', zeros(nlog, 1), ...         % log.SNR.nchaeff_noRF: 40.15373166349986,
		'snr_noRF', zeros(nlog, 1), ...             % log.SNR.log.SNR.snr_noRF: 336.6147228758923,
		'snr_DIETRICH1', zeros(nlog, 1), ...        % log.SNR.snr_DIETRICH1: 46.38108748927534,
		'snr_DIETRICH1_neff', zeros(nlog, 1), ...   % log.SNR.snr_DIETRICH1_neff: 46.37440198860684,
		'snr_CONSTANTINI', zeros(nlog, 1), ...      % log.SNR.snr_CONSTANTINI: 232.6461743685336,
		'snr_CONSTANTINI_neff', zeros(nlog, 1), ... % log.SNR.snr_CONSTANTINI_neff: 222.2452733805509,
		'snr_FRIEDMAN', zeros(nlog, 1), ...         % log.SNR.snr_FRIEDMAN: 178.9761725093214,
		'snr_DIETRICH2', zeros(nlog, 1), ...        % log.SNR.snr_DIETRICH2: 176.1695022642294
      	'perc_drift', zeros(nlog, 1), ...           % log.FBIRN.perc_drift: 1.9890024383606,
		'perc_fluct', zeros(nlog, 1), ...           % log.FBIRN.perc_fluct: 0.0516477692122591,
		'SFNR_voxel', zeros(nlog, 1), ...           % log.FBIRN.SFNR_voxel: 175.8265856552203,
		'rdc', zeros(nlog, 1)) ...                  % log.FBIRN.rdc: 11.25787672649888
	);
% special initialisation for fields that are cell arrays (not sure why but
% if directly initialised above, DATA.toPlot.ACQ is a [42x1 struct] intead
% of a [1x1 struct])...
DATA.toPlot.ACQ.date = initcellfield;
DATA.toPlot.ACQ.scanner = initcellfield;
DATA.toPlot.ACQ.coils = initcellfield;


%==========================================================================
% retrieve data from the list 
%==========================================================================
% NB: this is not elegant, but I'll live with it...
for clog = 1:nlog
    cfile = DATA.logList(clog,:);
    log = spm_jsonread(cfile);
    try DATA.toPlot.ACQ.date{clog} = log.PARAMS.date; 
    catch
        fprintf(1,'\nWARNING: no field PARAMS.date in %s\n', cfile);
    end
    try 
        DATA.toPlot.ACQ.scanner{clog} = log.PARAMS.scanner; 
        % to filter data according to the scanner of interest, if any
        %disp(log.PARAMS.scanner);
        if ~isempty(SCANNERfilter)
            DATA.toPlot.ACQ.scannerid(clog) = ~isempty(strfind(lower(log.PARAMS.scanner),lower(SCANNERfilter)));
        else
            DATA.toPlot.ACQ.scannerid(clog) = true;
        end

    catch
        fprintf(1,'\nWARNING: no field PARAMS.scanner in %s\n', cfile);
    end
    try DATA.toPlot.ACQ.ref_ampl(clog) = log.PARAMS.ref_ampl; 
    catch
        fprintf(1,'\nWARNING: no field PARAMS.ref_ampl in %s\n', cfile);
    end
    try DATA.toPlot.ACQ.rf_freq(clog) = log.PARAMS.rf_freq; 
    catch
        fprintf(1,'\nWARNING: no field PARAMS.rf_freq in %s\n', cfile);
    end
    try DATA.toPlot.ACQ.field_strength(clog) = log.PARAMS.field_strength; 
    catch
        fprintf(1,'\nWARNING: no field PARAMS.field_strength in %s\n', cfile);
    end
    try DATA.toPlot.ACQ.SAR(clog) = log.PARAMS.SAR; 
    catch
        fprintf(1,'\nWARNING: no field PARAMS.SAR in %s\n', cfile);
    end
    try 
        DATA.toPlot.ACQ.coils{clog} = log.PARAMS.coils;  
        % to filter data according to the coil of interest, if any
        %disp(log.PARAMS.coils);
        if ~isempty(COILfilter)
            DATA.toPlot.ACQ.coilsid(clog) = ~isempty(strfind(lower(log.PARAMS.coils),lower(COILfilter)));
        else
            DATA.toPlot.ACQ.coilsid(clog) = true;
        end
    catch
        fprintf(1,'\nWARNING: no field PARAMS.coils in %s\n', cfile);
    end
    try DATA.toPlot.ACQ.ncha(clog) = log.PARAMS.ncha; 
    catch
        fprintf(1,'\nWARNING: no field PARAMS.ncha in %s\n', cfile);
    end
    try DATA.toPlot.RES.sigma_noRF(clog) = log.RES.SNR.sigma_noRF; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.sigma_noRF in %s\n', cfile);
    end
    try DATA.toPlot.RES.nchaeff_noRF(clog) = log.RES.SNR.nchaeff_noRF; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.nchaeff_noRF in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_noRF(clog) = log.RES.SNR.snr_noRF; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_noRF in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_DIETRICH1(clog) = log.RES.SNR.snr_DIETRICH1; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_DIETRICH1 in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_DIETRICH1_neff(clog) = log.RES.SNR.snr_DIETRICH1_neff; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_DIETRICH1_neff in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_CONSTANTINI(clog) = log.RES.SNR.snr_CONSTANTINI; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_CONSTANTINI in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_CONSTANTINI_neff(clog) = log.RES.SNR.snr_CONSTANTINI_neff; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_CONSTANTINI_neff in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_FRIEDMAN(clog) = log.RES.SNR.snr_FRIEDMAN; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_FRIEDMAN in %s\n', cfile);
    end
    try DATA.toPlot.RES.snr_DIETRICH2(clog) = log.RES.SNR.snr_DIETRICH2; 
    catch
        fprintf(1,'\nWARNING: no field RES.SNR.snr_DIETRICH2 in %s\n', cfile);
    end
    try DATA.toPlot.RES.perc_drift(clog) = log.RES.FBIRN.perc_drift; 
    catch
        fprintf(1,'\nWARNING: no field RES.FBIRN.perc_drift in %s\n', cfile);
    end
    try DATA.toPlot.RES.perc_fluct(clog) = log.RES.FBIRN.perc_fluct; 
    catch
        fprintf(1,'\nWARNING: no field RES.FBIRN.perc_fluct in %s\n', cfile);
    end
    try DATA.toPlot.RES.SFNR_voxel(clog) = log.RES.FBIRN.SFNR_voxel; 
    catch
        fprintf(1,'\nWARNING: no field RES.FBIRN.SFNR_voxel in %s\n', cfile);
    end
    try DATA.toPlot.RES.rdc(clog) = log.RES.FBIRN.rdc; 
    catch
        fprintf(1,'\nWARNING: no field RES.FBIRN.rdc in %s\n', cfile);
    end
end

%==========================================================================
% save all results and parameters
%==========================================================================
% [todayy, todaym, todayd, nowh, nowm, nows] = fix(clock);
% fsave = sprintf('summary_report_%0.4d%0.2d%0.2d_%0.2d%0.2d%0.2d', todayy, todaym, todayd, nowh, nowm, nows);
currenttime = now;
savetag = REStag;
if ~isempty(REStag)
    savetag(strfind(savetag,' ')) = '_';
    savetag = [savetag '_'];
end
fsave = sprintf('summary_report_%s%s', savetag, datestr(currenttime,'yyyymmdd_HHMMSS'));
spm_jsonwrite(fullfile(PATHS.output, [fsave '.json']), DATA, struct('indent','\t'));

%==========================================================================
% save the job:
%==========================================================================
clear matlabbatch;
matlabbatch{1}.spm.tools.mriq.epi.epiqareport = job;
save(fullfile(PATHS.output, [fsave '_batch.mat']), 'matlabbatch');

%==========================================================================
% display summary figure
%==========================================================================
% [TODO: Here could be implemented some refinement and options to display only
% part of the collected data. For now, the only option is to display
% basically everything] eb20210506
    
% get SUMMARY REPORT window or create if not existing yet:
fRES = figure;
% fRES = spm_figure('GetWin','SUMMARY REPORT');
% clear it
% spm_figure('Clear',fRES);
% define display default params
figpos = [1 12 27 14]; 
set(fRES,'units','centimeters','position',figpos,'color',[1 1 1]);
def_fontname = 'Times';

% Plotting figures with multiple Y-axes is tricky. Built-in functions to
% handle multiple Y-axes have been implemented into Matlab recently, but
% not available in older versions, so using them might become an issue on
% computers running older versions of Matlab. For now, I just gather on a
% single plot all values that are spanning a similar range of values, and
% label them properly, and create new plot for different dynamic ranges...

% PLOT 0: values in the [0-400 range] non-snr-related: 
PLOT0 = [DATA.toPlot.ACQ.ref_ampl';
         DATA.toPlot.ACQ.rf_freq'*1.0e-06; % MHz
         DATA.toPlot.ACQ.ncha';
		 DATA.toPlot.RES.nchaeff_noRF';
    	 DATA.toPlot.RES.sigma_noRF';
		 %DATA.toPlot.RES.snr_noRF';
		 %DATA.toPlot.RES.snr_DIETRICH1';
		 %DATA.toPlot.RES.snr_DIETRICH1_neff';
		 %DATA.toPlot.RES.snr_CONSTANTINI';
		 %DATA.toPlot.RES.snr_CONSTANTINI_neff';
		 %DATA.toPlot.RES.snr_FRIEDMAN';
		 %DATA.toPlot.RES.snr_DIETRICH2';
		 %DATA.toPlot.RES.SFNR_voxel'
         ];
         
PLOT0legend = {'Reference amplitude [V]'
    'Frequency [MHz]'
    'Number of Rx channels'
    'Effective number of Rx channels'
    'Gaussian Noise SD'
    %'SNR (noRF)'
    %'SNR (DIETRICH1)'
    %'SNR (DIETRICH1 Neff)'
    %'SNR (CONSTANTINIDES)'
    %'SNR (CONSTANTINIDES Neff)'
    %'SNR (FRIEDMAN)'
    %'SNR (DIETRICH2)'
    %'tSNR'
    };

% PLOT 1: values in the [0-400 range], snr-related: 
PLOT1 = [%DATA.toPlot.ACQ.ref_ampl';
         %DATA.toPlot.ACQ.rf_freq'*1.0e-06; % MHz
         %DATA.toPlot.ACQ.ncha';
    	 %DATA.toPlot.RES.sigma_noRF';
		 %DATA.toPlot.RES.nchaeff_noRF';
		 DATA.toPlot.RES.snr_noRF';
		 DATA.toPlot.RES.snr_DIETRICH1';
		 %DATA.toPlot.RES.snr_DIETRICH1_neff';
		 DATA.toPlot.RES.snr_CONSTANTINI';
		 DATA.toPlot.RES.snr_CONSTANTINI_neff';
		 DATA.toPlot.RES.snr_FRIEDMAN';
		 DATA.toPlot.RES.snr_DIETRICH2';
		 DATA.toPlot.RES.SFNR_voxel'];
         
PLOT1legend = {%'Reference amplitude [V]'
    %'Frequency [MHz]'
    %'Number of Rx channels'
    %'Noise SD'
    %'Effective number of RX channels'
    'SNR (noRF)'
    'SNR (DIETRICH1)'
    %'SNR (DIETRICH1 Neff)'
    'SNR (CONSTANTINIDES)'
    'SNR (CONSTANTINIDES Neff)'
    'SNR (FRIEDMAN)'
    'SNR (DIETRICH2)'
    'tSNR'};
    
% PLOT 2: values in the [0-10 range]:
PLOT2 = [% DATA.toPlot.RES.perc_drift';
         % DATA.toPlot.RES.perc_fluct';
         DATA.toPlot.RES.rdc';
         % DATA.toPlot.ACQ.field_strength';
         % DATA.toPlot.ACQ.SAR'*100
         ]; %#ok<NBRAK>

PLOT2legend = { %'Intensity drift [%]'
    % 'Intensity fluctuation [%]'
    'RDC'
    % 'B_0 [T]'
    % 'SAR [%]'
    };

% PLOT 3: values in the [0-2 range]:
PLOT3 = [DATA.toPlot.RES.perc_drift';
         DATA.toPlot.RES.perc_fluct';
         % DATA.toPlot.RES.rdc';
         % DATA.toPlot.ACQ.field_strength';
         DATA.toPlot.ACQ.SAR'*100];

PLOT3legend = {'Intensity drift [%]'
    'Intensity fluctuation [%]'
    % 'RDC'
    % 'B_0 [T]'
    'SAR [%]'};

% Create x-axis data in numerical format:
xdatenum = datenum(DATA.toPlot.ACQ.date,'yyyymmdd');
[xdatesorted, sortidx] = sort(xdatenum);
PLOT0 = PLOT0(:,sortidx);
PLOT1 = PLOT1(:,sortidx);
PLOT2 = PLOT2(:,sortidx);
PLOT3 = PLOT3(:,sortidx);
filter = find(DATA.toPlot.ACQ.coilsid(sortidx) & DATA.toPlot.ACQ.scannerid(sortidx));

% Plot PLOT0
hax0 = subplot(4,1,1);
plot(xdatesorted(filter), PLOT0(:,filter),'marker','.');%,'linestyle','none');
datetick('x','dd/mm/yyyy','keepticks','keeplimits');
hax0.XAxis.FontSize = hax0.XAxis.FontSize*0.8;
hleg0 = legend(PLOT0legend,'Location','EastOutside');
set(hleg0,'FontSize', get(hleg0,'FontSize')*0.7);
% when run line by line, the position of axes and legend are settled when
% the refposax and refposleg are stored. It seems not to be the case
% otherwise, wrong value is ....  
pause(0.5); 
refposax = get(hax0,'Position');
%curposax = get(hax,'Position');
refposleg = get(hleg0,'Position');
%curposleg = get(hleg,'Position');
set(hax0,'Position', refposax);
set(hleg0,'Position', refposleg);


% Plot PLOT1
hax1 = subplot(4,1,2);
plot(xdatesorted(filter), PLOT1(:,filter),'marker','.');%,'linestyle','none');
datetick('x','dd/mm/yyyy','keepticks','keeplimits');
hax1.XAxis.FontSize = hax1.XAxis.FontSize*0.8;
hleg1 = legend(PLOT1legend,'Location','EastOutside');
set(hleg1,'FontSize', get(hleg1,'FontSize')*0.7);
curposax = get(hax1,'Position');
curposleg = get(hleg1,'Position');
set(hax1,'Position', [refposax(1) curposax(2) refposax(3) curposax(4)]);
set(hleg1,'Position', [refposleg(1) curposleg(2) refposleg(3) curposleg(4)]);

% Plot PLOT2
hax2 = subplot(4,1,3);
plot(xdatesorted(filter), PLOT2(:,filter),'marker','.');%,'linestyle','none');
datetick('x','dd/mm/yyyy','keepticks','keeplimits');
hax2.XAxis.FontSize = hax2.XAxis.FontSize*0.8;
hleg2 = legend(PLOT2legend,'Location','EastOutside');
set(hleg2,'FontSize', get(hleg2,'FontSize')*0.7);
curposax = get(hax2,'Position');
curposleg = get(hleg2,'Position');
set(hax2,'Position', [refposax(1) curposax(2) refposax(3) curposax(4)]);
set(hleg2,'Position', [refposleg(1) curposleg(2) refposleg(3) curposleg(4)]);

% Plot PLOT3
hax3 = subplot(4,1,4);
plot(xdatesorted(filter), PLOT3(:,filter),'marker','.');%,'linestyle','none');
datetick('x','dd/mm/yyyy','keepticks','keeplimits');
hax3.XAxis.FontSize = hax3.XAxis.FontSize*0.8;
hleg3 = legend(PLOT3legend,'Location','EastOutside');
set(hleg3,'FontSize', get(hleg3,'FontSize')*0.7);
curposax = get(hax3,'Position');
curposleg = get(hleg3,'Position');
set(hax3,'Position', [refposax(1) curposax(2) refposax(3) curposax(4)]);
set(hleg3,'Position', [refposleg(1) curposleg(2) refposleg(3) curposleg(4)]);

% general title for the figures
if ~isempty(REStag), REStag = [REStag ' - '];end
gentitl = sprintf('Summary report - %screated %s', REStag, datestr(currenttime,'dd/mm/yyyy at HH:MM:SS'));
A = annotation(fRES,'textbox',[0.0 0.95 1.0 0.05],'String',gentitl);
set(A,'LineStyle','none','HorizontalAlignment','center')
set(A,'FontName',def_fontname,'fontsize',12,'backgroundcolor',[1 1 1]);

% save summary figure
set(fRES,'PaperPositionMode','auto')
print(fRES,'-dpng',fullfile(PATHS.output, [fsave '.png']));
% saveas(gcf,fullfile(PATHS.output, [fsave '.fig']));

% define outputs - TO BE IMPLEMENTED
out = [];

end