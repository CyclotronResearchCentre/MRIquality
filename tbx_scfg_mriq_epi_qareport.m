function epiqareport = tbx_scfg_mriq_epi_qareport
% (Sub)configuration file for the MRIquality toolbox, partim EPI/QA Report.
% The report generates a series of figures summarizing the results
% collected over time from QA runs, in order to visualize any evolution,
% drift, outliers, etc...
%
% Warning and disclaimer: This software is for research use only.
% Do not use it for clinical or diagnostic purposes.
%
%==========================================================================
% Written by Evelyne Balteau, Cyclotron Research Centre, May 2021
%==========================================================================

%==========================================================================
% QA log files (JSON)
%==========================================================================
RESlist           = cfg_files;
RESlist.tag       = 'RESlist';
RESlist.name      = 'QA log files';
RESlist.help      = {'Select a list of QA log files.' ...
    ['The files have typically the following format (for example) ' ...
    '''Prisma_20190131_stud0042_ser0002.json'' ' ... 
    'and result from data acquired on the same scanner ' ...
    'and same EPI sequence. ']};
RESlist.dir       = mriq_get_defaults('path_input');
RESlist.filter    = '^*\.json$';
RESlist.ufilter   = '.*';
RESlist.num       = [0 Inf];
RESlist.val       = {''};

%==========================================================================
% Filters for scanner & coil (mixture usually not useful)
%==========================================================================
RESfilter_scanner         = cfg_entry;
RESfilter_scanner.tag     = 'RESfilter_scanner';
RESfilter_scanner.name    = 'Scanner';
RESfilter_scanner.help    = {['Enter a filter to select results from a ' ...
                              'specific scanner (e.g. Prisma, Allegra, ' ...
                              'Terra, ...). Leave empty if all scanners ' ...
                              '(if several) to be considered. ']};
RESfilter_scanner.strtype = 's';
RESfilter_scanner.num     = [0 Inf];
RESfilter_scanner.val     = {''};

RESfilter_coil         = cfg_entry;
RESfilter_coil.tag     = 'RESfilter_coil';
RESfilter_coil.name    = 'Coil';
RESfilter_coil.help    = {['Enter a filter to select results from a ' ...
                              'specific coil (e.g. HE for 20-channel HE/NE, ' ...
                              'HC for 64-channel HC/NC coil, etc). ' ...
                              'Leave empty if all coils (if several) ' ...
                              'to be considered. Check field ' ...
                              'log.PARAMS.coils in mriq_run_epi_qareport.m ' ...
                              'for a clue about the names used for the coils...']};
RESfilter_coil.strtype = 's';
RESfilter_coil.num     = [0 Inf];
RESfilter_coil.val     = {''};

REStag         = cfg_entry;
REStag.tag     = 'REStag';
REStag.name    = 'Tag';
REStag.help    = {['Enter a more explicit tag for the current result, ' ...
                              'agreeing with the selected coil and scanner ' ...
                              '(e.g. Prisma 64-channel Head-Neck coil). ' ...
                              'NOTE: avoid characters that cannot be used within a file name!']};
REStag.strtype = 's';
REStag.num     = [0 Inf];
REStag.val     = {''};

%==========================================================================
% Parameters of interest to be displayed
%==========================================================================
% % PARoI              = cfg_branch;
% % PARoI.tag       = 'PARoI';
% % PARoI.name      = 'Parameters';
% % PARoI.help      = {'Select the parameters of interest to be displayed.'};
% % PARoI.val       = {PARsnr};

% % PARacq: "PARAMS"
% %             "ref_ampl": 222.432723999,
% %             "rf_freq": 123254542,
% %             "SAR": 0.02550401926618,
% %     		"ncha": 44,
% % 		
% % PARsnr:	"RES": 
% % 			"SNR": 
% %             	"sigma_noRF": 35.81144175649234,
% %                 "nchaeff_noRF": 40.15373166349986,
% %                 "snr_noRF": 336.6147228758923,
% % 			
% % PARpercdrift:
% %         "RES":
% %             "FBIRN":
% %                 "perc_drift": 1.9890024383606,
% % 			"perc_fluct": 0.0516477692122591,
% % 			"SFNR_voxel": 175.8265856552203,
% % 			"rdc": 11.25787672649888
% %              PARpercdrift PARpercfluct PARrdc

%==========================================================================
% PARsnr
%==========================================================================
% % PARsnr        = cfg_menu;
% % PARsnr.tag    = 'PARsnr';
% % PARsnr.name   = 'SNR & tSNR';
% % PARsnr.help   = {'ON/OFF:', ...
% %     ['SNR and tSNR values are displayed according to various derivation methods. ']};
% % PARsnr.labels = {
% %                'ON'
% %                'OFF'}';
% % PARsnr.values = {1 0};
% % PARsnr.val    = {1};



%==========================================================================
% Tools for EPI quality report
%==========================================================================
epiqareport         = cfg_exbranch;
epiqareport.tag     = 'epiqareport';
epiqareport.name    = 'Quality report';
epiqareport.val     = {RESlist RESfilter_scanner RESfilter_coil REStag};
epiqareport.help    = {'To display the QA results and their evolution over time.'};
epiqareport.prog    = @mriq_run_epi_qareport;
epiqareport.vout    = @vout_epi_qareport; 

end


%% VOUT & OTHER SUBFUNCTIONS
% ========================================================================
% The RUN function:
% - out = hmri_run_create(job)
% is defined separately.
%_______________________________________________________________________

function dep = vout_epi_qareport(job)
% This depends on job contents, which may not be present when virtual
% outputs are calculated.

% All output from a given dataset have a common filename formed as follows:
% <Scanner>_<DateYYYYMMDD>_stud0000_ser0000 (see the asterisk * below). The
% file name is appended with suffix and the file extension varies according
% to the output. 
% For each run of the QA batch, the following output are produced:  
% - *.json file with the results
% - *.txt file with same information but a bit easier to read
% - *.png screenshot of the general results window
% - *_batch.mat job saved for the current run
% - *_MEAN.nii mean volume
% - *_NOISE_DISTRIB.png screenshot of the noise distribution figure
% - *_SD.nii standard deviation volume
% - *_SNR_noRF.nii SNR volume (calculated with noise image)
% - *_tSNR.nii tSNR volume

dep = [];
% k=1;
% cdep(1,5*numel(job.subj)) = cfg_dep;
% for i=1:numel(job.subj)
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('R1_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','R1','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('R2s_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','R2s','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('MT_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','MT','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('A_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','A','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('T1w_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','T1w','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('MTw_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','MTw','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
%     cdep(k)            = cfg_dep;
%     cdep(k).sname      = sprintf('PDw_subj%d',i);
%     cdep(k).src_output = substruct('.','subj','()',{i},'.','PDw','()',{':'});
%     cdep(k).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
%     
%     k=k+1;
%     
% end
% dep = cdep;
    
end