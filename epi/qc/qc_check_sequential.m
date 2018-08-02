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

%=========================================================================%
% This file is part of the MRI quality toolbox.
% Copyright (C) 2013-2018 - Cyclotron Research Centre
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

function out = qc_check_sequential(job)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: sequential display of EPI series to detect spikes, artefacts
% and head movements. Suitable for fMRI time series, resting-state,
% diffusion images.
%--------------------------------------------------------------------------
% Warning and disclaimer: This software is for research use only. 
% Do not use it for clinical or diagnostic purposes.
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - 2013-2018
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global handle P;

disp('----- Loading EPI series -----');

P = char(job.EPIimages);
nmeas = size(P,1);
Y = spm_read_vols(spm_vol(P(1:min(nmeas,10),:)));
maxmag = max(Y(:));
clear Y;

disp('----- Prepare display -----');

S0 = get(0,'ScreenSize');   %-Current screen size [1 1 width height]
WS = S0(1,3:4);             %-Window size [width height]          
WSF = 0.8;                  %-Window scaling factors
MAINSIZE = min(WS)*WSF;
F = figure('color',[1 1 1],'position',[50 50 MAINSIZE MAINSIZE]);
M = subplot(1,1,1);
set(M,'position',[0.06 0.06 0.92 0.92]);

handle.slider = uicontrol(F,'Style','Slider',...
    'BackgroundColor',[0.9 0.9 0.88],...
    'Max',nmeas,'Min',1,...
    'Tag','slider',...
    'Value',1,...    
    'SliderStep',[1 1]./nmeas,...
    'CallBack','manage_UI(''refresh'')',...
    'Position',MAINSIZE*[0.02 0.06 0.02 0.92]);
handle.slidercentre = uicontrol(F,'Style','Slider',...
    'BackgroundColor',[0.9 0.9 0.88],...
    'Max',maxmag,'Min',0,...
    'Tag','slider',...
    'TooltipString','Centre of the display window',...
    'Value',30,...    
    'SliderStep',[1 1]./maxmag,...
    'CallBack','manage_UI(''refresh'')',...
    'Position',MAINSIZE*[0.60 0.02 0.18 0.02]);
handle.sliderwidth = uicontrol(F,'Style','Slider',...
    'BackgroundColor',[0.9 0.9 0.88],...
    'Max',maxmag,'Min',1,...
    'Tag','slider',...
    'TooltipString','Width of the display window',...
    'Value',40,...    
    'SliderStep',[1 1]./(maxmag-1),...
    'CallBack','manage_UI(''refresh'')',...
    'Position',MAINSIZE*[0.80 0.02 0.18 0.02]);

manage_UI('refresh');

disp('----- Ready for inspection -----');

out = job;
end