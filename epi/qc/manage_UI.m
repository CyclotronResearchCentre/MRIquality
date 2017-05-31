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

function manage_UI(narg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: user interface displaying the series of mosaic images with
% adjustable windowing (centre and width) and scrolling back and forth in
% the series.
%--------------------------------------------------------------------------
% Warning and disclaimer: This software is for research use only. 
% Do not use it for clinical or diagnostic purposes.
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - 2013
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global handle P;
switch narg
    case 'refresh'
        cmeas = round(get(handle.slider,'Value'));
        centre = round(get(handle.slidercentre,'Value'));
        width = round(get(handle.sliderwidth,'Value'));
        contrast = centre + width*[-1 1];
        Ydis = get_mosaic(spm_read_vols(spm_vol(P(cmeas,:))));
        imshow(Ydis,contrast);
        dims = size(Ydis);
        [cpath cnam] = fileparts(P(cmeas,:));
        % text(dims(1)*0.03,dims(2)*0.97,[num2str(cmeas) ' [' num2str(contrast(1)) ' ' num2str(contrast(2)) ']'],'color',[1 1 0],'fontsize',15,'fontweight','bold');
        % text(dims(1)*0.03,dims(2)*0.97,['Volume #' num2str(cmeas) ' [C:' num2str(centre) '; W:' num2str(width) ']'],'color',[1 1 0],'fontsize',15,'fontweight','bold');
        text(dims(1)*0.03,dims(2)*0.97,[cnam ' [C:' num2str(centre) '; W:' num2str(width) ']'],'color',[1 0.8 0],'fontsize',15,'fontweight','bold');
        drawnow;
    case 'quit'
        close gcf;
end
end