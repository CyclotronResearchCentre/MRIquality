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

function Y = get_mosaic(Yin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: returns a 2D mosaic image from a 3D volume
%--------------------------------------------------------------------------
% Warning and disclaimer: This software is for research use only. 
% Do not use it for clinical or diagnostic purposes.
%--------------------------------------------------------------------------
% Written by Evelyne Balteau - 2013
% Cyclotron Research Centre, University of Liege, Belgium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dims = size(Yin);
nsqrt = ceil(sqrt(dims(3)));
% dimensions need to be inverted to have the same "look" as at the scanner
% (NB: not sure about the left-right orientation, might be inverted but I
% cannot be bothered to check!)
Y = zeros(nsqrt*dims(2), nsqrt*dims(1));
for i = 1:dims(3)
    %disp([floor((i-1)/nsqrt)*dims(2) (mod(i-1,nsqrt))*dims(1)]);
    Y(floor((i-1)/nsqrt)*dims(2)+[1:dims(2)],(mod(i-1,nsqrt))*dims(1)+[1:dims(1)]) = squeeze(Yin(end:-1:1,end:-1:1,i))';
end
end