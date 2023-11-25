%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================


function[] = run_CellSegm(file_dir, percentage, pixel_size)

a = file_dir;
p = str2num(percentage);
z = str2num(pixel_size);
minCell = 10*1000/z*p/100;
maxCell = 40*1000/z*p/100;
addpath(a);

b = pwd;
[parentDir, ~, ~] = fileparts(b);
addpath(parentDir);


% load the data
imnucl = double(read(Tiff('nucleus.tif')));
imsegm = double(read(Tiff('membrane.tif')));


% Smoothing
prm.smooothim.method = 'dirced';

% No ridge filtering
prm.filterridges = 0;

% threshold for nucleus markers
prm.getminima.nucleus.segmct.thrs.th = 0.50;

% edge enhancing diffusion with a suitable threshold
prm.getminima.nucleus.segmct.smoothim.method = 'eed';
prm.getminima.nucleus.segmct.smoothim.eed.kappa = 0.05;

% method for markers
prm.getminima.method = 'nucleus';

% Subtract the nucleus channel from the surface staining to reduce the
% cross talk effect.
imsegm1 = imsegm;
filt = fspecial('gaussian',3,2);
imsegm = imfilter(imsegm1,filt) - imfilter(imnucl,filt);

[cellbw,wat,imsegmout,minima,minimacell,info] = ...
    cellsegm.segmsurf(imsegm,minCell,maxCell,'imnucleus',imnucl,'prm',prm);



imwrite(uint8(cellbw), strcat(a, filesep, 'mask_CellSegm.png'));
end
