function raw = readmeas(name,nread,necho,ncoil,ccoil,MDH_END)
%---------------------------------------------------------------------
% raw = readmeas(name,nread,ncoil,ccoil)
% extracts raw data from meas.out or meas.dat files
% version for phased array data
% - name is the name of the file
% - nread is the number of complexe values per echo 
%   (including oversampling)
% - necho is the number of lines in one volume
% - ncoil is the total number of coils
% - ccoil is the current coil from which data are extracted
% - raw contains the raw data for ccoil as stream of double
%   (format  real/imaginary/real/imaginary etc.)
%   raw data is in chronological order, without echo reversals
%---------------------------------------------------------------------
% tips about the meas.out/meas.dat file structure (VA35):
% - leading dummy global header of 32 bytes (meas.out) or more (meas.dat,
%   read the first uint32 of the file to find out where the first 
%   block starts)
% - blocks of measurement data (one echo acquired by one channel=one block) 
%   arranged with a 128 bytes long MDH followed by measurement itself. 
%   Each block is therefore 128 + nread*2*4 bytes long.
% - data represented as 'float32' (4 bytes) real, imag, real, 
%   imag... (factor 2 above for each data point)
% - dummy scan at the end of each measurement and carrying 
%   MDH_ACQEND flag (256 bytes per channel)
%---------------------------------------------------------------------
% Adapted from code by Ralf Deichmann
% Last modif: July 2008 - Evelyne Balteau
% Last modif: October 2013 - Evelyne Balteau (VD13)
%---------------------------------------------------------------------

% open file
ftp=fopen(name,'r');

% some cst values
ACQEND_SIZE = 256;
MDH_SIZE = 128;
FLOAT32_SIZE = 4;
MDH_COIL = 8*FLOAT32_SIZE; % 32
MDH_ECHO = 48*FLOAT32_SIZE; % 192
MDH_END = MDH_END*FLOAT32_SIZE; % 448

% % memory of a hacker to work out what should be MDH_END
% % go to 1024 floats before the end of file
% fseek(ftp,-4096,1)
% precis = sprintf('%d*float32=>double',nread*2);
% skipn = 0;
% % read the last 1024 floats
% rawend = fread(ftp,nread*2*necho,precis,skipn);
% % display to visually determine how many floats after the last echo
% figure;plot(abs(rawend),'marker','.');

% get file size
% fseek(ftp,0,1); % go to eof
% size=ftell(ftp);
% fseek(ftp,0,-1); % go to bof

% lets reach the beginning of the first block starting from eof:
fseek(ftp, - MDH_END - necho*(MDH_ECHO + ncoil*(MDH_COIL + nread*2*FLOAT32_SIZE)) + MDH_ECHO,1);

% then reaches the first data point of the block corresponding to ccoil:
fseek(ftp,(FLOAT32_SIZE*nread*2 + MDH_COIL)*(ccoil-1) + MDH_COIL,0);

% now reads the pelines*nslices = nread*necho/ncoil data 
% corresponding to ccoil. This is done by:
% - reading a total of nread*2*necho float
% - in blocks of nread*2 (complex), 
% - read as float and saved as double
% - skip MDH_SIZE*ncoil + FLOAT32_SIZE*nread*(ncoil-1)
%   between blocks.
precis = sprintf('%d*float32=>double',nread*2);
skipn = MDH_ECHO + (MDH_COIL + FLOAT32_SIZE*nread*2)*(ncoil-1) + MDH_COIL;
raw = fread(ftp,nread*2*necho,precis,skipn);

% close file
fclose(ftp);
