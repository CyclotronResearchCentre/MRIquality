function data=make_sense_nosmooth(varargin)

data=varargin{1};
vox=varargin{2};

% Make sensitivity maps from image volume
% FORMAT data=make_sense(data,V);
%
% Data should be reconstructed and complex data with dimensions:
% nm,np,ncoils,ns
% The function does the following steps
% 1) Make a sqrt sum of squares voulme from all channels
% 2) Create a mask of the above using a threshold of 80% of mean 
dim=size(data);
ncoils=size(data,3);
vox;

ssqvol=squeeze(sqrt(sum(abs(data).^2,3)));
meanthresh=0.8; 
maskvol=ssqvol>(mean(ssqvol(:)*meanthresh));

sensedata=zeros(size(data));
rsmoothvol=zeros(size(squeeze(data(:,:,1,:))));
ismoothvol=zeros(size(squeeze(data(:,:,1,:))));
nvox=12./vox;
for nc=1:ncoils
    sensedata(:,:,nc,:)=(squeeze(data(:,:,nc,:)).*maskvol)./ssqvol;    
end
data=sensedata;
