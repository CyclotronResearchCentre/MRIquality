function eb_save_NIFTI(vol,hdr,fname)
% version for 2D acquisition

% Orientation information
%-----------------------------------------------------------------------
% Axial Analyze voxel co-ordinate system:
% x increases     right to left
% y increases posterior to anterior
% z increases  inferior to superior

% DICOM patient co-ordinate system:
% x increases     right to left
% y increases  anterior to posterior
% z increases  inferior to superior

% T&T co-ordinate system:
% x increases      left to right
% y increases posterior to anterior
% z increases  inferior to superior

% Retrieve parameters from header 
SpacingBetweenSlices = hdr.asc{3}.sSliceArray.asSlice(1).dThickness * (1+hdr.asc{3}.sGroupArray.asGroup(1).dDistFact);
PixelSpacing = hdr.asc{3}.sSliceArray.asSlice(1).dReadoutFOV/hdr.Meas{2}.DICOM.lBaseResolution*[1 1];

% Orientation
tmp = hdr.asc{3}.sSliceArray.asSlice(1).sNormal;
if isfield(tmp,'dCor'); NormCor = tmp.dCor; else NormCor = 0; end
if isfield(tmp,'dSag'); NormSag = tmp.dSag; else NormSag = 0; end
if isfield(tmp,'dTra'); NormTra = tmp.dTra; else NormTra = 0; end

ImageOrientationPatient = [1     0     0     0     NormTra  -NormCor];
ImageOrientationPatient = reshape(ImageOrientationPatient,[3 2]);
ImageOrientationPatient(:,3) = null(ImageOrientationPatient');
if det(ImageOrientationPatient)<0, ImageOrientationPatient(:,3) = -ImageOrientationPatient(:,3); end;

% Position of the central point of the volume
if isfield(hdr.asc{3}.sSliceArray.asSlice(1),'sPosition') % if isocentre, sPosition is not specified
    tmpfirst = hdr.asc{3}.sSliceArray.asSlice(1).sPosition;
else
    tmpfirst = [];
end
if isfield(tmpfirst,'dCor'); PosCor = tmpfirst.dCor; else PosCor = 0; end
if isfield(tmpfirst,'dSag'); PosSag = tmpfirst.dSag; else PosSag = 0; end
if isfield(tmpfirst,'dTra'); PosTra = tmpfirst.dTra; else PosTra = 0; end
if isfield(hdr.asc{3}.sSliceArray.asSlice(end),'sPosition') % if isocentre, sPosition is not specified
    tmplast = hdr.asc{3}.sSliceArray.asSlice(end).sPosition;
else
    tmplast = [];
end
if isfield(tmplast,'dCor'); PosCor = (PosCor+tmplast.dCor)/2; else PosCor = PosCor/2; end
if isfield(tmplast,'dSag'); PosSag = (PosSag+tmplast.dSag)/2; else PosSag = PosSag/2; end
if isfield(tmplast,'dTra'); PosTra = (PosTra+tmplast.dTra)/2; else PosTra = PosTra/2; end

% The ImagePositionPatient is the position of the corner of the first voxel
% of the first slice -> starts from centre point and add diagonal distance
% to corner of first slice (middle of slice thickness) * rotation given by the orientation
ImagePositionPatient = [PosSag PosCor PosTra] ...
    - (ImageOrientationPatient* [hdr.asc{3}.sSliceArray.asSlice(1).dReadoutFOV/2; 
                                     hdr.asc{3}.sSliceArray.asSlice(1).dPhaseFOV/2;
                                     ((hdr.asc{3}.sGroupArray.asGroup(1).nSize-1)/2)*SpacingBetweenSlices])';

% Orientation & spatial parameters for NIFTII 
dim    = size(vol);
analyze_to_dicom = [diag([1 -1 1]) [0 (dim(2)-1) 0]'; 0 0 0 1]*[eye(4,3) [-1 -1 -1 1]'];
vox    = [PixelSpacing SpacingBetweenSlices];
pos    = ImagePositionPatient';
orient = ImageOrientationPatient;

% The image position vector is not correct. In dicom this vector points to
% the upper left corner of the image. Perhaps it is unlucky that this is
% calculated in the syngo software from the vector pointing to the center of
% the slice (keep in mind: upper left slice) with the enlarged FoV.
dicom_to_patient = [orient*diag(vox) pos ; 0 0 0 1];
patient_to_tal   = diag([-1 -1 1 1]);
mat = patient_to_tal*dicom_to_patient*analyze_to_dicom;

% Save...
V = struct('fname',   fname,...
    'dim',     dim,...
    'dt',      [spm_type('int16') spm_platform('bigend')],...
    'mat',     mat,...
    'descrip', fname);

spm_write_vol(V,vol);

return;