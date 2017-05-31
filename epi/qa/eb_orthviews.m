function A = eb_orthviews(fnam, i, j, k)

% returns a diplayable matrix showing 3 orthogonal views from a nifti image
% whose file name is fnam.
% Evelyne Balteau - August 2016

V = spm_vol(fnam);
Y = spm_read_vols(V);
dim = V.dim;

A = zeros(dim(2)+dim(3),dim(1)+dim(2));
A(1:dim(3),1:dim(1)) = squeeze(Y(end:-1:1,j,end:-1:1))';
A(dim(3)+(1:dim(2)),1:dim(1)) = squeeze(Y(end:-1:1,end:-1:1,k))';
A(1:dim(3),dim(1)+(1:dim(2))) = squeeze(Y(i,end:-1:1,end:-1:1))';
