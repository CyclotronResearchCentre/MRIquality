function psii = eb_calc_psi(noise_array)
% calculate the noise matrix psi from a noise measurements
% (acquisition without RF)...
% follows Weiger et al 2001 MRM 45:495-504

[resol, npe, ncoil, nslices] = size(noise_array);
psii = eye(ncoil);
for i=1:ncoil
    for j=1:ncoil
        b = noise_array(:,:,i,:);
        c = conj(noise_array(:,:,j,:));
        psii(i,j) = mean(b(:).*c(:));
    end
end

% normalize (sans importance pour les calculs vu que psi et inv(psi) interviennent...
% c'est juste plus joli pour le isplay et pour la comparaison avec une matrice unite)
psii = psii*ncoil./(sum(sum(abs(psii))));

% figure('Position',[50 50 600 600]);
% subplot(1,1,1);
% set(gca,'position',[0.1 0.1 0.8 0.8])
% imshow(abs(psii),[]);colorbar;
% f = getframe(gcf);
% imwrite(f.cdata,['psi.bmp'],'bmp');
% 
% figure('Position',[50 50 600 600]);
% subplot(1,1,1);
% set(gca,'position',[0.1 0.1 0.8 0.8])
% imshow(abs(psii),[0 1]);colorbar;
% f = getframe(gcf);
% imwrite(f.cdata,['psi_scale_1.bmp'],'bmp');
