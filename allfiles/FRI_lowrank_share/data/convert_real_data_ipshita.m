load('/nfs/s-iibi52/projects/cbig/CBIG_DATA/SHARE/Ipshita_share/Yasir_data/AX_T2_FID138092.mat');
%AX_T1_FID138096.mat
%AX_T1_FID138100.mat
%AX_T2_FID138092.mat
%AX_T2_FID138098.mat
%%
slice = squeeze(final(61,:,:));
imagesc(abs(slice)); colormap('gray'); colorbar;
%%
sliceedit = slice;
sliceedit(:,1:200) = 0;
sliceedit(:,700:end) = 0;
slice = sliceedit;
%%
fslice = fft2(slice);
fslice = fslice(:,1:2:end);
down_slice = ifft2(fslice);
%%
data = fftshift(down_slice,2);
%data = circshift(data,[0,70]);
%data = flipud(data);
data = data/max(abs(data(:)));
%%
imagesc(abs(data)); colormap('gray'); colorbar;
%%
res = size(data);
m = fft2(data);
indx = [0:((res(1)/2)-1), -(res(1)/2):-1];
indy = [0:((res(2)/2)-1), -(res(2)/2):-1];
[kx,ky] = meshgrid(indx,indy);
clear k;
k(1,:) = kx(:);
k(2,:) = ky(:);

save('MR_realbrain_cor1.mat','m','k');