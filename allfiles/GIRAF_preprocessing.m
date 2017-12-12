function [X0,b0,N1,N2,N3,res] = GIRAF_preprocessing(data_reconstructed)

% GIRAF reconstructs on a grid on odd length.
% This indexing sets up the grid
[N1,N2,N3] = size(data_reconstructed);

res=[N1 N2];
indx = [0:((res(1)/2)-1), -(res(1)/2):-1];
indy = [0:((res(2)/2)-1), -(res(2)/2):-1];
[kx,ky] = meshgrid(indx,indy);kx=kx';ky=ky';
k(1,:) = kx(:);
k(2,:) = ky(:);

% GIRAF reconstructs the 2D time domain data, this step
% converts 2D spectral data to 2D time

b0 = ifft2(data_reconstructed);

%reindex input data
ind_trim = find((abs(k(1,:)) <= (res(1)-1)/2 ) & (abs(k(2,:)) <= (res(2)-1)/2));
b0 = reshape(b0(ind_trim),res-[1 1]);

% Fully sampled data on the new grid
X0 = fft2(b0);

