res = [256, 256];
%I = imread('~/Dropbox/Presentations/sampta2015/bird2.png');
I=phantom(256);
I = I(:,:,1)/255;
ind = mod((1:256)-(129),size(I,1))+1;
tmp = fft2(I);
data = ifft2(ifftshift(tmp(ind,ind)));

%%
data = data/max(abs(data(:)));
res = size(data);
m = fft2(data);
indx = [0:((res(1)/2)-1), -(res(1)/2):-1];
indy = [0:((res(2)/2)-1), -(res(2)/2):-1];
[kx,ky] = meshgrid(indx,indy);
clear k;
k(1,:) = kx(:);
k(2,:) = ky(:);

save('MR_sampta2.mat','m','k');
%%
m = zeros(256);
m(ind_samples) = X(:).';
res = size(m);
m0 = m;
K = 10;
knew = load('data/MR_isbi_logo.mat','k');
k = knew.k;