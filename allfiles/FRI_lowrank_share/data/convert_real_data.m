load thigh;

data = I(20:219,267:466);
imshow(abs(data),[]);
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

save('MR_realbrain_singlecoil.mat','m','k');