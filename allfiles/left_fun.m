function [X] = left_fun(x,overres,mu,wind_fun,lambda,mask_pad)

x = reshape(x,overres);

X1= ifft2(conj(wind_fun).*mu.*wind_fun.*fft2(x));


X2=mask_pad.*x;

X=X1+lambda*X2;
X=X(:);