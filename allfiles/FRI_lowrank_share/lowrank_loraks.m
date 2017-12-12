%% Weighted LORAKS regularization
%Solves min_x 0.5*||Ax-b||_2^2 + lambda*phi(Qx)
%where phi(Y) := min_{Z: rank Z <= r} ||Y-Z||_2^2 
%(i.e. the l2-norm of the residual of best rank r approximation),
%A is a Fourier undersampling mask, and
%Qx is the block Toeplitz matrix built from k-space data x
%Uses a LORAKS-style non-convex hard thresholding of singular values
%% Load data
clear all;
addpath('./algorithms','./data','./etc');
load MR_SL_convert;
%load MR_brain_convert;
%load MR_SL_linear;

%% Set model order
res = [63,63];  %output resolution
filter_siz = [15,15]; %filter dimensions

ind_trim = find((abs(k(1,:)) <= (res(2)-1)/2 ) & (abs(k(2,:)) <= (res(1)-1)/2));
m = reshape(m(ind_trim),res);
k = k(:,ind_trim);
%% Generate fourier undersampling operator
usf = 0.5; %undersampling factor
pdf = genPDF(res,3,usf); %density compensation
mask = im2double(genSampling(pdf,2,5));
mask(1) = 1; %need DC component

figure(1); imshow(fftshift(mask));
ind_samples = find(mask~=0);
[A,At] = defAAt(ind_samples,res);
b0 = A(m);
sig = 0; %std deviation 
b = b0+sig*(randn(size(b0))+ 1i*randn(size(b0)));
fprintf('Sample SNR = %5.2f dB\n',20*log10(norm(A(m))/norm(A(m)-b0)));
%% Build Toeplitz operators
dz{1} = reshape((-1i*pi*(k(1,:))).',res)/(res(1)/4);
dz{2} = reshape((-1i*pi*(k(2,:))).',res)/(res(2)/4);
clear Q;
Q = ForwardModelLowRankThin(dz,res,filter_siz,[1,1]);
QtQ = Q'*(Q*ones(res));
disp(size(Q*zeros(res))) %size of QX matrix
%% run algorithm
lambda = 10;
r = 200; %rank threshold

Atb = At(b);
AtA = mask;
x = Atb;

Qx = Q*x;
Z = zeros(size(Qx)); %aux variable

cost = [];
constraint = [];
iter = 100;
for i=1:iter
    Qx = Q*x;
    [U, S, V] = svd(Qx,'econ');
    s = diag(S);
    figure(4); plot(s); title('sing vals');
    ranknorm = 0.5*norm(s((r+1):end))^2;    
    s((r+1):end) = 0;
    Z = U*diag(s)*V';
    
    quad = AtA + lambda*QtQ;
    quadinv = 1./quad;      
    x = quadinv.*(Atb + lambda*(Q'*Z));
    
    diff = A(x)-b;
    thiscost = 0.5*norm(diff(:)).^2 + lambda*ranknorm;
    cost = [cost,thiscost];        

    figure(11); imshow(abs(ifft2(x)),[]);
    figure(12); plot(cost); title('cost');    
    figure(13); imagesc(fftshift(log(1+abs(x)))); colorbar;
end
toc
%% Plot Original Data
x0 = reshape(m,res);
figure(2); imshow(abs(ifft2(x0)),[]); title('original');
figure(3); imagesc(fftshift(log(1+abs(x0)))); colorbar;

