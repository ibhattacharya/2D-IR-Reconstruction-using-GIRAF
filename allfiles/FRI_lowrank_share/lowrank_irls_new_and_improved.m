%% IRLS Nuclear Norm Minimization
%Solves min_X (1/p)*||T(X)||_p^p + (lambda/2)*||AX-b||_2^2
%where A is a Fourier undersampling mask, and
%QX is the block Toeplitz matrix built from k-space data X
%Uses a two-stage IRLS MM approach, which alternates between updating
%an annihilating filter/edge mask and a least squares annihilation
%% Load data
clear all;
addpath('./algorithms','./data','./etc');
load MR_SL_convert;
%load MR_brain_convert;
%load MR_realbrain_sag;

%% Set model order
res = [255,255];  %output resolution (use odd nums)
filter_siz = [33,33]; %filter dimensions (use odd nums)
filter_siz2 = 2*filter_siz - [1,1]; %squared filter dimensions
overres = res + 2*filter_siz; %over-resolved reconstruction grid

%reindex input data
ind_trim = find((abs(k(1,:)) <= (res(2)-1)/2 ) & (abs(k(2,:)) <= (res(1)-1)/2));
m = reshape(m(ind_trim),res);
X0 = ifft2(m);
figure(31); imshow(abs(X0),[]);
figure(32); imshow(fftshift(abs(fft2(X0))).^(1/4),[]);
%define k-space indices for over-resolved grid
clear k;
k = zeros(2,prod(overres));
kx = ifftshift(repmat(-(overres(2)-1)/2:(overres(2)-1)/2,overres(1),1));
ky = ifftshift(repmat((-(overres(1)-1)/2:(overres(1)-1)/2).',1,overres(2)));
k(1,:) = kx(:);
k(2,:) = ky(:);
%
ind_full = find((abs(k(1,:)) <= (res(2)-1)/2) & (abs(k(2,:)) <= (res(1)-1)/2));
ind_filter = find((abs(k(1,:)) <= (filter_siz(2)-1)/2 ) & ( abs(k(2,:)) <= (filter_siz(1)-1)/2));
ind_filter2 = find((abs(k(1,:)) <= (filter_siz(2)-1) & (abs(k(2,:)) <= (filter_siz(1)-1))));
%% Generate fourier undersampling operator
usf = 0.33; %undersampling factor
%pdf = genPDF(res,6,usf); %density compensation
%mask0 = genSampling(pdf,2,5);
randind = randperm(prod(res)); %uniform random samples
mask0 = zeros(res);
numind = floor(usf*prod(res));
mask0(randind(1:numind)) = 1;
mask0(1) = 1; %need DC component
figure(1); imshow(fftshift(mask0));
ind_samples = find(mask0~=0);
[A,At] = defAAt(ind_samples,res);
b0 = A(m);
sig = 50; %added noise std deviation
b = b0+sig*(randn(size(b0))+ 1i*randn(size(b0)));

noiseSNR = 20*log10(norm(b)/norm(b-b0));
fprintf('Sample SNR = %5.2f dB\n',noiseSNR);
x1 = zeros(res);
x1(ind_samples) = b;
x1 = reshape(x1,res);
X_IFFT = ifft2(x1);
figure; imshow(abs(X_IFFT),[]);
stats.IFFT = imstats(X0,X_IFFT);
%% init variables and operators
%LORAKS spatial support penalty
% clear dz;
% dz(:,:,1) = ones(overres);

%1st order derivatives
clear dz;
dz(:,:,1) = reshape((-1i*2*pi*(k(1,:))).',overres)/(res(1));
dz(:,:,2) = reshape((-1i*2*pi*(k(2,:))).',overres)/(res(1));
 
%2nd order derivatives
% clear dz;
% dz(:,:,1) = reshape((-4*pi^2*(k(1,:).^2)).',overres)/(res(1)^2);
% dz(:,:,2) = reshape((-4*pi^2*(k(2,:).^2)).',overres)/(res(2)^2);
% dz(:,:,3) = reshape((-8*pi^2*(k(1,:).*k(2,:))).',overres)/(res(1)*res(2));

M = @(z) repmat(z,[1,1,size(dz,3)]).*dz;
Mt = @(Z) sum(Z.*conj(dz),3);
MtMmask = Mt(M(ones(overres)));

mask_pad = zeros(overres);
mask_pad(ind_full) = mask0;
ind_samples_pad = find(mask_pad~=0);
AtA = @(x) proj(x,ind_samples_pad);

Atb = At(b);
Atb_pad = zeros(overres);
Atb_pad(ind_full) = Atb;

x = Atb;
x_pad = zeros(overres);
x_pad(ind_full) = x;
X_old = X_IFFT;
%% run alg.
lambda0 = 1e-5; %regularization parameter (high value enforces low-rank constraint)
lambda = 1/lambda0;
eps = 1;         %epsilson initialization; typically 0.01*max(s)
eta = 1.3;       %epsilon update is eps = eps/eta; (1.3--1.5 work well)
epsmin = 1e-5;   %minimum possible epsilon value
 
p = 0;           %Shatten-p norm value
q = 1-(p/2);      

y = lambda*Atb_pad;
y = y(:);

cost = [];
epsrecord = [];
SNR = [];                                                           
relerr = [];

iter = 10; %typically 5-15 iterations. 
for i=1:iter
%Stage 1: Compute annihilating mask
gradx = x_pad;%M(x_pad);
gradx_ifft = ifft2(gradx);
sos = fft2(sum(conj(gradx_ifft).*gradx_ifft,3));
sos2 = fftshift(reshape(sos(ind_filter2),filter_siz2));
R = im2colstep(real(sos2),filter_siz,[1,1]) + 1i*im2colstep(imag(sos2),filter_siz,[1,1]);
R = rot90(R,-1);
[U,S] = eig(R+eps*eye(size(R)));
s = diag(S); %note: s = sing. values squared.
figure(11); plot(0.5*log10(s)); title('sing vals, log scale'); drawnow;

%calculate cost
if p == 0
    shatten = 0.5*sum(log(s-eps));
else
    shatten = (1/p)*sum((s-eps).^(p/2));
end
diff = A(x)-b;
thiscost = lambda0*shatten + 0.5*norm(diff(:)).^2; %objective function
cost = [cost,thiscost];
figure(12); plot(cost); title('cost'); drawnow;

%update epsilon
eps = max(eps/eta,epsmin);

%compute sos-polynomial
mu = zeros(overres);
for j=1:length(s)
    filter = zeros(overres);
    filter(ind_filter) = ifftshift(reshape(U(:,j),filter_siz));
    mu = mu + ((1/s(j))^q)*(abs(ifft2(filter)).^2);    
end

sqrtmu = sqrt(abs(mu));
figure(14); imagesc(sqrtmu.^(1/2)); colorbar; colormap jet; title('mask');    
drawnow    

%ADMM solution of WL2 problem
gam = max(mu(:))/10;  %tenth of sos-mask max value
L = zeros(size(dz));
for l = 1:150
    % Y subprob 
    Z = gam*(M(x_pad)-L);
    muinv = repmat((mu + gam).^(-1),[1,1,size(dz,3)]);
    Y = fft2(muinv.*ifft2(Z));

    % x subprob
    W = Atb_pad + (gam/lambda)*Mt(Y+L);
    x_pad = W./(mask_pad + (gam/lambda)*MtMmask);

    % L update
    L = L + Y - M(x_pad);
end

S_all(:,i)=s;
Mu_all(:,:,i)=mu;
XP_all(:,:,i)=x_pad;

x = reshape(x_pad(ind_full),res);
X = ifft2(x);
figure(15); imshow(abs(X),[]); title('current solution'); drawnow;
figure(16); imagesc(fftshift(log(1+abs(x)))); colorbar; title('current solution, k-space'); drawnow;

thesestats = imstats(X0,X);
SNR = [SNR,thesestats.SNR];
relerr = [relerr,norm(X-X_old,Inf)/norm(X,Inf)];
figure(17); plot(SNR); title('SNR'); drawnow;
figure(18); plot(relerr); title('relative error between iterates');
X_old = X;
end
X_LRG = X;
stats.LRG = imstats(X0,X_LRG);
stats.LRG
%% Comparison with TV
foptions.mask = mask0;
foptions.siz = res;
operators = defAAt_fourier(foptions);
options.siz = res;
options.Nouter = 100;   %ADMM iterations
options.Ninner = 1;
options.beta = 10;  %ADMM parameter
x0 = zeros(res);
x0(logical(foptions.mask)) = b;
x0 = ifft2(x0);
btv = operators.A(x0);
% run TV minimization
lambda_tv = 1e-3;
[X_TV, cost] = OpTV_AL(btv,operators,lambda_tv,options);
stats.TV = imstats(X0,X_TV);
stats.TV

figure(7); imshow(abs(X_TV),[]); title('standard TV');