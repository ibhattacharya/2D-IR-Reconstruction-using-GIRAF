%% IRLS Nuclear Norm Minimization
%Solves min_X 0.5*||AX-b||_2^2 + lambda*||QX||_*
%where A is a Fourier undersampling mask, and
%QX is the block Toeplitz matrix built from k-space data X
%Uses a two-stage IRLS MM approach, which alternates between updating
%an annihilating filter/edge mask and a least squares annihilation
%% Load data
clear all;
addpath('./algorithms','./data','./etc');
load MR_SL_convert;
%load MR_brain_convert;
%load MR_SL_linear;

%% Set model order
res = [91,91];  %output resolution (use odd nums)
filter_siz = [25,33]; %filter dimensions (use odd nums)
overres = res + 2*filter_siz;

ind_trim = find((abs(k(1,:)) <= (overres(2)-1)/2 ) & (abs(k(2,:)) <= (overres(1)-1)/2));
m = reshape(m(ind_trim),overres);
k = k(:,ind_trim);

ind_full = find((abs(k(1,:)) <= (res(2)-1)/2) & (abs(k(2,:)) <= (res(1)-1)/2));
ind_filter = find((abs(k(1,:)) <= (filter_siz(2)-1)/2 ) & ( abs(k(2,:)) <= (filter_siz(1)-1)/2));
%% Generate fourier undersampling operator
usf = 0.25; %undersampling factor
pdf = genPDF(res,7,usf); %density compensation
mask0 = im2double(genSampling(pdf,2,5));
mask0(1) = 1; %need DC component
mask = zeros(overres);
mask(ind_full) = mask0;
%ind_box = find((abs(k(1,:)) <= 12) & (abs(k(2,:)) <= 12));
%mask(ind_box) = 1;

figure(1); imshow(fftshift(reshape(mask(ind_full),res)));
ind_samples = find(mask~=0);
[A,At] = defAAt(ind_samples,overres);
b0 = A(m);
sig = 0; %std deviation 
b = b0+sig*(randn(size(b0))+ 1i*randn(size(b0)));
fprintf('Sample SNR = %5.2f dB\n',20*log10(norm(b)/norm(b-b0)));
%% init variables and operators 
dz{1} = reshape((-1i*pi*(k(1,:))).',overres)/(overres(1)/4);
dz{2} = reshape((-1i*pi*(k(2,:))).',overres)/(overres(2)/4);

Q = ForwardModelLowRankThin(dz,overres,filter_siz,[1,1]);
QhQ = @(x,mu) defQhQ( x, mu, dz, overres, ind_full, 0 );
AtA = @(x) proj(x,ind_samples);

Atb = At(b);
x = Atb;

%% run alg.
lambda = 1e5;  %regularization parameter
%parameters for epsilon update: eps = min(eps,gamma*s(r+1));
eps = 0.1;       %epsilson initialization
gamma = 10;     %epsilon multiplier for update
r = 500;       %rank cut-off value

r2 = 0;        %singular value cut-off for mask recon
q = 2;         %q=1 is nuclear norm, q=2 is Shatten p=1/2

lambdainv = 1/lambda;
y = lambdainv*Atb;
y = y(:);
%%
cost = [];
epsrecord = [];
for i=1:10
% Stage 1: Compute annihilating mask
% Run SVD -- could possibly replace with lansvd to speed up 
%(see /etc/lansvd) folder
[U,S,~] = svd(Q*x,'econ');
s = diag(S);
figure(11); plot(s); title('sing vals');

%calculate cost
nuclear = norm(s,1);
diff = A(x)-b;
thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
cost = [cost,thiscost];
figure(12); plot(cost); title('cost');

%update epsilon
eps = min(eps,gamma*s(r+1)); 
epsrecord(end+1) = eps;
figure(13); plot(epsrecord); title('eps');

mu = zeros(overres);
for j=(r2+1):length(s)
    filter = zeros(overres);
    filter(ind_filter) = ifftshift(reshape(U(:,j),filter_siz));
    mu = mu + ((1/max(s(j),eps))^q)*(abs(ifft2(filter)).^2);
end

mask = sqrt(abs(mu));
figure(14); imagesc(mask); colorbar; colormap jet; title('mask');    
    
% Stage 2: Weighted L2 problem    
R = @(x) QhQ(x,mu) + lambdainv*AtA(x); %update gradient with new mask mu
[x,flag,relres] = pcg(R,y,1e-9,2000,[],[],x(:)); 
%note: sometimes needs >1000 iterations for good convergence. 
%We should try and think about a preconditioner?

x = reshape(x,overres);
xtrim = reshape(x(ind_full),res);
X = ifft2(xtrim);
figure(15); imagesc(abs(X)); colorbar; title('current solution')
figure(16); imagesc(fftshift(log(1+abs(xtrim)))); colorbar; title('current solution, k-space')
end 

%% Plot Original Data
x0 = reshape(m(ind_full),res);
figure(2); imagesc(abs(ifft2(x0))); colorbar; title('original');
figure(3); imagesc(fftshift(log(1+abs(x0)))); colorbar;
%% Plot Original Data
x1 = reshape(Atb(ind_full),res);
figure(4); imagesc(abs(ifft2(x1))); colorbar; title('original');
figure(5); imagesc(fftshift(log(1+abs(x1)))); colorbar;

