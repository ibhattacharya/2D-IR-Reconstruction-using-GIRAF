%% Singular value thresholding algorithm with UV trick
%Solves min_X 0.5*||AX-b||_2^2 + lambda||QX||_*
%where A is a Fourier undersampling mask, and
%QX is the block Toeplitz matrix built from k-space data X
%Uses the UV trick: ||Y||_* = min_{Y=UV^*} ||U||_F^2 + ||V||_F^2
%Code can be easily modified to run on GPU (uncomment lines with %GPU)
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
fprintf('Sample SNR = %5.2f dB\n',20*log10(norm(b)/norm(b-b0)));
%% Build Toeplitz operators
step_siz = [1,1]; 
disp((res-filter_siz)./step_siz); %must be integers!
%QX decimation option; e.g. if set to [2,2] skips every other annihilation 
%equation which trims down the row dimension of QX. Doesn't seem to affect 
%results if set to small integers. 

dz{1} = reshape((-1i*pi*(k(1,:))).',res)/(res(1)/4);
dz{2} = reshape((-1i*pi*(k(2,:))).',res)/(res(2)/4);
clear Q;
Q = ForwardModelLowRankThin(dz,res,filter_siz,step_siz);
QtQ = Q'*(Q*ones(res));
disp(size(Q*zeros(res))) %size of QX matrix
%% Run algorithm
lambda = 10;    %reg param
r = 200;       %rank of U,V
beta = 1e-2;   %ADMM quad param--seems to work best at 1e-2 when max(QtQ) ~= 10-50.

Atb = At(b);
x = Atb;
[mQ, nQ] = size(Q*x);
Z = zeros(mQ,nQ);
G = zeros(mQ,nQ); %lagrange multiplier
AtA = mask;       %ones(res);
U = randn(mQ,r); %random init U and V
V = randn(nQ,r);
eyer = eye(r,r);
%U = gpuArray(U); %GPU
%V = gpuArray(V); %GPU
            
cost = [];
constraint = [];
tic
iter = 100;
for i=1:iter
    %C = gpuArray(beta*(Q*X+G)); %GPU
    C = beta*(Q*x+G);    
    U = (C*V)*inv(eyer + beta*V'*V);
    V = (C'*U)*inv(eyer + beta*U'*U);     
    clear C;
    Z = U*V';
    %Z = gather(Z); %GPU
    
    x = (Atb + lambda*beta*(Q'*(Z-G)))./(AtA + lambda*beta*QtQ);
            
    %Lagrange multiplier update
    G = G + (Q*x - Z);
    
    %Calculate error
    nuclear = 0.5*(norm(U,'fro').^2 + norm(V,'fro').^2);
    diff = A(x)-b;
    %lagrange = norm(QX-Z+G,'fro');
    thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; % + 0.5*beta*lagrange.^2; %objective function
    cost = [cost,thiscost];
    thisconstraint = norm(Q*x-U*V','fro')/norm(Q*x,'fro');
    constraint = [constraint, thisconstraint];
    
    %Plots
    figure(11); imshow(abs(ifft2(x)),[]);
    figure(12); plot(cost); title('cost');
    figure(13); imagesc(fftshift(log(1+abs(x)))); colorbar;
    figure(14); plot(constraint); title('constraint');    
end
toc
%% Plot Original Data
x0 = reshape(m,res);
figure(2); imshow(abs(ifft2(x0)),[]); title('original');
figure(3); imagesc(fftshift(log(1+abs(x0)))); colorbar;

