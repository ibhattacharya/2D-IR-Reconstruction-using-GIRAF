%% Singular value thresholding ADMM algorithm
%Solves min_x 0.5*||Ax-b||_2^2 + lambda*||Qx||_*
%where A is a Fourier undersampling mask, and
%Qx is the block Toeplitz matrix built from k-space data x
%This code is meant for benchmarking only
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
usf = 0.2; %undersampling factor
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
%% Run SVT algorithm
lambda = 10; %regularization parameter
beta = 1e-2; %ADMM parameter

Atb = At(b);
x = Atb; %initialize X
Qx = Q*x;
G = zeros(size(Qx)); %lagrange multiplier
Z = zeros(size(Qx)); %aux variable
AtA = mask;

cost = [];
constraint = [];
for i=1:20     
    %X subproblem  
    quad = AtA + lambda*beta*QtQ;
    quadinv = 1./quad;      
    x = quadinv.*(Atb + lambda*beta*(Q'*(Z-G)));

    figure(11); imshow(abs(ifft2(x)),[]); title('current solution');

    %Calculate error--switch off to save comp time          
    Qx = Q*x;
    [~, S, ~] = svd(Qx,'econ');    
    nuclear = norm(diag(S),1);
    diff = A(x)-b;
    thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
    cost = [cost,thiscost];

    figure(12); plot(cost); title('cost');
    figure(13); imagesc(fftshift(log(1+abs(x)))); colorbar;
    figure(14); plot(diag(S)); title('sing vals');
        
    
    %Shrinkage step
    Z = Qx + G;

    [U, S, V] = svd(Z,'econ');
    s = diag(S);
    s0 = shrink1(s,1/beta);
    Z = U*diag(s0)*V';     
      
    G = G + Qx - Z;
end
%% Plot Original Data
x0 = reshape(m,res);
figure(2); imshow(abs(ifft2(x0)),[]); title('original');
figure(3); imagesc(fftshift(log(1+abs(x0)))); colorbar;

