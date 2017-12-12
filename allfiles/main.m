%% GIRAF reconstruction for 2D IR Spectroscopy
% Author : Ipshita Bhattacharya, Jonathan Humston

% Step 1: Load data, intialize parameters
% Step 2: Solve GIRAF based reconstruction
% Step 3: Post processing : Calculate centerline slope of the main peak

%GIRAF solves min_X (1/p)*||T(X)||_p^p + (lambda/2)*||AX-b||_2^2
%where A is a Fourier undersampling mask in time domain, and
%T(X) is the block Toeplitz matrix built from 2D time domain data X
%Uses a two-stage IRLS MM approach, which alternates between updating
%an annihilating filter mask and a least squares annihilation


%% Add paths for suppoting files

clear all;close all;clc

addpath('./FRI_lowrank_share/');
addpath('./im2col/');  % this flder has mex files which might need to be recompiled. Refer to MATLAB mex commands

%% Load and preprocess data

load('2D_IRSpec_GIRAF_data.mat');
% Replace Nans and infinities in the data by 0
data(isnan(data)) = 0; data(isinf(data)) = 0;

% Convert from nonlinear w3 to linearized w3 data
run('CreateAxesAll.m');
run('LinearGriddedData');


%Truncate in spectral domain to remove noise from ends. Result must be an even length vector
opts.w3Trunc_low=274;opts.w3Trunc_high=775;
opts.w1Trunc_low=1;opts.w1Trunc_high=1024;

%Choose a scaling factor for the data
opts.magScale = 1e4;
opts.Nf2=1024;

% Cosine filtering 

[data_reconstructed,opts.Nt3,opts.N,opts.Nt] = cosine_filtering(data,opts);
% data_reconstructed : data after fitering

%% Set zoom in region for spectra visualization and analysis

opts.xlow=120;opts.xhigh=380;
opts.ylow =175; opts.yhigh = 350;

%% GIRAF Reconstruction : Set up data for reconstruction

[X0,b0,opts.N1,opts.N2,opts.N3] = GIRAF_preprocessing(data_reconstructed);
% X0 : true fully sampled data
% b0 ; fully sampled time domain data
% N1,N2,N3 : data dimensions

% Choose image scale for visulization
% im_high : scale for displaying
Xp = X0(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh);  opts.im_high = max(abs(Xp(:)));

%% GIRAF Reconstruction : Initialize parameters


% Set optimization parameters
opts.lambda0 =50;  %Larger lambda0 puts more weight on the model and less on the data
opts.iter = 25; %iterations - can increase for better convergence
opts.US = 5; % under samling factor

% load undersampling mask

load(sprintf('undersamplingMask_US%d.mat',opts.US));

% Create new sampling mask for variable US (change value of US & uncomment the next section)

% samp =  Hybrid_DownsamplingMASK(N*2, N1-1, N3, US); % undersampling by US
% s1 = permute(samp,[2,1,3]);
% s2= s1(:,N+1:end,:);
% s2 = [s2,ones(N1-1,Nf2-1-N,N3)]; % the zero padded dimension t1 has ones
% for points N through Nf2 in the mask

% s2(:,1)=1;  % the first point in t1 axis is always collected
% mask0=s2;

%% undersample data

b1 = mask0.*b0;
ind_samples = find(mask0~=0);
b = b1(ind_samples);
opts.res = [opts.N1-1 opts.N2-1];

% forward and backward opeator for undersampling
[A,At] = defAAt(ind_samples,opts.res);

% Fourier reconstructionof undersampled data
data_US_fourier = fft2(At(b));

% RMSE calculation
[a1,b1]=RMSE(data_US_fourier,X0);

figure(1);

suptitle('Fully sampled data and Fourier reconstruction of undersampled data')
subplot(2,1,1);contour(real(X0(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh)));title('Fully sampled')
subplot(2,1,2);contour(real(data_US_fourier(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh)));title('Undersampled Fourier Transform Recon'); xlabel(sprintf('RMSE = %d',a1));




% Set model order
opts.filter_siz = [21,21]; %filter dimensions (use odd nums)
opts.filter_siz2 = 2*opts.filter_siz - [1,1]; %squared filter dimensions
opts.overres = opts.res + 2*opts.filter_siz; %over-resolved reconstruction grid

clear k;
k = zeros(2,prod(opts.overres));
kx = ifftshift(repmat(-(opts.overres(2)-1)/2:(opts.overres(2)-1)/2,opts.overres(1),1));
ky = ifftshift(repmat((-(opts.overres(1)-1)/2:(opts.overres(1)-1)/2).',1,opts.overres(2)));
k(1,:) = kx(:);
k(2,:) = ky(:);
%
opts.ind_full = find((abs(k(1,:)) <= (opts.res(2)-1)/2) & (abs(k(2,:)) <= (opts.res(1)-1)/2));
opts.ind_filter = find((abs(k(1,:)) <= (opts.filter_siz(2)-1)/2 ) & ( abs(k(2,:)) <= (opts.filter_siz(1)-1)/2));
opts.ind_filter2 = find((abs(k(1,:)) <= (opts.filter_siz(2)-1) & (abs(k(2,:)) <= (opts.filter_siz(1)-1))));



%% initialize variables and operators

dz = ones(opts.overres);

M = @(z) repmat(z,[1,1,size(dz,3)]).*dz;
Mt = @(Z) sum(Z.*conj(dz),3);
MtMmask = Mt(M(ones(opts.overres)));


opts.mask_pad = zeros(opts.overres);
opts.mask_pad(opts.ind_full) = mask0;
ind_samples_pad = find(opts.mask_pad~=0);
AtA = @(x) proj(x,ind_samples_pad);

Atb = At(b);
Atb_pad = zeros(opts.overres);
Atb_pad(opts.ind_full) = Atb;


opts.wind_fun = ones(size(Atb_pad));

x = Atb;
x_pad = zeros(opts.overres);
x_pad(opts.ind_full) = x;
%% Run Algorithm .

    % Zoom into spectral region for better visualization
   opts.xl_trunc = 125; opts.xh_trunc = 375;
    opts.yl_trunc = 235; opts.yh_trunc = 290;

 [data_US_GIRAF] = runGIRAF(x_pad,x,Atb_pad,M,A,At,b,data_US_fourier,X0,opts);

%% Final Reconstructed data processing
data_final = real(data_US_GIRAF(:,1:opts.Nf2/2));

figure;
[C,H]=contour(real(data_final(opts.xl_trunc:opts.xh_trunc,opts.yl_trunc:opts.yh_trunc)),12);title('Final GIRAF Reconstruction'); 
set (H, 'LineWidth', 3);
%% Calculate the centerline slope

% zoom in to axis for CLS analysis
opts.w3_low = 155; opts.w3_high=340;
opts.w1_low = 252; opts.w1_high=280;
opts.n = n;
opts.w3_axisL=w3_axisL;

cls = cls_analysis(data_final,opts);

%% Output and Save
Spec = 'Recon_';
USstring = num2str(US);
currentLambda = num2str(opts.lambda0);
lambda0 = opts.lambda0;
iter = opts.iter;
Specname = strcat(Spec,'US',USstring,'_lambda0_',currentLambda,'.mat');

save(Specname,'data_final','cls','lambda0','US','iter');





