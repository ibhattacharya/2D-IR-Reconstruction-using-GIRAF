%% IRLS Nuclear Norm Minimization

% min_X 0.5*||AX-b||_2^2 + lambda*||ToepX||_*
%% Load data

clear all; addpath('./algorithms','./data','./etc'); clc; close all;
res = [64 64]; img = phantom(res(1));
%% Set model order

filtSiz = [12,12]; p = prod(filtSiz); N = prod(res);
%filter dim (odd nums) %overres = res + 2*filter_siz;
%% Generate fourier undersampling operator A

u = 0.25; pdf = genPDF(res,7,u); mask = im2double(genSampling(pdf,2,5));
mask(1) = 1; %figure; imshow(fftshift(mask)); title('Undersampling mask');
ind = find(mask~=0); clear usf mask pdf

[A,At] = defAAt(ind,res);
b = A(img);
%% init variables and operators

AtA = @(x) proj(x,ind);
Atb = At(b); x = Atb; %........in k-space (1)

%% parameters

lambda = 1; lambdainv = 1/lambda;  %regularization parameter
y = lambdainv*Atb; clear Atb
y = y(:); Niter = 2; cgIter = 2000;
%% 2 stage algorithm

for iter=1:Niter
    % 1 filter update
    W = zeros(p);
    for i=1:p
        W(i,i:p) = cconv(conj(fliplr(x(i:N-p+i))),(x(i:N-p+i)),p-i+1);
    end
     figure;imagesc(abs(W));
    W = W+W.'+diag(-diag(W)); 
    
    W = inv(W); W = sqrtm(W); W = sqrtm(W); %........(2)
   W = ifftshift(W);
    % 2 Weighted L2 problem
    T = @(x) getTms(x,filtSiz,res,W);
    
    %cost = [];
    %     % calculate cost
    %     [~,S,~] = svd(toeplitz(x,p),'econ'); s = diag(S); %%% bad
    %     figure; plot(s); title('Sing vals');
    %     nuclear = norm(s,1); diff = A(x)-b;
    %     thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
    %     cost = [cost,thiscost]; clear s S diff nuclear
    %     figure; plot(cost); title('Cost');
    
    R = @(x) T(x) + lambdainv*AtA(x); %.........(4)
    [x,flag,relres] = pcg(R,y,1e-9,cgIter,[],[],x(:));
    
    figure(2+iter); imagesc(abs(reshape(x,res))); colorbar; title('Current solution');axis off;
end

%% Plot Original Data

%figure(Niter+3); imagesc(img); colorbar; title('Original');