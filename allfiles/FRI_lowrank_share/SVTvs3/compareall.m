clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
k = 5; r = 5; Ns = 10;

U = zeros(128,r);
index = randperm(128);

U(index(1:k),:) = randn(k,r);

X = U*randn(r,Ns); % original matrix
Xhat = fft(X,[],1); % in k space

figure(1); imagesc(abs(X));title('Original'); 

sampset = randperm(length(Xhat(:)));
sampset = sampset(1:(128+Ns*k+Ns-2*k)*k); % 6kr
b = Xhat(sampset);

% %% cvx
% lambda1 = 1e-6; lambda2 = 1e-5;
% % b1 = zeros(128,Ns); b1(sampset) = b;
% cvx_begin
%     variable xc(128,Ns)
%     xch = dftmtx(128)*xc;
%     minimize( lambda1*norm(xc,1) +  lambda2*norm_nuc(xc) )
%     subject to 
%           xch(sampset) == b;
%     %xc == ifft(b1,[],1); %l21nrm_cvx check % xch = fft(xc,[],1); xch(sampset) == b;
% cvx_end
% figure(2); imagesc(abs((xc))); title('Rank+Sparsity CVX'); 
%  return
% %% nuclear only
% mu1 = 1e-7; mu2 = 1e-8;
% X1 = main_sajan(Xhat,b,sampset,mu1,mu2);
% %figure(3); imagesc(abs((X1))); title('Rank only'); %clear X1
% X2 = flipud(X1); X2 = circshift(X2,[1 0]); 
% figure(4); imagesc(abs((X2))); 
% % return
% % %% l2l1 only
% % mu1 = 0; mu2 = 1e-3;
% % [e2,X1] = main_sajan(Xhat,b,sampset,mu1,mu2);
% % figure(4); imagesc(abs(ifft(X1,[],1))); title('Sparsity only'); clear X1
% % %% nuclear + l2l1
% % mu1 = 0.06; mu2 = 0.005;
% % [e3,X1] = main_sajan(Xhat,b,sampset,mu1,mu2);
% % figure(5); imagesc(abs(ifft(X1,[],1))); title('Rank+Sparsity'); clear X1

%% svt approach
Np = k+1;
lambda = 0.5; %regularization parameter
beta = 1e-2; %ADMM parameter

Q = @(x) giveConcatHankMtx(x,Np);
Qt = @(x) giveAvgConcatHank(x,Np);
A = @(x) x(sampset);
At = @(b) Atranspose(b,sampset,128,Ns);

Atb = At(b);
x = Atb;   figure(12);imagesc(abs(ifft(x,[],1)));% init recon
Qx = Q(x);
G = zeros(size(Qx)); %lagrange multiplier
Z = zeros(size(Qx)); %aux variable
AtA = zeros(size(Xhat)); AtA(sampset)=1;
QtQ = Qt(Q(ones(size(Xhat)))); % why ones?

cost = [];
constraint = [];
err = [];

quad = AtA + lambda*beta*QtQ;
quadinv = 1./quad;
    
for i=1:100
    %X subproblem
    x = quadinv.*(Atb + lambda*beta*(Qt(Z-G))); % recon in k space
    
    figure(6); imagesc(abs(ifft(x,[],1))); title('SVT solution');pause(0.1);
    
    %Calculate error--switch off to save comp time
    Qx = Q(x);
    %[~, S, ~] = svd(Qx,'econ');
    %nuclear = norm(diag(S),1);     diff = A(x)-b;
    %thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
    %cost = [cost,thiscost];
    err = [err, norm(x-Xhat,'fro')/norm(x,'fro')];
    %figure(12); plot(cost); title('cost');
    %figure(13); imagesc(fftshift(log(1+abs(x)))); colorbar;
    %figure(14); plot(diag(S)); title('sing vals');
    %figure(15); plot(err); title('error');
    
    %Shrinkage step
    Z = Qx + G;
    [U, S, V] = svd(Z,'econ');
    s = diag(S);
    s0 = shrink1(s,1/beta);
    Z = U*diag(s0)*V';
    
    G = G + Qx - Z;
end
%-mag2db([e1 e2 e3 err(end)])