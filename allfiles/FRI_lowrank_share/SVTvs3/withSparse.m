clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
r = 4; Ns = 24; m = 128*128;
load data_2D_withoutmotion_complex.mat

X = reshape(data_complex,m,24);
% take dwt2 here Np init
Xhat = fft(X,[],1);
[a, b, c, d] = dwt2(reshape(X(:,1),128,128),'Haar');
ee = [a,b;c,d];
Np = length(find(ee > max(ee(:))*0.8))+1;
figure(2); imagesc(abs(X));


%%
sampset = randperm(length(Xhat(:)));
sampset = sampset(1:2*(Np-1)*r); % 6kr

b = Xhat(sampset);
%%
% Run SVT algorithm

lambda = 0.5; %regularization parameter
beta = 1e-2; %ADMM parameter

Q = @(x) giveConcatHankMtx(x,Np);
Qt = @(x) giveAvgConcatHank(x,Np);
A = @(x) Aforward(x,sampset);
At = @(b) Abackward(b,sampset,m,Ns);

Atb = At(b);
x = Atb; %space domain
Qx = Q(x);
G = zeros(size(Qx)); %lagrange multiplier
Z = zeros(size(Qx)); %aux variable
AtA = zeros(size(Xhat)); AtA(sampset)=1;
QtQ = Qt(Q(ones(size(Xhat)))); % why ones?

cost = [];
constraint = [];
err = [];

for i=1:100     
    %X subproblem  
    quad = AtA + lambda*beta*QtQ;
    quadinv = 1./quad;   % same for all iters, define above   
    x = quadinv.*(Atb + lambda*beta*(Qt(Z-G)));

    figure(11); imagesc(abs((x))); title('current solution');

    %Calculate error--switch off to save comp time          
    Qx = Q(x);
%     [~, S, ~] = svd(Qx,'econ');    
%     nuclear = norm(diag(S),1);
%     diff = A(x)-b;
%     thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
%     cost = [cost,thiscost];
%     err = [err, norm(x-Xhat,'fro')];
%     figure(12); plot(cost); title('cost');
%     %figure(13); imagesc(fftshift(log(1+abs(x)))); colorbar;
%     figure(14); plot(diag(S)); title('sing vals');
%     figure(15); plot(err); title('error');

    
    %Shrinkage step
    Z = Qx + G;

    [U, S, V] = svd(Z,'econ');
    s = diag(S);
    
    s0 = shrink1(s,1/beta);
    Z = U*diag(s0)*V';     
      
    G = G + Qx - Z;
end
