clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
r = 4; Ns = 24; 
load data_2D_withoutmotion_complex.mat
c1 = 30; c2 = 17;
d = data_complex(c1:end-c1,c2:end-c2,:); m = size(d,1)*size(d,2);

X = reshape(d,m,24); Xhat = zeros(m,Ns);
for i = 1:Ns
   Xhat(:,i) =  reshape(fft2(d(:,:,i)),m,1);
end

Np = 90+1;
B = giveConcatHankMtx(Xhat,Np);

figure(2); imagesc(abs(B));
figure(3); imagesc(abs(Xhat),[0 50]);
%tt = zeros(1,Ns);
for i = 1:Ns
%tt(i) = rank(B(:,(i-1)*Np+1:i*Np));
[~,cc,~]=svd(B(:,(i-1)*Np+1:i*Np),'econ');figure(8);plot(diag(cc));
pause(0.2);
end
%tt
rank(B)
[~,a,~]=svd(X,'econ'); figure(5);plot(diag(a));
[~,b,~]=svd(Xhat,'econ');figure(6);plot(diag(b));
[~,c,~]=svd(B,'econ');figure(7);plot(diag(c));
return
% [u,S,v] = svd(B,'econ');
% s = diag(S)'
% 
% figure(2); imagesc(abs(X));


%%
sampset = randperm(length(Xhat(:)));
sampset = sampset(1:6*k*r); % 6kr

b = Xhat(sampset);
%%
% Run SVT algorithm

lambda = 0.5; %regularization parameter
beta = 1e-2; %ADMM parameter

Q = @(x) giveConcatHankMtx(x,Np);
Qt = @(x) giveAvgConcatHank(x,Np);
A = @(x) x(sampset);
At = @(b) Atranspose(b,sampset,m,Ns);

Atb = At(b);
x = Atb; 
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

    figure(11); imagesc(abs(ifft(x,[],1))); title('current solution');

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
