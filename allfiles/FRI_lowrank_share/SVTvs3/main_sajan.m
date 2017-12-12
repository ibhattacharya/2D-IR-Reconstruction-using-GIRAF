function [ Recon ] = main_sajan(X,b,sampset,mu1,mu2) %(X,Y1,A1,mu1,mu2) old

m = size(X,1); n = size(X,2); 

A = @(z) afwd(z,sampset);
At = @(b) abwd(b,sampset,m,n);
%A = @(x) x(sampset);
%At = @(b) Atranspose(b,sampset,128,n);
%A = @(z)A1*z(:); % The forward Fourier sampling operator old
%At= @(z) reshape(A1'*z(:),m,n); % The backward Fourier sampling operator old

step_size = [1,1,1];
D = @(z) z;
Dt = @(z) z;

% First guess, direct IFFT
x_init  = At(b);
 %figure(11);imagesc(abs(x_init));
%b = Y1(:); old
%x_init = At(b); old

% Call k-t SLR using augmented Lagrangian with continuation (refer: S.G.Lingala et al, ISBI 2011)
%%
%mu1 =1e1; % Regularization parameter for the schatten p-norm
%mu2 =5e0; % Regularization parameter for the L1-L2 norm
opts.mu1 = mu1; 
opts.mu2 = mu2;
opts.p=0.1; % The value of p in Schatten p-norm; p = 0.1- non-convex, p = 1-convex nuclear norm. 
opts.beta1=1e-7; % The continuation parameter for the Schatten p-norm, initialize it
opts.beta2=1e-7; % The continuation parameter for the TV norm, initialize it
opts.beta1rate = 25; % Continuation parmeter increment for the Schatten p-norm
opts.beta2rate = 25; % Continuation parameter increment for the TV norm
opts.outer_iter =30; % Number of outer iterations
opts.inner_iter = 30; % Number of inner iterations


[Recon,cost,opts] = minSNandL1L2(A,At,D,Dt,x_init,b,opts); 
%err = norm(Recon-ifft(X,[],1),'fro')/norm(Recon,'fro');

end

