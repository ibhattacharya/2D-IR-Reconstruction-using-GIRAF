clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');

sp = [6 ]; ra = [ 3]; index = randperm(128);
error1 = zeros(length(sp));error2 = error1; error3 = error1; error4 = error1; error5 = error1;
for i=1:length(sp)
    for j=1:length(ra)
        k = sp(i); r = ra(j); Ns = 45;
        U = zeros(128,r);
        U(index(1:k),:) = randn(k,r);
        X = U*randn(r,Ns); % original matrix
        Xhat = fft(X,[],1); % in k space
        figure(1); imagesc(abs(X));title('Original');
        
        sampset = randperm(length(Xhat(:)));
        sampset = sampset(1:r*(k+Ns-r)); % 2kr
        b = Xhat(sampset);
        
%         %% cvx
%         lambda1 = 1e-6; lambda2 = 1e-5;
%         
%         cvx_begin
%             variable xc(128,Ns)
%             xch = dftmtx(128)*xc;
%             minimize( lambda1*norm(xc,1) +  lambda2*norm_nuc(xc) )
%             subject to
%             xch(sampset) == b;
%         cvx_end
%         e1 = norm(xc-X,'fro')/norm(xc,'fro');
%         figure(2); imagesc(abs((xc))); title('Rank+Sparsity CVX');
%         %% nuclear only
%         mu1 = 1e-7; mu2 = 0;
%         X1 = main_sajan(Xhat,b,sampset,mu1,mu2);
%         X2 = flipud(X1); X2 = circshift(X2,[1 0]);
%         e2 = norm(X2-X,'fro')/norm(X2,'fro');
%         figure(3); imagesc(abs((X2))); title('Rank only'); clear X1 X2
%         %% l2l1 only
%         mu1 = 0; mu2 = 1e-7;
%         X1 = main_sajan(Xhat,b,sampset,mu1,mu2);
%         X2 = flipud(X1); X2 = circshift(X2,[1 0]);
%         e3 = norm(X2-X,'fro')/norm(X2,'fro');
%         figure(4); imagesc(abs((X2))); title('Sparsity only'); clear X1 X2
%         %% nuclear + l2l1
%         mu1 = 1e-7; mu2 = 1e-8;
%         X1 = main_sajan(Xhat,b,sampset,mu1,mu2);
%         X2 = flipud(X1); X2 = circshift(X2,[1 0]);
%         e4 = norm(X2-X,'fro')/norm(X2,'fro');
%         figure(5); imagesc(abs((X2))); title('Rank+Sparsity'); clear X1 X2
        
        %% svt approach
        
        Np = k+1;
        lambda = 0.5; %regularization parameter
        beta = 1e-2; %ADMM parameter
        
        Q = @(x) giveConcatHankMtx(x,Np);
        Qt = @(x) giveAvgConcatHank(x,Np);
        A = @(x) x(sampset);
        At = @(b) Atranspose(b,sampset,128,Ns);
        
        Atb = At(b);
        x = Atb;  % init recon
        Qx = Q(x);
        G = zeros(size(Qx)); %lagrange multiplier
        Z = zeros(size(Qx)); %aux variable
        AtA = zeros(size(Xhat)); AtA(sampset)=1;
        QtQ = Qt(Q(ones(size(Xhat)))); % why ones?
        
        cost = [];        constraint = [];       
        
        quad = AtA + lambda*beta*QtQ;        quadinv = 1./quad;
        
        for ii=1:100
            %X subproblem
            x = quadinv.*(Atb + lambda*beta*(Qt(Z-G))); % recon in k space
                        
            %Calculate error--switch off to save comp time
            Qx = Q(x);
            %[~, S, ~] = svd(Qx,'econ'); nuclear = norm(diag(S),1);     diff = A(x)-b;
            %thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
            %cost = [cost,thiscost]; figure(12); plot(cost); title('cost');
            %figure(13); imagesc(fftshift(log(1+abs(x)))); colorbar;
            %figure(14); plot(diag(S)); title('sing vals'); %figure(15); plot(err); title('error');
            
            %Shrinkage step
            Z = Qx + G;   [U, S, V] = svd(Z,'econ');
            s = diag(S);    s0 = shrink1(s,1/beta);
            Z = U*diag(s0)*V';    G = G + Qx - Z;
        end
        e5 = norm(x-Xhat,'fro')/norm(x,'fro');
        figure(6); imagesc(abs(ifft(x,[],1))); title('SVT solution');%pause(0.1);
   %     error1(i,j) = -mag2db(e1);
%         error2(i,j) = -mag2db(e2);
%         error3(i,j) = -mag2db(e3);
%         error4(i,j) = -mag2db(e4);
        error5(i,j) = -mag2db(e5);
    end
end
ttt = [error1;error5];
figure(7);imagesc([error1;error5]);
save('kr650.mat','ttt');
e5 = norm(ifft(x,[],1)-X,'fro')/norm(ifft(x,[],1),'fro');-(e5)
e6 = norm(x-Xhat,'fro')/norm(x,'fro');-(e6)