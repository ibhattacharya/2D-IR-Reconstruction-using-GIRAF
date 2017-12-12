function [U,cost,opts] = minSNandL1L2(A,At,D,Dt, x_init,b,opts);

U=x_init; [m n d] = size(U);

%% Lam1, Lam2, Lam3 and Lam4 are respectively the Lagrange multipliers. 
% Lam4 is for low rank norm while Lam1, Lam2 and Lam3 are for the x, y and t gradients
% refer SG Lingala et al, ISBI 2011 for details
Lam4 = zeros(m,n,d);Lam1=Lam4;

%%
Dfwd=D(U);

o=0;cost=[];
for out = single(1:opts.outer_iter),
    o=o+1;
    for in = single(1:opts.inner_iter)
        
     
 % ================================
        %  Begin Alternating Minimization
        % ----------------
        %   w-subprolem (L1-L2 shrinkage)
        % ----------------
     
    Z1 = Dfwd + Lam1/opts.beta2;
    V = sqrt(sum(abs(Z1).^2,2));
    V(V==0) = 1/opts.beta2;
    V = max(V - 1/opts.beta2, 0)./V;
    W = Z1.*repmat(V,[1,n]);
      % --------------------
        %  Lambda - subproblem (singular value shrnkage)
        % --------------------
         [u,sigma,v] = givefastSVD(reshape((U+Lam4/opts.beta1),m*n,d));
         s=diag(sigma);
      
       
        thres=(1/opts.beta1).*(s.^(opts.p-1));
        
        s=(s-thres);
        s = s.*(s>0);
        Lambda=u*(diag(s))*v';
        Lambda=reshape(Lambda,m,n,d);
        
        % ----------------
        %  Solve for the recon: Conjugate gradient update
        % ----------------
      

 [U,earray1] = CG_solver(b(:),A, At,D,Dt,Lambda,Lam1,Lam4,W,opts, U, 1e-7,30);
% U = real(U); U = U.*(U>0);
%% plot a frame of the recon while it iterates
  %figure(3); imagesc(abs((reshape(U,128,10))));

 %% cost calculations
 e = A(U) - b;     
 [uq,sigmaq,vq] = givefastSVD(reshape((U),m*n,d));
 V1 = sum( sqrt(sum(abs(Dfwd).^2,2)),1);
      
 cost = [cost, sum(abs(e(:)).^2)  +  sum(abs(sigmaq(:)).^(opts.p)./(opts.p))*opts.mu1 + V1*opts.mu2 ];
 %figure(30);plot(double(cost)); hold on; pause(0.01); hold off; title('Cost');
 
 if in>1
        if abs(cost(end) - cost(end-1))/abs(cost(end)) < 1e-3
            break;
  end
        end

Dfwd=D(U);

%%        Update rules for the Lagrange multipliers
    Lam4 = Lam4 - 1.618*opts.beta1*(Lambda - U);
    Lam1 = Lam1 - 1.618*opts.beta2*(W - Dfwd);
    
   
    
    
    
    end 
    
    %% increment beta 1 and beta2
    opts.beta1=opts.beta1*opts.beta1rate;

    opts.beta2=opts.beta2*opts.beta2rate;
end
end 