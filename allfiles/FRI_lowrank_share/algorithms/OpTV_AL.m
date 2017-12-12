function [X, cost] = OpTV_AL(b,operators,lambda,options)
% OPTV: Solves TV regularized inverse problems with an alternating 
% minimization algorithm. Returns X, the recovered image, and earray,
% the values of the cost function for all iterations.
%
% Based on the fTVd implementation: 
% http://www.caam.rice.edu/~optimization/L1/ftvd/
A = operators.A;
At = operators.At;
Nouter = options.Nouter;
Ninner = options.Ninner;
beta = options.beta;
siz = options.siz;


%Define AtA
p_image = zeros(siz,'double'); p_image(1,1) = 1;
AtA = fft2(At(A(p_image)));

%Define derivative operators D, Dt, and DtD
[D,Dt] = defDDt_TV(siz);
DtD = fft2(Dt(D(p_image)));

%Initialize X
Atb = At(b);
X = Atb;
DX = D(X);
G = zeros([siz,2]);

% Begin alternating minimization alg.
cost = [];
for i=1:Nouter 
    for ii=1:Ninner                       
        %Shrinkage step
        Z = DX + G;
        Z1 = Z(:,:,1);
        Z2 = Z(:,:,2);
        AZ = sqrt(abs(Z1).^2 + abs(Z2).^2);
        shrinkZ = shrink(AZ,1/beta); %shrinkage of gradient mag.
        Z(:,:,1) = shrinkZ.*Z1; 
        Z(:,:,2) = shrinkZ.*Z2;

        %Inversion step 
        F1 = fft2(2*Atb + lambda*beta*Dt(Z-G));
        F2 = lambda*beta*DtD + 2*AtA;
        X = ifft2(F1./F2);
        
        %Calculate error        
        DX = D(X);
        DX1 = DX(:,:,1); 
        DX2 = DX(:,:,2);
        ADX = sqrt(abs(DX1).^2 + abs(DX2).^2); %isotropic TV           
        diff = A(X)-b;
        thiscost = norm(diff(:)).^2 + lambda*sum(ADX(:)); %objective function
        cost = [cost,thiscost];

        %Convergence test
        if ii > 2
            if abs((cost(end)-cost(end-1)))<tol
                break;
            end
        end                  
    end
    %Lagrange multiplier update
    G = G + DX - Z; 
    
    %Increase continuation parameter
    %beta = beta*bfactor;
        
end
