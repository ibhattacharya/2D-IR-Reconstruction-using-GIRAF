%--------------------------------------------------------------------------
% CG solution to region limited problem
% [X,Potential,ErrorArray,ErrorIndex,efinal] = xupdateHOTV(A,b,baseline,mask,kappa,lambda,mu,Niter,Potential)
% Solves {X*} = arg min_{X} ||Af-b||^2 + mu ||Rf||_{l_1}
%--------------------------------------------------------------------------

function [X,earray1] = CG_solver(b,A, At,D,Dt,w,Lam1,Lam4,Y1,C, X, THRESHOLD,Niter)

oldcost = 0;
earray1 = [];
lam1 = 0.5*C.mu1*C.beta1;
%lam1 = 10000;
lam2 = 0.5*C.mu2*C.beta2;

LamTV = Lam1; 
Lam4=Lam4+1i*1e-18;
eTV=double(0);eNN=double(0);
for i=1:Niter,
    
    resY = (A(X) - b');
    eY = sum(abs(resY(:)).^2);
    
    resw = X-w; 
    eNN = lam1*sum(abs(resw(:)).^2);
    eNN = eNN + abs(C.mu1*(Lam4(:)'*resw(:)));
    
    resTV = X-Y1; 

    LamTV=LamTV+1i*1e-18;
    resTV1=resTV; 
    
    eTV = lam2*(abs(resTV(:)'*resTV(:)).^2);
    
    eTV = eTV + C.mu2 * abs(conj(LamTV(:)'*resTV(:))); 
    
    
    cost1 = eY + eNN + eTV;
    
    earray1 = [earray1,cost1];
    
    if(abs(cost1-oldcost)/abs(cost1) < THRESHOLD)
       % i
        break;
    end
    oldcost = cost1;
    
  %  conjugate gradient direction
   % ------------------------------
    
    % gradient: gn
    
    gn = At(A(X)-b') + lam1*(X-w) + lam2*Dt(resTV);
    gn = 2*gn+C.mu1*Lam4 + C.mu2*Dt(LamTV);
    
    % search direction: sn  
    if(i==1)
        sn = gn;                                          
        oldgn = gn;
    else
        gamma = abs(sum(gn(:)'*gn(:))/sum(oldgn(:)'*oldgn(:)));
        sn = gn + gamma*sn; 
        oldgn = gn;
    end
    
    % line search
    %-------------
    Asn = A(sn);  
    Dsn = D(sn);
    
    numer = Asn(:)'*resY(:) + lam1*sn(:)'*resw(:) + 0.5* C.mu1*sn(:)'*Lam4(:);
    numer = numer + lam2*sum(sum(sum(conj(resTV).*Dsn)))   +  0.5*C.mu2*sum(sum(sum(conj(LamTV).*Dsn))); 
    
    denom = Asn(:)'*Asn(:) + lam1*sn(:)'*sn(:); 
    denom = denom + lam2*sum(sum(sum(conj(Dsn).*Dsn))); 
    if(denom < 1e-18)
        break;
    end
    alpha = -real(numer)/real(denom);
   
    % updating
    %-------------
    
    X = (X + alpha*sn);
end

    
