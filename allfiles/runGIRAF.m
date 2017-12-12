function [data_US_GIRAF] = runGIRAF(x_pad,x,Atb_pad,M,A,At,b,data_US_fourier,X0,opts)%x_pad,x,Atb_pad,M,A,At,b,data_US_fourier,X0,lambda0,iter,res,overres,ind_full,ind_filter,ind_filter2,wind_fun,mask_pad,filter_siz,filter_siz2,im_high,xl_trunc,xh_trunc,yl_trunc,yh_trunc,xlow,xhigh,ylow,yhigh)

lambda = 1/opts.lambda0;

% calculate initial value of eps-
% This is a heuristic and may need manual tuning

gradx = M(x_pad);
gradx_fft = fft2(gradx);
sos = ifft2((conj(gradx_fft).*gradx_fft));
sos2 = fftshift(reshape(sos(opts.ind_filter2),opts.filter_siz2));
B = im2colstep(real(sos2),opts.filter_siz,[1,1]) + 1i*im2colstep(imag(sos2),opts.filter_siz,[1,1]);
B = rot90(B,-1);
[U,S] = eig(B);
epscurve = abs(diag(S));
eps = mean(epscurve(8:12));


eta = 1.3;       %epsilon update is eps = eps/eta; (1.3--1.5 work well)
epsmin = 1e-2;   %minimum possible epsilon value

p = 0;           %Shatten-p norm value
q = 1-(p/2);

y = lambda*Atb_pad;
y = y(:);

cost = [];cost_rank=[];cost_data=[];
RMSE_val=[];

for i=1:opts.iter
    
    %Stage 1: Compute annihilating mask
    gradx = M(x_pad);
    gradx_fft = fft2(gradx);
    sos = ifft2((conj(gradx_fft).*gradx_fft));
    sos2 = fftshift(reshape(sos(opts.ind_filter2),opts.filter_siz2));
    R = im2colstep(real(sos2),opts.filter_siz,[1,1]) + 1i*im2colstep(imag(sos2),opts.filter_siz,[1,1]);
    R = rot90(R,-1);
    [U,S] = eig(R+eps*eye(size(R)));
    s = abs(diag(S));
    
    %note: s = sing. values squared.
    % figure(11); plot(0.5*log10(s)); title('sing vals, log scale'); drawnow;
    
    %calculate cost
    if p == 0
        shatten = 0.5*sum(log(s-eps));
    else
        shatten = (1/p)*sum((s-eps).^(p/2));
    end
    diff = A(x)-b;
    thiscost = opts.lambda0*shatten + 0.5*norm(diff(:)).^2; %objective function
    cost = [cost,thiscost];
    c_rank = opts.lambda0*shatten;
    cost_rank = [cost_rank, c_rank];
    c_data = 0.5*norm(diff(:)).^2;
    cost_data =[cost_data,c_data];
    
    % plot cost convergence curves
    figure(12);
    subplot(1,3,1);plot(cost); title('cost');
    subplot(1,3,2);plot(cost_data); title('cost data consistency');
    subplot(1,3,3);plot(cost_rank); title('cost rank');
    drawnow;
    
    %update epsilon
    eps = max(eps/eta,epsmin);
    
    %compute sos-polynomial
    mu = zeros(opts.overres);
    for j=1:length(s)
        filter = zeros(opts.overres);
        filter(opts.ind_filter) = ifftshift(reshape(U(:,j),opts.filter_siz));
        mu = mu + ((1/s(j))^q)*(abs(fft2(filter)).^2);
    end
    
    sqrtmu = sqrt(abs(mu));
    
    % Generate annihilating mask. Should have low values in the
    % in the spectral support region
    figure(14); imagesc(sqrtmu.^(1/2)); colorbar; colormap jet; title('mask');
    drawnow
    
    
    % save singular values and masks at all iterations
    %     S_all(:,i,kk)=s;
    %     Mu_all(:,:,i,kk)=mu;
    %     XP_all(:,:,i,kk)=x_pad;
    
    % Solve the least square problem using the annihilating mask
    
    lhs=@(x)left_fun(x,opts.overres,mu,opts.wind_fun,lambda,opts.mask_pad);
    rhs= lambda*Atb_pad; rhs=rhs(:);
    
    [XX,FLAG,RELRES,ITER,RESVEC] = pcg(lhs,rhs,1e-9,2000,[],[],x_pad(:));
    x_pad = reshape(XX,opts.overres);
    
    
    % reshape GIRAF reconstruction
    x = reshape(x_pad(opts.ind_full),opts.res);
    X = fft2(x);
    
    data_US_GIRAF = X;
    
    % plot results at every iteration
    figure(15);
    subplot(2,3,2);imagesc(abs(data_US_GIRAF(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh)),[0 opts.im_high]); title('GIRAF current solution'); drawnow;
    subplot(2,3,3);imagesc(abs(data_US_fourier(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh)),[0 opts.im_high]);title('Fourier Recon'); drawnow;
    subplot(2,3,1);imagesc(abs(X0(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh)),[0 opts.im_high]);title('Fully sampled');drawnow;
    
     
    TruncX0=real(X0(opts.xl_trunc:opts.xh_trunc,opts.yl_trunc:opts.yh_trunc));
    TruncX=real(data_US_GIRAF(opts.xl_trunc:opts.xh_trunc,opts.yl_trunc:opts.yh_trunc));
    TruncX_fft=real(data_US_fourier(opts.xl_trunc:opts.xh_trunc,opts.yl_trunc:opts.yh_trunc));
    
    % set scale
    Xtest = real(X0(opts.xlow:opts.xhigh,opts.ylow:opts.yhigh));
    maxVal = max(Xtest(:));
    minVal = min(Xtest(:));
    
    
    [a1,b1]=RMSE(TruncX_fft,TruncX0);
    [a2,b2]=RMSE(TruncX,TruncX0);
    
    subplot(2,3,5);
    [C,H]=contour(real(TruncX),12);xlabel(sprintf('RMSE=%f',a2));title('contour plot of current solution');drawnow;
    caxis([minVal,maxVal]);
    set (H, 'LineWidth', 3);
    subplot(2,3,6);
    [C,H]=contour(real(TruncX_fft),12);xlabel(sprintf('RMSE=%f',a1));title('contour plot of fourier recon');drawnow;
    caxis([minVal,maxVal]);
    set (H, 'LineWidth', 3);
    subplot(2,3,4);
    [C,H]=contour(real(TruncX0),12);title('contour plot of true data');drawnow;
    caxis([minVal,maxVal]);
    set (H, 'LineWidth', 3);
    %
    
    %figure(16); imagesc(fftshift(log(1+abs(x)))); colorbar; title('current solution, time-domain'); drawnow;
    
    RMSE_val = [RMSE_val,a2];
    % plot RMSE
    figure(21);plot(RMSE_val); title('RMSE'); drawnow;
    
    
    % Save recon at all iterations
    
    X_rec(:,:,i) = X;
    
    
    if (i>1 && (cost(i)-cost(i-1))>0)
        break;
    end
    
end