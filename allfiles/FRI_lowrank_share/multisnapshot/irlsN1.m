clear all; clc; %close all;
addpath('./algorithms','./etc','./data','./multisnapshot'); load MR_SL_convert; load maskind;
set(0,'DefaultFigureWindowStyle','docked');
res = [91,91];  %output resolution (use odd nums)
filter_siz = [25,33]; %filter dimensions (use odd nums)
overres = res + 2*filter_siz;

ind_trim = find((abs(k(1,:)) <= (overres(2)-1)/2 ) & (abs(k(2,:)) <= (overres(1)-1)/2));
m = reshape(m(ind_trim),overres);
k = k(:,ind_trim);

ind_full = find((abs(k(1,:)) <= (res(2)-1)/2) & (abs(k(2,:)) <= (res(1)-1)/2));
ind_filter = find((abs(k(1,:)) <= (filter_siz(2)-1)/2 ) & ( abs(k(2,:)) <= (filter_siz(1)-1)/2));
% usf = 0.25; %undersampling factor
% pdf = genPDF(res,7,usf); %density compensation
% mask0 = im2double(genSampling(pdf,2,5)); mask0(1) = 1; %need DC component
% mask = zeros(overres); mask(ind_full) = mask0;
 scale = 1; snap = 1; %%%change
% %figure(1); imshow(fftshift(reshape(mask(ind_full),res)));
% ind_samples = find(mask~=0);

[A,At] = defAAt(ind_samples,overres);
%b = zeros(size(m,1),size(m,2)*snap);

b1=A(m); %b2=A(m.*scale);          %%%change
%dz = ones(overres);
dz{1} = reshape((-1i*pi*(k(1,:))).',overres)/(overres(1)/4);
dz{2} = reshape((-1i*pi*(k(2,:))).',overres)/(overres(2)/4);

Q = ForwardModelLowRankThinsb(dz,overres,filter_siz,[1,1]);
QhQ = @(x,f) defQhQsb( overres, snap, x, f,dz);
AtA = @(x) projN(x,ind_samples,snap,overres);

Atb1 = At(b1); %Atb2 = At(b2);      %%%change
x1 = Atb1; %x2 = Atb2;     x = [x1(:), x2(:)];         %%%change
x = x1(:);
%% run alg.
lambda = 1e5;  %regularization parameter %parameters for epsilon update: eps = min(eps,gamma*s(r+1));
eps = 0.1;       %epsilson initialization
gamma = 10;     %epsilon multiplier for update
r = 500;       %rank cut-off value
r2 = 0;        %singular value cut-off for mask recon
q = 2;         %q=1 is nuclear norm, q=2 is Shatten p=1/2
lambdainv = 1/lambda;
%y = [lambdainv*Atb1;lambdainv*Atb2]; y = y(:);  %%%change
y = lambdainv*Atb1; y = y(:);
%%
Niter = 10; cost = []; epsrecord = [];
for i=1:Niter
    i
    % Stage 1: Compute annihilating mask
   %[U,S,~] = svd([Q*x1;Q*x2],'econ');  %%%change
    [U,S,~] = svd(Q*x1,'econ');
    s = diag(S);   figure(2); plot(s); title('sing vals');
    
    %calculate cost
    nuclear = norm(s,1);    diff =A(x1)-b1; %
    %diff =[A(x1)-b1; A(x2)-b2];
    thiscost = 0.5*norm(diff(:)).^2 + lambda*nuclear; %objective function
    cost = [cost,thiscost]; clear S nuclear diff thiscost
    figure(3); plot(cost); title('cost');
    
    %update epsilon
    eps = min(eps,gamma*s(r+1));
    epsrecord(end+1) = eps;
    figure(4); plot(epsrecord); title('eps');
    
    %%% filters
    f = zeros(prod(overres),snap*snap); loop = 0; 
    for i3=1:snap
        for i4=1:snap
            loop = loop+1; mu = zeros(overres);
            for i2=1:prod(filter_siz)*snap
                filt = zeros(overres); %filtC = filt; 
                
                filt(ind_filter)  =  (ifftshift(reshape(U(1+(i3-1)*prod(filter_siz):i3*prod(filter_siz),i2),filter_siz)));
                %filtC(ind_filter) = ((ifftshift(reshape(U(1+(i4-1)*prod(filter_siz):i4*prod(filter_siz),i2),filter_siz))));
                
                mu = mu + ((1/max(s(i2),eps))^q)*( abs(ifft2(filt)).*abs(ifft2(filt)) );
                %figure(18);imagesc(abs(reshape(f(:,1),overres)));
            end
            f(:,loop) = reshape(((mu)),prod(overres),1);
        end
    end
    figure(5);imagesc(sqrt(abs(mu)));
%     figure;imagesc(abs(reshape(f(:,2),overres)));
%     figure;imagesc(abs(reshape(f(:,3),overres)));
%     figure;imagesc(abs(reshape(f(:,4),overres)));
    clear U loop i2 i3 i4
    % Stage 2: Weighted L2 problem
    R = @(x) QhQ(x,f) + lambdainv*AtA(x); % nature of x passed and returned check
    [x,flag,relres] = pcg(R,y,1e-9,2000,[],[],x(:));
    
    x1 = reshape(x(1:end),overres); %% /2
    xtrim1 = reshape(x1(ind_full),res);
    X1 = ifft2(xtrim1);
    
%     x2 = reshape(x(end/2+1:end),overres);
%     xtrim2 = reshape(x2(ind_full),res);
%     X2 = ifft2(xtrim2);
    figure(6); imagesc(abs([X1])); colorbar; title('current solution'); %%
    clear xtrim1 xtrim2 X1 X2
end

%% Plot Original Data
m1 = reshape(m(ind_full),res); m2 = reshape(m(ind_full).*scale,res);
%figure(7); imagesc([abs(ifft2(m1))]); colorbar; title('original'); %% ,abs(ifft2(m1.*scale))