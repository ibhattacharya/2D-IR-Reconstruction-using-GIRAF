clear all; clc; close all;
addpath('./algorithms','./etc','./data','./multisnapshot'); load MR_SL_convert;
set(0,'DefaultFigureWindowStyle','docked');
res = [71,71];  %output resolution (use odd nums)
fsz = [19,19]; %filter dimensions (use odd nums)
overres = res + 2*fsz;

ind_trim = find((abs(k(1,:)) <= (overres(2)-1)/2 ) & (abs(k(2,:)) <= (overres(1)-1)/2));
m = reshape(m(ind_trim),overres);
k = k(:,ind_trim);

ind_full = find((abs(k(1,:)) <= (res(2)-1)/2) & (abs(k(2,:)) <= (res(1)-1)/2));
ind_filter = find((abs(k(1,:)) <= (fsz(2)-1)/2 ) & ( abs(k(2,:)) <= (fsz(1)-1)/2));
usf = 0.25; %undersampling factor
pdf = genPDF(res,7,usf); %density compensation
mask0 = im2double(genSampling(pdf,2,5)); mask0(1) = 1; %need DC component
mask = zeros(overres); mask(ind_full) = mask0;
scale = 1; snap = 2; %%%change
figure(1); imshow(fftshift(reshape(mask(ind_full),res)));
ind_samples = find(mask~=0);
[A,At] = defAAt(ind_samples,overres);
%b = zeros(size(m,1),size(m,2)*snap);

b1=A(m); b2=A(m.*scale);          %%%change
%dz = ones(overres);
dz{1} = reshape((-1i*pi*(k(1,:))).',overres)/(overres(1)/4);
dz{2} = reshape((-1i*pi*(k(2,:))).',overres)/(overres(2)/4);

Atb1 = At(b1); Atb2 = At(b2);      %%%change
x1 = Atb1; x2 = Atb2;     x = [x1(:), x2(:)];         %%%change

Q = ForwardModelLowRankThinsb(dz,overres,fsz,[1,1]);
QhQ = @(x,f) defQhQsb( overres, snap, x, f,dz);
AtA = @(x) projN(x,ind_samples,snap,overres);

%% run alg.
lambda = 1e6;  %regularization parameter %parameters for epsilon update: eps = min(eps,gamma*s(r+1));
eps = 0.1;       %epsilson initialization
gamma = 10;     %epsilon multiplier for update
r = 400;       %rank cut-off value
r2 = 0;        %singular value cut-off for mask recon
q = 2;         %q=1 is nuclear norm, q=2 is Shatten p=1/2
lambdainv = 1/lambda;
y = [lambdainv*Atb1 ; lambdainv*Atb2]; y = y(:);  %%%change

%%
Niter = 10; cost = []; epsrecord = [];
for i=1:Niter
    i
    % Stage 1: Compute annihilating mask
    [U,S,~] = svd([Q*x1;Q*x2],'econ');  %%%change
    s = diag(S); figure(2); plot(s); title('sing vals');
    
    %calculate cost
    nuclear = norm(s,1);
    diff =[(A(x1)-b1);(A(x2)-b2)];
    thiscost = 0.5*(norm(diff(:)).^2 )+ lambda*nuclear; %objective function
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
            loop = loop+1; %gg = zeros(overres);
            for i2=1:prod(fsz)*snap
%                   filt = zeros(overres); filtC = filt;
%                   filt(ind_filter)  = ifftshift(reshape(U(1+(i3-1)*prod(fsz):i3*prod(fsz),i2),fsz));
%                   filtC(ind_filter) = (ifftshift(reshape(U(1+(i4-1)*prod(fsz):i4*prod(fsz),i2),fsz)));
                 filt = padarray((reshape(U(1+(i3-1)*prod(fsz):i3*prod(fsz),i2),fsz)),[(overres(1)-fsz(1))/2 (overres(1)-fsz(2))/2],'both');
                filtC = padarray((reshape(U(1+(i4-1)*prod(fsz):i4*prod(fsz),i2),fsz)),[(overres(1)-fsz(1))/2 (overres(1)-fsz(2))/2],'both');
                
                filter = ((1/max(s(i2),eps))^q)*(conj(ifft2((filt))).*(ifft2((filtC))));
                %figure(11);imagesc([ifft2(filt) conj(ifft2(filtC))]);
                f(:,loop) = f(:,loop) + reshape(((filter)), prod(overres),1);
            end
            %f(:,loop) = gg;
        end
    end
    figure(5);imagesc(abs(reshape(f(:,1),overres)));
    figure(6);imagesc(abs(reshape(f(:,2),overres)));
    figure(7);imagesc(abs(reshape(f(:,3),overres)));
    figure(8);imagesc(abs(reshape(f(:,4),overres)));
    clear U loop i2 i3 i4
    % Stage 2: Weighted L2 problem
    R = @(x) QhQ(x,f) + lambdainv*AtA(x);
    [x,flag,relres] = pcg(R,y,1e-9,2000,[],[],x(:));
    
    x1 = reshape(x(1:end/2),overres);
    xtrim1 = reshape(x1(ind_full),res);
    X1 = ifft2(xtrim1);
    
    x2 = reshape(x(end/2+1:end),overres);
    xtrim2 = reshape(x2(ind_full),res);
    X2 = ifft2(xtrim2);
    figure(9); imagesc(abs([X1,X2])); colorbar; title('current solution');
    clear xtrim1 xtrim2 X1 X2
end

%% Plot Original Data
m1 = reshape(m(ind_full),res); m2 = reshape(m(ind_full).*scale,res);
figure(10); imagesc([abs(ifft2(m1)),abs(ifft2(m1.*scale))]); colorbar; title('original');