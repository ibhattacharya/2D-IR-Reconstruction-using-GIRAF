clear all; close all; clc;
set(0,'DefaultFigureWindowStyle','docked');
k = 4; r = 2; Ns = 11;

X1 = randn(128,r)*randn(r,Ns); X = zeros(128,Ns);
for i =1:Ns
    index = randperm(128);
    X (index(1:k),i) = X1(index(1:k),i);
end

Xhat = fft(X,[],1); % in k space
%figure;imagesc([X,X1,abs(Xhat)]);
[rank(X1),rank(X),rank(Xhat)]

% close all; clear all; clc;
% k = 5; N = 11; m = 14; n = 0:m-1; r = 1; w = [pi/3, pi/4, pi/6,pi/5,pi/8]; 
% X = zeros(m,N); M = randn(k,r)*randn(r,N); 
% for j = 1:N
%     for i = 1:length(w)
%         X(:,j) = X(:,j) + (M(i,j)*exp(-w(i)*n*1i)).';
%     end
% end
T = giveConcatHankMtx(X,k+1);
rank(T)

[~,c1,~]=svd(X,'econ');figure(4);plot(diag(c1));
[~,c2,~]=svd(Xhat,'econ');figure(5);plot(diag(c2));
[~,cc,~]=svd(T,'econ');figure(6);plot(diag(cc));
return
[rank(X)-r, rank(T)-k, rank(M)-r]
nx = null(X);

% from null(X) or vij kron I
nV = zeros(size(nx,1)*(k+1),size(nx,2)*(k+1));
for i=1:size(nx,1)
    for j=1:size(nx,2)
        nV((i-1)*(k+1)+1:i*(k+1),(j-1)*(k+1)+1:j*(k+1)) = nx(i,j)*eye(k+1);
    end
end
max(max(T*nV))
size(nx)
size(T*nV)
size(null(T))
r*(k+1)-k
t1n = null(T(:,1:(k+1)));

qq = [t1n;zeros((k+1)*(N-1),1)];
 
for i=1:N
   Q(:,i) =  circshift(qq,[(i-1)*(k+1),1]);
end
% max(max(T*[Q,nV]))
% figure;imagesc(abs(T*[Q,nV]));
[size(Q);size(nV); size([Q,nV])]
[rank(Q)-N,rank(nV)-((N-r)*(k+1)),rank([Q,nV]), size(null([Q,nV]),2)-(N-r)]
return
% qnv = [Q,nV];
% for i = 1:size(qnv,1)
%     for j = 1:size(qnv,1)
% q1 = (qnv(33,:));k1 = q1((q1~=0));
% q2 = (qnv(34,:));k2 = q2((q2~=0));
% qn = null(qnv); 
% cc = ([ k1.',k2.',qn(q1~=0,1), qn(q2~=0,1)]);
% ca = real(cc);
% [cc(:,1).'*cc(:,3),
% cc(:,2).'*cc(:,4),
% ca(:,1).'*ca(:,3),
% ca(:,2).'*ca(:,4)]
%     end
% end
ttt = (null([Q,nV]));

% for i = 1:N
%     a1 = randn(N-r,1);a2 = repmat(a1,k+1,1); a3 = sort(a1,'ascend');
%     a4 = nx(i,:)*a1;
%     a5 = a4*(-1/t1n(1));
%     
figure;imagesc(ttt);