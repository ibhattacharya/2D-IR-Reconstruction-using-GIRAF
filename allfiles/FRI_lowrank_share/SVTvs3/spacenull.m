
clear all; close all; clc;
k = 9; r = 6; N = 50;

U = zeros(128,r);
index = randperm(128);

U(index(1:k),:) = randn(k,r);

X = U*randn(r,N); 
Xhat = fft(X,[],1); 
%figure;imagesc(abs(Xhat));

Np = k+1;
B = giveConcatHankMtx(Xhat,Np);

nx = null(Xhat);

% from null(X) or vij kron I
nV = zeros(size(nx,1)*(k+1),size(nx,2)*(k+1));
for i=1:size(nx,1)
    for j=1:size(nx,2)
        nV((i-1)*(k+1)+1:i*(k+1),(j-1)*(k+1)+1:j*(k+1)) = nx(i,j)*eye(k+1);
    end
end
[rank(Xhat), rank(B)]
max(max(B*nV))
size(nx)
size(B*nV)
r*(k+1)-k
size(null(B))

t1n = null(B(:,1:(k+1)));
qq = [t1n;zeros((k+1)*(N-1),1)];

for i=1:N
   Q(:,i) =  circshift(qq,[(i-1)*(k+1),1]);
end
[size(Q);size(nV); size([Q,nV])]
[rank(Q),rank(nV),rank([Q,nV]), size(null([Q,nV]))]
figure;imagesc(abs([Q,nV]));
return


%% vary r and k
close all; clear all; clc;
k = 3; N = 10; m = 15; n = 0:m-1; r = 1; w = [pi/3, pi/4, pi/6]; X = zeros(m,N);
M = randn(k,r)*randn(r,N); 
% fc = randn(1,N);  be = 0.8;
% M = [fc*0.2; fc; randn(1,N); fc*be; randn(1,N)];
for j = 1:N
    for i = 1:length(w)
        X(:,j) = X(:,j) + (M(i,j)*exp(-w(i)*n*1i)).';
    end
    
end
T = giveConcatHankMtx(X,k+1);
[rank(X), rank(T), rank(M)]
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
3*(k+1)-k
t1n = null(T(:,1:(k+1)));
qq = [t1n;zeros((k+1)*(N-1),1)];
 

for i=1:N
   Q(:,i) =  circshift(qq,[(i-1)*(k+1),1]);
end
[size(Q);size(nV); size([Q,nV])]
[rank(Q),rank(nV),rank([Q,nV]), size(null([Q,nV]))]
figure;imagesc(abs([Q,nV]));
%figure(1);imagesc(abs(n20));figure(2);imagesc(abs(T*n20));

% % from full nulling filter
% c1 = exp(-w(1)*1i); c2 = exp(-w(2)*1i); c3 = exp(-w(3)*1i); c4 = exp(-w(4)*1i);
% h = [1, -(c1+c2+c3+c4), c1*c2+c2*c3+c3*c1+c1*c4+c2*c4+c3*c4, -(c1*c2*c3+c1*c2*c4+c2*c3*c4+c3*c1*c4), c1*c2*c3*c4];
% for i=1:N
%     y = filter(h,1,X(:,i));
%     figure(3); stem(n,y,'filled','sr')
%     axis([0 50 -1.1 1.1]); grid;
% end
















return
%% 3 expo 2 snaps k = 3 N = 2
close all; clear all; clc;
m = 50; n = 0:m-1; k = 3; N = 40;
M = randn(k,N); w = [pi/3, pi/4, pi/6]; X = zeros(m,N);
for j = 1:N
    for i = 1:length(w)
        X(:,j) = X(:,j) + (M(i,j)*exp(-w(i)*n*1i)).';
    end
end

c1 = exp(-w(1)*1i); c2 = exp(-w(2)*1i); c3 = exp(-w(3)*1i);
h = [1, -(c1+c2+c3), c1*c2+c2*c3+c3*c1, -c1*c2*c3];
h1 = [1, -exp(-w(1)*1i)]; h2 = [1, -exp(-w(2)*1i)]; h3 = [1, -exp(-w(3)*1i)];

y = zeros(1,m);
for i=1:N
    t = circshift(1:N,[2 -1]);
    if(i==1)
        ff = flipud(M(2,t(i)));
    else
        ff = flipud(-M(2,t(i)));
    end
    y = y + filter(h1,1,X(:,i)*ff);
    figure(1); stem(n,y,'filled','sr')
    axis([0 50 -1.1 1.1]); grid;
end
clear t i
y = zeros(1,m);
for i=1:N
    t = circshift(1:N,[2 -1]);
    if(i==1)
        ff = flipud(M(3,t(i)));
    else
        ff = flipud(-M(3,t(i)));
    end
    y = y + filter(h2*ff,1,X(:,i));
    figure(2); stem(n,y,'filled','sr')
    axis([0 50 -1.1 1.1]); grid;
end
clear t i
y = zeros(1,m);
for i=1:N
    t = circshift(1:N,[2 -1]);
    if(i==1)
        ff = flipud(M(1,t(i)));
    else
        ff = flipud(-M(1,t(i)));
    end
    y = y + filter(h3*ff,1,X(:,i));
    figure(3); stem(n,y,'filled','sr')
    axis([0 50 -1.1 1.1]); grid;
end
clear t y i
for i=1:N
    y = filter(h,1,X(:,i));
    figure(4); stem(n,y,'filled','sr')
    axis([0 50 -1.1 1.1]); grid;
end

return











%% 3 expo 2 snaps k = 3 N = 2
close all; clear all; clc;
n = 0:80; k = 3; a = 2; b = -3; c = -4; d = 1; e = 1; f = -2;
M = [a b c;d e f]
w1 = pi/3; w2 = pi/4; w3 = pi/6;
x1 = a*exp(-w1*n*1i) + b*exp(-w2*n*1i) + c*exp(-w3*n*1i);
x2 = d*exp(-w1*n*1i) + e*exp(-w2*n*1i) + f*exp(-w3*n*1i);
c1 = exp(-w1*1i); c2 = exp(-w2*1i); c3 = exp(-w3*1i);
h = [1, -(c1+c2+c3), c1*c2+c2*c3+c3*c1, -c1*c2*c3];
h1 = [1, -exp(-w1*1i)]; h2 = [1, -exp(-w2*1i)]; h3 = [1, -exp(-w3*1i)];

y11 = filter(h1*f,1,x1);
y12 = filter(h1*c,1,x2);

y21 = filter(h2*d,1,x1);
y22 = filter(h2*a,1,x2);

y31 = filter(h3*e,1,x1);
y32 = filter(h3*b,1,x2);
y = filter(h,1,x1);

figure(1); hold on;
%stem(n,x1,'filled','b')
%stem(n,y,'filled','sr')
stem(n,y11-y12,'filled','g')
stem(n,y21-y22,'filled','*m')
stem(n,y31-y32,'filled','sk')
axis([0 50 -1.1 1.1]); grid;
return
%% 3 expo 2 snaps k = 2 N = 2
close all; clear all; clc;
n = 0:80; a = 1; b = 2; c = -2; d = 4; M = [a,b;c,d]'; k = 2;
w0 = pi/3; w1 = pi/4;
x1 = a*exp(-w0*n*1i) + b*exp(-w1*n*1i) ;
x2 = c*exp(-w0*n*1i) + d*exp(-w1*n*1i) ;

h = [1, -exp(-w0*1i)-exp(-w1*1i), exp(-(w0+w1)*1i)];
h1 = [1, -exp(-w0*1i)]; h2 = [1, -exp(-w1*1i)];

y11 = filter(h1*d,1,x1);
y12 = filter(h1*b,1,x2);

y21 = filter(h2*c,1,x1);
y22 = filter(h2*a,1,x2);

y = filter(h,1,x2);

figure(1); hold on;
stem(n,x1,'filled','b')
stem(n,y,'filled','sr')
stem(n,y11-y12,'filled','g')
stem(n,y21-y22,'filled','*m')
axis([0 50 -1.1 1.1]); grid;
return
M
[x1.',x2.']*null(M)
yv = [hankmtx(x1,k),hankmtx(x2,k)]*(kron(eye(k),null(M)));

