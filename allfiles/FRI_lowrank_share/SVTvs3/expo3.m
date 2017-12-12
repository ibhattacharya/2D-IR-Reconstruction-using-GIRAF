close all; clear all; clc;
n = 0:80; k = 3; a = 2; b = -3; c = -2; d = 1; e = 1; f = -1;
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
M
[x1.',x2.']*null(M)
yv = [hankmtx(x1,k),hankmtx(x2,k)]*(kron(eye(k),null(M)));
