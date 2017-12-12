function x = Abackward(b,sampset,m,n)

x = zeros(m,n);
x(sampset) = b;
%x = ifft(x,[],1);
for i=1:n
t = ifft(reshape(x(:,i),128,128)); t = reshape(t,m,1);
x(:,i) = t;
end
end
