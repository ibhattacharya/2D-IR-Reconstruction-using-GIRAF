function [x]= abwd(b,sampset,m,n)
b1 = zeros(m,n);
b1(sampset) = b;
x = ifft(b1,[],1);
end
