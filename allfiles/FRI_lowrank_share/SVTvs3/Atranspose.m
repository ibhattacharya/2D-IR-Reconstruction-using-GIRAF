
function x = Atranspose(b,sampset,m,n)

x = zeros(m,n);
x(sampset) = b;
