function [b] = afwd(x,sampset)
xh = fft(x,[],1);
b = xh(sampset);
end