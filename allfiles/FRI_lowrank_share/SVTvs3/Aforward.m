function xh = Aforward(x,sampset)
xh = fft(x,[],1);
xh = xh(sampset);
end