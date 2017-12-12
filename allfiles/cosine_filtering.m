function [data_reconstructed,Nt3,N,Nt] = cosine_filtering(data,opts)

%Cosine window filtering and zero padding

dd = data(opts.w3Trunc_low:opts.w3Trunc_high,:);
dd_ifft = ifft(dd,[],1);
[Nt3,N,Nt] = size(dd_ifft);% dimensions of data
cosD = ones(size(dd_ifft));
temp = cos((floor(N/2)+1:N)*pi/N - pi/2);
temp= repmat(temp,[Nt3,1,Nt]);
cosD(:,floor(N/2)+1:N,:) = temp;
data_cos = zeros(Nt3,opts.Nf2,Nt);
data_cos(:,1:N,:) = dd_ifft.*cosD;
data_cos = opts.magScale*data_cos; %Scale the data to an appropriate factor
data_reconstructed = fft2(data_cos);
