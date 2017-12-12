function out = defQhQsb( res, snap, x, f,dz)
x = reshape(x,prod(res),snap);
G = zeros(prod(res),snap);

% for jj=1:snap
%     filtered = zeros(res);
%     for kk = 1:length(dz)
%          g = zeros(res);
%         for ii=1:snap
%             fpad = (reshape(f(:,snap*(jj-1)+ii),res)); 
%             g = g + fpad.*ifft2(dz{kk}.*(reshape(x(:,ii),res))); 
%         end
%         filtered = filtered + (conj(dz{kk}).*fft2(g)); 
%     end
%     G(:,jj) = filtered(:); 
% end
% out = G(:); %out = reshape(G,prod(res)*snap,1); 

%%
for jj=1:snap
    g = zeros(res);
    for ii=1:snap
        filtered = zeros(res);
        for kk = 1:length(dz)
            if(jj==ii)
                fpad = (reshape(f(:,snap*(jj-1)+ii),res));
                filtered = filtered + conj(dz{kk}).*(fft2(fpad.*ifft2(dz{kk}.*(reshape(x(:,ii),res)))));
            else
                fpad = (reshape(f(:,snap*(jj-1)+ii),res));
                filtered = filtered + conj(dz{kk}).*(fft2(fpad.*ifft2(dz{kk}.*(reshape(x(:,ii),res))))+(fft2(fpad.*ifft2(dz{kk}.*(reshape(x(:,ii),res)))))');
            end
        end
        g = g + filtered;
    end
    G(:,jj) = g(:);
end
out = G(:); %out = reshape(G,prod(res)*snap,1);



end

