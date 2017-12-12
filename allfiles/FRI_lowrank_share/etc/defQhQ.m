function out = defQhQ( x, mu, dz, res, ind_full, truncate_flag )
    function y = truncate(x)
        LPFILTER = zeros(res);
        LPFILTER(ind_full) = 1;
        y = ifft2(fft2(x).*LPFILTER);
    end           

    x = reshape(x,res);
    filtered = zeros(res);
    if truncate_flag
        for i=1:length(dz)
            dx = ifft2((dz{i}.*x));%
            filtered = filtered + conj(dz).*fft2(conj(mu).*truncate(mu.*dx));%
        end
    else
        for i=1:length(dz)
            filtered = filtered + conj(dz{i}).*fft2(mu.*ifft2((dz{i}.*x)));%figure(5);imagesc(abs(mu.*ifft2((dz{i}.*x))));
        end            
    end
    figure(1);imagesc(abs(filtered));
    %out = filtered(ind_full);
    out = filtered(:); 
end

