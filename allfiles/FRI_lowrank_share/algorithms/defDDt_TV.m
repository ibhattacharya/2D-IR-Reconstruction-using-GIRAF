function [D,Dt] = defDDt_TV(siz)
    % defines derivative operator D
    % and its transpose operator Dt

    D = @(X) D_TV(X,siz);
    Dt = @(Y) Dt_TV(Y,siz);
    
    %define finite difference filters
    dx = zeros(siz,'double'); dx(1,1) = 1; dx(1,end) = -1;
    dy = zeros(siz,'double'); dy(1,1) = 1; dy(end,1) = -1;
    DX = fft2(dx);
    DY = fft2(dy);

    function D = D_TV(X,siz)
        % Forward finite difference operator
        D = zeros([siz,2],'double');        
        D(:,:,1) = ifftn(DX.*fft2(X));
        D(:,:,2) = ifftn(DY.*fft2(X));    
    end

    function Dt = Dt_TV(Y,siz)
        % Transpose of the forward finite difference operator
        Dt = zeros(siz,'double');
        D1 = ifftn(conj(DX).*fft2(Y(:,:,1)));
        D2 = ifftn(conj(DY).*fft2(Y(:,:,2)));          
        Dt = D1+D2;
    end
        
end