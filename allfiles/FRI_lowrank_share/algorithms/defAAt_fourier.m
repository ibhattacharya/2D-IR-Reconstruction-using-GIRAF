%Defines Fourier Undersampling operators A, At for Compressed Sensing MRI
%Inputs: acc = acceleration factor, siz = image size
%Outputs: A = measurement operator, At = measurement operator transpose
function operators = defAAt_fourier(options)
    mask = options.mask;
    siz = options.siz;
    
    S = find(mask~=0);
    %tim=para.siz(1)*para.siz(2)/length(S); %Undersampling factor

    A = @(z) A_fhp(z,S,siz);
    At = @(z) At_fhp(z,S,siz);

    function A = A_fhp(z, S, siz)
        A = zeros(siz,'double');
        p = 1/sqrt(siz(1)*siz(2))*fft2(z);
        A(S) = p(S);
        A = A(:);
    end

    function At = At_fhp(z, S, siz)
        z = reshape(z,siz);
        p = zeros(siz,'double');
        p(S) = z(S);
        At = sqrt(siz(1)*siz(2))*ifft2(p);
    end

    operators.A = A;
    operators.At = At;
    if(isfield(options,'mu'))
        operators.mu = options.mu;
    end
end