function out = getTms(x,f,res,W)

% res is the image size ; [res(1) res(2)] size of each snapshot
% prod(f) is the flter size ; [f(1) f(2)] filter
% N snapshots
% x is the full signal ; prod(res) x N
% W is the big weight martix ; prod(f).N x prod(f)
N = 1;
out = zeros(size(x));xh = zeros(res(1));
for jj=1:N
    Wh = zeros(res(1));
    for ii = 1:prod(f)
        Wh = Wh + (ifft2(padarray(reshape(W((jj-1)*prod(f)+1:jj*prod(f),ii),f(1),f(2)),[(res(1)-f(1))/2 (res(1)-f(2))/2],'both'))).^2;
    end
    xh = xh + (fft2(reshape(x(:,jj),res(1),res(2))).*Wh);
    out(:,jj) = reshape(ifft2(xh),prod(res),1);
end

end

