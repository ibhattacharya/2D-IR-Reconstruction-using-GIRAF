function s = imstats(x0,x1)
%x0 ground truth
%x1 comparison image
%assumes images scaled to [0,1]
    s.MSE = (1/length(x0(:)))*norm(x1(:)-x0(:)).^2;
    s.SNR = 20*log10(norm(x0(:))/norm(x1(:)-x0(:)));
    s.PSNR = -10*log10(s.MSE);
    %s.SSIM = ssim(im2uint8(abs(x0)),im2uint8(abs(x1)));
    s.relerr_inf = norm(x1(:)-x0(:),Inf)/norm(x0(:),Inf);
    s.relerr_2 = norm(x1(:)-x0(:))/norm(x0(:));
end