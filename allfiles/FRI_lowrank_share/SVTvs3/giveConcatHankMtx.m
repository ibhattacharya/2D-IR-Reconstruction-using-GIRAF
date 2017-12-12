function B = giveConcatHankMtx(Xhat,Np)

Ns = size(Xhat,2); B = [];
%%%
% for i = 1:Ns
%     [a, b, c, d] = dwt2(reshape(Xhat(:,i),128,128),'Haar');
%     ee = reshape(fft2([a,b;c,d]),128*128,1);
%     B = [B,hankmtx(ee,Np)];
% end

%%%
for i = 1:Ns,
    B = [B,hankmtx(Xhat(:,i),Np)];
end
end
