function [anewhat,QtQ] = giveAvgConcatHank (Bnew,Np)

NsNew = size(Bnew,2)/Np;
Mnew = size(Bnew,1)+Np-1;
anewhat = zeros(Mnew,NsNew);
if(nargout > 1)
    QtQ = anewhat;
end

for i=1:NsNew,
    if(nargout > 1)
        [anewhat(:,i),QtQ(:,i)] = avghank(Bnew(:,(i-1)*Np+1:i*Np));
    else
        %%%
%         t = avghank(Bnew(:,(i-1)*Np+1:i*Np));
%         t1 = (ifft2(reshape(t,128,128)));
%         t2 = idwt2(t1(1:end/2,1:end/2),t1(1:end/2,1+end/2:end),t1(1+end/2:end,1:end/2),t1(1+end/2:end,1+end/2:end),'Haar');
%         anewhat(:,i) = reshape(t2,128*128,1);
        %%%
        anewhat(:,i) = avghank(Bnew(:,(i-1)*Np+1:i*Np));
    end
end