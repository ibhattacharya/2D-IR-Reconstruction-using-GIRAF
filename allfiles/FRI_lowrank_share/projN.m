function out = projN(x,ind,snap,res)
x = reshape(x,prod(res),snap);
out = zeros(prod(res),snap);
for ii=1:snap
    y = zeros(res);
    y(ind) = x(ind,ii);
    out(:,ii) = y(:);
end
out = out(:);


% %%
% x = reshape(x,prod(res),snap); 
% out = zeros(prod(res),snap);
% for ii=1:snap
%     y = zeros(res);
%     y(ind) = x(ind,ii);
%     out(:,ii) = y(:);
% end
% out = sum(out,2); out = out(:);


end