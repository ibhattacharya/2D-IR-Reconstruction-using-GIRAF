
function l = l21nrm_cvx(Z)

l = 0;
for ii = 1: size(Z,1)
    l = l + norm(Z(ii,:));
end
end


