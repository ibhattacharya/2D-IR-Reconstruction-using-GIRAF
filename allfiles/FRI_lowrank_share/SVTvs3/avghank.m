function [a,Npoints] = avghank(B)
   [m,n] = size(B);
   M = m+n-1; 
   a = zeros(M,1);
   Npoints = zeros(size(a));
   for i=1:n,
       a(i:i+m-1) = a(i:i+m-1) + B(:,i);
       Npoints(i:i+m-1) = Npoints(i:i+m-1) + 1;
   end
   
   if(nargout == 1)
        a = a./Npoints;
   end
end