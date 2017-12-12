function out = hankmtx(a,N)
   a = a(:);
   m = length(a)-N+1;
   out = zeros(m,N);
   for i=1:m,
       out(i,:) = (a(i:i+N-1));
   end