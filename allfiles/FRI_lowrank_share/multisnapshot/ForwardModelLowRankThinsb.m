classdef ForwardModelLowRankThinsb
    properties
        dz; %derivative masks
        adjoint = 0;
        res;
        ndx;
        filter_siz;
        conv_siz;
        step_siz;
        scale;
    end
    
    methods
        
          function obj = ForwardModelLowRankThinsb(dz,res,filter_siz,step_siz)
             obj.dz = dz; % array of derivatives
             obj.ndx = 2;
             obj.res = res;
             obj.filter_siz = filter_siz;
             obj.conv_siz = [prod(filter_siz),prod(((res-filter_siz)./step_siz)+[1,1])];
             obj.step_siz = step_siz;
             obj.scale = sqrt(prod(step_siz)/prod(filter_siz));
          end
          
          % setting the adjoint flag
          function c = ctranspose(obj)
              obj.adjoint = 1;
              c = obj;
          end
              
          % multiplication by adjoint
          function out = mtimes(obj,x)
              if(obj.adjoint==1)
                    m = zeros(obj.res);
                    x = reshape(x,[obj.conv_siz,obj.ndx]);
                    for i=1:obj.ndx
                        T = x(:,:,i);
                        I = col2imstep(real(T),obj.res,obj.filter_siz,obj.step_siz);
                        I = I + 1i*col2imstep(imag(T),obj.res,obj.filter_siz,obj.step_siz);
                        I = ifftshift(I).*conj(obj.dz{i});
                        m = m + I;                        
                    end
                    out = obj.scale*m;
          % multiplication
              else
                  T = zeros([obj.conv_siz,obj.ndx]);
                  for i=1:obj.ndx
                      D = x.*obj.dz{i}; % {}
                      I = fftshift(D); %
%                                             T = im2colstep(real(I),obj.filter_siz,obj.step_siz);
%                                             T = T + 1i*im2colstep(imag(I),obj.filter_siz,obj.step_siz);
                      T(:,:,i) = im2colstep(real(I),obj.filter_siz,obj.step_siz);
                      T(:,:,i) = T(:,:,i) + 1i*im2colstep(imag(I),obj.filter_siz,obj.step_siz);
                      
                  end
                  T = reshape(T,[obj.conv_siz(1),obj.conv_siz(2)*obj.ndx]);
                  out = obj.scale*T;
              end
          end
    end
    
    
end

