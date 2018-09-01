function fluid = getOGSatResiduals(fluid, SWOF,SGOF )
  [Swr, Sor] = getSatResiduals(SWOF);
  [Sgr, Sorg] = getSatResiduals(SGOF);  
  fluid.sr = [Swr, Sor, Sgr, Sorg]; 
end

function [s2r, s3r] = getSatResiduals(relPerm)
%{
% Copyright 2017 Zhonglin Yang, China University of Petroleum, Beijing
%  Email: Yzhonglin_edu@163.com
%
% SYNOPSIS:
%     [s2r, s3r] = getSatResiduals(relPerm)
%
% RETURNS:
%     s2r  -  The first non-zero saturation of colume 2 in relative
%               permeability table
%     s3r  -  The first non-zero saturation of colume 3 in relative
%               permeability table
%}
   s = relPerm(:, 1);
   s2r = s(findIndex(relPerm(:, 2)));
   s3r = s(findIndex(relPerm(:, 3)));
end

function index = findIndex(vector)
       assert(all(vector>=0), 'Saturation cannot be negative');
      for i = 2 : numel(vector)
          if  (vector(i-1) == 0 & vector(i) > 0) + (vector(i-1) > 0 & vector(i) == 0) == 1
              if vector(i-1) == 0
                  index = i-1;
              else
                  index =i;
              end
              return
          end
      end     
end