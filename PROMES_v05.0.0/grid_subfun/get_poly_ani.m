function [XX_poly] = get_poly_ani(xx,ord_vec,poly_type,poly_elem)
% [XX_poly] = get_poly_ani(xx,ord_vec,poly_type,poly_elem)
%
% Constructs an anisotropic multivariate polynomial XX_poly
%
% INPUT: 
%       xx          = matrix (mm x nn) with mm gridpoints in nn dimensions
%       ord_vec     = maximum order in dimension nn (1 x nn)
%       poly_type   = polynomial type (either 'cheb' or 'mono')
%       poly_elem   = matrix (pp x nn) to construct anisotropic polynomial (see poly_elem_ani)
%
% OUTPUT:
%       XX_poly          = matrix (mm x pp) with anisotropic polynomial in nn dimensions

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.


[mm,nn]     = size(xx);

max_ord     = max(ord_vec);
if size(ord_vec,2) == 1
    ord_vec = ord_vec*ones(1,nn);    
end

%% CONSTRUCT UNIVARIATE BASES:
PHI = NaN(mm,max_ord,nn);
for in = 1:nn
   PHI(:,1:ord_vec(in),in) =  constr_univar_basis(xx(:,in),ord_vec(in),poly_type);    
end

%% CONSTRUCT POLYBASE:
pp = size(poly_elem,1);
XX_poly = NaN(mm,pp);
XX_poly(:,1) = 1;
for ip = 2:pp
    X_tmp = ones(mm,1);    
    for in = 1:nn
        if poly_elem(ip,in) > 0%Only multiply if ord > 0
            X_tmp = X_tmp .* PHI(:,poly_elem(ip,in),in);            
        end            
    end    
    XX_poly(:,ip) = X_tmp;
end

end