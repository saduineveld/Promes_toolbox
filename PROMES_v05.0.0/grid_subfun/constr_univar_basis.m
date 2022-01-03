function [PHI] = constr_univar_basis(xx,order,poly_type)
%[PHI] = constr_univar_basis(xx,order,poly_type)
% 
% Constructs univariate polynomial, excluding degree 0 (ie. constant)
%            (as in Judd, 1998, p. 239)
%
%   Input:
%           xx          = mm x 1 matrix with mm datapoints
%           order       = order of polynomial
%           poly_type 	= 'cheb' for Chebyshev, or 'mono' for monomial
%
%   Output:
%           PHI = matrix of polynomial (excl. degree 0). Dimensions:
%                   - rows:     mm
%                   - columns:  order

% Sijmen Duineveld, December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

if size(xx,2) > 1
    error('xx should be column vector')
end
mm = size(xx,1);

PHI = NaN(mm,order);
PHI(:,1) = xx;
if order > 1
    %Add quadratic term:
    if strcmp(poly_type,'mono')
        PHI(:,2) = xx.^2;
	elseif strcmp(poly_type,'cheb')
        PHI(:,2) = 2*xx.^2 - 1;
    end
    %Add orders 3 and higher
    if order > 2
        for io = 3:order
            if strcmp(poly_type,'mono')
                PHI(:,io) = PHI(:,io-1).*xx;
            elseif strcmp(poly_type,'cheb')
                PHI(:,io) = 2*PHI(:,io-1).*xx - PHI(:,io-2);
            end
        end
    end

end

end