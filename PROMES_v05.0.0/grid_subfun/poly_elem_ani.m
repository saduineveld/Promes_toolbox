function [poly_elem,poly_iso] = poly_elem_ani(nn,ord_vec)
% [poly_elem,poly_iso] = poly_elem_ani(nn,ord_vec)
%
% Constructs a matrix poly_elem which is used to construct multivariate 
% complete polynomials (see function get_poly_ani)
%
% INPUT:
%   nn          = number of dimensions (or state variables)
%   ord_vec     = order of the polynomial in each dimension (1 x nn vector)
%
% OUTPUT:
%  poly_elem  = matrix (pp x nn) to construct multivariate complete polynomials:
%                  - each column refers to a variable
%                  - index number is the degree of the polynomial (either 
%                       Chebyshev or monomial)
%
% Example:
%       poly_elem = [0,2,1;3,2,4]; (3 dimensions, 2 polynomial terms)
%
% With monomials this poly_elem constructs the polynomial terms:
%  Phi = [x1^0*x2^2*x3^1,x1^3*x2^2*x3^4]; dim: mm x 2 (each xi is mm x 1)


% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

max_order = max(ord_vec);

ss=zeros(max_order+1,1);
poly_iso = [];
for j=0:max_order
    if j==0
        Q_n = zeros(1,nn);   
    else
        Q = Q_n;
        Q_n = [];
        for q=nn:-1:1%loop over state variables
            rr       = comb_pr(q+j-2,q-1);
            tt=[zeros(rr,nn-q) ones(rr,1) zeros(rr,q-1)]+[zeros(rr,nn-q) Q(ss(j)-rr+1:ss(j),nn-q+1:nn)];
            Q_n = cat(1,Q_n,tt);
        end        
    end
    ss(j+1) = size(Q_n,1);
    poly_iso = cat(1,poly_iso,Q_n);   
end  

poly_elem = poly_iso;
ind_drop=[];
if size(ord_vec,2) > 1%if ord_vec is scalar then use isotropic polynomial
    if ~(all(ord_vec(:)==ord_vec(1)))%is ord_vec anisotropic?
        %Drop all rows which have ord_ind(:,in) > ord_vec(in)
        for in = 1:nn
            ind_drop =[ind_drop poly_iso(:,in) > ord_vec(in)];
            
            %poly_elem(ind_drop,:) = [];%deletes row
        end 
        ind_drop=sum(ind_drop,2);
        ind_drop=ind_drop>0;
        poly_elem(ind_drop,:) = [];%deletes rows
    end 
end    
    
end

function [a]=comb_pr(b,c)
%function [a]=comb(b,c);

a=factorial(b)/(factorial(c)*factorial(b-c));

end
