function [gridVecs] = constr_vecs(qq,nod_type,scale_type,lb,ub)
% function [gridVecs] = constr_vecs(qq,nod_type,scale_type,lb,ub)
%
% Construct cell array (1 x nn) of vectors of grid in each dimension
% (each vector is 1 x qq(ii))
%
% Input:
%       qq          = vector of gridpoints in each dimension       (1 x nn vector)%
%      
%       nod_type    = 'cheb' for Chebyshev or 'equi' for equidistant nodes
%            
%       scale_type  = (optional) either 'dw' (default for cheb') or 
%                       'up' (default for 'equi')
%       
%       Required for scale_type 'up': 
%       lb = vector of lower bounds in each dimension       (1 x nn vector)
%       ub = vector of higher bounds in each dimension  	(1 x nn vector)   
%
% Output:
%       gridvecs  = cell array (1 x nn) containing row vectors x(i) of length qq(i)
%
% Uses:
%       For method 'cheb':
%       - chebnodes, which generates the Cheybshev nodes in the interval
%           [-1,1]

% Sijmen Duineveld, Updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

nn=size(qq,2);

%check if input vectors are of same length
if nargin > 4
    if size(lb,2)~=size(ub,2);
        error('Vector lb and ub should be of same length')
    elseif length(lb)~=size(qq,2);
        error('Vector lb and qq should be of same length')
    end
end

if strcmp(scale_type,'up') && nargin < 5 
    error('For scale_type up both lb and ub need to be assigned')
end

gridVecs = cell(1,nn);
for ii=1:nn
    if nargin == 5
      if lb(ii)>ub(ii)
        error(['lb(',int2str(ii),') is larger than ub(',int2str(ii),')'])
      end
    end
    
    if strcmp(nod_type,'equi')    
        if strcmp(scale_type,'dw')
            xi = linspace(-1,1,qq(ii));
        elseif strcmp(scale_type,'up')
            xi = linspace(lb(ii),ub(ii),qq(ii));
        end
    elseif strcmp(nod_type,'cheb')
        xi = (chebnodes(qq(ii)))';
        if strcmp(scale_type,'up')
            xi = sc_cheb_up(lb(ii),ub(ii),xi);
        end            
    else
        error('Invalid type chosen');
    end   
    
    gridVecs{1,ii} = xi;
end

end