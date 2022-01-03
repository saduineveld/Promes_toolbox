function [GRID] = gridstruct(nn,qq,lb,ub,grid_type,ord_vec)
%[GRID] = grid_structure(nn,qq,lb,ub,grid_type,ord_vec)
%
%   Set all necessary properties of GRID, given grid type (grid_type)
%
%   INPUT:
%           nn          = scalar of number of state variables
%           qq_vec      = vector of number of gridpoints (1 x nn)
%           lb          = vector of lower bounds (1 x nn)
%         	ub          = vector of upper bounds (1 x nn)          
%
%           grid_type 	= 'mono', 'cheb' or 'spline'
%
%           ord_vec 	= vector of maximum order (per dimension) of polynomial 
%                               for grid_type 'cheb' and 'mono'
%
%   OUTPUT:
%           GRID        = structure with fields:
%                           - nn,qq,lb,ub
%                           - mm: total number of nodes
%                           - grid_type ('cheb','mono' or 'spline'), gridVecs, xx
%                           - poly_elem (only for 'cheb','mono')
%                           For grid_type 'mono': 
%                               - XX_poly 
%                           For grid_type 'cheb':
%                               - gridVecs_dw, xx_dw, XX_poly_dw

% Sijmen Duineveld, updated for Promes_v05 December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.


%% TEST IF INPUTS ARE CORRECTLY SPECIFIED
if size(qq,2) ~= nn || size(qq,1) ~= 1
    error('qq should be 1 x nn vector');
end
if (size(ord_vec,2) ~= nn || size(ord_vec,1) ~= 1) && (strcmp(grid_type,'mono') || strcmp(grid_type,'cheb'))
    error('ord_vec should be 1 x nn vector');
end
if mod(nn,1) ~= 0
    error('nn should be an integer');
end
for in = 1:nn
    if mod(qq(in),1) ~= 0
        error('qq should consist of integers');
    end
end
if size(lb,2) ~= nn || size(lb,1) ~= 1
    error('lb should be 1 x nn vector');
end
if size(ub,2) ~= nn || size(ub,1) ~= 1
    error('ub should be 1 x nn vector');
end
if ~strcmp(grid_type,'mono') && ~strcmp(grid_type,'cheb') && ~strcmp(grid_type,'spline')
    error('Invalid grid_type');  
end

if nargin < 6 || isempty(ord_vec) && (strcmp(grid_type,'mono') || strcmp(grid_type,'cheb'))
    error('ord_vec is unspecified');
end

%Check if enough nodes:
if strcmp(grid_type,'cheb')
    for in = 1:nn
        if qq(in) <= ord_vec(in)
        error('Minimum number of grid points (in each dimension) has to be the order + 1.');
        % If this condition is not met the complete polynomial will contain a
        % column with zeros, because q-nodes are the roots of the qth order
        % polynomial. Therefore the number of nodes needs to be larger than
        % the order of the approximation (in every dimension).
        end
    end
end

%% Construct gridVecs (and gridVecs_dw for 'cheb')
if strcmp(grid_type,'spline') || strcmp(grid_type,'mono')
    %equidistant ('equi') and 'up':    
    GRID.gridVecs = constr_vecs(qq,'equi','up',lb,ub);
elseif strcmp(grid_type,'cheb')
    %chebyshev ('cheb') and 'dw':
    GRID.gridVecs_dw = constr_vecs(qq,'cheb','dw',lb,ub);
    % scaled up version of Chebyshev nodes:
    GRID.gridVecs = constr_vecs(qq,'cheb','up',lb,ub);
end

%% Construct full grid xx (and xx_dw for 'cheb')
%Grid of scaled up variables:
GRID.xx = constr_grid(GRID.gridVecs);
if strcmp(grid_type,'cheb')
    %Grid of scaled down variables:
    GRID.xx_dw = constr_grid(GRID.gridVecs_dw);
end

%% Construct polynomial for grid_type 'cheb' and 'mono'
if strcmp(grid_type,'cheb') || strcmp(grid_type,'mono')
    GRID.poly_elem = poly_elem_ani(nn,ord_vec);

    if strcmp(grid_type,'cheb')
        GRID.XX_poly_dw  = get_poly_ani(GRID.xx_dw,ord_vec,'cheb',GRID.poly_elem);
    elseif strcmp(grid_type,'mono')
        GRID.XX_poly     = get_poly_ani(GRID.xx,ord_vec,'mono',GRID.poly_elem);
    end
end

%% Assign properties of input fields:
GRID.nn = nn;
GRID.qq = qq;
GRID.mm = prod(qq);
GRID.lb = lb;
GRID.ub = ub;
GRID.grid_type = grid_type;
if strcmp(grid_type,'cheb') || strcmp(grid_type,'mono')
    GRID.ord_vec = ord_vec;    
end

end

