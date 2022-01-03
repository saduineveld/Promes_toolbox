function [GRID] = prepgrid(nn,lb,ub,algo,algo_spec)
% [GRID] = prepgrid(nn,lb,ub,algo,algo_spec)
% Set all necessary properties of GRID
%
%   INPUT:
%       nn          = scalar of number of state variables
%       lb          = vector of lower bounds (1 x nn)
%       ub          = vector of upper bounds (1 x nn)
%
%       algo    = algorithm. Options:
%                       'cheb_gal','cheb_tmi','cheb_mse' (grid_type='cheb')
%                       'spl_tmi','spl_dir' (grid_type = 'spline')  
%                       'smol_tmi','smol_dir' (grid_type = 'smolyak')
%                       'mono_mse' (grid_type = 'mono')
%
%       algo_spec   = (optional) structure with (algorithm specific) fields
%                           (default values per dimension in brackets):
%                           grid_type = 'cheb': qq (6), ord_vec (5)
%                           grid_type = 'spline': qq (5)
%                           grid_type = 'smolyak': mu_vec (2)
%                           grid_type = 'mono': qq (4), ord_vec (3)
%
%   OUTPUT:
%       GRID        = structure with fields:
%                       - nn,lb,ub
%                       - mm: total number of nodes
%                       - grid_type, xx
%                   	- qq, gridVecs (all algo except 'smolyak')        
%                       - ord_vec, poly_elem (only for 'cheb','mono')
%                       For grid_type 'mono': 
%                           - XX_poly 
%                       For grid_type 'cheb'
%                           - gridVecs_dw, xx_dw, XX_poly_dw
%                       For grid_type 'smolyak':
%                           - xx_dw, XX_poly_dw, inv_XX_poly_dw, mu_vec, mu_max, 
%                               smol_el_ani 

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

% For the Smolyak grid this code uses:
%   Rafa Valero (2021). Smolyak Anisotropic Grid (https://www.mathworks.com/matlabcentral/fileexchange/50963-smolyak-anisotropic-grid), 
%       MATLAB Central File Exchange. Retrieved December 9, 2021. 
%
% Relevant original code and license agreement can be found in the folder
%   "smolyak_subfun".
%
% The accompanying article is:
%   "Smolyak method for solving dynamic economic models: 
%   Lagrange interpolation, anisotropic grid and adaptive domain" 
%   by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero (2014). 
%   Journal of Economic Dynamics and Control. Volume 44, July 2014, Pages 92–123. 

if nargin < 5
    algo_spec = [];
end

% Assign default values is algo_spec is not set
if strcmp(algo,'cheb_mse') || strcmp(algo,'cheb_gal') || strcmp(algo,'cheb_tmi')
    grid_type = 'cheb';
    if isempty(algo_spec)
        algo_spec.ord_vec   = 5*ones(1,nn);
        algo_spec.qq        = algo_spec.ord_vec + 1;
    end
elseif strcmp(algo,'mono_mse') 
    grid_type = 'mono';    
    if isempty(algo_spec)
        algo_spec.ord_vec   = 3*ones(1,nn);
        algo_spec.qq        = algo_spec.ord_vec + 1;
    end
elseif strcmp(algo,'spl_tmi') || strcmp(algo,'spl_dir')
    grid_type = 'spline';
    if isempty(algo_spec)
        algo_spec.qq        = 5*ones(1,nn);
    end
elseif strcmp(algo,'smol_tmi') || strcmp(algo,'smol_dir')
    grid_type = 'smolyak';
    if isempty(algo_spec)
        algo_spec.mu_vec = 2*ones(1,nn);
    end
else
    error('Invalid algorithm');
end


%% Construct grid
if strcmp(grid_type,'smolyak')
    GRID = gridstruct_smolyak(nn,lb,ub,algo_spec.mu_vec);
elseif strcmp(grid_type,'mono') || strcmp(grid_type,'cheb')
    GRID = gridstruct(nn,algo_spec.qq,lb,ub,grid_type,algo_spec.ord_vec);
elseif strcmp(grid_type,'spline')
    GRID = gridstruct(nn,algo_spec.qq,lb,ub,grid_type,[]);
end


end