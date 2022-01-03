function [GRID] = prepgrid_fpa(nn,lb,ub,algo,algo_spec)
% [GRID] = prepgrid_fpa(nn,lb,ub,algo,algo_spec)
% Set all necessary properties of GRID


if nargin < 5
    algo_spec = [];
end

% Assign default values is algo_spec is not set
if strcmp(algo,'cheb_mse') || strcmp(algo,'cheb_gal') || strcmp(algo,'cheb_tmi') || strcmp(algo,'cheb_fpa')
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
elseif strcmp(algo,'spl_tmi') || strcmp(algo,'spl_dir') || strcmp(algo,'spl_fpa')
    grid_type = 'spline';
    if isempty(algo_spec)
        algo_spec.qq        = 5*ones(1,nn);
    end
elseif strcmp(algo,'smol_tmi') || strcmp(algo,'smol_dir') || strcmp(algo,'smol_fpa')
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