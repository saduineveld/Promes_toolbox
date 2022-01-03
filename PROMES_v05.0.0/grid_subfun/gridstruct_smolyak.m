function [GRID] = gridstruct_smolyak(nn,lb,ub,mu_vec)
%[GRID] = grid_struct_smolyak(nn,lb,ub,mu_vec)
%
% Set all necessary properties of GRID for Smolyak gridtype
%   
% This function is completely based on:
%       Rafa Valero (2021). Smolyak Anisotropic Grid (https://www.mathworks.com/matlabcentral/fileexchange/50963-smolyak-anisotropic-grid), 
%           MATLAB Central File Exchange. Retrieved December 9, 2021. 
%
% Relevant original code and license agreement can be found in the folder
% "smolyak_subfun".
%
% The accompanying article is:
%   "Smolyak method for solving dynamic economic models: 
%   Lagrange interpolation, anisotropic grid and adaptive domain" 
%   by Kenneth L. Judd, Lilia Maliar, Serguei Maliar and Rafael Valero (2014). 
%   Journal of Economic Dynamics and Control. Volume 44, July 2014, Pages 92–123. 
%
%
%   INPUT:
%           nn          = scalar of number of state variables
%           lb          = vector of lower bounds (1 x nn)
%         	ub          = vector of upper bounds (1 x nn)          
%           mu_vec      = vector of maximum 'accuracy' (mu) (1 x nn)
%
%   OUTPUT:
%           GRID        = structure with fields
%                           - nn, lb, ub
%                           - grid_type ('smolyak')
%                           - mm: total number of nodes
%                           - xx,  xx_dw, XX_poly_dw, inv_XX_poly_dw
%                           - mu_vec, mu_max, smol_el_ani


% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.


%% TEST IF INPUTS ARE CORRECTLY SPECIFIED
if mod(nn,1) ~= 0
    error('nn should be an integer');
end
if size(mu_vec,2) ~= nn || size(mu_vec,1) ~= 1
    error('mu_vec should be 1 x nn vector');
end
for in = 1:nn
    if mod(mu_vec(in),1) ~= 0
        error('mu_vec should consist of integers');
    end
end
if size(lb,2) ~= nn || size(lb,1) ~= 1
    error('lb should be 1 x nn vector');
end
if size(ub,2) ~= nn || size(ub,1) ~= 1
    error('ub should be 1 x nn vector');
end

GRID.grid_type  = 'smolyak';
GRID.nn         = nn;
GRID.lb         = lb;
GRID.ub         = ub;  
GRID.mu_vec     = mu_vec;

GRID.mu_max  = max(mu_vec);%max of level of approximation

% Construct the matrix of indices of multidimesional Smolyak elements 
% (grid points and polynomial basis functions) that satisfy the usual 
% isotropic Smolyak rule for the approximation level equal to "mu_max":
smol_el_iso = Smolyak_Elem_Isotrop(nn,GRID.mu_max);

% Select from the matrix of isotropic indices "Smol elem_iso" a subset of 
% indices that correspond to the given anisotropic "vector_mus_dimensions":
GRID.smol_el_ani = Smolyak_Elem_Anisotrop(smol_el_iso,mu_vec);

% Construct the Smolyak grid for the given subindices of anisotropic  
% Smolyak elements "Smol_elem_ani":
GRID.xx_dw  = Smolyak_Grid(nn,GRID.mu_max,GRID.smol_el_ani);

GRID.xx     = scal_mat_up(GRID.xx_dw,GRID.lb,GRID.ub);

GRID.mm     = size(GRID.xx_dw,1);  

% Matrix of the polynomial basis functions evaluated in the grid points:
GRID.XX_poly_dw = Smolyak_Polynomial(GRID.xx_dw,nn,GRID.mu_max,GRID.smol_el_ani);

GRID.inv_XX_poly_dw = GRID.XX_poly_dw\eye(GRID.mm);

end

