function [AA,BB] = constr_mat_csd(MOD,DD,Rho)
%[AA,BB] = constr_mat_csd(MOD,DD,Rho)
%
% Construct numerical matrices

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner, 
% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

na = MOD.nx + MOD.nz + MOD.ny;

% Set up the linearized model
BB    = [-DD(:,1+na:2*na);[zeros(MOD.nz,MOD.nx), Rho, zeros(MOD.nz,MOD.ny)]];
AA    = [DD(:,1:na);[zeros(MOD.nz,MOD.nx),eye(MOD.nz),zeros(MOD.nz,MOD.ny)]];

end