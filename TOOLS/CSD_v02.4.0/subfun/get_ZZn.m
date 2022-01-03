function [ZZn] = get_ZZn(SOL,ZZ_dev,epsilon)
% [ZZn] = get_ZZn(SOL,ZZ_dev,epsilon)
%
% Simulate exogenous state variable
%
% INPUT: 
%   - ZZ_dev: nz by cols matrix with rows:
%           [ZZ(1,:)- ZZ_ss(1);...,ZZ(nz,:)-ZZ_ss(nz)]
%               with each state variable a row vector and cols the
%               datapoints
%	- epsilon: nz x cols matrix of shocks
%
% OUTPUT:
%   - ZZn: nz x cols matrix

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner.
% Sijmen Duineveld, updated April 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

ZZn = SOL.Rho*ZZ_dev + SOL.Omega*epsilon + SOL.ZZ_ss';

end