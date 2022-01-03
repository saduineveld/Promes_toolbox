function [xx_dw] = scal_mat_dw(xx,lb,ub)
% function [xx_dw] = scal_mat_dw(xx,lb,ub)
%
% Purpose: Construct scaled down matrix from xx to xx_dw 
%
% Input:
%       xx      = mm x nn matrix: mm datapoints, nn variables
%
%       lb      = lower bounds of grid (1 x nn vector)
%
%       ub      = upper bounds of grid (1 x nn vector)
%
% Ouptut:
%       xx_dw   = xx mapped linearly from [lb,ub] into [-1,1] 
%
% Uses: 
%       sc_cheb_dw (maps linearly from [lb,ub] into [-1,1])

% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

nn = size(xx,2);
mm = size(xx,1);

xx_dw = NaN(mm,nn);
for in = 1:nn     
    xx_dw(:,in) = sc_cheb_dw(lb(in),ub(in),xx(:,in));
end
    
end