function [xx] = scal_mat_up(xx_dw,lb,ub)
% function [xx] = scal_mat_up(xx_dw,lb,ub)
%
% Purpose: Construct scaled up matrix from xx_dw to xx 
%
% Input:
%       xx_dw   = mm x nn matrix: mm datapoints, nn variables
%
%       lb      = lower bounds of grid (1 x nn vector)
%
%       ub      = upper bounds of grid (1 x nn vector)
%
% Ouptut:
%       xx      = xx_dw mapped linearly from [-1,1] to [lb,ub]
%
% Uses: 
%       sc_cheb_up (maps linearly from [lb,ub] to [-1,1])

% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

nn = size(xx_dw,2);
mm = size(xx_dw,1);

xx = NaN(mm,nn);
for in = 1:nn
    xx(:,in) = sc_cheb_up(lb(in),ub(in),xx_dw(:,in));
end 
    
end