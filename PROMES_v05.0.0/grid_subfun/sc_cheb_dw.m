function [x_dw] = sc_cheb_dw(lb,ub,x)
% function [x_dw] = sc_cheb_dw(lb,ub,x)
%
% Purpose: maps column vector (x) in interval [lb,ub] into [-1,1] (x_dw);
%
% Input:
%       x       = vector of input data with bounds [lb,ub]
%       lb      = scalar of lower bound (for x)
%       ub      = scaler of upper bound (for x)
%
% Output:
%       x_dw    = scaled down input vector, which is mapped from interval
%       [lb,ub] into [-1,1] using linear transformation. 

% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

x_dw = 2*x/(ub-lb) - (lb+ub)/(ub-lb);

end

