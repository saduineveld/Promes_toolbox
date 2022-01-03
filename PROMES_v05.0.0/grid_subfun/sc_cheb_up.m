function [x] = sc_cheb_up(lb,ub,x_dw)
% function [x] = sc_cheb_up(lb,ub,x_dw)
%
% Purpose: maps column vector (x_dw) in interval [-1,1] into [lb,ub] (x);
%
% Input:
%       x_dw    = vector of input data with bounds [-1,1]. 
%       lb      = scalar of lower bound (for x)
%       ub      = scaler of upper bound (for x)
%
% Output:
%       x    = scaled up input vector, which is mapped from interval
%                   [-1,1] into [lb,ub] using linear transformation.


% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

x = (x_dw+1)*(ub-lb)/2 + lb;

end

