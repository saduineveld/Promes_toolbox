function [YY] = get_YY(SOL,WW_dev,ord,WW2_dev,WW3_dev)
%[YY] = get_YY(SOL,WW_dev,ord,WW2_dev,WW3_dev)
%
% Get control variables, given state variables
%
% INPUT: see eval_sol_csd
%
% OUTPUT:
%   - YY: ny x cols matrix

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner.
% Sijmen Duineveld, updated April 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

nz = size(SOL.ZZ_ss,2);

YY = SOL.YY_ss' + SOL.Hy_w*WW_dev;
if ord > 1.5
    if nz > 0.5
       Hy_ss_pol = 0.5*SOL.Hy_ss;
    else
       Hy_ss_pol = zeros(size(YY));
    end
    YY = YY + Hy_ss_pol + 0.5*SOL.vec.Hy2*WW2_dev;
    %OLD: YY  = YY    + Hy_ss_pol + 0.5*kron(eye(ny),WW') * SOL.Hy_ww*WW;
    if ord > 2.5 
        if nz > 0.5
            Hy_s3rd = (1/6)*SOL.Hy_sss+0.5*SOL.vec.Hyssw*WW_dev;    
        else
            Hy_s3rd = zeros(size(YY));
        end
        YY = YY   + Hy_s3rd + (1/6)*SOL.vec.Hy3*WW3_dev;   
    end
end

end
