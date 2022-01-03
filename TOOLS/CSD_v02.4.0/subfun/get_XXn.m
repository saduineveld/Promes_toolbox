function [XXn] = get_XXn(SOL,WW_dev,ord,WW2_dev,WW3_dev)
%[XXn] = get_XXn(SOL,WW_dev,ord,WW2_dev,WW3_dev)
%
% Get endogenous state variables in t+1, given current state variables
%
% INPUT: see eval_sol_csd
%
% OUTPUT:
%   - XXn: nx x cols matrix

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner.
% Sijmen Duineveld, updated April 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

nz = size(SOL.ZZ_ss,2);

XXn = SOL.XX_ss' + SOL.Hx_w*WW_dev;

if ord > 1.5
    if nz > 0.5
       Hx_ss_pol = 0.5*SOL.Hx_ss;
    else
       Hx_ss_pol = zeros(size(XXn));
    end
    XXn = XXn + Hx_ss_pol + 0.5*SOL.vec.Hx2*WW2_dev;
    if ord > 2.5          
        if nz > 0.5            
            Hx_s3rd = (1/6)*SOL.Hx_sss+0.5*SOL.vec.Hxssw*WW_dev;            
        else
            Hx_s3rd = zeros(size(XXn));
        end
        XXn = XXn + Hx_s3rd + (1/6)*SOL.vec.Hx3*WW3_dev;        
    end
end

end
