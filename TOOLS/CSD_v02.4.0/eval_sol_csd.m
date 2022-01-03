function [XXn,YY,ZZn] = eval_sol_csd(SOL,XX,ZZ,epsilon,order)
% [XXn,YY,ZZn] = eval_sol_csd(SOL,XX,ZZ,epsilon,order)
%    
% INPUT:
%       - SOL: structure containing policy functions (Hi_qqq's), XX_ss,
%       ZZ_ss, YY_ss (1 by nx, nz, ny vectors, respectively)
%
%       - XX: nx by cols matrix with each endo. state variable a row:
%               [XX(1,:);...;XX(nx,:)]
%
%       - ZZ: nz by cols matrix with each exo. state variable a row:
%               [ZZ(1,:),...,ZZ(nz,:)]
%
%       - epsilon: nz by cols matrix of shocks
%
%       - order: order of simulation
%
% OUTPUT:
%       - XXn = XX(t+1), nx by cols  matrix
%       - YY  = YY(t), ny by cols matrix
%       - ZZn = ZZ(t+1), nz by cols matrix

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner, 
% Sijmen Duineveld, Updated May 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.


nx = size(SOL.XX_ss,2);
ny = size(SOL.YY_ss,2);
nz = size(SOL.ZZ_ss,2);

% Check consistency:
if size(XX,1) ~= nx
    error('XX does not have nx rows')
elseif size(ZZ,1) ~= nz
    error('ZZ does not have nz rows')

elseif order ~= 1 && order ~= 2 && order ~= 3
    error('Order should be 1, 2 or 3')
end

XX_dev = XX-SOL.XX_ss';
if nz > 0
    ZZ_dev = ZZ-SOL.ZZ_ss';
else
    ZZ_dev = [];
end
WW_dev = [XX_dev;ZZ_dev];

% Construct higher order terms of WW:
if order == 2 || order == 3
    if order == 2
        WW2_dev = high_order_WW(WW_dev);
        WW3_dev = [];
        
    elseif order == 3        
        [WW2_dev,WW3_dev] = high_order_WW(WW_dev);
        
    end
else
    WW2_dev = [];
    WW3_dev = [];
end

% Calculate variables:
[XXn]   = get_XXn(SOL,WW_dev,order,WW2_dev,WW3_dev);
if nargout > 1
    [YY]    = get_YY(SOL,WW_dev,order,WW2_dev,WW3_dev);
end
if nargout > 2 
    if nz > 0
        if size(epsilon,1) ~= nz
            error('epsilon does not have nz rows')
        end
        [ZZn]   = get_ZZn(SOL,ZZ_dev,epsilon);
    else
       ZZn = [];
    end
end

end

%% Construct higher order terms of WW:
function [WW2_dev,WW3_dev] = high_order_WW(WW_dev)

nw    = size(WW_dev,1);
cols    = size(WW_dev,2);
WW2_dev     = NaN(nw^2,cols);
WW3_dev     = NaN(nw^3,cols);

cnt2 = 0;
cnt3 = 0;
for i1 = 1:nw
    for i2 = 1:nw
        cnt2 = cnt2+1;
        WW2_dev(cnt2,:) = WW_dev(i1,:).*WW_dev(i2,:);
        if nargout > 1
            for i3 = 1:nw
                cnt3 = cnt3+1;    
                WW3_dev(cnt3,:) = WW2_dev(cnt2,:).*WW_dev(i3,:);    
            end       
        end
    end
end

end
