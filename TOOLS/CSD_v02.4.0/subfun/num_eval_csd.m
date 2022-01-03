function [NUM,SOL] = num_eval_csd(MOD)
%[NUM,SOL] = num_eval_csd(MOD)
%
% Numerically evaluate derivatives

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner.
% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

for is = 1:size(MOD.par_nms,2)
    eval([MOD.par_nms{1,is},' = MOD.par_val(is);']);
end
clear is;


for iv = 1:size(MOD.var_bs_nms,2)
    eval([MOD.var_bs_nms{1,iv},'_t = MOD.SS_vec(iv);']);
    eval([MOD.var_bs_nms{1,iv},'_n = MOD.SS_vec(iv);']);
end

% Get X_ss, Z_ss, Y_ss
SOL.XX_ss = NaN(1,MOD.nx);
SOL.ZZ_ss = NaN(1,MOD.nz);
SOL.YY_ss = NaN(1,MOD.ny);
cnt_x = 0; cnt_z = 0; cnt_y = 0;
for in = 1:MOD.nx + MOD.nz + MOD.ny
    idx = strcmp(MOD.qq_str{1,in},MOD.var_bs_nms);
    if in < MOD.nx + 0.5
       cnt_x = cnt_x + 1;
       SOL.XX_ss(1,cnt_x) = MOD.SS_vec(1,idx);
    elseif in > MOD.nx +0.5 && in < MOD.nx + MOD.nz + 0.5
        cnt_z = cnt_z + 1;     
        SOL.ZZ_ss(1,cnt_z) = MOD.SS_vec(1,idx);
    else
    	cnt_y = cnt_y + 1;     
        SOL.YY_ss(1,cnt_y) = MOD.SS_vec(1,idx);
    end        
end

NUM.DD = eval(MOD.DS);

if MOD.ord > 1.5
    NUM.HH = eval(MOD.HS);
    if MOD.ord > 2.5
        NUM.TT = eval(MOD.TS);
    end
end 

end
