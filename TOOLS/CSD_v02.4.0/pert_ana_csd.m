function [SOL,NUM,MOD,stab] = pert_ana_csd(MOD,Rho,order,Omega,excl_der)
%[SOL,NUM,MOD,stab] = pert_ana_csd(MOD,Rho,order,Omega,excl_der)
%
% Solves model with pertubation methods
% 
% INPUT:
%       MOD: structure, should contain the following fields (see
%       documentation for details and examples):
%           - FS: symbolic model with endogenous equations (not exogenous
%           variables variables): E_t g(xt+1, zt+1, yt+1, xt , zt , yt) = 0
%           IMPORTANT: variables names have to end in '_n' & '_t' for next
%           period and current period, respectively. So if model contains
%           a variable 'c', then model equations should contain 'c_n' and 'c_t';
%           The base name in 'var_bs_nms' then has to be 'c';
%
%           - XX, XXn, ZZ, ZZn, YY, YYn vectors of symbols of endogenous state (X), ...
%           stochastic state (Z) & control variables (Y)
%           , where 'n' indicates next period
%
%           - var_bs_nms: cell array with base names of variables,
%               (dropping '_n' or '_t', see under 'FS'); Ordering as in SS_vec;
% 
%           - par_nms: names of parameters (in symbolic model)
%
%       Rho:        Auto correlation matrix of exog. variables:
%                       ZZ_t+1 = Rho*ZZ_t + Omega*eta_t+1;
%                   Set to "[]" for non-stochastic model
%
%       order:      order of pertubation solution. Maximum is 3
%
%       Omega:      variance-covariance matrix of exogenous shocks: 
%                       ZZ_t+1 = Rho*ZZ_t + Omega*eta_t+1;
%                   Set to "[]" for non-stochastic model
%
%      excl_der:    (optional) set to 1 if MOD already includes the derivatives

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner. 
% Sijmen Duineveld, updated June 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.


if nargin < 4.5
    excl_der = 0;
end

Skew = 0;%setting Skewness to zero
Syl  = 0;

if order > 3.5 
    error('Maximum order of approximation is 3')
end

if excl_der ~= 1
    MOD             = get_deriv_csd(MOD,order);
end

[NUM,SOL]       = num_eval_csd(MOD);
SOL.Rho         = Rho;
SOL.Omega       = Omega;

[NUM.AA,NUM.BB]     = constr_mat_csd(MOD,NUM.DD,SOL.Rho);

[SOL.Hx_w,SOL.Hy_w,SOL.solab.PP,SOL.solab.FF,stab] = solab_adj_stab(NUM.AA,NUM.BB,MOD.nx,MOD.nz);

if stab == 1
if order > 1.5
    
    SOL.Omega   = Omega;
    [SOL.Hx_ww,SOL.Hy_ww,SOL.Hx_ss,SOL.Hy_ss,NUM.aux] = quad_csd(SOL.Hx_w,SOL.Hy_w,NUM.DD,NUM.HH,SOL.Omega,SOL.Rho,MOD.nx,MOD.nz,MOD.ny,Syl);
    
    %Transformations for vectorization of simulation:
    SOL.vec.Hx2 = vec_Hv_ww(SOL.Hx_ww,MOD.nx,MOD.nx,MOD.nz);
    SOL.vec.Hy2 = vec_Hv_ww(SOL.Hy_ww,MOD.ny,MOD.nx,MOD.nz);
    
    if order > 2.5
        NUM.Skew    = Skew;
        [SOL.Hx_www,SOL.Hy_www,SOL.Hx_ssw,SOL.Hy_ssw,SOL.Hx_sss,SOL.Hy_sss] = ...
            cubic_csd(SOL.Hx_w,SOL.Hy_w,NUM.DD,NUM.HH,SOL.Omega,SOL.Rho,MOD.nx,MOD.nz,MOD.ny,Syl,...
            NUM.TT,NUM.Skew,SOL.Hx_ww,SOL.Hy_ww,SOL.Hx_ss,SOL.Hy_ss,NUM.aux);
        
        %Transformations for vectorization of simulation:
        SOL.vec.Hx3 = vec_Hv_www(SOL.Hx_www,MOD.nx,MOD.nx,MOD.nz);
        SOL.vec.Hy3 = vec_Hv_www(SOL.Hy_www,MOD.ny,MOD.nx,MOD.nz);
        
        if MOD.nz > 0.5
            SOL.vec.Hxssw = vec_Hv_ssw(SOL.Hx_ssw,MOD.nx,MOD.nx,MOD.nz);
            SOL.vec.Hyssw = vec_Hv_ssw(SOL.Hy_ssw,MOD.ny,MOD.nx,MOD.nz);
        end
    end
end
end


end

%% TRANSFORMATIONS OF POLICIES FOR VECTORIZATION OF SIMULATION
% Second order policy transformation for Hx_ww and Hyww:
function [Hv2_rows] = vec_Hv_ww(Hv_ww,nv,nx,nz)

Hv2_rows = NaN(nv,(nx+nz)^2);

cnt_r = 0;
for iv = 1:nv    
    st_r = cnt_r + 1;
    cnt_r = cnt_r + (nx+nz);
    
    Hv2_rows(iv,:) = reshape(Hv_ww(st_r:cnt_r,:),1,[]);   
end

end

% Third order policy transformation for Hx_www and HY_www
function [Hv3_rows] = vec_Hv_www(Hv_www,nv,nx,nz)

Hv3_rows = NaN(nv,(nx+nz)^3);

cnt_r = 0;
for iv = 1:nv    
    st_r = cnt_r + 1;
    cnt_r = cnt_r + (nx+nz)^2;
    
    Hv3_rows(iv,:) = reshape(Hv_www(st_r:cnt_r,:),1,[]);   
end

end

% Third order policy transformation for Hx_ssw and HY_ssw
function [Hvssw_rows] = vec_Hv_ssw(Hv_ssw,nv,nx,nz)

Hvssw_rows = NaN(nv,nx+nz);

cnt_r = 0;
for iv = 1:nv    
    st_r = cnt_r + 1;
    cnt_r = cnt_r + (nx+nz);
    
    Hvssw_rows(iv,:) = reshape(Hv_ssw(st_r:cnt_r,:),1,[]);   
end

end