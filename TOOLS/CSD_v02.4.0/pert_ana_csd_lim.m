function [SOL,NUM,MOD,GAL] = pert_ana_csd_lim(MOD,Rho,order,Omega,excl_der)
% [SOL,NUM,MOD,GAL] = pert_ana_csd_lim(MOD,Rho,ord,Omega)
%
% Solves pertubation solution, with possible limit cycles
%
% Output SOL:
%           - if one candidate solution: SOL is structure with solution
%           - if multiple candidate solutions exist: SOL is cell array with
%               each cell a candidate solution
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
%       order:      order of pertubation solution. Maximum is 3
%
%       Rho:        Auto correlation matrix of exog. variables:
%                       ZZ_t+1 = Rho*ZZ_t + Omega*eta_t+1;
%                   Set to "[]" for non-stochastic model
%
%       Omega:      variance-covariance matrix of exogenous shocks: 
%                       ZZ_t+1 = Rho*ZZ_t + Omega*eta_t+1;
%                   Set to "[]" for non-stochastic model
%
%      excl_der:    (optional) set to 1 if MOD already includes the derivatives
%
% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner,
% and integrated the function InvSubGen by Dana Galizia
%
% Sijmen Duineveld, Updated June 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

if nargin < 4.5
    excl_der = 0;
end

if order > 3.5 
    error('Maximum order of approximation is 3')
end

if excl_der ~= 1
    MOD             = get_deriv_csd(MOD,order);
end

[NUM,SOL] 	= num_eval_csd(MOD);

NUM.Skew    = 0;%setting Skewness to zero
NUM.Syl     = 0;

SOL.Rho     = Rho;
SOL.Omega   = Omega;

% Construct matrices for InvSubGen:
na      = MOD.nx + MOD.nz + MOD.ny;
GAL.Gx1 = [NUM.DD(:,1:MOD.nx),NUM.DD(:,MOD.nx+MOD.nz+1:na)];
GAL.Gx0 = [NUM.DD(:,na+1:na+MOD.nx),NUM.DD(:,na+MOD.nx+MOD.nz+1:2*na)];
GAL.Gt1 = NUM.DD(:,MOD.nx+1:MOD.nx+MOD.nz);    
GAL.Gt0 = NUM.DD(:,na+MOD.nx+1:na+MOD.nx+MOD.nz);

%Du_v is pol for 'u', given 'v'
[GAL.nW,GAL.othstb,GAL.Dx_x,GAL.Dy_x,GAL.Dx_z,GAL.Dy_z] = InvSubGen(GAL.Gx0,GAL.Gx1,GAL.Gt0,GAL.Gt1,Rho,MOD.nx);

if GAL.nW == 1
    SOL = compile_solution(SOL,MOD,NUM,order,GAL.Dx_x,GAL.Dx_z,GAL.Dy_x,GAL.Dy_z);
else
    for in = 1:GAL.nw
        SOL{1,in} = compile_solution(SOL,MOD,NUM,order,GAL.Dx_x(:,:,in),GAL.Dx_z(:,:,in),GAL.Dy_x(:,:,in),GAL.Dy_z(:,:,in));
    end
end

end

%% Compile solution (reorder solution, and solve higher order:
function [SOL_i] = compile_solution(SOL_i,MOD,NUM,order,Dx_xi,Dx_zi,Dy_xi,Dy_zi)

[SOL_i.Hx_w,SOL_i.Hy_w] = reorder_maus_mat(Dx_xi,Dx_zi,Dy_xi,Dy_zi);

if order > 1.5
	
	[SOL_i.Hx_ww,SOL_i.Hy_ww,SOL_i.Hx_ss,SOL_i.Hy_ss,NUM.aux] = quad_csd(SOL_i.Hx_w,SOL_i.Hy_w,NUM.DD,NUM.HH,SOL_i.Omega,SOL_i.Rho,MOD.nx,MOD.nz,MOD.ny,NUM.Syl);

 	%Transformations for vectorization of simulation:
 	SOL_i.vec.Hx2 = vec_Hv_ww(SOL_i.Hx_ww,MOD.nx,MOD.nx,MOD.nz);
	SOL_i.vec.Hy2 = vec_Hv_ww(SOL_i.Hy_ww,MOD.ny,MOD.nx,MOD.nz);
        
    if order > 2.5
        
        [SOL_i.Hx_www,SOL_i.Hy_www,SOL_i.Hx_ssw,SOL_i.Hy_ssw,SOL_i.Hx_sss,SOL_i.Hy_sss] = ...
                cubic_csd(SOL_i.Hx_w,SOL_i.Hy_w,NUM.DD,NUM.HH,SOL_i.Omega,SOL_i.Rho,MOD.nx,MOD.nz,MOD.ny,NUM.Syl,...
            NUM.TT,NUM.Skew,SOL_i.Hx_ww,SOL_i.Hy_ww,SOL_i.Hx_ss,SOL_i.Hy_ss,NUM.aux);
        
        %Transformations for vectorization of simulation:
        SOL_i.vec.Hx3 = vec_Hv_www(SOL_i.Hx_www,MOD.nx,MOD.nx,MOD.nz);
        SOL_i.vec.Hy3 = vec_Hv_www(SOL_i.Hy_www,MOD.ny,MOD.nx,MOD.nz);
        
        if MOD.nz > 0.5
            SOL_i.vec.Hxssw = vec_Hv_ssw(SOL_i.Hx_ssw,MOD.nx,MOD.nx,MOD.nz);
            SOL_i.vec.Hyssw = vec_Hv_ssw(SOL_i.Hy_ssw,MOD.ny,MOD.nx,MOD.nz);
        end
    end
end

end


%% Reorder system in Maussner's style
function [Hx_w,Hy_w] = reorder_maus_mat(Dx_x,Dx_z,Dy_x,Dy_z)

Hx_w = [Dx_x,Dx_z];
Hy_w = [Dy_x,Dy_z];

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