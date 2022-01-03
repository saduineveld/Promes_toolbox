function [RES] = STND_RBC_proj(par,GRID,POL)
% [RES] = stnd_rbc_proj(par,GRID,POL)
%
% Calculates Euler residuals for standard RBC model, for Promes toolbox
%
% INPUTS:
%		par     = parameters of model
%
%       GRID     = structure with necessary fields (all assigned by prepgrid):
%
%       POL     = structure with necessary fields assigned by solve_proj:
%                   - algo
%                   - theta (polynomial) or pp_y (spline)
%                   - for 'tmi': theta_old (polynomial) or pp_y_old (spline)
%
% OUTPUT:
%       RES     = column vector Euler residuals (mm x 1)
%
%
% USES:
%       get_pol_var (in PromeS toolbox), stnd_gr_aux (model specific file)

% Sijmen Duineveld, updated for Promes_v05.0.0 December 2021, s.a.duineveld@outlook.com

LK = GRID.xx(:,1);%first state variable, log(K_t)
LZ = GRID.xx(:,2);%second state variable, log(Z_t)

%policy variable, log(C_t):
if ~(strcmp(POL.algo,'spl_tmi') || ...
    strcmp(POL.algo,'smol_tmi') ||...
    strcmp(POL.algo,'cheb_tmi') )
    
    %use initial grid for polynomials 
    % (or ignore for spl_dir)
  	spec_opt  = 'ini_grid';    
    LC  = get_pol_var(POL,[LK,LZ],GRID,[],spec_opt);    

else%for 'spl_tmi','cheb_tmi','smol_tmi':
    
    %LC is set directly for Time Iteration
    LC  = POL.YY;
end

%Capital in next period:
LK_n    = stnd_rbc_aux(par,LK,LZ,LC);

if strcmp(POL.algo,'spl_tmi') || ...
	strcmp(POL.algo,'smol_tmi') ||...
	strcmp(POL.algo,'cheb_tmi') 

    % use old policy function in t+1 
    % for Time Iteration 
    % (pp_y_old or theta_old)
    spec_opt_next = 'old_pol';
else
    spec_opt_next = [];
end
       
%Allocate empty matrix for RHS of Euler equation:
rhs_l   = NaN(size(LK,1),par.her.gh_nod);
for ll = 1:par.her.gh_nod
    % Shock to TFP (using Gauss-Hermite nodes):
    EPS_n       = sqrt(2)*par.her.xi(ll);
    
    %log(Z_t+1):
    LZ_n        = par.rho_z*LZ + par.sigma_z*EPS_n;
    
    %log(C_t+1)
    LC_n        = get_pol_var(POL,[LK_n,LZ_n],GRID,[],spec_opt_next);

    %log(MPK_t+1) (marginal prod. of capital)
    [~,LMPK_n]    = stnd_rbc_aux(par,LK_n,LZ_n,LC_n);
    
    % RHS of Euler equation, 
    % weighted by Gauss-Hermite weights
    rhs_l(:,ll) = par.her.wi(ll)/sqrt(pi) * par.beta * exp(-par.nu*LC_n) .* (exp(LMPK_n) + 1 - par.delta);
end

%Right hand side of Euler equation:
RHS = sum(rhs_l,2);

% Euler residuals (scaled by C^-nu):
RES     =  RHS./exp(-par.nu *LC) - 1;
  
end
