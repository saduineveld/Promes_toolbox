function [LC_new] = STND_RBC_fpa(par,GRID,POL)
% [LC_new] = stnd_rbc_fpa(par,GRID,POL)
%
% Calculates new consumption policy, given old policy for standard RBC model, for Promes toolbox

LK = GRID.xx(:,1);%first state variable, log(K_t)
LZ = GRID.xx(:,2);%second state variable, log(Z_t)

%policy variable, log(C_t):
% Evaluate the policy function log(C):
if strcmp(POL.algo,'spl_fpa') || ...
    strcmp(POL.algo,'smol_fpa') || ...
    strcmp(POL.algo,'cheb_fpa') 

    % use old policy for RHS of Euler:
    LC1     = POL.YY_old; 
end

%Capital in next period:
LK_n    = stnd_rbc_aux(par,LK,LZ,LC1);

if strcmp(POL.algo,'spl_fpa') || ...
    strcmp(POL.algo,'smol_fpa') || ...
    strcmp(POL.algo,'cheb_fpa') 

    % use old policy function in t+1 
    % (pp_y_old or theta_old)
    spec_opt_next = 'old_pol';
end
       
%Allocate empty matrix for RHS of Euler equation:
rhs_l   = NaN(size(LK,1),par.her.gh_nod);
for ll = 1:par.her.gh_nod
    % Shock to TFP (using Gauss-Hermite nodes):
    EPS_n       = sqrt(2)*par.her.xi(ll);
    
    %log(Z_t+1):
    LZ_n        = par.rho_z*LZ + par.sigma_z*EPS_n;
    
    %log(C_t+1)
    LC_n        = get_pol_var_fpa(POL,[LK_n,LZ_n],GRID,[],spec_opt_next);

    %log(MPK_t+1) (marginal prod. of capital)
    [~,LMPK_n]    = stnd_rbc_aux(par,LK_n,LZ_n,LC_n);
    
    % RHS of Euler equation, 
    % weighted by Gauss-Hermite weights
    rhs_l(:,ll) = par.her.wi(ll)/sqrt(pi) * par.beta * exp(-par.nu*LC_n) .* (exp(LMPK_n) + 1 - par.delta);
end

%Right hand side of Euler equation:
RHS = sum(rhs_l,2);

% Compute LC_new:
LC_new =  -1/par.nu*log(RHS);
  
end
