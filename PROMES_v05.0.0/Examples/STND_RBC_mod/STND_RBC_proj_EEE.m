%% PROJECTION MODEL
function [ERR] = STND_RBC_proj_EEE(par,xx,GRID,POL)
% [ERR] = STND_RBC_EEE(par,xx,GRID,POL)
%
% Calculate normalized Euler Equation Errors

LK = xx(:,1);%first state variable, log(K_t)
LZ = xx(:,2);%second state variable, log(Z_t)

%policy variable, log(C_t):
LC  = get_pol_var(POL,[LK,LZ],GRID);    

%Capital in next period:
LK_n    = stnd_rbc_aux(par,LK,LZ,LC);
      
%Allocate empty matrix for RHS of Euler equation:
lng_h = size(par.her.xi,1);
rhs_j   = NaN(size(LK,1),lng_h);
for jj = 1:lng_h
    % Shock to TFP (using Gauss-Hermite nodes):
    EPS_n       = sqrt(2)*par.her.xi(jj);
    
    %log(Z_t+1):
    LZ_n        = par.rho_z*LZ + par.sigma_z*EPS_n;
    
    %log(C_t+1)
    LC_n        = get_pol_var(POL,[LK_n,LZ_n],GRID);

    %log(MPK_t+1) (marginal prod. of capital)
    [~,LMPK_n]    = stnd_rbc_aux(par,LK_n,LZ_n,LC_n);
    
    % RHS of Euler equation, 
    % weighted by Gauss-Hermite weights
    rhs_j(:,jj) = par.her.wi(jj)/sqrt(pi) * par.beta * exp(-par.nu*LC_n) .* (exp(LMPK_n) + 1 - par.delta);
end

RHS = sum(rhs_j,2);

ERR    = (RHS.^(-1/par.nu) - exp(LC))./exp(LC);
 
end