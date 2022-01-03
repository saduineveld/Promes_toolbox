function [SOL] = solve_stnd_rbc_csd(par,SS,order)
% [SOL] = solve_stnd_rbc_csd(par,SS,order)
%
% Solves standard RBC model with CSD toolbox

% Sijmen Duineveld, September 2021, s.a.duineveld@outlook.com

MOD         = stnd_rbc_csd;%get symbolic model;

% Numerical values of steady state:
MOD.SS_vec  = [log(SS.Kss),log(SS.Zss),log(SS.Css),log(SS.Hss),log(SS.Yss)];

% Numerical values of parameters:
MOD.par_val = [par.alpha,par.beta,...
        par.delta,par.eta,par.nu,par.chi];

% Get solution:
SOL    = pert_ana_csd(MOD,par.rho_z,order,par.sigma_z);

end