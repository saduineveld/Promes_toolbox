% Compares the perturbation solution of Dynare (4.6.4) and CSD toolbox.
% Model: standard RBC model.

% Sijmen Duineveld, September 2021, s.a.duineveld@outlook.com

clear all;
clc;
close all;

dbstop if error;

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

addpath stnd_rbc;
addpath TOOLS;
addpath ('..')
addpath ('..\subfun')
addpath C:\dynare\4.6.4\matlab


%% Set parameters and get steady state
par.alpha       = 0.36;		%capital share income 
par.beta        = 0.985;	% discount factor
par.delta       = 0.025;	% depreciation of capital 
par.nu          = 2; 	% risk aversion 
par.eta         = 4; 	% Frisch elasticity of lab. supply
par.chi         = 1; 	% scalar for s.s hours worked

par.rho_z       = 0.95;	% autocorrelation coefficient TFP
par.sigma_z     = 0.01;	% standard dev. shocks in TFP 


%% Solve steady state
SS              = stnd_rbc_ss(par,1);
% the 1 is steady state TFP

% Settings for simulation
par.sim.TT      = 1000;%Number of periods per simulation
par.sim.cols    = 25;% Number of simulations 
par.sim.T_ini   = 10;% Initial number of periods at steady state

par.opt.order   = 3;% Order of the approximation (in simulation)


%% Solve with dynare:
run_dyn = 0;%Set to 0 to load saved Dynare solution; Set to 1 to recompute Dynare solution

if run_dyn == 1;
    % Add the folder with your dynare installation:
    addpath C:\dynare\4.6.4\matlab;

    save ('stnd_rbc_dyn\params','par','SS');

    cd stnd_rbc_dyn;
    
    dynare stnd_rbc_dyn noclearall;
    clearvars -except par SS;
    
    cd ..\;
end

load stnd_rbc_dyn\stnd_rbc_dyn_results oo_ M_ options_;


%% Simulate with dynare eval_dyn_pol_beta:
par.sim.sol_typ = 'dyn';
DR2              = oo_.dr;
DR2.state_var    = M_.state_var;
SIM_d = stnd_rbc_sim_pert(par,SS,DR2,par.sim,par.opt.order);


%% Solve with CSD:
[SOL] = solve_stnd_rbc_csd(par,SS,par.opt.order);

par.sim.sol_typ = 'csd';
SIM_c = stnd_rbc_sim_pert(par,SS,SOL,par.sim,par.opt.order);


%% Compare differences with eval_dyn
dZ = abs(SIM_c.LZ - SIM_d.LZ);
dZ_max = max(dZ(:));

dK = abs(SIM_c.LK - SIM_d.LK);
dK_max = max(dK(:));

dC = abs(SIM_c.LC - SIM_d.LC);
dC_max = max(dC(:));

fprintf('\n\nMaximum differences between Dynare and CSD Toolbox in simulation\nlog(Z): %.2e; log(K): %.2e; log(C): %.3e\n',dZ_max,dK_max,dC_max);

