function main_repl_bgp_aer2020()
% main_repl_bgp_aer202()
%
% Replicates the non-linear risk premium model of BGP (2020, AER, p. 35)
%
% Written for CSD version 2.4.0

% Sijmen Duineveld, December 2021, s.a.duineveld@outlook.com

%% STEP 0:  Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

%Add folder with CSD and subfun:
addpath ('..\');
addpath ('..\subfun');

%% LOAD PARAMETERS
load 'Params_2_50_bgp_nlrp';

par.alpha1  = -MPsv.al1ht;
par.alpha2  = MPsv.al2ht;
par.alpha3  = MPsv.al3ht;
par.alpha4  = MPsv.tau/MPsv.kap;

par.delta   = MPsv.del;
par.psi     = MPsv.psi;
par.vrho1   = MPsv.r1;
par.vrho2   = MPsv.r2/2;
par.vrho3   = MPsv.r3/6;

par.rho     = MPsv.rho;
par.sigma   = MPsv.sight;

%% GET MODEL
[MOD] = BGP_AER2020_CSD();

% Parameters:
MOD.par_val = [par.alpha1,par.alpha2,par.alpha3,par.alpha4,par.delta,par.psi,...
                par.vrho1,par.vrho2,par.vrho3,par.rho,par.sigma];
 
% Steady state:
MOD.SS_vec = [0,0,0,0];

par.opt.order = 3;
[SOL,NUM]  = pert_ana_csd_lim(MOD,par.rho,par.opt.order,par.sigma);

%% SIMULATE MODEL
par.sim.TT      = 350;
par.sim.cols    = 1;
par.sim.T_ini   = 5;

par.sim.sim_type = 'sim';
SIM             = bgp_aer2020_sim_csd(par,SOL,par.sim,par.opt.order);

par.sim.sim_type = 'nsh';
NSH            = bgp_aer2020_sim_csd(par,SOL,par.sim,par.opt.order);

%% PLOTS
% Plot stochastic sim:
f1 = figure;
plot(SIM.Le(5:273),'LineWidth',1.5);
xlim([0 269]);
ylabel('Hours (log deviation)');
title('Stochastic simulation');

% Plot non-stochastic sim:
f2 = figure;
plot(NSH.Le(47:315),'LineWidth',1.5);
xlim([0 269]);
ylabel('Hours (log deviation)');
title('Non-stochastic simulation');

end

function [MOD] = BGP_AER2020_CSD()

%% BLOCK 1: Define parameters and variables
%Parameters:
syms alf1 alf2 alf3 alf4 del ps vrho1 vrho2 vrho3 rho sgm;

%Variables:
syms LX_n LX_t Le1_t Le1_n Le_t Le_n mu_n mu_t;

% Auxiliary:
r_t = vrho1*Le_t + vrho2*Le_t^2 + vrho3*Le_t^3;

%% BLOCK 2: MODEL EQUATIONS (excl. stochastic process)
f1 = (1-del)*LX_t + ps*Le_t - LX_n;

f2 = alf1*LX_t + alf2*Le1_t + alf3*Le_n - alf4*r_t + alf4*mu_t - Le_t;

f3 = Le_t - Le1_n;


%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2;f3];

%Endogenous state variables:
MOD.XX = [LX_t, Le1_t];
%Exogenous state variables:
MOD.ZZ = [mu_t];
%Control variables:
MOD.YY = Le_t;

MOD.XXn = [LX_n, Le1_n];
MOD.ZZn = [mu_n];
MOD.YYn = Le_n;

% Variable names (strings):
MOD.var_bs_nms  = {'LX','Le1','mu','Le'};
% Parameter names (strings):
MOD.par_nms     = {'alf1','alf2','alf3','alf4','del','ps',...
                    'vrho1','vrho2','vrho3','rho','sgm'};

end

%% SIMULATION given perturbation solution (Dynare or CSD toolbox)
function [SIM] = bgp_aer2020_sim_csd(par,SOL,opt_sim,order)
%[SIM] = bgp_aer2020_sim_csd(par,SOL,opt_sim,order)
%
%       opt_sim     = structure with settings for simulation. Contains:
%                       TT          - number of periods in simulation (scalar)
%                       cls         - number of columns in simulation (scalar)
%                       T_ini       - initial periods at s.s.
%                       sim_type    - 'nsh' for no shock; 'sim' for
%                                       standard
%
%       order       = order of solution to be used
%
% OUTPUT:
%       SIM         = structure with series in logs

% Sijmen Duineveld, December 2021, s.a.duineveld@outlook.com

cols = opt_sim.cols;


% Initialize variables (LX Le mu)
nms                 = {'LX','Le','mu'}; 
for in = 1:size(nms,2)
   eval([nms{in},' = NaN(opt_sim.T_ini+opt_sim.TT,cols);']); 
   eval([nms{in},'(1:opt_sim.T_ini,:) = 0;']);  
end   
% pre-allocate one extra row for endo. state variable
LX = [LX;NaN(1,cols)];


% initialize shocks:

if strcmp(opt_sim.sim_type,'sim')
    rng(406430064);
    epsilon = [zeros(cols,opt_sim.T_ini),randn(cols,opt_sim.TT+1)]';
    
    Le(opt_sim.T_ini,:) = 0.032;
    
elseif strcmp(opt_sim.sim_type,'nsh')
    
    epsilon = [zeros(cols,opt_sim.T_ini),zeros(cols,opt_sim.TT+1)]';
    
    Le(opt_sim.T_ini,:) = 0.01;    
end

% Set LX_t+1
LX(opt_sim.T_ini+1,:)  = 0;

% Set mu_t+1 = rho*mu(t)+sigma*eps(t+1)
mu(opt_sim.T_ini+1,:) = get_ZZn(SOL,zeros(1,cols),epsilon(opt_sim.T_ini+1,:));

% loop over time (from T_ini +1 till T_ini + TT):
for it = opt_sim.T_ini + 1: opt_sim.T_ini + opt_sim.TT   
    
    % Get XXn, YY, ZZn using eval_sol_csd:
    [XXn,Le(it,:),mu(it+1,:)] = ...
            eval_sol_csd(SOL,[LX(it,:)';Le(it-1,:)'],mu(it,:),epsilon(it+1,:),order);
        LX(it+1,:) = XXn(1,:);                 
end

% remove last row:
LX = LX(1:end-1,:);
mu = mu(1:end-1,:);


% Put variables in structure SIM:
for in = 1:size(nms,2)
   eval(['SIM.',nms{in},' = ',nms{in},';']); 
end  

SIM.epsilon = epsilon;

end