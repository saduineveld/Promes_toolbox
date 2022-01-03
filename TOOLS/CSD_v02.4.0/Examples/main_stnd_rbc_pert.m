function [SOL,par,SS,SIM] = main_stnd_rbc_pert()
% [SOL,par,SS,SIM] = main_stnd_rbc_pert()

% Sijmen Duineveld, updated May 2021, s.a.duineveld@outlook.com

%% STEP 0:  Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

%Add folder with CSD and subfun:
addpath ('..');
addpath ('..\subfun');

addpath stnd_rbc;

%% STEP 1: GET SYMBOLIC MODEL & DIFFERENTIATE
MOD         = STND_RBC_mod;

% Differentiate symbolic model:
par.opt.order = 3;
[MOD] = get_deriv_csd(MOD,par.opt.order);


%% STEP 2: SET PARAMETERS
par.alpha       = 0.36;		%capital share income 
par.beta        = 0.985;	% discount factor
par.delta       = 0.025;	% depreciation of capital 
par.nu          = 2; 	% risk aversion 
par.eta         = 4; 	% Frisch elasticity of lab. supply
par.chi         = 1; 	% scalar for s.s hours worked

par.rho_z       = 0.95;	% autocorrelation coefficient TFP
par.sigma_z     = 0.01;	% standard dev. shocks in TFP 

% Assign the numerical values to MOD structure:
MOD.par_val = [par.alpha,par.beta,...
        par.delta,par.eta,par.nu,par.chi];

%% STEP 3: SET STEADY STATE VALUES
SS = stnd_rbc_ss(par,1);
% the 1 is steady state TFP

% Assign S.S. values to the MOD structure :
MOD.SS_vec      = [log(SS.Kss),log(SS.Zss),log(SS.Css),log(SS.Hss)];


%% STEP 4: SOLVE MODEL NUMERICALLY
%exclude calculation of derivates (see STEP 1)
excl_der = 1;

% Solve model numerically:
SOL  = pert_ana_csd(MOD,par.rho_z,par.opt.order,par.sigma_z,excl_der);


%% STEP 5: SIMULATION
par.sim.TT      = 1000;
par.sim.cols    = 2;
par.sim.T_ini   = 10;

SIM = stnd_rbc_sim_pert(par,SS,SOL,par.sim,par.opt.order);


%% Plot two simulated series of output
figure;
par.sim.T_plot = 200;%Plot last 50 years in simulation
plot([1:par.sim.T_plot]/4,SIM.LY(end-par.sim.T_plot+1:end,1)-log(SS.Yss),'LineWidth',1.5);
hold all;
plot([1:par.sim.T_plot]/4,SIM.LY(end-par.sim.T_plot+1:end,2)-log(SS.Yss),'LineWidth',1.5);
xlabel('Time (years)');
ylabel('Output (% dev.)');
legend({'Simulation 1','Simulation 2'},'Location','best');


end

%% Perturbation model of Standard RBC model
function [MOD] = STND_RBC_mod()

% Note: Law of Motion for exogenous states 
% is defined outside model file:
% ZZ_t = Rho*ZZ_t-1 + Omega*epsilon

%% BLOCK 1: Define parameters and variables
%Parameters:
syms aa bb dd ee nu ch;

%Variables:
syms LC_t LC_n LK_t LK_n LZ_t LZ_n LH_t LH_n;

%% BLOCK 2: MODEL EQUATIONS (excl. stochastic process)
%Auxiliary variables:
LY_t      = LZ_t + aa*LK_t + (1-aa)*LH_t;
LY_n      = LZ_n + aa*LK_n + (1-aa)*LH_n;

%Euler equation:
f1 = bb *exp(-nu*LC_n)*( aa*exp(LY_n-LK_n) + (1-dd) ) - exp(-nu*LC_t);

% Capital accumulation:
f2 = exp(LY_t) + (1-dd)*exp(LK_t) - exp(LC_t) - exp(LK_n);

% Labour market:
f3 = log(1-aa) + LZ_t + aa*(LK_t-LH_t) + -nu*LC_t - log(ch) - 1/ee*LH_t;

%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2;f3];

%Endogenous state variables:
MOD.XX = LK_t;
%Exogenous state variables:
MOD.ZZ = LZ_t;
%Control variables:
MOD.YY = [LC_t,LH_t];

MOD.XXn = LK_n;
MOD.ZZn = LZ_n;
MOD.YYn = [LC_n,LH_n];

% Variable names (strings):
MOD.var_bs_nms  = {'LK','LZ','LC','LH'};
% Parameter names (strings):
MOD.par_nms     = {'aa','bb','dd','ee','nu','ch'};

end

%% Analytical steady state
function [SS] = stnd_rbc_ss(par,Zss)
%[SS] = stnd_rbc_ss(par,Zss)
%
% Calculation of steady state of Standard RBC model 
%
% INPUTS:
%		par     = parameters of model
%
%       Zss     = steady state level of TFP

% Sijmen Duineveld, updated March 2021, s.a.duineveld@outlook.com

Omega = (1-par.beta*(1-par.delta))/ (par.alpha*par.beta*Zss);

Kss = ( ((1-par.alpha)/par.chi*Zss.*(Zss.*Omega-par.delta).^-par.nu ).^par.eta .* Omega.^((1+par.eta*par.alpha)/(par.alpha-1)) ).^(1/(1+par.eta*par.nu));

Hss = Omega.^(1/(1-par.alpha)) .* Kss;

Css = (Zss.*Omega - par.delta) .* Kss;

Yss = Zss *Kss^par.alpha*Hss^(1-par.alpha);

nms                 = {'K','C','H','Z','Y'}; 
for in = 1:size(nms,2)
   eval(['SS.',nms{in},'ss = ',nms{in},'ss;']);  
end

end

%% SIMULATION given perturbation solution
function [SIM] = stnd_rbc_sim_pert(par,SS,SOL,opt_sim,order)
%[SIM] = stnd_rbc_sim_pert(par,SS,SOL,opt_sim,order)
%
%       opt_sim     = structure with settings for simulation. Contains:
%                       TT      - number of periods in simulation (scalar)
%                       cls     - number of columns in simulation (scalar)
%                       T_ini   - initial periods at s.s.
%
%       order       = order of solution to be used
%
% OUTPUT:
%       SIM         = structure with series in logs: LK (capital), 
%                       LC (consumption), LH (hours worked), LZ (TFP)

cols = opt_sim.cols;
nx = size(SOL.XX_ss,2);
nz = size(SOL.ZZ_ss,2);
ny = size(SOL.YY_ss,2);

% Initialize variables (LK, LC, LH, LZ)
nms                 = {'K','C','H','Z'}; 
for in = 1:size(nms,2)
   eval(['L',nms{in},' = NaN(opt_sim.T_ini+opt_sim.TT,cols);']); 
   eval(['L',nms{in},'(1:opt_sim.T_ini,:) = log(SS.',nms{in},'ss);']);  
end   
% pre-allocate one extra row for LK_t+1 and LZ_t+1
LK = [LK;NaN(1,cols)];
LZ = [LZ;NaN(1,cols)];

% initialize shocks:
rng('default');
rng(1);
epsilon = [zeros(cols,opt_sim.T_ini),randn(cols,opt_sim.TT+1)]';

% Set k_t+1 = k_ss 
LK(opt_sim.T_ini+1,:)  = log(SS.Kss);

% Set z_t+1 = lz(t)+eps(t+1)
LZ(opt_sim.T_ini+1,:) = get_ZZn(SOL,zeros(1,cols),epsilon(opt_sim.T_ini+1,:));

% loop over time (from T_ini +1 till T_ini + TT):
for it = opt_sim.T_ini + 1: opt_sim.T_ini + opt_sim.TT   

   % Get XXn, YY, ZZn using eval_sol_csd:
    [LK(it+1,:),YY_t,LZ(it+1,:)] = ...
        eval_sol_csd(SOL,LK(it,:),LZ(it,:),epsilon(it+1,:),order);
    LC(it,:) = YY_t(1,:);%Cons. is first control variable
    LH(it,:) = YY_t(2,:);%Hours is second control variable    
end

% remove last k_t+1 and z_t+1 (ie. adjust size to (T_ini+TT) x cols )
LK = LK(1:end-1,:);
LZ = LZ(1:end-1,:);

% Add output:
SIM.LY = LZ + par.alpha*LK + (1-par.alpha)*LH;

% Put variables in structure SIM:
for in = 1:size(nms,2)
   eval(['SIM.L',nms{in},' = L',nms{in},';']); 
end  

SIM.epsilon = epsilon;

end