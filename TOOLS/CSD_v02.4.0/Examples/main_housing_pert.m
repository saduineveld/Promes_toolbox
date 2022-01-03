function [SOL,par,MOD,SS,SIM] = main_housing_pert()
% [POL,par,GRID,SS] = main_housing_pert()
% 
% Solves an RBC model with housing

% Sijmen Duineveld, Updated May 2021, s.a.duineveld@outlook.com

%% STEP 0:  Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

%Add folder with CSD and subfun:
addpath ('..');
addpath ('..\subfun');

%% STEP 1: SET PARAMETERS
par.alpha       = 0.36;     % coefficient on capital in production function
par.beta        = 0.985;    % discount rate

par.delta_k     = 0.025;    % depreciation capital
par.delta_d     = 0.01;     % depreciation housing  

par.eta         = 5;        % risk aversion for housing
par.nu          = 2;      % risk aversion for consumption
  
par.varrho      = 0.1;      % weight of housing utility

par.rho_z       = 0.95;	% autocorrelation coefficient TFP
par.sigma_z     = 0.01;	% standard dev. shocks in TFP 


%% STEP 2: SOLVE STEADY STATE (numerically)
[SS] = housing_ss(par,1);


%% STEP 3: GET SYMBOLIC MODEL
MOD         = HOUSING_pert;


%% STEP 4: SET NUMERICAL VALUES
% Values of parameters:
MOD.par_val = [par.alpha,par.beta,...
        par.delta_k,par.delta_d,par.eta,par.nu,par.varrho,par.rho_z,par.sigma_z];

% Vector of steady states:
MOD.SS_vec  = [log(SS.Kss),log(SS.Dss),log(SS.Zss),log(SS.Css)];
    

%% STEP 5: SOLVE MODEL NUMERICALLY
order = 3;
SOL         = pert_ana_csd(MOD,par.rho_z,order,...
        par.sigma_z);
    
    
%% STEP 6: SIMULATION
par.sim.TT      = 1000;
par.sim.cols    = 100;
par.sim.T_ini   = 10;

SIM = housing_sim_pert(par,SS,SOL,par.sim,3);


end

%% Housing model for CSD toolbox
function [MOD] = HOUSING_pert()
% Note: Law of Motion for exogenous states 
% is defined outside model file:
% ZZ_t = Rho*ZZ_t-1 + Omega*epsilon

%% BLOCK 1: Define parameters and variables
%Parameters:
syms aa bb dk dh ee nu varrho rhoz sgmz;
%Variables:
syms LK_t LK_n LD_t LD_n  LZ_t LZ_n LC_t LC_n;

%% BLOCK 2: MODEL EQUATIONS (excl. stochastic process)
%Auxiliary variables:
Lambda_t = exp(-nu*LC_t);
Lambda_n = exp(-nu*LC_n);

%Marginal prod. capital in t+1:
MPKn  = aa*exp(LZ_n + (aa-1)*LK_n );
%Output in t:
YY    = exp(LZ_t + aa*LK_t);


%Euler 1 (capital):
f1 = bb *Lambda_n*(MPKn + 1-dk) - Lambda_t;

% Euler 2 (housing):
f2 = bb *( varrho*exp(-ee*LD_n) + Lambda_n*(1-dh) ) - Lambda_t;

% Accumulation of assets:
f3 = YY + (1-dk)*exp(LK_t) - exp(LK_n)  + (1-dh)*exp(LD_t) - exp(LD_n) - exp(LC_t);

%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2;f3];

%Endogenous state variables:
MOD.XX = [LK_t,LD_t];
%Exogenous state variables:
MOD.ZZ = [LZ_t];
%Control variables:
MOD.YY = LC_t;

MOD.XXn = [LK_n,LD_n];
MOD.ZZn = [LZ_n];
MOD.YYn = LC_n;

% Variable names (strings):
MOD.var_bs_nms  = {'LK','LD','LZ','LC'};
% Parameter names (strings):
MOD.par_nms     = {'aa','bb','dk','dh','ee','nu','varrho','rhoz','sgmz'};

end

%% Steady state
function [SS] = housing_ss(par,Zss)

SS.Zss  = Zss;

SS.Kss  = (SS.Zss*par.alpha*par.beta/(1-par.beta*(1-par.delta_k)))^(1/(1-par.alpha));

% Solve Css and Dss numerically:
X0      = 0.05*SS.Kss*ones(2,1);
[Xss,~,ex_fl]     = fsolve(@(Xss)res_housing_ss(par,SS.Zss,SS.Kss,Xss),X0,...
    optimoptions(@fsolve,'Display','off','FiniteDifferenceType','central',...
                    'FunctionTolerance',1e-12,'StepTolerance',1e-12,'OptimalityTolerance',1e-12));
if ex_fl < 1
    error('No proper solution found');
end

SS.Css     = Xss(1);
SS.Dss     = Xss(2);

end

%% Residual function for steady state
function [RES] = res_housing_ss(par,Zss,Kss,Xss)
        
Css     = Xss(1);
Dss     = Xss(2);

Yss  = output_MPK(par,Kss,Zss);
        
RES(1) = Yss - par.delta_k*Kss - Css - par.delta_d*Dss; 

RES(2) = par.varrho*Dss^-par.eta - Css^-par.nu*(1-par.beta*(1-par.delta_d))/par.beta;

end

%% Output (Y) and Marginal Productivity of Capital (MPK)
function [YY,MPK] = output_MPK(par,KK,ZZ)

YY  = ZZ.*KK.^par.alpha;

MPK = par.alpha*ZZ.*KK.^(par.alpha-1);

end

%% SIMULATION HOUSING MODEL, using CSD solution
function [SIM] = housing_sim_pert(par,SS,SOL,opt_sim,order)
%
% opt_sim     = structure with settings for simulation. Contains:
%                       TT      - number of periods in simulation (scalar)
%                       cols     - number of columns in simulation (scalar)
%                       T_ini   - initial periods at s.s.
%
% OUTPUT:
%       SIM         = structure with series in logs: LK (capital), 
%                       LD (housing), LZ (TFP), LC (cons.), LY (output)

cols = opt_sim.cols;
nx = size(SOL.XX_ss,2);
nz = size(SOL.ZZ_ss,2);
ny = size(SOL.YY_ss,2);

% Initialize variables (LK, LD, LC, LZ)
nms                 = {'K','D','C','Z'}; 
for in = 1:size(nms,2)
   eval(['L',nms{in},' = NaN(opt_sim.T_ini+opt_sim.TT,cols);']); 
   eval(['L',nms{in},'(1:opt_sim.T_ini,:) = log(SS.',nms{in},'ss);']);  
end   
% pre-allocate one extra column 
% for LK_t+1,LD_t+1 and LZ_t+1
LK = [LK;NaN(1,cols)];
LD = [LD;NaN(1,cols)];
LZ = [LZ;NaN(1,cols)];

% initialize shocks:
rng('default');
rng(1);
epsilon = [zeros(cols,opt_sim.T_ini),randn(cols,opt_sim.TT+1)]';

% Set K_t+1 = k_ss  and D_t+1 = d_ss
LK(opt_sim.T_ini+1,:)  = log(SS.Kss);
LD(opt_sim.T_ini+1,:)  = log(SS.Dss);

% Set LZ_t+1 = LZ(t)+eps(t+1)
LZ(opt_sim.T_ini+1,:) = get_ZZn(SOL,zeros(1,cols),epsilon(opt_sim.T_ini+1,:));

% loop over time (from T_ini +1 till T_ini + TT):
for it = opt_sim.T_ini + 1: opt_sim.T_ini + opt_sim.TT
    [XXn,LC(it,:),LZ(it+1,:)] = eval_sol_csd(SOL,[LK(it,:);LD(it,:)],LZ(it,:),epsilon(it+1,:),order);
   
    LK(it+1,:) = XXn(1,:);
    LD(it+1,:) = XXn(2,:);    
end

% remove last k_t+1 and z_t+1 (ie. adjust size to (T_ini+TT) x cols )
LK = LK(1:end-1,:);
LD = LD(1:end-1,:);
LZ = LZ(1:end-1,:);

% Add output:
SIM.LY = LZ + par.alpha*LK;

% Put variables in structure SIM:
for in = 1:size(nms,2)
   eval(['SIM.L',nms{in},' = L',nms{in},';']); 
end  

SIM.epsilon = epsilon;

end