% Sijmen Duineveld, updated May 2021, s.a.duineveld@outlook.com

clearvars;
close all;
clc;
dbstop if error;

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

%Add folder with CSD and subfun:
addpath ('..');
addpath ('..\subfun');

%% STEP 1: GET SYMBOLIC MODEL
[MOD]       = SGU_mod;


%% STEP 2: SET PARAMETERS
par.alpha        = 0.3;         % coefficient on capital in production function
par.beta         = 0.95;        % discount rate
par.delta        = 1;       % depreciation
par.nu           = 2;           % risk aversion

par.rho_a       = 0;%0.95;%0;
par.sigma_a     = 1;

par.opt.order = 3;%order of solution

% Assign the numerical values to MOD structure:
MOD.par_val = [par.alpha,par.beta,par.delta,par.nu];


%% STEP 3: SET STEADY STATE VALUE
kappa 	= (1-par.beta*(1-par.delta))/(par.alpha*par.beta);
SS.Kss	= kappa ^ (1/(par.alpha-1));        
SS.Css	= SS.Kss^par.alpha - par.delta * SS.Kss;
clear kappa;

SS.Ass	= 1;
SS.Yss  = SS.Ass*SS.Kss^par.alpha;

% Assign S.S. values to the MOD structure :
MOD.SS_vec      = [log(SS.Kss),log(SS.Ass),log(SS.Css)];


%% STEP 4: SOLVE MODEL NUMERICALLY
[SOL,NUM,MOD] = pert_ana_csd(MOD,par.rho_a,par.opt.order,par.sigma_a);


%% STEP 5: EVALUATE POLICY FUNCTION: 1st order IRF
%IRF: 1 standard deviation shock in period 3:
ini_T   = 3;%shock in period 3
TT      = 15;
LK = NaN(ini_T+TT+1,1);
LA = NaN(ini_T+TT+1,1);
LC = NaN(ini_T+TT,1);

%Starting values:
LK(1:3,1)   = log(SS.Kss);
LC(1:2,1)   = log(SS.Css);
LA(1:2,1)   = log(SS.Ass);
%Unexpected shock in period 3:
LA(3,1)     = get_ZZn(SOL,0,1);

for it = ini_T:ini_T+TT
    [LK(it+1,:),LC(it,:),LA(it+1,:)] = eval_sol_csd(SOL,LK(it,:),LA(it,:),0,1);    
end
%Remove last row for LK and LA:
LK = LK(1:end-1,:);
LA = LA(1:end-1,:);

%Output:
LY = LA+par.alpha*LK;

%Plot IRF:
f1 = figure;
plot(LY - log(SS.Yss),'LineWidth',1.5);
xlabel('Time')
ylabel('Output (% dev.)')
xlim([1 ini_T+TT]);
title('Impuls Response Function');

%print(f1,'SGU_irf_y', '-depsc')


