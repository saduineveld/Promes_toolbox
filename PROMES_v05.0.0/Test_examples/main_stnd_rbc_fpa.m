% Solves standard RBC model with projection
% for a single variable policy function
% 
%  Uses: 
%    - Promes toolbox (v05.0)
%    - Optimization toolbox
%    - hernodes (in folder TOOLS) 
%    - models functions: STND_RBC_proj, stnd_rbc_ss, stnd_rbc_aux, 
%       stnd_rbc_sim, plot_pol_stnd_rbc, STND_RBC_proj_EEE
%       (all in subfolder STND_RBC_mod)
% 
%    For perturbation solution (optional): 
%    - CSD_vXX.Y.Z (in TOOLS folder), 
%    - STND_RBC_pert (in subfolder STND_RBC_mod)

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

%% STEP 0: Matlab settings
clearvars;
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% Add relevant folders of Promes toolbox:
%addpath ('..');
addpath ('toolbox_fpa');
addpath ('..\grid_subfun');
addpath ('..\smolyak_subfun');

% Add folder TOOLS
addpath ('..\..\TOOLS');
addpath ('..\..\TOOLS\CSD_v02.4.0');
addpath ('..\..\TOOLS\CSD_v02.4.0\subfun');

% Add folder with model files:
addpath ('STND_RBC_mod_fpa');
addpath ('..\Examples\STND_RBC_mod');

tic;

%% STEP 1: INITIAL BLOCK 
% STEP 1.A: Set parameters of the model
par.alpha       = 0.36;		%capital share income 
par.beta        = 0.985;	% discount factor
par.delta       = 0.025;	% deprec. of capital 
par.nu          = 2; 	% risk aversion 
par.eta         = 4; 	% el. of lab. supply
par.chi         = 1; 	% scalar disut. work

par.rho_z       = 0.95;	% autocorr. coeff. TFP
par.sigma_z     = 0.01;	% standard dev. shocks in TFP 

par.her.gh_nod  = 5;% number of Gauss-Hermite nodes
[par.her.xi,par.her.wi] = hernodes(par.her.gh_nod);
% xi are roots, wi are weights

% STEP 1.B: Solve steady state
SS              = stnd_rbc_ss(par,1);
% the 1 is steady state TFP


%% STEP 2: Construct the grid
%Step 2.A: Set grid parameters
gin.nn      = 2;%number of state variables

%Set lower and upper bound for capital:
gin.lK_dev  = 0.1275;% deviation from Kss
gin.lb(1)   = -gin.lK_dev + log(SS.Kss); 
gin.ub(1)   =  gin.lK_dev + log(SS.Kss);

% Set lower and upper bound for log(Z) 
gin.lZ_fac  = 2.6;%in multiple of stnd. deviation
gin.lb(2)   = -gin.lZ_fac*sqrt( par.sigma_z^2 / (1-par.rho_z^2) ); 
gin.ub(2)   =  gin.lZ_fac*sqrt( par.sigma_z^2 / (1-par.rho_z^2) );


% STEP 2.B: Set algorithm: 
% 'cheb_fpa','spl_fpa'                      
% 'smol_fpa'
POL.algo    = 'cheb_fpa';

% (OPTIONAL) Algorithm specific parameters
%CHEBYSHEV:
if strncmp(POL.algo,'cheb',4)
    %order of polyn. in each dim.
    algo_spec.ord_vec   = 5*ones(1,gin.nn);
    %# nodes in each dim.:
    algo_spec.qq        = algo_spec.ord_vec+1;
    
%SPLINE:
elseif strncmp(POL.algo,'spl',3)
    %# nodes in each dim.:
    algo_spec.qq        = 7*ones(1,gin.nn);
    
%SMOLYAK:
elseif strncmp(POL.algo,'smol',4)
    %accuracy param. in each dim.:
    algo_spec.mu_vec    = 3*ones(1,gin.nn);
    
%MONOMIALS:
elseif strncmp(POL.algo,'mono',4)
    %order of polyn. in each dim.:
    algo_spec.ord_vec   = 3*ones(1,gin.nn);
    %# nodes in each dim.:
    algo_spec.qq        = algo_spec.ord_vec+1;
    
else
    error('Invalid algo');
end
     
% STEP 2.C: Construct the grid
[GRID] = prepgrid_fpa(gin.nn,gin.lb,gin.ub,POL.algo,algo_spec);
clear algo_spec;


%% STEP 3: Handle for objective function 
% (ie. the model file) 
fun_new_pol     = @(POL)STND_RBC_fpa(par,GRID,POL);    


%% STEP 4: Initial guess for policy function for log(C) (Y0)
par.opt.get_pert_sol = 0;
% 0: use pre-determined linear perturbation solution (for given parameters);
% 1: calculate linear perturbation solution

if par.opt.get_pert_sol == 0    
    %For the given parameters the perturbation solution is:
    %PERT.Hy_w = [0.3456,0.3525];
    
    %Poor estimation:
    PERT.Hy_w = [0.25,0.25];
else% Solve model with perturbation     
    % (CSD perturbation toolbox needs to be on the searchpath)
    
    % Symbolic model file: 
    MOD         = STND_RBC_pert;

    % Vector of parameters:    
    MOD.par_val = [par.alpha,par.beta,...
        par.delta,par.eta,par.nu,par.chi];
    
    % Vector of steady states:
    MOD.SS_vec  = [log(SS.Kss),log(SS.Zss),log(SS.Css)];
    
    % Get solution:
    PERT         = pert_ana_csd(MOD,par.rho_z,1,...
        par.sigma_z);
    
    clear MOD;
end

PERT.Hy_w

%Initial policy function (for given grid):
Y0 = log(SS.Css) + ...
    PERT.Hy_w(1,1)*(GRID.xx(:,1)-log(SS.Kss)) + ...
    PERT.Hy_w(1,2)*(GRID.xx(:,2)-log(SS.Zss));


%% STEP 5: Solve the model
% Solve model
%POL.mem_Y   = 0.75;
POL         = solve_proj_fpa(GRID,POL,fun_new_pol,Y0);
toc; 
clear fun_res Y0;

% Plot policy function (2D and 3D):
plot_pol_stnd_rbc_fpa(par,SS,GRID,POL)


%% Step 6: Evaluate policy function (in simulation)
%Simulate the model
opt_sim.TT      = 1000; %number of periods in simulation
opt_sim.T_ini   = 10;   %initial periods at deterministic steady state
opt_sim.rws     = 2;   %number of simulated series

[SIM]           = stnd_rbc_sim_fpa(par,SS,POL,GRID,opt_sim);

% Plot two simulated series of output
figure;
opt_sim.T_plot = 200;%Plot last 50 years in simulatio
plot([1:opt_sim.T_plot]/4,SIM.LY(1,end-opt_sim.T_plot+1:end)-log(SS.Yss),'LineWidth',1.5);
hold all;
plot([1:opt_sim.T_plot]/4,SIM.LY(2,end-opt_sim.T_plot+1:end)-log(SS.Yss),'LineWidth',1.5);
%plot([1:opt_sim.T_plot]/4,repmat(log(SS.Yss),1,opt_sim.T_plot),':','LineWidth',1.5);
xlabel('Time (years)');
ylabel('Output (% dev.)');
legend({'Simulation 1','Simulation 2'},'Location','best');

%% Compute normalized Euler Equation Error (EEE)
% on initial grid
% and on equidistant grid with 1e6 nodes

% Compute EEE:
[ERR_ig_bs] = STND_RBC_proj_EEE_fpa(par,GRID.xx,GRID,POL);
ERR_ig = log10(abs(ERR_ig_bs));
fprintf('\nMax. norm. Euler Eq. Error on ini. grid is: %1.2f (in log10) \n',max(ERR_ig(:)))

% Construct equidistant grid:
G_EQ.qq         = 1000*ones(1,GRID.nn);
G_EQ.gridVecs   = constr_vecs(G_EQ.qq,'equi','up',GRID.lb,GRID.ub);
G_EQ.xx         = constr_grid(G_EQ.gridVecs);

% Compute EEE:
[ERR_bs] = STND_RBC_proj_EEE_fpa(par,G_EQ.xx,GRID,POL);

%Average, min and max (in log10 & abs)
ERR_avg = log10(mean(abs(ERR_bs)));
ERR = log10(abs(ERR_bs));
ERR_min = min(ERR(:));
ERR_max = max(ERR(:));

fprintf('\nMax. norm. Euler Eq. Error on equid. grid is: %1.2f (in log10) \n',ERR_max)
fprintf('Avg. norm. Euler Eq. Error on equid. grid is: %1.2f (in log10)\n',ERR_avg)
