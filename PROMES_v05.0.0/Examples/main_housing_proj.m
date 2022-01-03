function [POL,par,GRID,SS] = main_housing_proj()
% [POL,par,GRID,SS] = main_housing_proj()
% 
% Solves an RBC model with housing
%
% Uses:
%       - Promes toolbox (v05.0.0)
%       - Optimization toolbox
%       - hernodes (in folder TOOLS)
%       - CSD_v02.4.0 (optional, in folder TOOLS)
%
% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com


%% STEP 0: Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% Add relevant folders of Promes toolbox:
addpath ('..');
addpath ('..\grid_subfun');
addpath ('..\smolyak_subfun');

% Add folder TOOLS
addpath ('..\..\TOOLS');
addpath ('..\..\TOOLS\CSD_v02.4.0');
addpath ('..\..\TOOLS\CSD_v02.4.0\subfun');

tic;

%% STEP 1: Initial block
% STEP 1.A: Set parameters of the model
par.alpha       = 0.36;     % coefficient on capital in production function
par.beta        = 0.985;    % discount rate

par.delta_k     = 0.025;    % depreciation capital
par.delta_d     = 0.01;     % depreciation housing  

par.eta         = 5;        % risk aversion for housing
par.nu          = 2;      % risk aversion for consumption
  
par.varrho      = 0.1;      % weight of housing utility

par.rho_z       = 0.95;	% autocorrelation coefficient TFP
par.sigma_z     = 0.01;	% standard dev. shocks in TFP 

par.her.gh_nod  = 5;% number of Gauss-Hermite nodes
[par.her.xi,par.her.wi] = hernodes(par.her.gh_nod);
% xi are roots, wi are weights

% STEP 1.B: Solve steady state
[SS] = housing_ss(par,1);


%% STEP 2: Construct the grid
% STEP 2.A: Set parameters & bounds of grid, in log(capital)
gin.nn          = 3;% # of state variables

gin.lb(1)       = -0.2 + log(SS.Kss);% lower bound capital
gin.ub(1)       = 0.2 + log(SS.Kss);% upper bound capital

gin.lb(2)       = -0.07 + log(SS.Dss);% lower bound housing
gin.ub(2)       = 0.07 + log(SS.Dss);% upper bound housing

gin.lZ_fac  = 3;% 
gin.lb(3)   = -gin.lZ_fac*sqrt( par.sigma_z^2 / (1-par.rho_z^2) ); 
gin.ub(3)   =  gin.lZ_fac*sqrt( par.sigma_z^2 / (1-par.rho_z^2) );

% STEP 2.B: Set algorithm: 
% 'cheb_gal','cheb_tmi','cheb_mse'
% 'spl_tmi','spl_dir'                      
% 'smol_tmi','smol_dir'
% 'mono_mse'
POL.algo    = 'spl_tmi';

% STEP 2.C: Construct the grid
% use default settings
GRID    = prepgrid(gin.nn,gin.lb,gin.ub,POL.algo);
clear gin;


%% STEP 3: Handle for objective function 
% (ie. the model file) 
fun_res     = @(POL)housing_proj(par,GRID,POL); 


%% STEP 4: Initial guess policy function of log(K') and log(D')
% Use first order perturbation as initial guess:
Y0 = set_Y0_housing(par,SS,GRID,0);


%% STEP 5: Solve the model

POL         = solve_proj(GRID,POL,fun_res,Y0);
toc;

end

function [RES] = housing_proj(par,GRID,POL)
% [RES] = housing_proj(par,GRID,POL)
%
% Calculates dd=2 sets of Euler residuals for housig model
%
% First policy variable: K_t+1
% Second policy variable: C_t

LK = GRID.xx(:,1);%first state variable,    log(K_t)
LD = GRID.xx(:,2);%second state variable,   log(D_t)
LZ = GRID.xx(:,3);%third state variable,    log(Z_t)

%policy variables, log(K_t+1) and  log(C_t):
if ~(strcmp(POL.algo,'spl_tmi') || ...
    strcmp(POL.algo,'smol_tmi') ||...
    strcmp(POL.algo,'cheb_tmi') )

    %use initial grid for polynomials 
    % (or ignore for spl_dir)
    spec_opt    = 'ini_grid';     
    % fourth entry is index for policy variable (i_pol)
    LK_n  	= get_pol_var(POL,[LK,LD,LZ],GRID,1,spec_opt);
    LC 	= get_pol_var(POL,[LK,LD,LZ],GRID,2,spec_opt);

else%for 'spl_tmi','cheb_tmi','smol_tmi':
    
 	%Policy variables in columns of POL.YY for 'tmi':
    LK_n    = POL.YY(:,1);
    LC      = POL.YY(:,2);
end

%Housing in next period, marg. ut. of C, and marg. ut. housing:
[LD_n,~,~,lambda,dVdD_n]    = housing_aux(par,LK,LD,LZ,LK_n,LC);

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
rhs_j1   = NaN(size(LK,1),par.her.gh_nod);
rhs_j2   = NaN(size(LK,1),par.her.gh_nod);
for jj = 1:par.her.gh_nod
    % Shock to TFP (using Gauss-Hermite nodes):
    EPS_n       = sqrt(2)*par.her.xi(jj);   
    %log(Z_t+1):
    LZ_n        = par.rho_z*LZ + par.sigma_z*EPS_n;
    
    LK_n2        = get_pol_var(POL,[LK_n,LD_n,LZ_n],GRID,1,spec_opt_next);
    LC_n        = get_pol_var(POL,[LK_n,LD_n,LZ_n],GRID,2,spec_opt_next);
    
    [~,LMPK_n,~,lambda_n]    = housing_aux(par,LK_n,LD_n,LZ_n,LK_n2,LC_n);    
   
    % RHS of Euler equation, 
    % weighted by Gauss-Hermite weights
    rhs_j1(:,jj) = par.her.wi(jj)/sqrt(pi) * par.beta * lambda_n .* (exp(LMPK_n) + 1 - par.delta_k);
    
    rhs_j2(:,jj) = par.her.wi(jj)/sqrt(pi) * par.beta * (dVdD_n + lambda_n*(1 - par.delta_d));
end

% Euler residuals 
RES1    = sum(rhs_j1,2)./lambda - 1;
RES2 	= sum(rhs_j2,2)./lambda - 1;

% RES1 and RES2 are mm by 1 vectors, stacked vertically:
RES = [RES1;RES2];
 
end

%% Auxiliary variables housing model
function [LD_n,LMPK,LY,lambda,dVdD_n] = housing_aux(par,LK,LD,LZ,LK_n,LC)

[YY,MKP] = output(par,exp(LK),exp(LZ));

II      = exp(LK_n) - (1-par.delta_k)*exp(LK);  %capital investment
HH      = YY - II - exp(LC);%Housing investment
LD_n    =  log( (1-par.delta_d)*exp(LD) + HH );

if nargout > 1
    LMPK = log(MKP);
end
if nargout > 2
    LY = log(YY);
end
if nargout > 3
    lambda  = exp(-par.nu*LC);
end
if nargout > 4
    dVdD_n   = par.varrho*exp(-par.eta*LD_n);
end

end

%% Output
function [YY,MPK] = output(par,KK,ZZ)

YY  = ZZ.*KK.^par.alpha;

MPK = par.alpha*ZZ.*KK.^(par.alpha-1);

end

%% Housing model for use with CSD_PERT_SOLVER
function [MOD] = HOUSING_pert()

syms aa bb dk dh ee nu varrho rhoz sgmz;

syms LK_t LK_n LD_t LD_n  LZ_t LZ_n LC_t LC_n;

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

% Budget / law of motion assets:
f3 = YY + (1-dk)*exp(LK_t) - exp(LK_n)  + (1-dh)*exp(LD_t) - exp(LD_n) - exp(LC_t);

MOD.FS = [f1;f2;f3];

MOD.XX = [LK_t,LD_t];
MOD.ZZ = [LZ_t];
MOD.YY = LC_t;

MOD.XXn = [LK_n,LD_n];
MOD.ZZn = [LZ_n];
MOD.YYn = LC_n;

MOD.var_bs_nms  = {'LK','LD','LZ','LC'};
MOD.par_nms     = {'aa','bb','dk','dh','ee','nu','varrho','rhoz','sgmz'};

end

%% Steady state
function [SS] = housing_ss(par,Zss)

SS.Zss  = Zss;

SS.Kss  = (SS.Zss*par.alpha*par.beta/(1-par.beta*(1-par.delta_k)))^(1/(1-par.alpha));

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

Yss  = output(par,Kss,Zss);
        
RES(1) = Yss - par.delta_k*Kss - Css - par.delta_d*Dss; 
RES(2) = par.varrho*Dss^-par.eta - Css^-par.nu*(1-par.beta*(1-par.delta_d))/par.beta;

end

%% Solve the model with perturbation methods
function [PERT] = get_pert_sol(par,SS)
% (CSD perturbation solver needs to be on the searchpath)

MOD         = HOUSING_pert;%get symbolic model;

% Numerical values of parameters:
MOD.par_val = [par.alpha,par.beta,...
        par.delta_k,par.delta_d,par.eta,par.nu,par.varrho,par.rho_z,par.sigma_z];

% Vector of steady states:
%MOD.var_bs_nms  = {'LK','LD','LZ','LC'};
MOD.SS_vec  = [log(SS.Kss),log(SS.Dss),log(SS.Zss),log(SS.Css)];
    
% Get solution:
PERT         = pert_ana_csd(MOD,par.rho_z,1,...
        par.sigma_z);
    
end

%% Set initial policy function
function [Y0] = set_Y0_housing(par,SS,GRID,opt_run_CSD)

if opt_run_CSD == 1
    % Get perturbation solution:
    [PERT] = get_pert_sol(par,SS);
    
else
    % Solution for given parameters:
    PERT.Hx_w(1,1) = 0.9608;
    PERT.Hx_w(1,2) = 0.0540;
    PERT.Hx_w(1,3) = 0.0829;
    
    PERT.Hy_w(1,1) = 0.4722;
    PERT.Hy_w(1,2) = 0.0266;
    PERT.Hy_w(1,3) = 0.3865;
end

%pre-allocate dimensions:
Y0 = NaN(GRID.mm,2);

Y0(:,1) = log(SS.Kss) + PERT.Hx_w(1,1)*(GRID.xx(:,1)-log(SS.Kss)) + ...
                      PERT.Hx_w(1,2)*(GRID.xx(:,2)-log(SS.Dss)) + ...
                      PERT.Hx_w(1,3)*(GRID.xx(:,3)-log(SS.Zss));
                    
Y0(:,2) = log(SS.Css) + PERT.Hy_w(1,1)*(GRID.xx(:,1)-log(SS.Kss)) + ...
                      PERT.Hy_w(1,2)*(GRID.xx(:,2)-log(SS.Dss)) + ...
                      PERT.Hy_w(1,3)*(GRID.xx(:,3)-log(SS.Zss));

end


