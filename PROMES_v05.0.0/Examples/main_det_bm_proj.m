function [POL,par,GRID,SS] = main_det_bm_proj()
% [POL,par,GRID,SS] = main_det_bm_proj()
%
% Function gives projection solution for deterministic version of 
% Brock-Mirman model, using default settings of Promes toolbox
%
% (For higher accuracy: set options argument in solve_proj, and for time
% iteration adjust POL.diff_tol & POL.res_tol)
%
%  Uses:
%       - Promes toolbox (v05.0.0)
%       - Optimization toolbox

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


%% STEP 1: Initial block
% STEP 1.A: Set parameters of the model
par.alpha       = 0.33;      % production: K^alpha
par.beta        = 0.96;      % discount factor

% STEP 1.B: Solve steady state
SS              = det_bm_ss(par);


%% STEP 2: Construct the grid
% STEP 2.A: Set parameters & bounds of grid, in log(capital)
gin.nn          = 1;% number of state variables

% Boundaries of grid at steady state +/- 20%:
gin.lb(1)       = -0.2 + log(SS.Kss);% lower bound
gin.ub(1)       = 0.2 + log(SS.Kss);% upper bound

% STEP 2.B: Assign algorithm to POL:
% 'cheb_gal','cheb_tmi','cheb_mse'
% 'spl_tmi','spl_dir'                      
% 'smol_tmi','smol_dir'
% 'mono_mse'
POL.algo    = 'spl_dir';

% STEP 2.C: Construct the grid
% using default settings:
GRID            = prepgrid(gin.nn,gin.lb,gin.ub,POL.algo);
clear gin;


%% STEP 3: Handle for objective function 
% (ie. the model file) 
fun_res     = @(POL)det_bm_res(par,GRID,POL);


%% STEP 4: Initial guess policy function of log(c)
% steady state consumption 
% + small linear term in [log(K)-log(Kss)]:
Y0      = log(SS.Css) + 0.01*(GRID.xx-log(SS.Kss));


%% STEP 5: Solve the model
POL         = solve_proj(GRID,POL,fun_res,Y0);
clear Y0;


%% STEP 6: Evaluate & plot policy function:
LK  = GRID.xx; % = initial grid
%Sort LK for smolyak grid:
if strncmp(POL.algo,'smol',4);
    LK = sort(LK);
end
LC  = get_pol_var(POL,LK,GRID);

% Plot policy function
figure
plot(LK,LC,'LineWidth',1.5)
xlabel('Capital (log)')
ylabel('Consumption (log)')
title('Policy function')

% Compare the numerical result with analytical result
LC_ana = log(1-par.alpha*par.beta) + par.alpha*LK;
df_LC  = LC-LC_ana;

fprintf('The maximum absolute error in C(k) is %1.2e\n',max(abs(df_LC)));

figure
plot(LK,df_LC,'LineWidth',1.5);
xlabel('Capital (log)')
ylabel('Error in C (log)')
title('Errors')

end


%% Deterministic B-M model file, or residual function:
function [RES] = det_bm_res(par,GRID,POL)

% Initial grid is stored in GRID.xx:
LK = GRID.xx;%log(K_t)

% Evaluate the policy function log(C):
if ~(strcmp(POL.algo,'spl_tmi') || ...
    strcmp(POL.algo,'smol_tmi') || ...
    strcmp(POL.algo,'cheb_tmi')) 
    % standard format:
    LC     = get_pol_var(POL,LK,GRID);
else
    %tmi: solver directly sets log(C) to min. residuals
    LC    = POL.YY;
end

% Capital in next period (log): 
LK_n = log( exp(par.alpha*LK) - exp(LC) );

% log(C_t+1) from policy function, given log(K_t+1): 
if ~(strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'cheb_tmi')) 
    LC_n  	= get_pol_var(POL,LK_n,GRID);
    
else    
    % use old policy function in t+1 
    % for Time Iteration 
    % (pp_y_old or theta_old)
	spec_opt_next = 'old_pol'; 
	LC_n  	= get_pol_var(POL,LK_n,GRID,[],spec_opt_next); 
end 

% log(RR_t+1): marginal prod. of capital (logs)
LR_n = log(par.alpha) + (par.alpha-1)*LK_n;

% RHS of Euler equation:
RHS         = par.beta * exp(-LC_n) .* exp(LR_n) ;
   
% Euler equation (multiplied by C):
RES         =  exp(LC).*RHS - 1;

end


%% Analytical steady state:
function [SS] = det_bm_ss(par)

SS.Kss  	= (par.alpha*par.beta)^(1/(1-par.alpha));
SS.Css  	= SS.Kss^par.alpha - SS.Kss;

end
