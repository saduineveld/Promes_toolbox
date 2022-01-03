function [POL,par,GRID,SS] = main_det_bm_fp_algo()
% [POL,par,GRID,SS] = main_det_bm_proj()
%
% Function gives projection solution for deterministic version of 
% Brock-Mirman model, using fixed point algorithm



%% STEP 0: Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% Add relevant folders of Promes toolbox:
%addpath ('..');
addpath ('..\grid_subfun');
addpath ('..\smolyak_subfun');

addpath ('toolbox_fpa');


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
% 'cheb_fpa','spl_fpa'                      
% 'smol_fpa'
POL.algo    = 'smol_fpa';

% STEP 2.C: Construct the grid
% using default settings:
GRID            = prepgrid_fpa(gin.nn,gin.lb,gin.ub,POL.algo);
clear gin;


%% STEP 3: Handle for objective function 
% (ie. the model file) 
fun_y_new     = @(POL)det_bm_pol_new(par,GRID,POL);


%% STEP 4: Initial guess policy function of log(c)
% steady state consumption 
% + small linear term in [log(K)-log(Kss)]:
Y0      = log(SS.Css) + 0.2*(GRID.xx-log(SS.Kss));


%% STEP 5: Solve the model
POL.mem_Y = 0.5;
[POL,ex_fl] = solve_proj_fpa(GRID,POL,fun_y_new,Y0);



%% STEP 6: Evaluate & plot policy function:
LK  = GRID.xx; % = initial grid
%Sort LK for smolyak grid:
if strncmp(POL.algo,'smol',4);
    LK = sort(LK);
end
LC  = get_pol_var_fpa(POL,LK,GRID);

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
function [LC_new] = det_bm_pol_new(par,GRID,POL)

% Initial grid is stored in GRID.xx:
LK = GRID.xx;%log(K_t)

% Evaluate the policy function log(C):
if strcmp(POL.algo,'spl_fpa') || ...
    strcmp(POL.algo,'smol_fpa') || ...
    strcmp(POL.algo,'cheb_fpa') 

    % use old policy for RHS of Euler:
    LC1     = POL.YY_old;
end


% Capital in next period (log): 
LK_n = log( exp(par.alpha*LK) - exp(LC1) );

% log(C_t+1) from policy function, given log(K_t+1): 
if strcmp(POL.algo,'spl_fpa') || ...
    strcmp(POL.algo,'smol_fpa') || ...
    strcmp(POL.algo,'cheb_fpa') 

    % use old policy function pp_y_old  (theta_old) 
    LC_n  	= get_pol_var_fpa(POL,LK_n,GRID,[],'old_pol');
end 

% log(RR_t+1): marginal prod. of capital (logs)
LR_n = log(par.alpha) + (par.alpha-1)*LK_n;

% RHS of Euler equation:
RHS         = par.beta * exp(-LC_n) .* exp(LR_n) ;
   
% Euler equation (multiplied by C):
%RES         =  exp(LC).*RHS - 1;

LC_new = -log(RHS); 

end


%% Analytical steady state:
function [SS] = det_bm_ss(par)

SS.Kss  	= (par.alpha*par.beta)^(1/(1-par.alpha));
SS.Css  	= SS.Kss^par.alpha - SS.Kss;

end
