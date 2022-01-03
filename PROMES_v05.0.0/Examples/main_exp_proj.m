function [POL,par,GRID] = main_exp_proj()
% Solve a Simple Life Cycle model with C(x)=exp(x) as solution.
% Approximation of the solution over the interval [0,3];
%
%  Uses:
%       - Promes toolbox (v05.0.0)
%       - Optimization toolbox

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

%% STEP 0: Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;
restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% Add relevant folders of Promes toolbox:
addpath ('..');
addpath ('..\grid_subfun');
addpath ('..\smolyak_subfun');


%% STEP 1: Initial block
% Set parameter:
par.nu      = 2;%risk aversion


%% STEP 2: Construct the grid
gin.nn = 1;%number of state variables
gin.lb = 0;%lower bound
gin.ub = 3;%upper bound
gin.qq = 4;%number of nodes

POL.algo    = 'cheb_gal'; 
% 'cheb_gal','cheb_tmi','cheb_mse'
% 'spl_tmi','spl_dir'                      
% 'smol_tmi','smol_dir'
% 'mono_mse'

if strncmp(POL.algo,'cheb',4) || strncmp(POL.algo,'spl',3)
    grid_spec.qq = 4;%set number of nodes
    if strncmp(POL.algo,'cheb',4)
        grid_spec.ord_vec = 3;%set order of approximation
    end
else
    grid_spec = [];
end

% Construct grid with default settings:
GRID            = prepgrid(gin.nn,gin.lb,gin.ub,POL.algo,grid_spec);

%Print grid properties on screen:
print_grid_exp(POL,GRID);


%% STEP 3: Handle to model function 
fun_res     = @(POL)res_exp(par,GRID,POL);


%% STEP 4: Initial guess
% Use third order Taylor expansion of exp(x) around x_bar: 
% and evaluate at the GRD.xx (the initial grid)
x_bar   = (GRID.ub + GRID.lb)/2;
Y0      = exp(x_bar) + exp(x_bar)*(GRID.xx - x_bar) + exp(x_bar)*(GRID.xx - x_bar).^2/2 + exp(x_bar)*(GRID.xx - x_bar).^3/6;


%% STEP 5: Solve the model
POL         = solve_proj(GRID,POL,fun_res,Y0);

% Calculate maximum error using a fine grid:
xx_acc      = linspace(GRID.lb,GRID.ub,10001)';
max_error   = max(abs(exp(xx_acc) - get_pol_var(POL,xx_acc,GRID)));

fprintf('The maximum error in C(x) is %d\n',max_error)


%% PLOT APPROXIMATION, including initial guess:
f1 = figure;

LW = 2;
%Exact solution exp(x):
plot(xx_acc,exp(xx_acc),'LineWidth',LW);
hold all;

%Projection:
plot(xx_acc,get_pol_var(POL,xx_acc,GRID),'--','LineWidth',LW);

%Taylor series (second order):
Y0_acc      = exp(x_bar) + exp(x_bar)*(xx_acc - x_bar) + exp(x_bar)*(xx_acc - x_bar).^2/2;
plot(xx_acc,Y0_acc,':','Color',[0.9290, 0.6940, 0.1250],'LineWidth',LW);

%Mark policy at gridpoints:
C_i       = get_pol_var(POL,GRID.xx,GRID);
plot(GRID.xx,C_i,'xk','LineWidth',LW);

xlabel('x');
ylabel('C(x)')

if strncmp(POL.algo,'mono',4) 
    METH = ['Monomials (order ', num2str(GRID.ord_vec(1)),')'];
elseif strncmp(POL.algo,'cheb',4)
     METH = ['Chebyshev polyn. (order ', num2str(GRID.ord_vec(1)),')'];
elseif strncmp(POL.algo,'smol',4)
    METH = ['Smolyak polyn. (accuracy $\mu = ', num2str(GRID.mu_vec(1)),'$)'];
elseif strncmp(POL.algo,'spl',3)
    METH = ['Spline (', num2str(GRID.qq(1)),' nodes)'];
end
    
legend({'Exact solution ($e^x$)',METH,'Taylor series (order 3)','Gridpoints ($x_i$)'},'Location','best','Interpreter','latex');

%print(f1, [POL.algo], '-depsc')
%save (['POL_',POL.algo],'POL','GRID','par');

end

%% Residual function Simple Life Cycle model
function [RES]  = res_exp(par,GRID,POL)

%Initial grid of state variable:
xx = GRID.xx;

% Evaluate policy function,
% at the initial grid:
if ~(strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'cheb_tmi')) 
    % standard: log(C) from policy function 
    CC     = get_pol_var(POL,xx,GRID);
else
    %time iteration: solver directly sets C_i 
    % (at gridpoint x_i) 
    CC    = POL.YY;
end

% Budget constraint gives C2:
C2 = 2*exp(xx) - CC;

% Euler residuals:
RES = (CC./C2).^-par.nu - 1;

% Alternative formulation of Euler residuals, resulting in high errors for
% high levels of consumption with overidentified polynomials (p<q)
%RES = C2.^-par.nu - CC.^-par.nu;

end

%% Print the grid
function print_grid_exp(POL,GRID)

fprintf('The initial grid (xx) is:\n')
disp(GRID.xx)
if strncmp(POL.algo,'mono',4) 
 	fprintf('The complete polynomial (XX) of order %d is:\n',GRID.ord_vec(1))
	disp(GRID.XX_poly)    
elseif strncmp(POL.algo,'cheb',4)
    fprintf('The grid of Chebyshev nodes in [-1,1] (xx_dw) is:\n')
	disp(GRID.xx_dw)
    fprintf('The complete Chebyshev polynomial (XX_dw) of order %d is:\n',GRID.ord_vec(1))
  	disp(GRID.XX_poly_dw)  
elseif strncmp(POL.algo,'smol',4)
    fprintf('The grid of Smolyak nodes in [-1,1] (xx_dw) is:\n')
	disp(GRID.xx_dw)
    fprintf('The Smolyak polynomial (XX_dw) with accuray mu %d is:\n',GRID.mu_vec(1))
  	disp(GRID.XX_poly_dw)  
end

end



