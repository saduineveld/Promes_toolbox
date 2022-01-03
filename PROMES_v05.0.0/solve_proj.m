function [POL,ex_fl] = solve_proj(GRID,POL,fun_res,Y0,options)
%[POL,ex_fl] = solve_proj(GRID,POL,fun_res,Y0,options)
%
% General solver using projection, for multiple policy functions
%
% INPUT:
%       GRID     = structure with required properties (as set in
%                   'prepgrid')   
%
% 
%       POL 	= structure for policy function. Should contain fields:
%                   - algo: algorithm (see documentation). Options:
%                       Chebyshev poly: 'cheb_gal','cheb_tmi', 'cheb_mse'
%                       Spline:         'spl_tmi', 'spl_dir'     
%                       Smolyak:        'smol_tmi','smol_dir'     
%                       Monomial:       'mono_mse'
%
%                   For algo with time iteration (see documentation):
%                   - diff_tol (optional): tolerance for difference between
%                       two iterations (YY & Y_old) (default is 1e-8)
%                   - res_tol (optional): tolerance for maxim residual of
%                       Euler equation (default is 1e-8)
%                   - max_iter (optional): maximum number of iterations
%                   (default 500)
%                   - step_acc (optional): all tolerances are scaled
%                       with step_acc (default 0.1) when the solver stalls
%                   - mem_Y (optional): dampening/memory in updating Y (default 0)
%          
%
%       fun_res = handle of model specific residual function
%
%       Y0      = mm x dd matrix, with each column the initial guess for policy function 
%
%       options = (optional) structure with options for solvers (lsqnonlin and fsolve)
%                   if options.override_all == 1  
%                       use options: opt_solver = options.optimoptions;
%                   else
%                       only specified fields in options will replace default options
%                   
%                      
%
% OUTPUT:
%       POL: structure with policy function. Added fields are:
%                   For algo 'cheb_xxx', 'smol_xxx', 'mono_xxx':
%                   - theta, matrix of polynomials coefficients (pp x dd)%                           
%                   For algo 'spl_xxx':
%                   - pp_y: a spline as created by griddedInterpolant
%                   For all solution types:
%                   - Y0: initial guess at gridpoints (mm x dd)
%                   - dd: number of policy variables
%                   - pp: length of policy function
%                           for 'spl_xxx':  pp = GRID.mm;
%                           for 'cheb_xxx': pp = size(GRID.XX_poly_dw,2)
%                           for 'smol_xxx': pp = size(GRID.XX_poly_dw,2)
%                           for 'mono_mse': pp = size(GRID.XX_poly,2)
%
% Uses:     
%       - lsqnonlin (optimization toolbox) for 'mono_mse' and 'cheb_mse'
%       - fsolve (optimization toolbox) for all other methods

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

if nargin < 5 || isempty(options)
    options.override_all = 0; 
elseif ~isfield(options,'override_all')
    options.override_all = 0;
end

if  options.override_all == 1
    opt_solver = options.optimoptions;
else
    % Set options for solver in optimoptions structure:
    opt_solver = set_default_opt_solver(options,POL.algo);
end

% Check if Y0 has correct format:
if GRID.mm ~= size(Y0,1) 
    error('Y0 has to be mm x dd matrix');
end
POL.dd = size(Y0,2);%number of policy variables;

POL.Y0 = Y0; clear Y0;

% Get length of parameter vector (pp) of policy function 
POL.pp = get_length_para(POL.algo,GRID);

% For splines only: set default interpolation method:
if strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'spl_dir') 
    if ~isfield(POL,'spl_meth')
        POL.spl_meth = 'spline';
    end
end

% For polynomials only (but not 'smol_tmi' or 'cheb_tmi'):
% Set initial policy functon (theta0) 
if strcmp(POL.algo,'mono_mse') || strcmp(POL.algo,'cheb_gal') ||...
        strcmp(POL.algo,'cheb_tmi') || strcmp(POL.algo,'cheb_mse') ...
        || strcmp(POL.algo,'smol_dir')
    POL.theta0 = get_theta_ini_grid(GRID,POL,POL.Y0);
end


% Solve the policy functions for each type
%% Time iteration
if strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'cheb_tmi') 
    if ~isfield(POL,'diff_tol')
        POL.diff_tol = 1e-8;%default maximum difference between YY & Y_old
    end
    if ~isfield(POL,'res_tol')
        POL.res_tol = 1e-8;%default maximum tolerance of residuals 
    end
    if ~isfield(POL,'max_iter')
        POL.max_iter = 500;%default max. number of iterations;
        
    end
    if ~isfield(POL,'step_acc')
        POL.step_acc = 0.1;%default change in accuracy, when y stalls (ex_fl = 0)
    end
   	if ~isfield(POL,'mem_Y')
        POL.mem_Y = 0;%default memory for updating Y
    else
        if POL.mem_Y > 1 || POL.mem_Y < 0
            error('mem_Y should be between 0 and 1');
        end
    end
    
    max_diff    = inf;
    max_res     = inf;
    
    % Initial guess for policy function
    Y_old           = POL.Y0;
    
    % Initialize 'old' policy function
    % (for next period choices)
    if strcmp(POL.algo,'spl_tmi')
        POL.pp_y_old    = get_pp_y(GRID.gridVecs,Y_old,GRID.qq,POL.spl_meth);
    elseif strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'cheb_tmi') 
        POL.theta_old = get_theta_ini_grid(GRID,POL,Y_old);
    end
    
    % Use pattern in the Jacobian 
    
    if ~isfield(options,'JacobPattern')        
        PP = repmat(speye(GRID.mm),POL.dd);
        opt_solver = optimoptions(opt_solver,'JacobPattern',PP);
    end
         
    iter = 1;
    while (max_diff > POL.diff_tol || max_res > POL.res_tol) && iter < POL.max_iter 

        % Update policy, given old policy:        
        [POL.YY,max_diff,~,ex_fl] = fun_solve_tmi(Y_old,POL,fun_res,opt_solver);        
        
        if ex_fl == 99%solved at initial point: 
           %reduce tolerances in fsolve by factor POL.step_acc    
                
           
                acc_step = POL.step_acc*opt_solver.StepTolerance;
                acc_opt = POL.step_acc*opt_solver.OptimalityTolerance;
                acc_fun = POL.step_acc*opt_solver.FunctionTolerance;
                opt_solver = optimoptions(opt_solver,'FunctionTolerance',acc_fun,...
                    'OptimalityTolerance',acc_opt,'StepTolerance',acc_step);
        end        
        Y_old   = POL.YY;        
        
        %Update 'old' policy function 
        % (for next period choices) 
        if strcmp(POL.algo,'spl_tmi')
            POL.pp_y_old = get_pp_y(GRID.gridVecs,Y_old,GRID.qq,POL.spl_meth);
        elseif strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'cheb_tmi') 
            POL.theta_old = get_theta_ini_grid(GRID,POL,Y_old);            
        end
        
        % Check residuals, after updating policy:         
        % (only when max_diff < POL.diff_tol) 
        if max_diff < POL.diff_tol;
            RES         = fun_res(POL);
            max_res     = max(abs(RES));
        else
            max_res = [];
        end
        %fprintf('Iteration %d completed. Max. difference is %0.8f. Max. res. is %0.8f. \n',iter,max_diff,max_res)
        
        iter    = iter + 1;      
    end
    
    if iter == POL.max_iter  
        
        ex_fl = 0;%no solution found
    else         
        ex_fl = 1;     
    end
    
%% Spline, Direct Computation
elseif strcmp(POL.algo,'spl_dir') 
     
    % Handle to objective function: 
    f_obj               = @(YY_col)fun_res_spl_dir(YY_col,GRID,POL,fun_res);
    
    Y0_col              = reshape(POL.Y0,[],1);
    
    % Solve system of equations: 
    [YY_col,~,ex_fl]    = fsolve(f_obj,Y0_col,opt_solver);
    
    POL.YY              = reshape(YY_col,[],POL.dd);
    if ex_fl < 1
        error('No proper solution found');
    end
%% Smolyak, Direct Computation
elseif strcmp(POL.algo,'smol_dir') 
     
    % Handle to objective function: 
    f_obj               = @(theta_col)fun_res_smol_dir(theta_col,GRID,POL,fun_res);
    
    theta0_col          = reshape(POL.theta0,[],1);
    
    % Solve system of equations:     
    [theta_col,~,ex_fl]    = fsolve(f_obj,theta0_col,opt_solver);
    
    POL.theta = reshape(theta_col,[],POL.dd);
    
    if ex_fl < 1
        error('No proper solution found');
    end
end

if strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'spl_dir')% Construct splines (in cell array) from YY
    POL.pp_y = get_pp_y(GRID.gridVecs,POL.YY,GRID.qq,POL.spl_meth);
elseif strcmp(POL.algo,'cheb_tmi') || strcmp(POL.algo,'smol_tmi')% Regress to get theta from YY
    POL.theta      = get_theta_ini_grid(GRID,POL,POL.YY);
end

%% Minimized of squared errors (MSE and Mono)
if strcmp(POL.algo,'mono_mse') || strcmp(POL.algo,'cheb_mse') 
    % Handle to MSE residual function: 
    f_res_lsq  = @(theta_col)fun_res_lsq(POL,fun_res,theta_col);
    
    theta0_col = reshape(POL.theta0,[],1);
    
    % Minimize squared errors:
    [theta_col,~,~,ex_fl] = lsqnonlin(f_res_lsq,theta0_col,[],[],opt_solver);
    
    POL.theta = reshape(theta_col,[],POL.dd);
    
%% GALERKIN:
elseif strcmp(POL.algo,'cheb_gal')     

    % Handle to Galerkin residual function: 
    f_res_gal = @(theta_col)fun_res_gal(POL,fun_res,GRID.XX_poly_dw,theta_col);
    
    theta0_col = reshape(POL.theta0,[],1);
    
    % Solve residual function:   
    [theta_col,~,ex_fl]    = fsolve(f_res_gal,theta0_col,opt_solver);
    
    POL.theta = reshape(theta_col,[],POL.dd);
end

if ex_fl < 1
    
    warning('No proper solution found');
end


end

%% Update the policy function for Time Iteration;
function [YY,max_diff,RES,ex_fl,output] = fun_solve_tmi(Y_old,POL,fun_res,opt_solver)
% Solve equations (by setting YY), given next periods policy (in pp_y_old):

% Make column vector:
Y_old_col = reshape(Y_old,[],1);

% Handle to model file with YY_col as input:
f_res_tmi = @(YY_col)fun_res_tmi(POL,fun_res,YY_col);

% Solve equations:
[YY_col,RES,ex_fl,output] = fsolve(f_res_tmi,Y_old_col,opt_solver);
       
if output.iterations == 0
	% solved at initial point: 
	% rerun with better accuracy
	ex_fl = 99;
end
if ex_fl < 1;
   % output.message
	error('No proper solution found')
end

if POL.mem_Y ~= 0
    YY_new = POL.mem_Y*Y_old_col + (1-POL.mem_Y)*YY_col;
    
    %Reshape YY_col:
    YY      = reshape(YY_new,[],POL.dd);
else
    YY      = reshape(YY_col,[],POL.dd);
end

max_diff  = max(abs(YY_col-Y_old_col));

end

%% Residual function for time iteration 
function [RES] = fun_res_tmi(POL,fun_res,YY_col)
% YY_col: (dd*mm) x 1 vector

POL.YY      = reshape(YY_col,[],POL.dd);

% Residual vector from model file:
RES         = fun_res(POL);

end


%% Residual function for Direct Computation
function [RES] = fun_res_spl_dir(YY_col,GRID,POL,fun_res)

YY = reshape(YY_col,[],POL.dd);

if strcmp(POL.algo,'spl_dir')
	POL.pp_y = get_pp_y(GRID.gridVecs,YY,GRID.qq,POL.spl_meth); 
end

% Residual vector from model file:
RES         = fun_res(POL);

end

%% Residual function for Direct Computation
function [RES] = fun_res_smol_dir(theta_col,GRID,POL,fun_res)

POL.theta   = reshape(theta_col,[],POL.dd);

% Residual vector from model file:
RES         = fun_res(POL);

end



%% Residual for Galerkin ('cheb_gal')
function [RES] = fun_res_gal(POL,fun_res,XX_poly_dw,theta_col)

% note: theta_col is column vector, 
% but POL.theta = pp*dd matrix;
POL.theta   = reshape(theta_col,[],POL.dd);

% Residual vector from model file:
res     = fun_res(POL);

mm = size(XX_poly_dw,1);

% Pre-allocate dimensions of RES:
RES     = NaN(POL.dd*POL.pp,1);
  

cnt_d = 0;%counter for policy residual vector
cnt_p = 0;%counter for policy vector
for id = 1:POL.dd%loop over policy variables 
	st_d  = cnt_d + 1;
	cnt_d = cnt_d + mm;
	for jj = 1:POL.pp
        cnt_p = cnt_p + 1;
        
        % only orthogonal to relevant residuals 
        % (not other Euler residuals) 
        R_int = res(st_d:cnt_d,1).* XX_poly_dw(1:mm,jj);          
        RES(cnt_p,1) = sum(R_int);
        clear R_int
	end
end


end

%% Residual for minimization of squared error using lsqnonlin (for method 'cheb_mse' and 'mono_mse')
function [RES] = fun_res_lsq(POL,fun_res,theta_col)

% note: theta_col is column vector, 
% but POL.theta = pp*dd matrix;
POL.theta   = reshape(theta_col,[],POL.dd);

% Residual vector from model file:
RES         = fun_res(POL);

end


%% Get pp_y cell array
function [pp_y] = get_pp_y(gridVecs,YY,qq,spl_meth)
% YY - mm x dd matrix of policy variables (each variable a column)

dd = size(YY,2);

pp_y  = cell(1,dd);
for id = 1:dd
    if size(gridVecs,2) == 1
        Yi_grid      = reshape(YY(:,id),qq,1);
    else
        Yi_grid      = reshape(YY(:,id),qq);
    end    
    % Fit a spline through all points:
	pp_y{1,id}  = griddedInterpolant(gridVecs,Yi_grid,spl_meth);
end
    
end

%% Initial value for theta
function [theta] = get_theta_ini_grid(GRID,POL,YY)
% YY: mm x dd policy variables

theta = NaN(POL.pp,POL.dd);
if strcmp(POL.algo,'mono_mse')
    for id = 1:POL.dd
        theta(:,id)         = GRID.XX_poly\YY(:,id);
    end
elseif strcmp(POL.algo,'cheb_tmi') || strcmp(POL.algo,'cheb_mse')...
        || strcmp(POL.algo,'cheb_gal') 
    for id = 1:POL.dd
        theta(:,id)         = GRID.XX_poly_dw\YY(:,id);
    end
elseif strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'smol_dir')
    for id = 1:POL.dd
        theta(:,id)         = GRID.inv_XX_poly_dw*YY(:,id);
    end
end

end


%% Get length of parameter vector (pp) of policy function 
function [pp] = get_length_para(algo,GRID)

%(for each pol. variable, so total number of parameters is dd*pp):
if strcmp(algo,'cheb_mse') || strcmp(algo,'cheb_gal') || ...
        strcmp(algo,'smol_tmi') || strcmp(algo,'smol_dir') ||  strcmp(algo,'cheb_tmi')
    pp  = size(GRID.XX_poly_dw,2);%number of parameters in theta(:,id)    
elseif strcmp(algo,'mono_mse')
    pp  = size(GRID.XX_poly,2);%number of parameters in theta(:,id)     
elseif strcmp(algo,'spl_tmi') || strcmp(algo,'spl_dir')
    pp  = GRID.mm;%number of gridpoints    
end
if strcmp(algo,'cheb_mse') || strcmp(algo,'cheb_gal') || ...
        strcmp(algo,'mono_mse') || strcmp(algo,'cheb_tmi') || ...
        strcmp(algo,'smol_tmi') || strcmp(algo,'smol_dir')
    if pp > GRID.mm
        error('System underidentified: increase number of gridpoints');
    end
end

end


%% Set default options for solver, and override only speficified fiels in options
function [opt_solver] = set_default_opt_solver(options,algo)

if isfield(options,'Algorithm')%If algorithm is speficied: create optimoptions with algorithm
  algo = options.Algorithm;
    if strcmp(algo,'spl_tmi') ||strcmp(algo,'spl_dir') || strcmp(algo,'smol_tmi') ||strcmp(algo,'smol_dir')  || strcmp(algo,'cheb_gal') || strcmp(algo,'cheb_tmi')
        opt_solver =  optimoptions(@fsolve,'Algorithm',algo,'Display','off'); 
    elseif strcmp(algo,'mono_mse') || strcmp(algo,'cheb_mse') 
        opt_solver =  optimoptions(@lsqnonlin,'Algorithm',algo,'Display','off');       
    end
else%Create optimoptions with preferred algorithms
    if strcmp(algo,'spl_tmi') || strcmp(algo,'smol_tmi') || strcmp(algo,'cheb_tmi')
        opt_solver =  optimoptions(@fsolve,'Algorithm','trust-region','Display','off'); 
    elseif strcmp(algo,'spl_dir') || strcmp(algo,'smol_dir') ||strcmp(algo,'cheb_gal')
        opt_solver =  optimoptions(@fsolve,'Display','off'); 
    elseif strcmp(algo,'mono_mse') || strcmp(algo,'cheb_mse') 
        opt_solver =  optimoptions(@lsqnonlin,'Display','off');       
    end 
end

% Override previously set optimoptions with values set in options
field_nms   = fields(options);
vals        = struct2cell(options);

for ir = 1:size(field_nms,1)
    if ~(strcmp(field_nms{ir,1},'Algorithm') || strcmp(field_nms{ir,1},'override_all'))        
        opt_solver  = optimoptions(opt_solver,field_nms{ir,1},vals{ir,1});        
    end    
end


end