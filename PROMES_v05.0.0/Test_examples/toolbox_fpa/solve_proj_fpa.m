function [POL,ex_fl] = solve_proj_fpa(GRID,POL,fun_y_new,Y0)
%[POL,ex_fl] = solve_proj_fpa(GRID,POL,fun_res,Y0,options)

% Check if Y0 has correct format:
if GRID.mm ~= size(Y0,1) 
    error('Y0 has to be mm x dd matrix');
end
POL.dd = size(Y0,2);%number of policy variables;

POL.Y0 = Y0; clear Y0;

% Get length of parameter vector (pp) of policy function 
POL.pp = get_length_para(POL.algo,GRID);

% For splines only: set default interpolation method:
if strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'spl_dir') || strcmp(POL.algo,'spl_fpa') 
    if ~isfield(POL,'spl_meth')
        POL.spl_meth = 'spline';
    end
end

% For polynomials only (but not 'smol_tmi' or 'cheb_tmi'):
% Set initial policy functon (theta0) 
if strcmp(POL.algo,'mono_mse') || strcmp(POL.algo,'cheb_gal') ||...
        strcmp(POL.algo,'cheb_tmi') || strcmp(POL.algo,'cheb_mse') || ...
        strcmp(POL.algo,'smol_dir')
    POL.theta0 = get_theta_ini_grid(GRID,POL,POL.Y0);
end


% Solve the policy functions for each type
%% Fixed Point Iteration
if strcmp(POL.algo,'spl_fpa') || strcmp(POL.algo,'smol_fpa') || strcmp(POL.algo,'cheb_fpa') 
    if ~isfield(POL,'diff_tol')
        POL.diff_tol = 1e-8;%default maximum difference between YY & Y_old
    end
    if ~isfield(POL,'max_iter')
        POL.max_iter = 500;%default max. number of iterations;
        
    end
    if ~isfield(POL,'mem_Y')
        POL.mem_Y = 0;%default memory for updating Y
    end
    
    max_diff    = inf;
    max_res     = inf;
    
    % Initial guess for policy function
    POL.YY_old           = POL.Y0;
    
    % Initialize 'old' policy function
    % (for next period choices)
    if strcmp(POL.algo,'spl_fpa')
        POL.pp_y_old    = get_pp_y(GRID.gridVecs,POL.YY_old,GRID.qq,POL.spl_meth);
    elseif strcmp(POL.algo,'smol_fpa') || strcmp(POL.algo,'cheb_fpa') 
        POL.theta_old = get_theta_ini_grid(GRID,POL,POL.YY_old);
        
        if strcmp(POL.algo,'cheb_fpa') 
            %For consistency: recompute YY_old using the polynomial:
            %otherwise YY_old is not spanned by polynomial        
            for id = 1:POL.dd
                POL.YY_old(:,id) = GRID.XX_poly_dw*POL.theta_old(:,id);
            end
        end
        
    end
         
    iter = 1;
    while max_diff > POL.diff_tol && iter < POL.max_iter
        
        %POL.YY_old'
              
        % Update policy, given old policy:   
        YY_new = fun_y_new(POL);
        
        max_diff  = max(abs(YY_new(:)-POL.YY_old(:)));
        
        POL.YY_old = POL.mem_Y*POL.YY_old + (1-POL.mem_Y)*YY_new;       
   
        
        %Update 'old' policy function 
        % (for next period choices) 
        if strcmp(POL.algo,'spl_fpa')
            POL.pp_y_old = get_pp_y(GRID.gridVecs,POL.YY_old,GRID.qq,POL.spl_meth);
        elseif strcmp(POL.algo,'smol_fpa') || strcmp(POL.algo,'cheb_fpa') 
            POL.theta_old = get_theta_ini_grid(GRID,POL,POL.YY_old);            
        end  
        
        if strcmp(POL.algo,'cheb_fpa') 
            %For consistency: recompute YY_old using the polynomial:
            %otherwise YY_old is not spanned by polynomial        
            for id = 1:POL.dd
                POL.YY_old(:,id) = GRID.XX_poly_dw*POL.theta_old(:,id);
            end
        end
                
        iter    = iter + 1;      
    end
    
    if strcmp(POL.algo,'spl_fpa')
        POL.pp_y = POL.pp_y_old;
  	elseif strcmp(POL.algo,'smol_fpa') || strcmp(POL.algo,'cheb_fpa') 
        POL.theta = POL.theta_old;            
 	end  
    
    if iter == POL.max_iter  
        
        ex_fl = 0;%no solution found
    else         
        ex_fl = 1;     
    end
    
end

if ex_fl < 1
    
    warning('No proper solution found');
end


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
        || strcmp(POL.algo,'cheb_gal') || strcmp(POL.algo,'cheb_fpa') 
    for id = 1:POL.dd
        theta(:,id)         = GRID.XX_poly_dw\YY(:,id);
    end
elseif strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'smol_dir') || strcmp(POL.algo,'smol_fpa') 
    for id = 1:POL.dd
        theta(:,id)         = GRID.inv_XX_poly_dw*YY(:,id);
    end
end

end


%% Get length of parameter vector (pp) of policy function 
function [pp] = get_length_para(algo,GRID)

%(for each pol. variable, so total number of parameters is dd*pp):
if strcmp(algo,'cheb_mse') || strcmp(algo,'cheb_gal') || ...
        strcmp(algo,'smol_tmi') || strcmp(algo,'smol_dir') ||  ...
        strcmp(algo,'cheb_tmi') || ...
        strcmp(algo,'cheb_fpa') || strcmp(algo,'smol_fpa')
    pp  = size(GRID.XX_poly_dw,2);%number of parameters in theta(:,id)    
elseif strcmp(algo,'mono_mse')
    pp  = size(GRID.XX_poly,2);%number of parameters in theta(:,id)     
elseif strcmp(algo,'spl_tmi') || strcmp(algo,'spl_dir') || strcmp(algo,'spl_fpa')
    pp  = GRID.mm;%number of gridpoints    
end
if strcmp(algo,'cheb_mse') || strcmp(algo,'cheb_gal') || ...
        strcmp(algo,'mono_mse') || strcmp(algo,'cheb_tmi') || ...
        strcmp(algo,'smol_tmi') || strcmp(algo,'smol_dir') || ...
        strcmp(algo,'cheb_fpa') || strcmp(algo,'smol_fpa') 
    if pp > GRID.mm
        error('System underidentified: increase number of gridpoints');
    end
end

end