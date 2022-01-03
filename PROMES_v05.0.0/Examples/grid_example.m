% grid_example is a script that replicates the procedures in
% 'prepgrid' and prints (intermediate) output on screen

% Sijmen Duineveld, updated for Promes_v05.0.0 December 2021, s.a.duineveld@outlook.com

%% Matlab settings
close all;
clc;
clearvars;

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

% Add relevant folders of Promes toolbox:
addpath ('..');
addpath ('..\grid_subfun');
addpath ('..\smolyak_subfun');


%% Initializtion of grid parameters:
gin.nn = 2;%number of state variables
gin.lb = [1,10];%lower bounds for [x1,x2]
gin.ub = [3,25];%upper bounds for [x1,x2]
   
% Set solution method
% 'cheb_gal','cheb_tmi','cheb_mse',
% 'spl_tmi','spl_tmi',
% 'smol_tmi','smol_tmi',
% 'mono_mse';
algo = 'cheb_gal';
    
if strncmp(algo,'cheb',4) || strncmp(algo,'mono',4) 
     meth_spec.ord_vec  = 2*ones(1,gin.nn);%order in each dim.
     meth_spec.qq       = [3,4];%number of nodes in each dim.
     
elseif strncmp(algo,'spl',3) 
    meth_spec.qq       = [3,4];%number of nodes in each dim.
    
elseif strncmp(algo,'smol',4) 
    meth_spec.mu_vec = [2,2];%accuracy in each dim.
end
    
% Construct structure with grid:
[GRID] = prepgrid(gin.nn,gin.lb,gin.ub,algo,meth_spec);
clear gin;


%% PRINT GRID VARIABLES:
if strncmp(algo,'spl',3) || strncmp(algo,'mono',4)
	fprintf('gridVecs: cell array with following gridpoints in each dimension:\n');
    fprintf('\n')
    for in = 1:GRID.nn
        disp(GRID.gridVecs{1,in})        
    end
    fprintf('\n')
    
	fprintf('The full grid of all possible combinations (xx) is:\n')
	disp(GRID.xx)
        
    if strcmp(algo,'mono_mse')        
        fprintf('\n')
        fprintf('The complete polynomial (XX) of order [%d,%d] is:\n',meth_spec.ord_vec(1),meth_spec.ord_vec(2))
        disp(GRID.XX_poly)    
    end
	        
elseif strncmp(algo,'cheb',4)
 	fprintf('gridVecs_dw is cell array with Chebyshev nodes (in [-1,1]) \n')
    fprintf('with following gridpoints in each dimension:\n');
    fprintf('\n')
  	for in = 1:GRID.nn
        disp(GRID.gridVecs_dw{1,in})        
    end
    fprintf('\n')
        
	fprintf('The full grid of all possible combinations of Chebyshev nodes in [-1,1] (xx_dw) is:\n')
	disp(GRID.xx_dw)
        
  	fprintf('The full grid of Chebysheb nodes when scaled up to [lb,ub] (xx) is:\n')
  	disp(GRID.xx)
        
   	fprintf('The complete Chebyshev polynomial (XX_poly_dw) of order [%d,%d] is:\n',meth_spec.ord_vec(1),meth_spec.ord_vec(2))
  	disp(GRID.XX_poly_dw)     
elseif strncmp(algo,'smol',4)
    fprintf('The Smolyak grid in [-1,1] (xx_dw) consists of the points:\n')
	disp(GRID.xx_dw)
    
    fprintf('The Smolyak grid scaled up to [lb,ub] (xx) is:\n')
  	disp(GRID.xx)
    
    fprintf('The Smolyak polynomial (XX_poly_dw) with accuracy mu=[%d,%d] is:\n',meth_spec.mu_vec(1),meth_spec.mu_vec(2))
  	disp(GRID.XX_poly_dw)     
    
else
    error('Invalid algo');
end

figure;
plot(GRID.xx(:,1),GRID.xx(:,2),'x','LineWidth',2);
xlabel('x_1');
ylabel('x_2');
title(['Initial grid (',GRID.grid_type,')']);
xlim([GRID.lb(1) GRID.ub(1)]);
ylim([GRID.lb(2) GRID.ub(2)]);

        