function [SIM] = stnd_rbc_sim(par,SS,POL,GRID,opt_sim)
% Simulation for standard rbc model
%
% INPUTS:
%		par         = structure with parameters of model
%
%       SS          = structure with steady state levels
%
%       POL         = structure with policy function
%
%       GRID        = structure with all necessary GRID properties
%
%       opt_sim     = structure with settings for simulation. Contains:
%                       TT      - number of periods in simulation (scalar)
%                       rws     - number of rows in simulation (scalar)
%                       T_ini   - initial periods at s.s.
%
% OUTPUT:
%       SIM         = structure with series in logs: lK (capital), 
%                       lC (consumption), lh (hours worked), lZ (TFP)
%
% USES:
%       get_pol_var = function in toolbox PROMES
%
%       stnd_rbc_aux = calculates auxiliary variables for stnd. rbc model     

% Sijmen Duineveld, updated April 2021, s.a.duineveld@outlook.com


% Initialize variables (LK, LC, LH, LZ)
nms                 = {'K','C','H','Z'}; 
for in = 1:size(nms,2)
   eval(['L',nms{in},' = NaN(opt_sim.rws,opt_sim.T_ini+opt_sim.TT);']); 
   eval(['L',nms{in},'(:,1:opt_sim.T_ini) = log(SS.',nms{in},'ss);']);  
end   
% pre-allocate one extra column for LK_t+1
LK = [LK,NaN(opt_sim.rws,1)];
% Set k_t+1 = k_ss:
LK(:,opt_sim.T_ini+1)  = log(SS.Kss);


% initialize shocks:
rng('default');
rng(1);
epsilon = [zeros(opt_sim.rws,opt_sim.T_ini),randn(opt_sim.rws,opt_sim.TT)];

% loop over time (from T_ini +1 till T_ini + TT):
for it = opt_sim.T_ini + 1: opt_sim.T_ini + opt_sim.TT
    % Calculate TFP by adding a stochastic shock:
    LZ(:,it) = par.rho_z * LZ(:,it-1) + par.sigma_z * epsilon(:,it);
    
    % Calculate the policy variable:
    %LC(:,it) = get_pol_var(POL,[LK(:,it),LZ(:,it)],GRID);
    LC(:,it) = get_pol_var(POL,[LK(:,it),LZ(:,it)],GRID);
    
    % Calculate the other variables of interest:
    [LK(:,it+1),~,LH(:,it)] = ...
        stnd_rbc_aux(par,LK(:,it),LZ(:,it),LC(:,it));    
end

% remove last k_t+1 (ie. adjust size to rws x (T_ini+TT) )
LK = LK(:,1:end-1);

% Add output:
SIM.LY = LZ + par.alpha*LK + (1-par.alpha)*LH;

% Put variables in structure SIM:
for in = 1:size(nms,2)
   eval(['SIM.L',nms{in},' = L',nms{in},';']); 
end  

SIM.epsilon = epsilon;

end

