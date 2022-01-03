%% SIMULATION given perturbation solution (Dynare or CSD toolbox)
function [SIM] = stnd_rbc_sim_pert(par,SS,DRSOL,opt_sim,order)
%[SIM] = stnd_rbc_sim_pert(par,SS,DRSOL,opt_sim,order)
%
%       opt_sim     = structure with settings for simulation. Contains:
%                       sol_typ - 'dyn', or 'csd'
%                       TT      - number of periods in simulation (scalar)
%                       cls     - number of columns in simulation (scalar)
%                       T_ini   - initial periods at s.s.
%
%       order       = order of solution to be used
%
% OUTPUT:
%       SIM         = structure with series in logs: LK (capital), 
%                       LC (consumption), LH (hours worked), LZ (TFP), LY
%                       (output)

% Sijmen Duineveld, September 2021, s.a.duineveld@outlook.com

cols = opt_sim.cols;


% Initialize variables (LK, LC, LH, LZ,LH)
nms                 = {'K','C','H','Z','Y'}; 
for in = 1:size(nms,2)
   eval(['L',nms{in},' = NaN(opt_sim.T_ini+opt_sim.TT,cols);']); 
   eval(['L',nms{in},'(1:opt_sim.T_ini,:) = log(SS.',nms{in},'ss);']);  
end   
% pre-allocate one extra row for LK_t+1 and LZ_t+1
LK = [LK;NaN(1,cols)];
LZ = [LZ;NaN(1,cols)];

% initialize shocks:
rng('default');
rng(1);
epsilon = [zeros(cols,opt_sim.T_ini),randn(cols,opt_sim.TT+1)]';

% Set k_t+1 = k_ss 
LK(opt_sim.T_ini+1,:)  = log(SS.Kss);

% Set z_t+1 = lz(t)+eps(t+1)
if strcmp(opt_sim.sol_typ,'csd')
    LZ(opt_sim.T_ini+1,:) = get_ZZn(DRSOL,zeros(1,cols),epsilon(opt_sim.T_ini+1,:));
end

% loop over time (from T_ini +1 till T_ini + TT):
for it = opt_sim.T_ini + 1: opt_sim.T_ini + opt_sim.TT   
    
    if strcmp(opt_sim.sol_typ,'csd')

        % Get XXn, YY, ZZn using eval_sol_csd:
        [LK(it+1,:),YY_t,LZ(it+1,:)] = ...
            eval_sol_csd(DRSOL,LK(it,:),LZ(it,:),epsilon(it+1,:),order);
        LC(it,:) = YY_t(1,:);%Cons. is first control variable
        LH(it,:) = YY_t(2,:);%Hours is second control variable  
        LY(it,:) = YY_t(3,:);%Output is third control variable
    elseif strcmp(opt_sim.sol_typ,'dyn')
        Xp = [LK(it,:)',LZ(it-1,:)'];
        EPSt = epsilon(it,:)';
        
        [Xt,Yt] = eval_dyn_pol_beta(DRSOL,Xp,EPSt,order);
          
        % Transpose results (each simulation is column here, but a row in
        % output of eval_dyn_pol
        
        LK(it+1,:)  = Xt(:,1)';%K is 1st state variable in DR order
        LZ(it,:)    = Xt(:,2)';%Z is 2nd state variable in DR order
        
        LY(it,:) = Yt(:,1)';%Output is 1st control variable in DR order
        LC(it,:) = Yt(:,2)';%Cons. is 2nd control variable in DR order
        LH(it,:) = Yt(:,3)';%Hours is 3rd control variable in DR order       
        
 	end       
        
end

% remove last k_t+1 and z_t+1 (ie. adjust size to (T_ini+TT) x cols )
LK = LK(1:end-1,:);
LZ = LZ(1:end-1,:);


% Put variables in structure SIM:
for in = 1:size(nms,2)
   eval(['SIM.L',nms{in},' = L',nms{in},';']); 
end  

SIM.epsilon = epsilon;

end