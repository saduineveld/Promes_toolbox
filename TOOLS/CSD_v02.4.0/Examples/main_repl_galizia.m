function main_repl_galizia()
% Replicates Galizia's Saddle Cycle model (BGPe)

% Sijmen Duineveld, June 2021, s.a.duineveld@outlook.com

%% STEP 0:  Matlab settings
close all; %close all figures
clc; %clear command prompt
dbstop if error;%acces workspace if error

restoredefaultpath;
clear RESTOREDEFAULTPATH_EXECUTED;

%Add folder with CSD and subfun:
addpath ('..');
addpath ('..\subfun');


%% STEP 1: LOAD parameters
shock_type = 'mu';%'mu' shock or 'zz' shock

if strcmp(shock_type,'mu')
    load pars_mu_galizia;
else%if strcmp(shock_type,'zz') || strcmp(shock_type,'det')
    load pars_z_galizia;
end
par = p;%use 'par' instead  of 'p'
clear p;

% Constant in interest rate rule:
par.betTh = par.e_^(-par.phie)/(1+(1-par.e_)*par.phi0*par.Phi0);


%% STEP 2: Get steady state:
[SS] = get_ss_bgpe(par);


%% STEP 3: Get symbolic model:
MOD =  BGPe_cds_mod(shock_type);


%% STEP 4: Set numerical values:
MOD.par_val =  [par.del,par.al,par.e_,par.om,par.gam,par.psi,par.phie,par.phi0,...
    par.Phi0,par.Phi2,par.Phi3,par.betTh];

if strcmp(shock_type,'det')
    MOD.SS_vec = [log(SS.XX_ss),log(SS.YY_ss),SS.Letr_ss,log(SS.lambda_ss)];
else    
    MOD.SS_vec = [log(SS.XX_ss),log(SS.YY_ss),0,SS.Letr_ss,log(SS.lambda_ss)];
end


%% STEP 5: Solve model
order = 3;
if strcmp(shock_type,'det')
    SOL = pert_ana_csd_lim(MOD,[],order,[]);
else
    if strcmp(shock_type,'mu')
        %SOL = pert_ana_csd_lim_old(MOD,par.rhomu,order,par.sigmu);
        
        SOL = pert_ana_csd_lim(MOD,par.rhomu,order,par.sigmu);
    elseif strcmp(shock_type,'zz')
        SOL = pert_ana_csd_lim(MOD,par.rhoz,order,par.sigz);
    end
end


%% STEP 6: SIMULATE MODEL

%Stochastic simulation:
[sim_set_stoch] = get_sim_set(shock_type,'STOCH',SS);
sim_set_stoch.sol_type = 'partial';%'partial' uses pert. pol. fun for 'e', and calculates other variables using non-linear equations (not perturbation solution)
SIM_stoch = bgpe_sim_pert(par,SS,SOL,sim_set_stoch,order,shock_type);

% Non stochastic simulation
[sim_set_nss] = get_sim_set(shock_type,'NSS',SS);
sim_set_nss.sol_type = 'partial';%'partial' uses pert. pol. fun for 'e', and calculates other variables using non-linear equations (not perturbation solution)
SIM_nss = bgpe_sim_pert(par,SS,SOL,sim_set_nss,order,shock_type);

%% Plot simulations
if strcmp(shock_type,'zz')
    t_nm = '$z$ shock';
elseif strcmp(shock_type,'mu')
    t_nm = '$\mu$ shock';
elseif strcmp(shock_type,'det')
    t_nm = 'Deterministic';
end


f1 = figure;
%Non stochastic (limit cycle):
st = sim_set_nss.TT - 250;
plot(SIM_nss.XX(st:end,1),SIM_nss.YY(st:end,1),'-','LineWidth',2);
hold all
en1 = 120;
plot(SIM_nss.XX(1:en1,1),SIM_nss.YY(1:en1,1),':','LineWidth',1.5);
en2 = 60;
plot(SIM_nss.XX(1:en2,2),SIM_nss.YY(1:en2,2),'--','LineWidth',1.5);
plot(SS.XX_ss,SS.YY_ss,'*')
title(t_nm,'Interpreter','Latex');
xlabel('$X$','Interpreter','Latex');
ylabel('$Y$','Interpreter','Latex');

%print(f1,['limit_cycle_noshock_',shock_type], '-depsc')

%Non stochastic (employment)
f2 = figure;
plot(0:270-1,SIM_nss.ee(1:270,1),'LineWidth',1.5);
xlim([0 270-1]);
title(t_nm,'Interpreter','Latex');
xlabel('$t$','Interpreter','Latex');
ylabel('Employment','Interpreter','Latex');
%print(f2,['employment_noshock_',shock_type], '-depsc')

if ~strcmp(shock_type,'det')
f3 = figure;
%Stochastic:
plot(0:sim_set_stoch.TT-1,SIM_stoch.ee,'LineWidth',1.5);
xlim([0 sim_set_stoch.TT-1]);
title(t_nm,'Interpreter','Latex');
xlabel('$t$','Interpreter','Latex');
ylabel('Employment','Interpreter','Latex');
%print(f3,['employment_stoch_',shock_type], '-depsc')
end


end

%% MODEL FUNCTION
function [MOD] = BGPe_cds_mod(shock_type)
% shock_type = 'det': no shock
% shock_type = 'mu': mu shocks
% shock_type = 'zz': zz shocks

syms del al e_ om gam phie phi0 Phi0 Phi2 Phi3 betTh;

psi = sym('psi');

syms Llambda_t Llambda_n Letr_t Letr_n Lxx_t Lxx_n Lyy_t Lyy_n;
if strcmp(shock_type,'mu')
    syms Lmu_t Lmu_n;
    Lzz_t = 0;
elseif strcmp(shock_type,'zz')
  	syms Lzz_t Lzz_n; 
    Lmu_t = 0;
else
    Lmu_t = 0;
    Lzz_t = 0;    
end

ee_t = 1/(1+exp(-Letr_t));
ee_n = 1/(1+exp(-Letr_n));

edev_t = ee_t-e_;
Phit = Phi0*exp( 0.5*Phi2/Phi0*edev_t^2 + 1/6*Phi3/Phi0*edev_t^3 );
Qt = betTh*(1+(1-ee_t)*phi0*Phit);

frac_x = (1-del-gam)/(1-del);
frac_y = (1-del-psi)/(1-del);

% Model:
f1 = Qt*ee_n^phie*exp(Llambda_n) - exp(Lmu_t)*exp(Llambda_t);
f2 = (1-del)*exp(Lxx_t) + psi*exp(Lyy_n) - exp(Lxx_n);
f3 = exp(Lzz_t)*ee_t^al - exp(Lyy_n);
f4 = (exp(Lyy_n) + frac_x*exp(Lxx_t) - frac_y*gam*exp(Lyy_t))^-om - exp(Llambda_t);



%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2;f3;f4];

%Endogenous state variables:
MOD.XX = [Lxx_t,Lyy_t];
%Exogenous state variables:
if strcmp(shock_type,'mu')
    MOD.ZZ = [Lmu_t];
elseif  strcmp(shock_type,'zz')
    MOD.ZZ = [Lzz_t];
else
    MOD.ZZ = [];
end
%Control variables:
MOD.YY = [Letr_t,Llambda_t];


MOD.XXn = [Lxx_n, Lyy_n];
if strcmp(shock_type,'mu')
    MOD.ZZn = [Lmu_n];
elseif  strcmp(shock_type,'zz')
    MOD.ZZn = [Lzz_n];
else
    MOD.ZZn = [];
end
MOD.YYn = [Letr_n,Llambda_n];


% Variable names (strings):
if strcmp(shock_type,'mu')
    MOD.var_bs_nms  = {'Lxx','Lyy','Lmu','Letr','Llambda'};
elseif  strcmp(shock_type,'zz')
    MOD.var_bs_nms  = {'Lxx','Lyy','Lzz','Letr','Llambda'};
else
    MOD.var_bs_nms  = {'Lxx','Lyy','Letr','Llambda'};
end

% Parameter names (strings):
MOD.par_nms     = {'del','al','e_','om','gam','psi','phie','phi0',...
    'Phi0','Phi2','Phi3','betTh'};


end

%% SET PARAMETERS FOR SIMULATION
function [sim_set] = get_sim_set(shock_type,sim_type,SS)


%% STOCHASTIC SIM:
sim_set.sim_type = sim_type;%NSS for Non Stochastic Simulation or 'STOCH' for stoch. sim
if strcmp(sim_type,'NSS')
    sim_set.TT = 1000;
else
   sim_set.TT = 270; 
end

%Set initial endogenous state variables:
if strcmp(shock_type,'zz') || strcmp(shock_type,'det')
    if strcmp(sim_type,'NSS')
        sim_set.ini_state = [7.56, 7.48; 0.96, 0.94];        
    else
        sim_set.ini_state = 0.999*[SS.XX_ss;SS.YY_ss];
    end
elseif strcmp(shock_type,'mu')
    if strcmp(sim_type,'NSS')
        sim_set.ini_state = [7.5, 7.45;    0.96, 0.94];      
    else
        sim_set.ini_state = 0.999*[SS.XX_ss;SS.YY_ss];
    end    
end
sim_set.cols    = size(sim_set.ini_state,2);


end

%% SIMULATION BGPe model, given perturbation solution
function [SIM] = bgpe_sim_pert(par,SS,SOL,sim_set,order,shock_type)


cols = sim_set.cols;
nx = size(SOL.XX_ss,2);
nz = size(SOL.ZZ_ss,2);
ny = size(SOL.YY_ss,2);

% Initialize variables
nms                 = {'XX','YY','Letr','lambda','ee','Lexo'}; 
 
for in = 1:size(nms,2)
   eval([nms{in},' = NaN(sim_set.TT+1,cols);']); 
end  

% Set initial values:
XX(1,:) = sim_set.ini_state(1,:);
YY(1,:) = sim_set.ini_state(2,:);

% initialize shocks:
if strcmp(shock_type,'det')
    epsilon = [];
else
    if strcmp(sim_set.sim_type,'NSS')
        epsilon = [zeros(cols,sim_set.TT)]';
        Lexo(1,:) = 0;
    else
        if cols ~= 1
            error('Should only simulate one serie')
        end
        
        %rng('default');
        %rng(1);
        sd = 1546032;
        rng(sd)
        
        nt = 2;%Galizia simulates two different stochastic processes (but uses only one at a time)
        if strcmp(shock_type,'zz')
            th0simnrm = randn(nt,1)./sqrt(1-par.rhoz.^2); % initial exogenous state
            ini_exo = par.sigz*th0simnrm(1);
        elseif strcmp(shock_type,'mu')
             th0simnrm = randn(nt,1)./sqrt(1-par.rhomu.^2); % initial exogenous state
             ini_exo = par.sigmu*th0simnrm(2);
        end
      
        epsilon_both = randn(2,sim_set.TT);
    
        if strcmp(shock_type,'zz')
            epsilon = [epsilon_both(1,:)]';
        elseif strcmp(shock_type,'mu')
            epsilon = [epsilon_both(2,:)]';   
        end        
        Lexo(1,:) = get_ZZn(SOL,ini_exo,epsilon(1,:));
    end
end

epsilon = [epsilon;zeros(1,cols)];%this extra row will not be used anyway

% loop over time (from T_ini +1 till T_ini + TT):
for it = 1: sim_set.TT
    
    if strcmp(sim_set.sol_type,'full')

        if strcmp(shock_type,'det')
          	[XXn,JUMP_t] = eval_sol_csd(SOL,[log(XX(it,:));log(YY(it,:))],[],[],order);
        else
            [XXn,JUMP_t,Lexo(it+1,:)] = eval_sol_csd(SOL,[log(XX(it,:));log(YY(it,:))],Lexo(it,:),epsilon(it+1,:),order);
        end
        Letr(it,:) = JUMP_t(1,:);
        
        ee(it,:) = 1./(1+exp(-Letr(it,:)));
    
        XX(it+1,:) = exp(XXn(1,:));
        YY(it+1,:) = exp(XXn(2,:)); 
 
    else
        if strcmp(shock_type,'det')
          	[~,JUMP_t] = eval_sol_csd(SOL,[log(XX(it,:));log(YY(it,:))],[],[],order);
        else
            [~,JUMP_t,Lexo(it+1,:)] = eval_sol_csd(SOL,[log(XX(it,:));log(YY(it,:))],Lexo(it,:),epsilon(it+1,:),order);
        end
        Letr(it,:) = JUMP_t(1,:);
        
        ee(it,:) = 1./(1+exp(-Letr(it,:)));
        
        % Calculate next period state variables:
        if strcmp(shock_type,'zz')
            YY(it+1,:) = exp(Lexo(it,:)).*ee(it,:).^par.al;
        else
            YY(it+1,:) = ee(it,:).^par.al;
        end
        XX(it+1,:) = (1-par.del)*XX(it,:) + par.psi*YY(it+1,:);
    end
end

% remove TT+1 data:
for in = 1:size(nms,2)
   eval([nms{in},' = ',nms{in},'(1:sim_set.TT,:);']); 
end  

if strcmp(shock_type,'zz')
    SIM.Lzz = Lexo;
elseif strcmp(shock_type,'mu')
    SIM.Lmu = Lexo;
end

% Put variables in structure SIM:
for in = 1:size(nms,2)
   eval(['SIM.',nms{in},' = ',nms{in},';']); 
   if ~strcmp(nms{in},'Letr') && ~strcmp(nms{in},'Lexo')
        eval(['SIM.dvl.',nms{in},' = log(',nms{in},') - log( SS.',nms{in},'_ss);']);
   end
end  

SIM.epsilon = epsilon;

end

%% Get steady state, given e_
function [SS] = get_ss_bgpe(par)

SS.zz_ss = 1;
SS.mu_ss = 1;

SS.YY_ss = SS.zz_ss*par.e_^par.al;
SS.XX_ss = 1/par.del*par.psi*SS.YY_ss;

SS.lambda_ss = ((1-par.del-par.gam)/(1-par.del) * SS.XX_ss +  ...
                (1-par.gam*(1-par.del-par.psi)/(1-par.del))*SS.YY_ss)^-par.om;
            
SS.Letr_ss = -log(1./par.e_ - 1);

SS.ee_ss = par.e_;

end