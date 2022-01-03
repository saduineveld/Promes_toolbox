function [MOD] = SGU_mod
%Schmidt-Grohe-Uribe (2004, JEDC) model
%
% Note: Law of Motion for exogenous states 
% is defined outside model file:
% ZZ_t = Rho*ZZ_t-1 + Omega*epsilon

% Sijmen Duineveld, updated May 2021, s.a.duineveld@outlook.com

%% BLOCK 1: Define parameters and variables
%Parameters:
syms aa bb dd nu rho_a sigma_a;

%Variables:
syms LC_t LC_n LK_t LK_n LA_t LA_n;

%% BLOCK 2: MODEL EQUATIONS (excl. stochastic process)
%Euler equation:
f1 = bb *exp(-nu*LC_n)*( aa*exp((aa-1)*LK_n + LA_n) + (1-dd) ) - exp(-nu*LC_t);

%Budget constraint:
f2 = exp(aa*LK_t + LA_t) + (1-dd)*exp(LK_t) - exp(LK_n) -exp(LC_t);

%% BLOCK 3: ASSIGNMENTS
%Model equations:
MOD.FS = [f1;f2];

%Endogenous state variables:
MOD.XX = [LK_t];
%Exogenous state variables:
MOD.ZZ = [LA_t];
%Control variables:
MOD.YY = LC_t;

MOD.XXn = [LK_n];
MOD.ZZn = [LA_n];
MOD.YYn = LC_n;

% Variable names (strings):
MOD.var_bs_nms  = {'LK','LA','LC'};
% Parameter names (strings):
MOD.par_nms     = {'aa','bb','dd','nu'};

end