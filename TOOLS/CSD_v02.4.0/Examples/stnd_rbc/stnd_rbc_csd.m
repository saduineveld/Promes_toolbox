function [MOD] = stnd_rbc_csd()
% Model file of standard RBC model for use with CSD_PERT_SOLVER

% Sijmen Duineveld, updated April 2021, s.a.duineveld@outlook.com

syms aa bb dd ee nu ch;

syms l_cc_t l_cc_n l_kk_t l_kk_n  l_zz_t l_zz_n l_hh_t l_hh_n l_yy_t l_yy_n;


%Euler equation:
f1 = bb *exp(-nu*l_cc_n)*( aa*exp(l_zz_n+(aa-1)*(l_kk_n-l_hh_n)) +(1-dd) ) - exp(-nu*l_cc_t);

% Law of Motion capital:
f2 = exp(l_yy_t) + (1-dd)*exp(l_kk_t) - exp(l_cc_t) - exp(l_kk_n) ;

% Static 1: labor supply
f3  =  (ee/(1+aa*ee))*( log(1-aa) + l_zz_t + aa*l_kk_t - nu*l_cc_t - log(ch) ) - l_hh_t;

% Static 2: Output
f4 = l_zz_t + aa*l_kk_t + (1-aa)*l_hh_t - l_yy_t;


MOD.FS = [f1;f2;f3;f4];

MOD.XX = [l_kk_t];
MOD.ZZ = [l_zz_t];
MOD.YY = [l_cc_t,l_hh_t,l_yy_t];

MOD.XXn = [l_kk_n];
MOD.ZZn = [l_zz_n];
MOD.YYn = [l_cc_n,l_hh_n,l_yy_n];

MOD.var_bs_nms  = {'l_kk','l_zz','l_cc','l_hh','l_yy'};
MOD.par_nms     = {'aa','bb','dd','ee','nu','ch'};

end