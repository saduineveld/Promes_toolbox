function [LK_n,LMPK,LH] = stnd_rbc_aux(par,LK,LZ,LC)
%function [LK_n,LMPK,LH] = stnd_rbc_aux(par,LK,LZ,LC)
%
% Get log(K_t+1), log(MPK_t) and log(H_t) for standard RBC model
% where MPK is marginal prod. of capital

% Sijmen Duineveld, updated April 2021, s.a.duineveld@outlook.com

LH = par.eta/(1+par.alpha*par.eta) * ...
    ( -log(par.chi) -par.nu*LC + ...
    log(1-par.alpha) + LZ + par.alpha*LK );

LMPK = log(par.alpha) + LZ + (par.alpha-1)*(LK-LH);

LK_n = log( exp(LZ + par.alpha*LK + (1-par.alpha)*LH)...
    - exp(LC) + (1-par.delta)*exp(LK) );

end