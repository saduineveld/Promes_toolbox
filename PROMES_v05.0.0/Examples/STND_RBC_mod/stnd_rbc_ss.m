function [SS] = stnd_rbc_ss(par,Zss)
%[SS] = stnd_rbc_ss(par,Zss)
%
% Calculation of steady state of Standard RBC model 
%
% INPUTS:
%		par     = parameters of model
%
%       Zss     = steady state level of TFP

% Sijmen Duineveld, updated March 2021, s.a.duineveld@outlook.com

Omega = (1-par.beta*(1-par.delta))/ (par.alpha*par.beta*Zss);

Kss = ( ((1-par.alpha)/par.chi*Zss.*(Zss.*Omega-par.delta).^-par.nu ).^par.eta .* Omega.^((1+par.eta*par.alpha)/(par.alpha-1)) ).^(1/(1+par.eta*par.nu));

Hss = Omega.^(1/(1-par.alpha)) .* Kss;

Css = (Zss.*Omega - par.delta) .* Kss;

Yss = Zss *Kss^par.alpha*Hss^(1-par.alpha);

nms                 = {'K','C','H','Z','Y'}; 
for in = 1:size(nms,2)
   eval(['SS.',nms{in},'ss = ',nms{in},'ss;']);  
end

end