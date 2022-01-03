% Standard RBC model, model file for Dynare

% Sijmen Duineveld, September 2021, s.a.duineveld@outlook.com

%Variables:
var LC, LK, LZ, LY, LH;    
 
% TFP shocks:
varexo epsilon;

% Parameters:
parameters alpha beta delta eta nu chi rho_z sigma_z
    Css Kss Zss Yss Hss;

load params;

set_param_value('alpha',    par.alpha);
set_param_value('beta',     par.beta);
set_param_value('delta',    par.delta);
set_param_value('eta',      par.eta);
set_param_value('nu',       par.nu);
set_param_value('chi',       par.chi);

set_param_value('rho_z',     par.rho_z);
set_param_value('sigma_z',   par.sigma_z);

set_param_value('Css',      SS.Css);
set_param_value('Kss',      SS.Kss);
set_param_value('Zss',      SS.Zss);
set_param_value('Yss',      SS.Yss);
set_param_value('Hss',      SS.Hss);


model;

% Dynamic eq. 1: Euler equation:
beta*exp(-nu*LC(1))*( alpha*exp(LZ(1) + (alpha-1)*(LK-LH(1))) +1-delta ) - exp(-nu*LC);

% Dynamic eq. 2: Law of Motion Capital:
exp(LC)+exp(LK)= exp(LY) + (1-delta)*exp(LK(-1));

% Dynamic eq.3: Law of Motion TFP:
LZ = rho_z*LZ(-1) + epsilon*sigma_z;

% Static eq. 2: Output
LY = LZ + alpha*LK(-1) + (1-alpha)*LH;

% Static eq. 1: Labor supply
LH = (eta/(1+alpha*eta))*( log(1-alpha) + LZ + alpha* LK(-1) - nu*LC - log(chi));

end;

shocks;
  var epsilon; stderr 1;
end;

steady_state_model;
LC = log(Css);
LK = log(Kss);
LZ = log(Zss);
LY = log(Yss);
LH = log(Hss);
end;

steady;

stoch_simul(order=3,nocorr,nomoments,irf=0,print);
