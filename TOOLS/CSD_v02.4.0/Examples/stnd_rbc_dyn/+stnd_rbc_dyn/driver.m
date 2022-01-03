%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'stnd_rbc_dyn';
M_.dynare_version = '4.6.4';
oo_.dynare_version = '4.6.4';
options_.dynare_version = '4.6.4';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('stnd_rbc_dyn.log');
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'epsilon'};
M_.exo_names_tex(1) = {'epsilon'};
M_.exo_names_long(1) = {'epsilon'};
M_.endo_names = cell(5,1);
M_.endo_names_tex = cell(5,1);
M_.endo_names_long = cell(5,1);
M_.endo_names(1) = {'LC'};
M_.endo_names_tex(1) = {'LC'};
M_.endo_names_long(1) = {'LC'};
M_.endo_names(2) = {'LK'};
M_.endo_names_tex(2) = {'LK'};
M_.endo_names_long(2) = {'LK'};
M_.endo_names(3) = {'LZ'};
M_.endo_names_tex(3) = {'LZ'};
M_.endo_names_long(3) = {'LZ'};
M_.endo_names(4) = {'LY'};
M_.endo_names_tex(4) = {'LY'};
M_.endo_names_long(4) = {'LY'};
M_.endo_names(5) = {'LH'};
M_.endo_names_tex(5) = {'LH'};
M_.endo_names_long(5) = {'LH'};
M_.endo_partitions = struct();
M_.param_names = cell(13,1);
M_.param_names_tex = cell(13,1);
M_.param_names_long = cell(13,1);
M_.param_names(1) = {'alpha'};
M_.param_names_tex(1) = {'alpha'};
M_.param_names_long(1) = {'alpha'};
M_.param_names(2) = {'beta'};
M_.param_names_tex(2) = {'beta'};
M_.param_names_long(2) = {'beta'};
M_.param_names(3) = {'delta'};
M_.param_names_tex(3) = {'delta'};
M_.param_names_long(3) = {'delta'};
M_.param_names(4) = {'eta'};
M_.param_names_tex(4) = {'eta'};
M_.param_names_long(4) = {'eta'};
M_.param_names(5) = {'nu'};
M_.param_names_tex(5) = {'nu'};
M_.param_names_long(5) = {'nu'};
M_.param_names(6) = {'chi'};
M_.param_names_tex(6) = {'chi'};
M_.param_names_long(6) = {'chi'};
M_.param_names(7) = {'rho_z'};
M_.param_names_tex(7) = {'rho\_z'};
M_.param_names_long(7) = {'rho_z'};
M_.param_names(8) = {'sigma_z'};
M_.param_names_tex(8) = {'sigma\_z'};
M_.param_names_long(8) = {'sigma_z'};
M_.param_names(9) = {'Css'};
M_.param_names_tex(9) = {'Css'};
M_.param_names_long(9) = {'Css'};
M_.param_names(10) = {'Kss'};
M_.param_names_tex(10) = {'Kss'};
M_.param_names_long(10) = {'Kss'};
M_.param_names(11) = {'Zss'};
M_.param_names_tex(11) = {'Zss'};
M_.param_names_long(11) = {'Zss'};
M_.param_names(12) = {'Yss'};
M_.param_names_tex(12) = {'Yss'};
M_.param_names_long(12) = {'Yss'};
M_.param_names(13) = {'Hss'};
M_.param_names_tex(13) = {'Hss'};
M_.param_names_long(13) = {'Hss'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 5;
M_.param_nbr = 13;
M_.orig_endo_nbr = 5;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.linear_decomposition = false;
M_.nonzero_hessian_eqs = [1 2];
M_.hessian_eq_zero = isempty(M_.nonzero_hessian_eqs);
M_.orig_eq_nbr = 5;
M_.eq_nbr = 5;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 0 3 8;
 1 4 0;
 2 5 9;
 0 6 0;
 0 7 10;]';
M_.nstatic = 1;
M_.nfwrd   = 2;
M_.npred   = 1;
M_.nboth   = 1;
M_.nsfwrd   = 3;
M_.nspred   = 2;
M_.ndynamic   = 4;
M_.dynamic_tmp_nbr = [2; 1; 0; 0; ];
M_.model_local_variables_dynamic_tt_idxs = {
};
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'LZ' ;
  4 , 'name' , 'LY' ;
  5 , 'name' , 'LH' ;
};
M_.mapping.LC.eqidx = [1 2 5 ];
M_.mapping.LK.eqidx = [1 2 4 5 ];
M_.mapping.LZ.eqidx = [1 3 4 5 ];
M_.mapping.LY.eqidx = [2 4 ];
M_.mapping.LH.eqidx = [1 4 5 ];
M_.mapping.epsilon.eqidx = [3 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.state_var = [2 3 ];
M_.exo_names_orig_ord = [1:1];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(5, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(13, 1);
M_.endo_trends = struct('deflator', cell(5, 1), 'log_deflator', cell(5, 1), 'growth_factor', cell(5, 1), 'log_growth_factor', cell(5, 1));
M_.NNZDerivatives = [20; 21; 25; ];
M_.static_tmp_nbr = [1; 0; 0; 0; ];
M_.model_local_variables_static_tt_idxs = {
};
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
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (1)^2;
steady;
options_.k_order_solver = true;
options_.irf = 0;
options_.nocorr = true;
options_.nomoments = true;
options_.noprint = false;
options_.order = 3;
var_list_ = {};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);
save('stnd_rbc_dyn_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('stnd_rbc_dyn_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('stnd_rbc_dyn_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('stnd_rbc_dyn_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('stnd_rbc_dyn_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('stnd_rbc_dyn_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('stnd_rbc_dyn_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
