function g3 = dynamic_g3(T, y, x, params, steady_state, it_, T_flag)
% function g3 = dynamic_g3(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g3
%

if T_flag
    T = stnd_rbc_dyn.dynamic_g3_tt(T, y, x, params, steady_state, it_);
end
v3 = zeros(25,3);
v3(1,1)=1;
v3(2,1)=1;
v3(3,1)=1;
v3(4,1)=1;
v3(5,1)=1;
v3(6,1)=1;
v3(7,1)=1;
v3(8,1)=1;
v3(9,1)=1;
v3(10,1)=1;
v3(11,1)=1;
v3(12,1)=1;
v3(13,1)=1;
v3(14,1)=1;
v3(15,1)=1;
v3(16,1)=1;
v3(17,1)=1;
v3(18,1)=1;
v3(19,1)=1;
v3(20,1)=1;
v3(21,1)=1;
v3(22,1)=2;
v3(23,1)=2;
v3(24,1)=2;
v3(25,1)=2;
v3(1,2)=267;
v3(2,2)=932;
v3(3,2)=928;
v3(4,2)=933;
v3(5,2)=934;
v3(6,2)=884;
v3(7,2)=889;
v3(8,2)=890;
v3(9,2)=944;
v3(10,2)=945;
v3(11,2)=956;
v3(12,2)=400;
v3(13,2)=405;
v3(14,2)=406;
v3(15,2)=460;
v3(16,2)=461;
v3(17,2)=472;
v3(18,2)=1065;
v3(19,2)=1066;
v3(20,2)=1077;
v3(21,2)=1198;
v3(22,2)=267;
v3(23,2)=1;
v3(24,2)=400;
v3(25,2)=666;
v3(1,3)=(-((-params(5))*(-params(5))*(-params(5))*exp((-params(5))*y(3))));
v3(2,3)=(1+params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))-params(3))*params(2)*(-params(5))*(-params(5))*(-params(5))*exp((-params(5))*y(8));
v3(3,3)=params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
v3(4,3)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
v3(5,3)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
v3(6,3)=T(3)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v3(7,3)=T(3)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v3(8,3)=T(3)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(9,3)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*T(3);
v3(10,3)=T(3)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(11,3)=T(3)*params(1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(12,3)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v3(13,3)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v3(14,3)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(15,3)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v3(16,3)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(17,3)=T(1)*params(1)*(params(1)-1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(18,3)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v3(19,3)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(20,3)=T(1)*params(1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(21,3)=T(1)*params(1)*(-(params(1)-1))*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v3(22,3)=exp(y(3));
v3(23,3)=(-((1-params(3))*exp(y(1))));
v3(24,3)=exp(y(4));
v3(25,3)=(-exp(y(6)));
g3 = sparse(v3(:,1),v3(:,2),v3(:,3),5,1331);
end
