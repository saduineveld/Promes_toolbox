function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = stnd_rbc_dyn.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
v2 = zeros(21,3);
v2(1,1)=1;
v2(2,1)=1;
v2(3,1)=1;
v2(4,1)=1;
v2(5,1)=1;
v2(6,1)=1;
v2(7,1)=1;
v2(8,1)=1;
v2(9,1)=1;
v2(10,1)=1;
v2(11,1)=1;
v2(12,1)=1;
v2(13,1)=1;
v2(14,1)=1;
v2(15,1)=1;
v2(16,1)=1;
v2(17,1)=1;
v2(18,1)=2;
v2(19,1)=2;
v2(20,1)=2;
v2(21,1)=2;
v2(1,2)=25;
v2(2,2)=85;
v2(3,2)=81;
v2(4,2)=41;
v2(5,2)=86;
v2(6,2)=96;
v2(7,2)=87;
v2(8,2)=107;
v2(9,2)=37;
v2(10,2)=42;
v2(11,2)=92;
v2(12,2)=43;
v2(13,2)=103;
v2(14,2)=97;
v2(15,2)=98;
v2(16,2)=108;
v2(17,2)=109;
v2(18,2)=25;
v2(19,2)=1;
v2(20,2)=37;
v2(21,2)=61;
v2(1,3)=(-((-params(5))*(-params(5))*exp((-params(5))*y(3))));
v2(2,3)=(1+params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))-params(3))*params(2)*(-params(5))*(-params(5))*exp((-params(5))*y(8));
v2(3,3)=T(3)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v2(4,3)=v2(3,3);
v2(5,3)=params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*T(3);
v2(6,3)=v2(5,3);
v2(7,3)=T(3)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v2(8,3)=v2(7,3);
v2(9,3)=T(1)*params(1)*(params(1)-1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v2(10,3)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v2(11,3)=v2(10,3);
v2(12,3)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v2(13,3)=v2(12,3);
v2(14,3)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
v2(15,3)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v2(16,3)=v2(15,3);
v2(17,3)=T(1)*params(1)*(-(params(1)-1))*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
v2(18,3)=exp(y(3));
v2(19,3)=(-((1-params(3))*exp(y(1))));
v2(20,3)=exp(y(4));
v2(21,3)=(-exp(y(6)));
g2 = sparse(v2(:,1),v2(:,2),v2(:,3),5,121);
end
