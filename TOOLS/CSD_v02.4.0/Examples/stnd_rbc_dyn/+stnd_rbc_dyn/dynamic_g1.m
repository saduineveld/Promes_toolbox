function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = stnd_rbc_dyn.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(5, 11);
g1(1,3)=(-((-params(5))*exp((-params(5))*y(3))));
g1(1,8)=(1+params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))-params(3))*T(3);
g1(1,4)=T(1)*params(1)*(params(1)-1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g1(1,9)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)));
g1(1,10)=T(1)*params(1)*exp(y(9)+(params(1)-1)*(y(4)-y(10)))*(-(params(1)-1));
g1(2,3)=exp(y(3));
g1(2,1)=(-((1-params(3))*exp(y(1))));
g1(2,4)=exp(y(4));
g1(2,6)=(-exp(y(6)));
g1(3,2)=(-params(7));
g1(3,5)=1;
g1(3,11)=(-params(8));
g1(4,1)=(-params(1));
g1(4,5)=(-1);
g1(4,6)=1;
g1(4,7)=(-(1-params(1)));
g1(5,3)=(-((-params(5))*T(2)));
g1(5,1)=(-(params(1)*T(2)));
g1(5,5)=(-T(2));
g1(5,7)=1;

end
