function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = stnd_rbc_dyn.static_resid_tt(T, y, x, params);
end
residual = zeros(5, 1);
residual(1) = params(2)*exp((-params(5))*y(1))*(1+params(1)*exp(y(3)+(params(1)-1)*(y(2)-y(5)))-params(3))-exp((-params(5))*y(1));
lhs = exp(y(1))+exp(y(2));
rhs = exp(y(4))+exp(y(2))*(1-params(3));
residual(2) = lhs - rhs;
lhs = y(3);
rhs = y(3)*params(7)+x(1)*params(8);
residual(3) = lhs - rhs;
lhs = y(4);
rhs = y(3)+params(1)*y(2)+y(5)*(1-params(1));
residual(4) = lhs - rhs;
lhs = y(5);
rhs = T(1)*(params(1)*y(2)+y(3)+log(1-params(1))-params(5)*y(1)-log(params(6)));
residual(5) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
