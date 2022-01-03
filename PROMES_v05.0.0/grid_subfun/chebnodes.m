function [x] = chebnodes(dd)
% function [x] = chebnodes(dd)
%
% Purpose: creats dd Chebyshev nodes, in [-1,1]
%
% Input:
%       dd       = scalar of number of nodes
%
% Output:
%       x       = dd x 1 vector of nodes

% Code copied from Joris de Wind
% (this version: Sijmen Duineveld, s.a.duineveld@outlook.com)


ii	= (1:dd)';
x	= -cos( (pi*(ii-0.5))/dd );

end

