% First order solution:
function [Hx_w,Hy_w,FF,PP,stab] = solab_adj_stab(AA,BB,nx,nz)
%
% Function: solab
%
% Written by Paul Klein
%
% Rewritten in November 2007 after a suggestion by Alexander Meyer-Gohde
%
% Format: [f,p] = solab(a,b,nw);
%
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations.
%
% Inputs: Two square matrices a and b and a natural number nw
%
% a and b are the coefficient matrices of the difference equation
%
% a*x(t+1) = b*x(t)
% 
% where x(t) is arranged so that the state variables come first, and
%
% nw is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)   = f*k(t) and
%
% k(t+1) = p*k(t).
%
% Calls: qz, ordqz


% Adjusted to include output 'stab' with stability property
% Adjusted Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

nw = nx + nz;

[s,t,q,z] = qz(AA,BB);                % upper triangular factorization of the matrix pencil b-za
[s,t,q,z] = ordqz(s,t,q,z,'udo');   % reordering of generalized eigenvalues with the block inside the unit circle in the upper left
z21 = z(nw+1:end,1:nw);
z11 = z(1:nw,1:nw);

if rank(z11)<nw;
	error('Invertibility condition violated')
end

z11i = z11\eye(nw);
s11 = s(1:nw,1:nw);
t11 = t(1:nw,1:nw);

if abs(t(nw,nw))>abs(s(nw,nw)) || abs(t(nw+1,nw+1))<abs(s(nw+1,nw+1))
    if abs(t(nw,nw))>abs(s(nw,nw))    
        stab = -1;%steady state locally explosive
    else
        stab = 0;%indeterminate       
    end
    FF = NaN(size(z21*z11i));
    dyn = s11\t11;
    PP = NaN(size(z11*dyn*z11i));

    Hx_w = PP(1:nx,:);
    Hy_w = FF;
    
   warning('Wrong number of stable eigenvalues.');
else
    stab = 1;
    
    dyn = s11\t11;

    FF = real(z21*z11i);
    PP = real(z11*dyn*z11i);

    Hx_w = PP(1:nx,:);
    Hy_w = FF;
end



end