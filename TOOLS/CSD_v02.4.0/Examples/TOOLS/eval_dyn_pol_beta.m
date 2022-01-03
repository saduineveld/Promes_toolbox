function [Xt,Yt] = eval_dyn_pol_beta(DR,Xp,EPSt,order)
% [Xt,Yt] = eval_dyn_pol_beta(DR,Xp,EPSt,order)
%
% Evaluate dynare policy function, upto order 3
%
% Definitions:
% rws: number of simulations (rows in Xp and EPSt)
%
% Nx: number of state variables
% Ny: number of control variables
% Ne: number of shocks
%
% DCL order: declartion order
% DR order: Decision Rule order
% 
% INPUT:
%       DR - dynare structure (use oo_.dr) with required fields: ys,
%            order_var, state_var, and the decision rules (see below)
%               order_var: vector mapping DCL to DR order, dim: (Nx+Ny) x1 
%               state_var: state variables, referring to DCL order, dim: Nx x 1
%               ys: vector of steady states, dim: (Nx + Ny) x 1, in
%                   DCL order
%               decision rules:
%                   order 1: ghx, ghu
%                   order 2: order 1 + ghxx, ghxu, ghuu, ghs2
%                   order 3:  g_0, g_1, g_2, g_3
%                     
%       Xp - State variables (period t-1), dim: rws x Nx. Columns in DR
%       order (as in DR.state_var)
%
%       EPSt - Stochastic shocks (period t), dim: rws x Ne. Columns in
%           declaration order.
%
%       order - order of the evaluation
%
% OUTPUT:
%       Xt - Matrix of state variables, dim. rws x Nx. Order as in
%           state_var
%       Yt - Matrix of control variables, dim. rws x Ny. DR order after
%               removing state variables

% Sijmen Duineveld, September 2021, s.a.duineveld@outlook.com

if nargin < 4 || isempty(order)
    %Set default order to 1 if not speficied
    order = 1;
else
    if ~(order == 1 || order == 2 || order == 3)
        error('Order should be 1, 2 or 3')
    end
end


rws = size(Xp,1);
if exist('EPSt','var') && ~isempty(EPSt)
    %Check size Xp and EPSt:
    if rws ~= size(EPSt,1)
        error('Xp and EPSt do not have same number of rows');
    end
end

Nx = size(DR.ghx,2);%Number of state variables
Ny = size(DR.ys,1) - Nx;%Number of control variables
%Ne = size(DR.ghu,2);%Number of shocks

% Change order of ys to DR order:
ys_DR = DR.ys(DR.order_var);

% Vector of steady states of state variables, in DR order:
xs_DR = DR.ys(DR.state_var);
% State variables in deviation from steady state:
Xp_dev = Xp-xs_DR';

%Pre-allocate size of XYt with zeros
XYt = zeros(rws,Nx+Ny);

if order == 1 || order == 2
 	% Add 1st order terms:
    XYt = XYt + ys_DR' + Xp_dev*DR.ghx' + EPSt*DR.ghu';

    if order > 1
        % Construct order 2 polynomials:
        [Xp2_dev]   = kron_poly2(Xp_dev,Xp_dev);    
        EPSt2       = kron_poly2(EPSt,EPSt);
        XE          = kron_poly2(Xp_dev,EPSt);

        % Add 2nd order terms:
        XYt         = XYt + 0.5*DR.ghs2' + 0.5*Xp2_dev*DR.ghxx' + 0.5*EPSt2*DR.ghuu' + XE*DR.ghxu';
    end
else%3rd order:
    % Construct order 2 and 3 folded polynomials:
    ZZ2_fold   = kron_poly2_fold([Xp_dev,EPSt]);
    ZZ3_fold   = kron_poly3_fold([Xp_dev,EPSt]);  
    
    XYt = XYt + ys_DR' + DR.g_0' + [Xp_dev,EPSt]*DR.g_1'+ ...
                ZZ2_fold*DR.g_2' + ZZ3_fold*DR.g_3';    
end
    

% SPLIT XYt in controls Yt and state variables Xt
% State variables:
Xt = XYt(:,DR.state_var);
if nargout > 1
    XYt(:,DR.state_var) = []; %remove state variables from XYt;
    % Control variables:
    Yt = XYt;
end

end


%% Construct folded second order polynomial (to multiply with g_2)
function [ZZ2_fold] = kron_poly2_fold(ZZ)
% ZZ = [Xp_dev,EPSt];

nz          = size(ZZ,2);
rws         = size(ZZ,1);
ZZ2_fold     = NaN(rws,nz*(nz+1)/2);

cnt = 0;
for i1 = 1:nz
    for i2 = 1:nz            
        if i1 <= i2
            cnt = cnt+1;
            if i1 ~= i2
                fac = 2;
            else
                fac = 1;
            end
            ZZ2_fold(:,cnt) = fac*ZZ(:,i1).*ZZ(:,i2)  ;          
        end       
    end
end

end


%% Construct folded third order polynomial (to multiply with g_3)
function [ZZ3_fold] = kron_poly3_fold(ZZ)
% ZZ = [Xp_dev,EPSt];

nz          = size(ZZ,2);
rws         = size(ZZ,1);
ZZ3_fold     = NaN(rws,nz*(nz+1)*(nz+2)/6);

cnt = 0;
for i1 = 1:nz
    for i2 = 1:nz 
        if i1 <= i2
        
            for i3 = 1:nz
                if i2 <= i3
                    cnt = cnt+1;
                    if i1 ~= i2 && i1 ~= i3 && i2 ~= i3
                        fac = 6;
                    elseif i1 == i2 && i1 == i3 && i2 == i3
                        fac = 1;
                    else
                        fac = 3;
                    end
                    ZZ3_fold(:,cnt) = fac*ZZ(:,i1).*ZZ(:,i2).*ZZ(:,i3) ;
                end
            end
        end       
    end
end

end


%% Construct second order polynomial
function [UV] = kron_poly2(UU,VV)

nu          = size(UU,2);
nv          = size(VV,2);
rws         = size(UU,1);
if rws ~= size(VV,1)
    error('Number of rows in U and V not same')
end
UV     = NaN(rws,nu*nv);

cnt = 0;
for i1 = 1:nu
    for i2 = 1:nv
        cnt = cnt+1;
        UV(:,cnt) = UU(:,i1).*VV(:,i2);
    end
end

end
