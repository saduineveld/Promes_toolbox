function [nW,othstb,pi_y,phi_y,pi_t,phi_t] = InvSubGen(Gx0,Gx1,Gt0,Gt1,Pi,ny,allsol,tol)
%
% [nW,othstb,pi_y,phi_y,pi_t,phi_t] = InvSubGen(Gx0,Gx1,Gt0,Gt1,Pi,ny,allsol,tol)
%
% Code to implement the solution method of Galizia (2020)
%   "Solving Rational Expectations Models Featuring Limit Cycles (or Chaos)
%     Using Perturbation Methods"
%
% This function computes first-order approximations to elements of
% S_{empty}, where the linearized system is given by
%
%   Gx1*E_t[x(t+1)] + Gx0*x(t) + (Gt0+Gt1*Pi)*theta(t) = 0   ,
%   E_t[theta(t+1)] = Pi*theta(t)            ,
%
% where x(t) = (y(t),z(t)), y(t) is an ny-vector of pre-determined
% variables, z(t) is an nz-vector of jump variables, and theta(t) is an
% nt-vector of exogenous stochastic variables. Note that we do not include
% the perturbation parameter zeta, since it has no effect on the linear
% solution.
%
% Returned approximations correspond to invariant subspaces that:
%   (1) contain either all or none of any given RGE (see Assumption (A3);
%       and
%   (2) do not contain RGEs corresponding to eigenvalues that are real and
%       greater than one (see Proposition 3(i));
%
% Note regarding screening for indeterminacy:
%   - This code will return an error if "linear indeterminacy" is detected
%       (i.e., Blanchard-Kahn condition is violated).
%   - To return only candidate solutions for which the "stable" RGEs are
%       always contained in the associated w, set allsol=0 (the default).
%       This will typically be faster, but will detect fewer cases of
%       indeterminacy.
%   - To return candidate solutions corresponding to all of W*_{empty}, set
%       allsol=1. If candidate solution j excludes any stable eigenvalues
%       (indicated by setting othstb(j)=true) but satisfies the TVC this 
%       would imply the presence of indeterminacy.
%
% Inputs:
%   Gx0:    nx-by-nx endog-endog matrix
%   Gx1:    nx-by-nx endog-endog matrix
%   Gt0:   nx-by-nt endog-exog matrix (set empty for non-stochastic case)
%   Gt1:   nx-by-nt endog-exog matrix
%   Pi:    nt-by-nt exog-exog matrix
%   ny: number of pre-determined variables
%   allsol: see note above regarding indeterminacy (default=0)
%   tol: tolerance for equality of eigenvalues (default=1e-8); see MATLAB
%       functions 'ismembertol' and 'uniquetol' for details
%
% Outputs:
%   nW: number of invariant subspaces returned
%   othstb: nW-vector whose j-th element is true if there are stable
%       eigenvalues corresponding to RGEs not in the j-th subspace
%   phi_y and phi_t: nz-by-ny-by-nW and nz-by-nt-by-nW arrays containing
%       first-order approximations to phi (where phi_y is w.r.t. y, and
%       phi_t is w.r.t. to theta)
%   pi_y and pi_t: ny-by-ny-by-nW and ny-by-nt-by-nW real arrays containing
%       first-order approximations to pi (where pi_y is w.r.t. y, and pi_t
%       is w.r.t. to theta)

% Written by:
%       Dana Galizia
%       Carleton University
%       Last edit: June 30, 2020

% set allsol to default if not passed in
if nargin < 7
    allsol = 0;
elseif isempty(allsol)
    allsol = 0;
end

% set tolerance for equality of eigenvalues if not passed in
if nargin < 8
    tol = 1e-8;
elseif isempty(tol)
    tol = 1e-8;
end


nx = size(Gx0,1);
if size(Gx0,2) ~= nx
    error('Gx0 must be square.')
end
if any(size(Gx1) ~= nx)
    error('Gx1 must be the same size as Gx0.')
end

nt = size(Gt0,2);
if nt==0
    stch = false;
    if numel(Pi) > 0
        warning('Non-stochastic case, but Pi isn''t empty. Setting Pi=[].')
        Pi = [];
        nt = 0;
    end
else
    stch = true;
    if any(size(Pi) ~= nt)
        error('Pi must be square with same number of columns as Gt0.')
    end
    if size(Gt0,1) ~= nx
        error('Gt0 must have same number of rows as Gx0.')
    end
    if any(size(Gt1) ~= [nx,nt])
        error('Gt1 must be the same size as Gt0.')
    end
end

if (~isreal(Gx0)) || (~isreal(Gx1)) || (~isreal(Gt0)) || (~isreal(Gt1)) || (~isreal(Pi))
    error('Input matrices must be real.')
end

% make sure we have a valid choice for ny
if ny > nx
    error('ny cannot exceed the dimension of Gx0')
elseif ny <=0
    error('ny must be strictly positive.')
end

nz = nx - ny;    % number of jump variables

B = Gx1;
A = -Gx0;
C = -(Gt0+Gt1*Pi);

%% Set up initial factorization

[T,S,Rp,U] = qz(A,B,'real');  % obtain real generalized Schur decomposition
Ex = ordeig(T,S); % get generalized eigenvalues
Ex(abs(Ex)>1e10) = Inf; % very large eigenvalues probably actually infinite

% set up clusters for re-ordering: 4=stable; 1=Inf; 2=real,>1; 3=all others
clst = 3*ones(nx,1);    
clst(abs(Ex)<1) = 4;
clst((imag(Ex)==0) & (Ex>1)) = 2;
clst(isinf(Ex)) = 1;

[T,S,Rp,U] = ordqz(T,S,Rp,U,clst); % re-order to put stable eigenvalues first, infinite last
Ex = ordeig(T,S); % get generalized eigenvalues
Ex(abs(Ex)>1e10) = Inf; % very large eigenvalues probably actually infinite

%% Analyze eigenvalue structure

% Logical indices of Ex that are inside the unit circle
lg_stb = (abs(Ex) < 1);
nstb = nnz(lg_stb);       % number of such elements

if nstb > ny
    error('Linear indeterminacy: More stable eigenvalues than pre-determined variables.')
end

% In this block, we set 'lg_alw' as a vector of logical indices for the
% eigenvalues whose RGEs should always be included in the invariant
% subspace. If allsol=0, this is the set of RGEs corresponding to stable
% eigenvalues and exogenous eigenvalues. If allsol=1, then it is only
% exogenous eigenvlaues
if allsol == 0
    lg_alw = lg_stb;
    nalw = nstb;
else
	lg_alw = false(nx,1);
    nalw = 0;
end

% Logical indices of elements of Ex that are real and greater than one, or
% infinite, which are required to never be a part of a subspace
lg_rlgr1 = ((imag(Ex)==0) & (Ex>1));
lg_inf = isinf(Ex);
lg_nev = lg_rlgr1 | lg_inf;

% Indices of all other elements of E
lg_oth = ~(lg_alw | lg_nev);        % logical indices
ln_oth = find(lg_oth);         % linear indices
Exoth = Ex(lg_oth);               % extract these other elements
noth = numel(Exoth);            % number of them

nrem = ny-nalw;      % number of additional elements we need to choose from
                    % Ex over and above the required ones
ninf = nnz(lg_inf);  % number of infinite eigenvalues

%% Determine set of eigenvalue configurations to get projections for

if nrem == 0    % if no additional elements needed, then the set of projections 
    nW = 1;     % is trivial: it includes only the required elements of E
    othstb = 0;
    
    selcmb = lg_alw;
    
elseif nrem > noth       % in this case, to construct any invariant subspace
                            % we would have to choose eigenvalues that are
                            % real and greater than one, which is not allowed
    nW = 0;
    othstb = [];
    phi_y = [];
    pi_y = [];
    phi_t = [];
    pi_t = [];
    return;
else    % the remaining case
    
    % In this block we look for clusters of eigenvalues in Eoth that are either the same
    % or are complex conjugate pairs (or both)
    Eoth2 = [real(Exoth),abs(imag(Exoth))];   % if two rows of this matrix are identical, then
                                            % they are part of the same cluster
    [~,EtoC,CtoE] = uniquetol(Eoth2,tol,'ByRows',true);
                % The length of EtoC here gives the number of clusters
                % Each cluster is then given a different number, and CtoE 
                % gives the cluster number of the corresponding element of Eoth
    nC = length(EtoC);  % number of clusters found
    
    % nEC gives the number of elements of Eoth in each cluster.
    nEC = zeros(nC,1);
    for j = 1:nC
        nEC(j) = nnz(CtoE==j);
    end
    
    cmbs = combnk(1:noth,nrem); % compute all the different ways to choose
                % nrem elements from Eoth; each of these ways is stored in
                % a row of cmbs as indices of Eoth (rather than actual
                % values of Eoth).
    ncmbs = size(cmbs,1);   % number of such different ways
    
    % This block eliminates combinations that feature only part of a
    % cluster
    incl = true(ncmbs,1);   % 'true' here will indicate that only full clusters
                            % are included in the corresponding combination
    for j = 1:ncmbs     % for each combination
        for k = 1:nrem      % and each element of Eoth in that combination
            njk = nnz(CtoE(cmbs(j,:))==CtoE(cmbs(j,k))); % count the number of times that element
                                                         % appears in the combination
            if njk ~= nEC(CtoE(cmbs(j,k)))	% if it appears less than it should
                incl(j) = false;	% indicate that this is not a valid combination
                break;              % and move on to the next combination
            end
        end
    end
    
    cmbschk = cmbs(incl,:);     % extract only the valid combinations
    nW = size(cmbschk,1);       % number of valid combinations
    
    othstb = false(nW,1);
    
    % This block computes the logical indices of Ex to be included for each
    % valid combination. These are collected as the columns of selcmb.
    selcmb = repmat(lg_alw,1,nW);    % start from the required elements of E in all cases
    for j = 1:nW                    % for each valid combination
        selcmb(ln_oth(cmbschk(j,:)),j) = true;  % add in the elements corresponding to 
                                        % the elements of Eoth found in the j-th row of cmbschk
        othstb(j) = min(abs(Ex(~selcmb(:,j)))<1);
    end
        
end

%% Compute projections

% arrays to hold projections
phi_y = zeros(nz,ny,nW);
phi_t = zeros(nz,nt,nW);

% array to hold transition matrix of state variables
pi_y = zeros(ny,ny,nW);
pi_t = zeros(ny,nt,nW);

in1 = [true(ny,1);false(nz,1)];
in2 = [false(ny,1);true(nz-ninf,1);false(ninf,1)];
in3 = [false(nx-ninf,1);true(ninf,1)];
inz = [false(ny,1);true(nz,1)];

for j = 1:nW            % for each invariant subspace j found
    
    
    % First, re-order qz decomp so that a basis for subspace j is given by
    % the first N columns of U, and infinite eigenvalues ordered last
    
    % set up clusters for re-ordering: 3=selcmb; 1=Inf; 2=all others
    clst = 2*ones(nx,1);
    clst(selcmb(:,j)) = 3;
    clst(lg_inf) = 1;

    [Tj,Sj,Rpj,Uj] = ordqz(T,S,Rp,U,clst); % re-order
    
    
    % Now get projection onto that basis by choice of z(t).
    
    Uyy = Uj(in1,in1);
    Uyz = Uj(in1,inz);
    Uzz = Uj(inz,inz);
    S11 = Sj(in1,in1);
    T11 = Tj(in1,in1);
    
    phi_y(:,:,j) = -(Uzz')\Uyz';
    pi_y(:,:,j) = Uyy*(S11\T11)/Uyy;
    
    if stch        
        S22 = Sj(in2,in2);
        Syz = Sj(in1,inz);
        T22 = Tj(in2,in2);
        Tyz = Tj(in1,inz);
        Rp1 = Rpj(in1,:);
        Rp2 = Rpj(in2,:);
        if ninf>0
            S23 = Sj(in2,in3);
            S33 = Sj(in3,in3);
            T23 = Tj(in2,in3);
            T33 = Tj(in3,in3);
            Rp3 = Rpj(in3,:);
            
            if max(abs(S33(:))) > 1e-12     % if S33 contains any non-zero entries
                error('S_33 is non-zero. Static relationships can''t be solved uniquely.')
            end
            
            Psi3 = -T33\(Rp3*C);
            Phi2 = -S22\(T23*Psi3-S23*Psi3*Pi+Rp2*C);
            Psi2 = sylvester(S22\T22,-Pi,Phi2);
            Psi = [Psi2;Psi3];
            Phi1 = -S11\(Tyz*Psi - Syz*Psi*Pi + Rp1*C);
            
        else
            Phi2 = -S22\(Rp2*C);
            Psi = sylvester(S22\T22,-Pi,Phi2);
            Phi1 = -S11\(Tyz*Psi - Syz*Psi*Pi + Rp1*C);
        end
        
        phi_t(:,:,j) = (Uzz')\Psi;
        pi_t(:,:,j) = Uyz*Psi*Pi-Uyy*Phi1-pi_y(:,:,j)*Uyz*Psi;
        
    end
end

