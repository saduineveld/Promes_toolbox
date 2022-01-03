% Get second order terms
function [Hx_ww,Hy_ww,Hx_ss,Hy_ss,aux] = quad_csd(Hx_w,Hy_w,D,H,Omega,Rho,nx,nz,ny,Syl)

Muss = 0;%adjustment for expected mean of log normal distribution
%{

 Copyright: Alfred Maußner
 Revisions:     06 December 2016, first version
                12 January 2017 (adapted to the DSGE class)
                11 May     2017 (non zero means of innovations)

 This script computes the matrices of a second-order approximation of the DSGE model
 described in the program SolveModelA3.

 From the calling function or script it needs to know

           nx: integer, the number of state variables (elements in x(t))
           nz: integer, the number of shocks
           ny: integer, the number of control or jump variables (elements in y(t))
          Rho: nz by nz matrix with the autoregressive coefficients of the VAR(1) defining the dynamics of the shocks
        Omega: nz by nz matrix, defines the correlation structure of innovations
         Muss: nz by 1 vector, the second-order derivatives of the mean function wrt the perturbation parameter
            D: (nx+ny) by 2*(nx+nz+ny) matrix, the Jacobian of the models equations computed at the stationary point.
            H: (nx+ny)*2*(nx+nz+ny) by 2*(nx+nz+ny) matrix, the (nx+ny) stacked Hessian matrices of the model's equations
               computed at the stationary point
          Hx_w: nx by (nx+nz) matrix, the linear part of the solution for the vector of states x(t+1)
          Hy_w: ny by (nx+nz) matrix, the linear part of the solution for the vector of costates y(t)
           Syl: flag, if true, the solution will be obtained via solving the Sylvester equation
%}

na = nx + ny + nz;

p1 = Hx_w;
p2 = [zeros(nz,nx) Rho];
p3 = [Hy_w(:,1:nx)*Hx_w(:,1:nx) (Hy_w(:,1:nx)*Hx_w(:,nx+1:nx+nz)+Hy_w(:,nx+1:nx+nz)*Rho)];
p4 = [eye(nx) zeros(nx,nz)];
p5 = [zeros(nz,nx) eye(nz)];
p6 = Hy_w;

fw=[p1;p2;p3;p4;p5;p6];

A1=kron(eye(nx+ny),fw')*H*fw;
B1=kron(D(:,1+2*na-ny:2*na),eye(nx+nz)); % gy
B2=kron(D(:,1:nx),eye(nx+nz)); % gx'
B3=kron(D(:,1+nx+nz:nx+nz+ny),eye(nx+nz)); % gy'
C2 = [Hx_w;[zeros(nz,nx) Rho]];
C1=kron(eye(ny),C2');
C3=kron(Hy_w(:,1:nx),eye(nx+nz));
if Syl;
    atilde=[B1, (B2+B3*C3)];
    btilde=-C2;
    ctilde=-A1;
    dtilde=[B3*C1,zeros((nx+ny)*(nx+nz),nx*(nx+nz))];
    etilde=eye(nx+nz);
    ftilde=zeros((nx+ny)*(nx+nz),nx+nz);
    rmat=Sylvester(atilde,btilde,ctilde,dtilde,etilde,ftilde);
    Hy_ww=rmat(1:ny*(nx+nz),:);
    Hx_ww=rmat(ny*(nx+nz)+1:(nx+ny)*(nx+nz),:);
else
    bigmat=kron(eye(nx+nz),B1)+kron(C2',B3*C1);
    bigmat=[bigmat,kron(eye(nx+nz),(B2+(B3*C3)))];
    vecEG=linsolve(bigmat,-A1(:));
    Hy_ww=reshape(vecEG(1:ny*(nx+nz)^2),[ny*(nx+nz),nx+nz]);
    Hx_ww=reshape(vecEG(ny*(nx+nz)^2+1:(ny+nx)*(nx+nz)^2),[nx*(nx+nz),nx+nz]);
end;

if nz > 0.5;
% Solve for Hx_ss and Hy_ss
smat=Omega*(Omega');
nmat=[zeros(nx,nz);eye(nz);Hy_w(:,1+nx:nx+nz);zeros(nx+nz+ny,nz)];
xvec=tracem(kron(eye(nx+ny),nmat')*H*(nmat*smat));
M11=zeros(ny,1);
for i=1:ny;
   M11(i,1)=trace(smat*Hy_ww((i-1)*(nx+nz)+nx+1:i*(nx+nz),1+nx:nx+nz));
end;
xvec=xvec+(D(:,1+nx+nz:na)*M11); % gy'
xvec=xvec+(D(:,1+nx:nx+nz)+D(:,1+nx+nz:na)*Hy_w(:,1+nx:nx+nz))*Muss;
bigmat=D(:,nx+nz+1:nx+nz+ny)+D(:,2*(nx+nz)+ny+1:2*na);
bigmat=[bigmat,D(:,1:nx)+D(:,1+nx+nz:nx+nz+ny)*Hy_w(:,1:nx)];
VecEG=linsolve(bigmat,-xvec);
Hy_ss=VecEG(1:ny);
Hx_ss=VecEG(ny+1:ny+nx);
else
    Hy_ss = [];
    Hx_ss = [];
    M11   = [];
    nmat  = [];
    smat  = [];
end

aux.fw      = fw;
aux.M11     = M11;
aux.nmat    = nmat;
aux.smat    = smat;

end