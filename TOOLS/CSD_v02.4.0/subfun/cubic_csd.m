%% Get 3ord terms
function [Hx_www,Hy_www,Hx_ssw,Hy_ssw,Hx_sss,Hy_sss] = cubic_csd(Hx_w,Hy_w,D,H,Omega,Rho,nx,nz,ny,Syl,T,Skew,Hx_ww,Hy_ww,Hx_ss,Hy_ss,aux)
%{
 
    Copyright: Alfred Maußner
    First version: 07 December 2016
    Revision     : 12 January 2017 (adapted to DSGE class)

    This script computes the matrices of a thrid-order approximation of the DSGE model
    described in the program SolveModel_T.

    From the calling function or script it needs to know

           nx: integer, the number of state variables (elements in x(t))
           nz: integer, the number of shocks
           ny: integer, the number of control or jump variables (elements in y(t))
          Rho: nz by nz matrix with the autoregressive coefficients of the VAR(1) defining the dynamics of the shocks
            D: (nx+ny) by 2*(nx+nz+ny) matrix, the Jacobian of the models equations computed at the stationary point.
            H: (nx+ny)*2*(nx+nz+ny) by 2*(nx+nz+ny) matrix, the (nx+ny) stacked Hessian matrices of the model's equations
               computed at the stationary point
            T: (nx+ny)*4*(nx+nz+ny)^2 by 2*(nx+ny+nz) matrix, the (nx+ny) stacked matrices of thrid-order partial
               derivatives of the model's equilibrium conditions computed at the stationary point
         Hx_w: nx by (nx+nz) matrix, the linear part of the solution for the vector of states x(t+1)
         Hy_w: ny by (nx+nz) matrix, the linear part of the solution for the vector of controls y(t)
         Hx_ww: nx*(nx+nz) by (nx+nz) matrix, the quadratic part of the solution for the vector of states x(t+1)
         Hy_ww: ny*(nx+nz) by (nx+nz) matrix, the quadratic part of the solution for the vector of controls y(t)
         Hx_ss: nx by 1 matrix, the quadratic part of the solution for the states with respect to the perturbation parameter
         Hy_ss: ny by 1 matrix, the quadratic part of the solution for the costates with respect to the 
                perturbation parameter
            fw: 2*(nx+ny+nz) by (nx+nz) matrix, the derivative of the inner functions wrt the state variables 
           M11: ny by 1 vector which is part of the matrix fss
      Syl, flag, if true, the solution will be obtained via solving the Sylvester equation

%}

% ---------------------------------------------- Hy_www and Hx_www --------------------------------------- %

% first, get the reordered matrices Hy_wwt and Hx_wwt: due to the symmetry of the matrices 
% this can be achieved via the reshape command. The loop below writes each colum j=1, ..., nw
% of block i of Hy_ww (or Hx_ww) as a row and places them next to each other.

nmat    = aux.nmat;
smat    = aux.smat;
M11     = aux.M11;
fw      = aux.fw;

nw = nx + nz;
na = nx + ny + nz;

Hy_wwt=zeros(ny,nw^2);
for i=1:ny;
    Hy_wwt(i,:)=reshape(Hy_ww((i-1)*nw + (1:nw),:),1,nw^2);
end;

Hx_wwt=zeros(nx,nw^2);
for i=1:nx;
    Hx_wwt(i,:)=reshape(Hx_ww((i-1)*nw + (1:nw),:),1,nw^2);
end;

% second, get uw, uwt, and fwt
uw=[Hx_w;[zeros(nz,nx) Rho]]; % C2 in script Quadratic
uwt=zeros(nw^2,nw^2);
for i=1:nw;    
    uwt(1+(i-1)*nw:i*nw,:)=kron(eye(nw),uw(i,:));
end;
i1=2*na;
fwt=zeros(i1*nw,nw^2);
for i=1:i1;
    fwt(1+(i-1)*nw:i*nw,:)=kron(eye(nw),fw(i,:));
end;

% third, get uww, uwwt, fww, and fwwt
uww=[Hx_ww; zeros(nz*nw,nw)];
uwwt=[Hx_wwt;zeros(nz,nw^2)];
fww=[Hx_ww; ...
    zeros(nz*nw,nw);...
    kron(eye(ny),uw')*Hy_ww*uw+kron(Hy_w(:,1:nx),eye(nw))*Hx_ww;...
    zeros(nx*nw,nw);...
    zeros(nz*nw,nw);...
    Hy_ww];
fwwt=[Hx_wwt;...
      zeros(nz,nw^2);...
      Hy_wwt*kron(uw,uw)+Hy_w(:,1:nx)*Hx_wwt;...
      zeros(nw,nw^2);...
      Hy_wwt];

Q1=kron(D(:,1:nx),eye(nw^2));
Q2=kron(D(:,1+nx+nz:na),eye(nw^2));
Q3=kron(D(:,1+na+nx+nz:2*na),eye(nw^2));
R2=[Hx_w; [zeros(nz,nx) Rho]];
R1=kron(kron(eye(ny),R2'),R2');
R3=kron(Hy_w(:,1:nx),eye(nw^2));
P1=kron(kron(eye(nx+ny),fw'),fw')*T*fw;
P2=kron(eye(nx+ny),fwwt')*H*fw;
P3=kron(kron(eye(nx+ny),fw'),eye(nw))*kron(H,eye(nw))*fww;
P4=kron(eye(nx+ny),fwt')*kron(H,eye(nw))*fww;
P5=kron(eye(ny),uwwt')*Hy_ww*uw + kron(kron(eye(ny),uw'),eye(nw))*kron(Hy_ww,eye(nw))*uww + ...
   kron(eye(ny),uwt')*kron(Hy_ww,eye(nw))*uww;
P5=Q2*P5;
if Syl;
   atilde=[Q3, (Q1+Q2*R3)];
   btilde=-R2;
   ctilde=-(P1+P2+P3+P4+P5);
   dtilde=[Q2*R1,zeros((nx+ny)*nw^2,nx*nw^2)];
   etilde=eye(nw);
   ftilde=zeros((nx+ny)*nw^2,nw);
   rmat=Sylvester(atilde,btilde,ctilde,dtilde,etilde,ftilde);
   Hy_www=rmat(1:ny*nw^2,:);
   Hx_www=rmat(ny*nw^2+1:(nx+ny)*nw^2,:); 
else
    xvec=-(P1+P2+P3+P4+P5);    
    bigmat=[kron(eye(nw),Q3)+kron(R2',Q2*R1), ...
            kron(eye(nw),(Q1+Q2*R3))];
    vecEG=linsolve(bigmat,xvec(:));
    Hy_www=reshape(vecEG(1:ny*nw^3),ny*nw^2,nw);    
    Hx_www=reshape(vecEG(1+ny*nw^3:(nx+ny)*nw^3),nx*nw^2,nw);    
end;

if nz > 0.5

% ------------------------------------------- Hy_ssw and Hx_ssw --------------------------------------- %
fss=[Hx_ss; zeros(nz,1); (M11+Hy_ss+Hy_w(:,1:nx)*Hx_ss);zeros(nw,1);Hy_ss];
A1=kron(eye(nx+ny),fw')*H*fss;
Hy_wzt=zeros(ny,nw*nz);
for i=1:ny;
    Hy_wzt(i,:)=reshape(Hy_ww((i-1)*nw + (nx+1:nw),:),1,nw*nz);
end;
fsw=[zeros(nw,nw*nz);(Hy_wzt*kron(uw,eye(nz)));zeros(nw+ny,nw*nz)];
A2=2*tracem(kron(eye(nx+ny),fsw')*H*nmat*smat);
A3=tracem(kron(kron(eye(nx+ny),fw'),nmat')*T*nmat*smat);
A4=kron(D(:,nw+1:nw+ny),eye(nw))*kron(eye(ny),uw')*Hy_ww(:,1:nx)*Hx_ss;
Hy_wzz=zeros(ny*nw*nz,nz);
nyw=ny*nw;
for i=1:nyw;
    Hy_wzz((i-1)*nz+(1:nz),:)=Hy_www((i-1)*nw+(nx+1:nw),nx+1:nw);
end;
A5=kron(D(:,nw+1:nw+ny),eye(nw))*tracem(kron(eye(nw*ny),smat)*kron(kron(eye(ny),uw'),eye(nz))*Hy_wzz);
%B11=kron(D(:,1+na+nw:2*na),eye(nw));
%B12=kron(D(:,nw+1:na),eye(nw))*kron(eye(ny),uw');
bigmat=kron(D(:,1+na+nw:2*na),eye(nw))+kron(D(:,nw+1:na),eye(nw))*kron(eye(ny),uw');
%B2=kron(D(:,1:nx),eye(nw))+kron(D(:,nw+1:na),eye(nw))*kron(Hy_w(:,1:nx),eye(nw));
bigmat=[bigmat,kron(D(:,1:nx),eye(nw))+kron(D(:,nw+1:na),eye(nw))*kron(Hy_w(:,1:nx),eye(nw))];
VecEG=linsolve(bigmat,-(A1+A2+A3+A4+A5));
Hy_ssw=VecEG(1:ny*nw);
Hx_ssw=VecEG(1+ny*nw:(nx+ny)*nw);

% ---------------------------------------- Hx_sss and Hy_sss ------------------------------------------ %
if Skew==0; Hx_sss=zeros(nx,1); Hy_sss=zeros(ny,1); return; end;
Hy_zzt=zeros(ny,nz^2);
Hy_zzz=zeros(ny*nz^2,nz);
for i=1:ny;
    Hy_zzt(i,:)=reshape(Hy_ww((i-1)*nw+1+nx:i*nw,1+nx:nw),1,[]);
    Hy_zzz((i-1)*nz*nz+1:i*nz*nz,:)=Hy_wzz((i-1)*nw*nz+1+nx*nz:i*nw*nz,:);
end;
A1=tracem(kron(kron(eye(nx+ny),nmat'),nmat')*T*nmat*Skew);
A2=3*tracem(kron(eye(nx+ny),[zeros(nx+nz,nz*nz);Hy_zzt;zeros(nw+ny,nz*nz)]')*H*nmat*Skew);
A3=tracem(kron(eye(ny),Skew)*Hy_zzz);
A3=D(:,nw+1:na)*A3;
bigmat=[(D(:,1+nw:na)+D(:,1+na+nw:2*na)), (D(:,1:nx)+D(:,1+nw:na)*Hy_w(:,1:nx))];
VecEG=linsolve(bigmat,-(A1+A2+A3));
Hy_sss=VecEG(1:ny);
Hx_sss=VecEG(1+ny:nx+ny);
else
    
    Hy_ssw = [];
    Hx_ssw = [];
    
    Hy_sss = [];
    Hx_sss = [];
end

end