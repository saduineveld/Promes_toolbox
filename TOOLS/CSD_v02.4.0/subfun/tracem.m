function x = tracem(X)
% Implements the matrix trace operator


%{
 
    Copyright: Alfred Maußner

    First version: 05 December 2016

    Purpose: the matrix trace operator

    The matrix trace operator: given a m*n by n matrix with blocks
    X=[X1; X2; ...; Xm] the function returns the vector x=[trace(X1); trace(X2); ..., trace(Xm)];

%}

% get dimensions
[mn,n]=size(X);

% check for proper intput
if mod(mn,n)>0;
    fprintf('X is not a m*n by n matrix. x=0 will be returned');
    x=0;
    return;
end;   
m=mn/n;
x=zeros(m,1);
for i=1:m;
    x(i)=trace(X(1+(i-1)*n:i*n,:));
end;
return;
end

