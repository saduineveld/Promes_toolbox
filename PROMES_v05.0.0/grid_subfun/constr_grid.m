function [xx] = constr_grid(gridvecs)
% function [xx] = constr_grid(gridvecs)
%
% Constructs a matrix of all gridpoints, with each variable represented in 
% a column vector. 
%
%
% Input:
%       gridvecs    = cell array (1,nn) with each cell containing the
%                       vector of gridpoints in that dimension
%
% Output:
%       xx          = matrix of all gridpoints (dimensions: mm x nn, with 
%                       mm = prod(qq), and qq a 1 x nn vector of gridpoints
%                       in each dimension)

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

nn = size(gridvecs,2);

mm = 1;
for in = 1:nn
    mm      = mm*length(gridvecs{1,in});
end
xx      = NaN(mm,nn);

if nn == 1
    xx = reshape(gridvecs{1,1},[],1);
elseif nn == 2;
    [x1,x2] = ndgrid(gridvecs{1,1},gridvecs{1,2});
    
    xx(:,1) = reshape(x1,[],1);
    xx(:,2) = reshape(x2,[],1);
elseif nn == 3;
    [x1,x2,x3] = ndgrid(gridvecs{1,1},gridvecs{1,2},gridvecs{1,3});    
  
    xx(:,1) = reshape(x1,[],1);
    xx(:,2) = reshape(x2,[],1);
    xx(:,3) = reshape(x3,[],1);
    
elseif nn == 4;
    [x1,x2,x3,x4] = ndgrid(gridvecs{1,1},gridvecs{1,2},gridvecs{1,3},...
                                gridvecs{1,4});    
  
    xx(:,1) = reshape(x1,[],1);
    xx(:,2) = reshape(x2,[],1);
    xx(:,3) = reshape(x3,[],1);
    xx(:,4) = reshape(x4,[],1);
    
elseif nn == 5;
    [x1,x2,x3,x4,x5] = ndgrid(gridvecs{1,1},gridvecs{1,2},gridvecs{1,3},...
                                gridvecs{1,4},gridvecs{1,5});
    
    xx(:,1) = reshape(x1,[],1);
    xx(:,2) = reshape(x2,[],1);
    xx(:,3) = reshape(x3,[],1);
    xx(:,4) = reshape(x4,[],1);
    xx(:,5) = reshape(x5,[],1);
    
elseif nn == 6;
    [x1,x2,x3,x4,x5,x6] = ndgrid(gridvecs{1,1},gridvecs{1,2},gridvecs{1,3},...
                                gridvecs{1,4},gridvecs{1,5},gridvecs{1,6});
    
    xx(:,1) = reshape(x1,[],1);
    xx(:,2) = reshape(x2,[],1);
    xx(:,3) = reshape(x3,[],1);
    xx(:,4) = reshape(x4,[],1);
    xx(:,5) = reshape(x5,[],1);
    xx(:,6) = reshape(x6,[],1);
else
    error('nn > 6 not implemented')
end
    
end