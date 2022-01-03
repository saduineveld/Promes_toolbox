function [MOD] = get_deriv_csd(MOD,ord)
%[MOD] = get_deriv_csd(MOD,ord)
%
% Get analytical derivatives & some auxiliary variables

% Code is adjusted version of "CORRAM-M" package (2018) by Alfred Maussner, 
% Sijmen Duineveld, July 2019, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the CSD Toolbox. The CSD Toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The CSD Toolbox is distributed without any warranty.

MOD.nx = size(MOD.XX,2);
MOD.nz = size(MOD.ZZ,2);
MOD.ny = size(MOD.YY,2);

qq = [MOD.XXn MOD.ZZn MOD.YYn MOD.XX MOD.ZZ MOD.YY];

MOD.qq_str = cell(1,MOD.nx+MOD.nz+MOD.ny);
for in = 1:MOD.nx + MOD.nz + MOD.ny
    str_tmp = char(qq(1,in));
    str_tmp = str_tmp(1:end-2);
    MOD.qq_str{1,in} = str_tmp;     
end

[DS,HS,TS] = ana_deriv_csd(MOD.FS,qq,MOD.nx,MOD.nz,MOD.ny,ord);

MOD.DS = DS;
MOD.HS = HS;
MOD.TS = TS;

MOD.qq = qq;

MOD.ord = ord;

end

%% Get analytical derivatives
function [DS,HS,TS] = ana_deriv_csd(FS,qq,nx,nz,ny,ord)

nu      = 0;%number of intratemporal equations
nw      = nx+nz;
nxy     = nx+ny;
na      = nxy+nz;


%Compute the first and second derivatives of f
DS = jacobian(FS,qq);

if ord > 1.5 
    HS1 = jacobian(DS(:),qq);

    n_m=2*na;
    imat=zeros(n_m,nxy);
    for i=1:n_m;
        imat(i,:)=[(i-1)*nxy+1:i*nxy];
    end;
    HS=HS1(imat(:),:);
    
    if ord > 2.5
            TS=jacobian(HS1(:),qq);
            m2=n_m^2;
            imat=zeros(m2,nxy);
            %clear imat;
            for i=1:m2;
                imat(i,:)=[(i-1)*nxy+1:i*nxy];
            end;
            TS=TS(imat(:),:);
    else
        TS = [];
    end

else
    HS = [];
    TS = [];
end 

end