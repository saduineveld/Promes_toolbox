function [YY] = get_pol_var(POL,xx,GRID,i_pol,spec_opt)
% [YY] = get_pol_var(POL,xx,GRID,i_pol,spec_opt)
%
% Purpose: get policy variables, using policy functions and state vectors as
% input
%
% INPUT:
%   POL     = structure with policy function, which contains at least fields
%               - algo
%               For algo 'cheb_xxx','smol_xxx', and 'mono_xxx'
%               - theta
%               For algo 'spl_xxx'
%               - pp_y                           
%               If spec_opt='old_pol' (for 'tmi' variants):
%               - pp_y_old if algo is 'spl_xxx'
%               - theta_old for 'cheb_tmi, or 'smol_tmi'
%
%   xx      = mm x nn matrix with datapoints (mm data points, nn state variables)
%
%   GRID    = structure with necessary fields (see prepgrid for a list)
%
%   i_pol   = index of the policy variable (ie. if i_pol = 2, the policy 
%               function theta(:,i_pol)/pp_y{1,i_pol} is used
%               if POL.dd = 1 and i_pol unspecified then i_pol is set to 1;
%
%   spec_opt=  (optional)
%               if 'old_pol': theta_old/pp_y_old is used to evaluate YY
%               if 'ini_grid' (only for methods using polynomials, 
%                   'cheb_xxx','smol_xxx','mono_xxx'): use initial 
%                   polynomial GRID.XX_poly or GRID.XX_poly_dw to 
%                   save computation time. NOTE: spec_opt='ini_grid' will 
%                   ignore input argument xx.
%
% OUTPUT:
%       YY          = policy variable(i_pol) at the gridpoints (mm x 1)
%
% 
%
% Uses:
%       For methods 'cheb_xxx'
%           - scal_mat_dw, get_poly_ani
%       For methods 'mono_xxx': 
%           - get_poly_ani
%       For methods 'smol_xxx':
%           - scal_mat_dw, Smolyak_Polynomial (coded by others, see 
%               that file for the license)

% Sijmen Duineveld, updated December 2021, s.a.duineveld@outlook.com

% Copyright 2019-2021 Sijmen Duineveld
% This file is part of the Promes toolbox. The Promes toolbox is free software 
% under the terms of the GNU General Public License version 3. 
% The Promes toolbox is distributed without any warranty.

if nargin < 5
    spec_opt = 0;
end
if (nargin < 4 || isempty(i_pol)) && POL.dd > 1
    error('i_pol is unassigned, but there are multiple policy variables');
elseif POL.dd == 1
    i_pol = 1;%there is only one policy variable
end

if strcmp(POL.algo,'spl_tmi') || strcmp(POL.algo,'spl_dir')
    if strcmp(spec_opt,'ini_grid') && strcmp(POL.algo,'spl_tmi')
        error('With spl_tmi the policy at the initial grid is set directly in POL.YY. See Chapter Model File in the manual'); 
        %NOTE: with spl_dir the option 'ini_grid' is simply ignored
    end
    
    if strcmp(spec_opt,'old_pol')%use the old policy function
        if strcmp(POL.algo,'spl_dir')
            error('Algorithm spl_dir should not use old policy.')
        else
            YY = get_YY(POL.pp_y_old{1,i_pol},xx,GRID.nn);
        end
    else
        YY = get_YY(POL.pp_y{1,i_pol},xx,GRID.nn);
    end
    
elseif strcmp(POL.algo,'smol_tmi') || strcmp(POL.algo,'smol_dir')
    if strcmp(spec_opt,'ini_grid')
        if strcmp(POL.algo,'smol_tmi')
            error('With smol_tmi the policy at the initial grid is set directly in POL.YY. See Chapter Model File in the manual'); 
        else
            YY = GRID.XX_poly_dw*POL.theta(:,i_pol);
        end
    else
        xx_dw   = scal_mat_dw(xx,GRID.lb,GRID.ub);
        XX_poly_dw   = Smolyak_Polynomial(xx_dw,GRID.nn,GRID.mu_max,GRID.smol_el_ani);    
        if strcmp(spec_opt,'old_pol')%use the old policy function
            if strcmp(POL.algo,'smol_dir')
                error('Algorithm smol_dir should not use old policy.')
            else
                YY = XX_poly_dw*POL.theta_old(:,i_pol);
            end
        else       
            YY = XX_poly_dw*POL.theta(:,i_pol);
        end  
    end
    
elseif strcmp(POL.algo,'cheb_tmi') 
    if strcmp(spec_opt,'ini_grid')
        error('With cheb_tmi the policy at the initial grid is set directly in POL.YY. See Chapter Model File in the manual'); 
        %YY      = GRID.XX_poly_dw*POL.theta(:,i_pol);
    else
        xx_dw   = scal_mat_dw(xx,GRID.lb,GRID.ub);
        XX_poly_dw   = get_poly_ani(xx_dw,GRID.ord_vec,GRID.grid_type,GRID.poly_elem);
        if strcmp(spec_opt,'old_pol')%use the old policy function
            YY = XX_poly_dw*POL.theta_old(:,i_pol);
        else       
            YY = XX_poly_dw*POL.theta(:,i_pol);
        end 
    end
    
elseif strcmp(POL.algo,'cheb_mse') || strcmp(POL.algo,'cheb_gal')
    if strcmp(spec_opt,'ini_grid')
        YY      = GRID.XX_poly_dw*POL.theta(:,i_pol);
    else
        if strcmp(spec_opt,'old_pol')
            error('Algorithms cheb_mse and cheb_gal should not use old policy.')
        else
            xx_dw   = scal_mat_dw(xx,GRID.lb,GRID.ub);
            XX_poly_dw   = get_poly_ani(xx_dw,GRID.ord_vec,GRID.grid_type,GRID.poly_elem);
        
            YY      = XX_poly_dw*POL.theta(:,i_pol);
        end
    end
    
elseif strcmp(POL.algo,'mono_mse')    
    if strcmp(spec_opt,'ini_grid')
        YY      = GRID.XX_poly*POL.theta(:,i_pol); 
    else
        if strcmp(spec_opt,'old_pol')
            error('Algorithm mono_mse should not use old policy.')
        else
        	XX_poly      = get_poly_ani(xx,GRID.ord_vec,GRID.grid_type,GRID.poly_elem);
            YY      = XX_poly*POL.theta(:,i_pol);   
        end
    end    
   
else
    error('Invalid solution type')
end

end

function YY = get_YY(pp_y,xx,nn)

if nn == 1
    YY      = pp_y(xx(:,1));   
elseif nn == 2
    YY      = pp_y(xx(:,1),xx(:,2));        
elseif nn == 3
    YY      = pp_y(xx(:,1),xx(:,2),xx(:,3));        
elseif nn == 4
    YY      = pp_y(xx(:,1),xx(:,2),xx(:,3),xx(:,4));    
elseif nn == 5
    YY      = pp_y(xx(:,1),xx(:,2),xx(:,3),xx(:,4),xx(:,5));
elseif nn == 6
    YY      = pp_y(xx(:,1),xx(:,2),xx(:,3),xx(:,4),xx(:,5),xx(:,6));    
else
	error('Add extra statements if nn > 6')    
end 
    
end