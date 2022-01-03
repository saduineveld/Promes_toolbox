function plot_pol_stnd_rbc(par,SS,GRID,POL)
%function plot_pol_stnd_rbc(par,SS,GRID,POL)
%
% Plots the policy function of standard RBC model

% Sijmen Duineveld, updated for Promes_v5.0.0 December 2021, s.a.duineveld@outlook.com

%Set lower and upper bound for capital in plot:
LK_dev1  = 0.4;% percentage deviation from kss
lb1(1)   = -LK_dev1 + log(SS.Kss); 
ub1(1)   =  LK_dev1 + log(SS.Kss);

% Set lower and upper bound for log(z)  in plot:
LZ_fac1  = 6;% 
lb1(2)   = -LZ_fac1*...
    sqrt( par.sigma_z^2 / (1-par.rho_z^2) ); 
ub1(2)   =  LZ_fac1*...
    sqrt( par.sigma_z^2 / (1-par.rho_z^2) );

nod = 101; %number of nodes in plots




%% 2D plots
nl = 3;% nl: number of lines in plot
if mod(nl,2) ~= 1
    error('nl should be odd');
else 
    mid = (nl +1)/2;
end

[PLT] = get_vars_plot(GRID,POL,nl,101,lb1,GRID.lb,ub1,GRID.ub);

[colors] = gen_color(nl,[0,0.5,0],[0,0,0.5]);%colors of winter colormap


lin_stl = {'-.','-','--','-.',':'};
for in = 1:size(PLT,1)
    figure
    if in == 1
        xlbl ='$K_{t}$ (logs)';     
        lgd = {['log(Z) = ',num2str(PLT{in,1}.LZ(1))],['log(Z) = ',num2str(PLT{in,mid}.LZ(1))],['log(Z) = ',num2str(PLT{in,end}.LZ(1))]};
    elseif in == 2
        xlbl = '$Z_{t}$ (logs)';
        lgd = {['log(K) = ',num2str(PLT{in,1}.LK(1))],'log(K) = log(K_{ss})',['log(K) = ',num2str(PLT{in,end}.LK(1))]};
    end
    min_LC = + 1e6;
    max_LC = -1e6;
    for ij = 1:size(PLT,2)
        plt = PLT{in,ij};
        if ij ~=1 && ij ~= size(PLT,2) && ij ~= mid%in case more than 3 lines are used
            %plot(plt.LX,plt.LC,lin_stl{1,ij},'LineWidth',1.5,'HandleVisibility','off')
        else
            plot(plt.LX,plt.LC,lin_stl{1,ij},'LineWidth',1.5,'Color',colors(ij,:));
        end
        hold all;
        
        min_LC = min(min_LC,min(plt.LC));
        max_LC = max(max_LC,max(plt.LC));
    end
    plot(repmat(GRID.lb(in),1,2),[min_LC,max_LC],':k');%,'LineWidth',1.5)
        
    plot(repmat(GRID.ub(in),1,2),[min_LC,max_LC],':k','HandleVisibility','off');%,'LineWidth',1.5)
    
    clear ij plt;
    xlabel(xlbl,'Interpreter','latex')
    ylabel('$C_{t}$ (logs)','Interpreter','latex')
    
    lgd = [lgd,'bounds of grid'];
    
    legend(lgd,'Location','NorthWest');
    colormap winter
    title('2D policy function for consumption')
    
    clear lgd;
end

%% 3D plot
skip_3D = 0;
if skip_3D ~= 1
f3D = figure;
%construct grid with lb1 and ub1:
vecs_tmp    = constr_vecs([nod,nod],'equi','up',lb1,ub1);
xx_tmp      = constr_grid(vecs_tmp);
LC_tmp      = get_pol_var(POL,xx_tmp,GRID);

x1 = reshape(xx_tmp(:,1),[nod,nod]);
x2 = reshape(xx_tmp(:,2),[nod,nod]);
LC_grid = reshape(LC_tmp,[nod,nod]);


surf(x1,x2,LC_grid);
xlabel('$K_{t}$ (logs)','Interpreter','latex');
ylabel('$Z_{t}$ (logs)','Interpreter','latex');
zlabel('$C_{t}$ (logs)','Interpreter','latex');
title('3D policy function consumption')
end

end


%% Plot 2D - policy functions
function [PLT] = get_vars_plot(GRID,POL,nl,nod,lb1,lb2,ub1,ub2)
% Inputs:
% nl = number of lines in plot
% nod = number of nodes
% lb1 = lower bounds if variable is plotted (1 x 2 vector)
% lb2 = lower bounds if variable is not plotted
% ub1 = upper bounds if variable is plotted (1 x 2 vector)
% ub2 = upper bounds if variable is not plotted

LK_short = linspace(lb2(1),ub2(1),nl);
LZ_short = linspace(lb2(2),ub2(2),nl);

PLT = cell(GRID.nn,nl);
for stvr = 1:2%loop over state variables (LK,lZ)
    for jj = 1:nl%loop over lines in graph (value other state variable)
        if stvr == 1%LK on x-axis
            LK = linspace(lb1(stvr),ub1(stvr),nod)'; 
            LZ = repmat(LZ_short(jj),nod,1);
            LX = LK;
        else%LZ on x-axis
            LK = repmat(LK_short(jj),nod,1);
            LZ = linspace(lb1(stvr),ub1(stvr),nod)';
            LX = LZ;
        end
        xx = [LK,LZ];
        LC = get_pol_var(POL,xx,GRID);
        
        PLT{stvr,jj}.LC = LC;
        PLT{stvr,jj}.LX = LX;
        PLT{stvr,jj}.LZ = LZ;
        PLT{stvr,jj}.LK = LK;
    end
end

end

function [colors] = gen_color(nodes,lb,ub)

d1 = linspace(lb(1),ub(1),nodes);
d2 = linspace(lb(2),ub(2),nodes);
d3 = linspace(lb(3),ub(3),nodes);

colors = NaN(nodes,3);
for in = 1:nodes
    colors(in,:) = [d1(in),d2(in),d3(in)];   
end

end
