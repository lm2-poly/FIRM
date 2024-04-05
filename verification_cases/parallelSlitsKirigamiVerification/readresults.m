%Reads the data to plot the reconfiguration of a ribbon kirigami sheet
%Compares with Marzin et al.'s data (Marzin et al., 2022)
%Author: Danick Lamoureux
%Creation date: 2023-05-21

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

clear all;
close all;
clc;

green = [75, 204, 40]/256;
blue = [10, 157, 255]/256;
red = [256,0,0]/256;
orange = [255, 201, 40]/256;
white = [220,220,220]/256;

cut_thickness = 0.01;

%% Parameters
% Geometry
cut_thickness = 1E-5;

t = 100E-6;
Ls = 17.2E-3;
dx = 4.4E-3/2;
dy = 3.8E-3;

Nx = 5;
Ny = 31;

L = (Ny+1)*dy;
W = Nx*(Ls+2*dx);

S = 0.15*0.15;

% Material properties
E = 4E9; % Young's modulus
nu = 0.3; % Poisson's coefficient

K = 20*E*t^3*dy*Nx/((Ls-2*dx)^3*2*Ny);

%Load definition
D = E*t^3/(12*(1-nu^2));
CD = 2.0;
rho = 997;



%% Non blocked
cd data

%% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

middle_nodes = [];
middle_nodes_ys = [];
for i = 1:length(coord)
    if coord(i,1) >= (0.98*L/2) && coord(i,1) <= (1.02*L/2)
        middle_nodes(end+1) = i;
        middle_nodes_ys(end+1) = coord(i,2);
    end
end

[~, indices] = sort(middle_nodes_ys);
middle_nodes = middle_nodes(indices);

list = dir;
Us = [];
u = {};
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 8
        if item(1:8) == "solution"
            fprintf(item+"\n");
            res = open(item);
            u{end+1} = res.displacements;
            startind = find('U' == item);
            Us(end+1) = str2double(item(startind+2:end-4));
        end
    end
end


cd ..

[sorted_Us, indices] = sort(Us);
sorted_CYs = rho*sorted_Us.^2*W/K;
sorted_u = u(indices);
CYs = [];
CYsB = [];
Rs = [];
ymaxLs = [];
middle_defs = zeros(length(middle_nodes), 2);

for i = 1:length(sorted_Us)
    u = sorted_u{:,i};
    u = matrix2vec(u);
    Rs(end+1) = get_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, 0);
    factor = S/(S-L*W);
    CYs(end+1) = rho*(sorted_Us(i)/factor)^2*W/K;
    
    [x, disp_norms] = get_displacements(u);
    ymaxLs(end+1) = max(abs(x(:,3)))/L;
    def_nodes = get_deformed_nodes(coord, x);
end

%% Blockage
cd variable_blocked

%% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

middle_nodes = [];
middle_nodes_ys = [];
for i = 1:length(coord)
    if coord(i,1) >= (0.98*L/2) && coord(i,1) <= (1.02*L/2)
        middle_nodes(end+1) = i;
        middle_nodes_ys(end+1) = coord(i,2);
    end
end

[~, indices] = sort(middle_nodes_ys);
middle_nodes = middle_nodes(indices);

list = dir;
Us = [];
u = {};
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 8
        if item(1:8) == "solution"
            fprintf(item+"\n");
            res = open(item);
            u{end+1} = res.displacements;
            startind = find('U' == item);
            Us(end+1) = str2double(item(startind+2:end-4));
        end
    end
end


cd ..

[sorted_Us, indices] = sort(Us);
sorted_CYs = rho*sorted_Us.^2*W/K;
sorted_u = u(indices);
blockedCYs = [];
blockedRs = [];
blockedymaxLs = [];
middle_defs = zeros(length(middle_nodes), 2);

for i = 1:length(sorted_Us)
    u = sorted_u{:,i};
    u = matrix2vec(u);
    blockedRs(end+1) = get_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, 0);
    blockedCYs(end+1) = rho*sorted_Us(i)^2*W/K;
    
    [x, disp_norms] = get_displacements(u);
    blockedymaxLs(end+1) = max(abs(x(:,3)))/L;
    def_nodes = get_deformed_nodes(coord, x);
    if mod(i,10) == 0
        for j = 1:length(middle_nodes)
            middle_defs(j,1) = def_nodes(middle_nodes(j), 2)/L;
            middle_defs(j,2) = def_nodes(middle_nodes(j), 3)/L;
        end
        f4 = figure(4);
        hold on;
        factor = i/length(sorted_CYs);
        color = factor*red + (1 - factor)*white;
        plot(middle_defs(:,1), middle_defs(:,2), '-', 'Color', color);
        axis("tight");
        xlabel("$x/L$");
        ylabel("$z/L$");
    end        
end

figure(4);
c = colorbar;
c.Location = "northoutside";
c.TickLabels = [0.001 100];
c.Ticks = [0 1];
c.Label.Interpreter = 'latex';
c.Label.String = '$C_Y$';

f4.Units = "centimeters";
f4.Position = [1 1 6, 4];

%% Displacement plot

model = table2array(readtable("MarzinEtAl.csv"));

f5 = figure(3);
hold on;
grid on;
plot(CYs, ymaxLs, '-', "color", red, "linewidth", 2)
plot(blockedCYs, blockedymaxLs, '-', "color", red, "linewidth", 2)
scatter(model(:,1), model(:,2), 'ok', "filled", "linewidth", 0.5);
xlabel('$C_Y$')
ylabel('$z/L$')
set(gca, 'XScale', 'log')
axis("tight")
xlim([10^(-5), 10^2])
f5.Units = "centimeters";
f5.Position = [1 1 20.5, 8];