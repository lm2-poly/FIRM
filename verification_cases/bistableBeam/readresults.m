%Reads the data to plot the reconfiguration during snapping of a bistable
%beam
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
red = orange;

map = [linspace(1, red(1), 150)', linspace(1, red(2), 150)', linspace(1, red(3), 150)';
        linspace(red(1), 0, 150)', linspace(red(2), 0, 150)', linspace(red(3), 0, 150)'];

cut_thickness = 0.01;

%% Parameters
% Geometry
Lx = 50; % Length in direction 1 (before applying symmetry, if any)
Ly = 10; % ... in direction 2 (before applying symmetry, if any)
h = 0.1; % Thickness

dh = 20;

% Material properties
E = 3E3; % Young's modulus
nu = 0.; % Poisson's coefficient

%Load definition
D = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 1.225*(10^(-3))^3;

%% Symmetry
cd half

% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));
L = Lx/2;
mid = 0;
for i = 1:length(coord)
    if coord(i,1) < L
        s = coord(i,1)/L;
        coord(i,3) = dh*(3*s-s^3)/2;
    else
        s = 1-(coord(i,1) - L)/L;
        coord(i,3) = dh*(3*s-s^3)/2;
    end
    if s == 1
        mid = i;
    end
end

list = dir;
CY = [];
u = {};
names = [];
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 4
        if item(end-3:end) == ".mat"
            fprintf(item+"\n");
            res = open(item);
            names(end+1) = str2double(item(1:end-4));
            u{end+1} = res.u;
            CY(end+1) = res.lambda;
        end
    end
end


cd ..

[~, indices] = sort(names);
sorted_CYs = CY(indices);
sorted_u = u(indices);
ymaxLshalf = []; xmaxLshalf = [];
CYshalf = [];
for i = 1:length(sorted_CYs)
    CYshalf(end+1) = sorted_CYs(i);
    u = sorted_u{i};
    [x, ~] = get_displacements(u);
    def_nodes = get_deformed_nodes(coord, x);
    ymaxLshalf(end+1) = def_nodes(mid,3)/L;
    xmaxLshalf(end+1) = def_nodes(mid,2)/L - 1;
    
end

%% Whole
cd whole

% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));
L = Lx/2;
mid = 0;
for i = 1:length(coord)
    if coord(i,1) < L
        s = coord(i,1)/L;
        coord(i,3) = dh*(3*s-s^3)/2;
    else
        s = 1-(coord(i,1) - L)/L;
        coord(i,3) = dh*(3*s-s^3)/2;
    end
    if s == 1
        mid = i;
    end
end

list = dir;
CY = [];
u = {};
names = [];
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 4
        if item(end-3:end) == ".mat"
            fprintf(item+"\n");
            res = open(item);
            names(end+1) = str2double(item(1:end-4));
            u{end+1} = res.u;
            CY(end+1) = res.lambda;
        end
    end
end


cd ..

[~, indices] = sort(names);
sorted_CYs = CY(indices);
sorted_u = u(indices);
ymaxLswhole = []; xmaxLswhole = [];
CYswhole = [];
for i = 1:length(sorted_CYs)
    CYswhole(end+1) = sorted_CYs(i);
    u = sorted_u{i};
    [x, ~] = get_displacements(u);
    def_nodes = get_deformed_nodes(coord, x);
    ymaxLswhole(end+1) = def_nodes(mid,3)/L;
    xmaxLswhole(end+1) = def_nodes(mid,2)/L - 1;
    
end


%% Displacement plot

f1 = figure(1);
hold on;
grid on;
plot(CYshalf, ymaxLshalf, '-', "color", blue, "linewidth", 2)
plot(CYswhole, ymaxLswhole, '--', "color", blue, "linewidth", 2)
xlabel('$C_Y$')
ylabel('$z/L$')
axis("tight")
xlim([-100, 200]);
f1.Units = "centimeters";
f1.Position = [1 1 21, 9];
figureSaver(f1, "BistableBeams.pdf");