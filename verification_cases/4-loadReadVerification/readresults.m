%Reads the data to plot the displacement of a plate
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

%% Load data
% Geometry
Lx = 1; %Length of the plate
Ly = 10; % Width of the plate
h = 10E-1; % Thickness

% Material properties
E = 1.2E3; % Young's modulus
nu = 0.; % Poisson's coefficient

%Defining the applied loads
I = Ly*h^3/12;
q0 = E*I/Lx^3;

% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

loadnodes = [];

for i = 1:length(coord)
    if coord(i,2) == Ly
        tip = i;
        loadnodes(end+1) = i;
    end
end

cd data;

list = dir;
Fs = {};
u = {};
names = [];
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 8
        if item(1:8) == "solution"
            fprintf(item+"\n");
            names(end+1) = str2num(item(10:end-4));
            res = open(item);
            u{end+1} = res.displacements;
            Fs{end+1} = res.forces;
        end
    end
end
cd ..;

[~, indices] = sort(names);
sorted_F = Fs(indices);
sorted_u = u(indices);
ws = [];
us = [];

Fs = [];
map1 = linspace(blue(1), 220/256, 100);
map2 = linspace(blue(2), 220/256, 100);
map3 = linspace(blue(3), 220/256, 100);
map = [map1', map2', map3'];

f2 = figure(2);
hold on;
axis('equal');
axis off
set(gcf,'color','white')
f2.Units = "centimeters";
f2.Position = [1 1 8.81, 7.01];
trisurf(connec, coord(:, 1), coord(:, 2), coord(:, 3), zeros(length(coord), 1), 'EdgeColor', 'None');

for i = 1:length(sorted_F)
    if mod(i,200) == 0
        disp = sorted_u{i};
        u1 = matrix2vec(disp);
        [disp, disp_norms] = get_displacements(u1);
        def_nodes = get_deformed_nodes(coord, disp);
        trisurf(connec, def_nodes(:, 1), def_nodes(:, 2), def_nodes(:, 3), disp_norms, 'EdgeColor', 'None');
    end
    
    x = sorted_u{:,i};
    u = matrix2vec(x);
    us(end+1) = -x(tip,2);
    ws(end+1) = x(tip,3);
    
    F = sorted_F{i};
    loads = F(loadnodes, 3);
    Fs(end+1) = sum(loads);   
end
view(-45, 30);
figureSaver(f2, "overlappeddeformation.png");



%%
%Comparison models
modelw = "modelw.csv";
modelw = table2array(readtable(modelw));

modelu = "modelu.csv";
modelu = table2array(readtable(modelu));

%Plotting the tip deflections (to compare with theory)
f1 = figure();
hold on;
plot(ws/Ly, Fs*Ly^3/(E*I), '-k');
scatter(modelw(:,1), modelw(:,2), 'ok');
plot(us/Ly, Fs*Ly^3/(E*I), '--k');
scatter(modelu(:,1), modelu(:,2), 'sk');
xlabel('$\frac{\Delta}{L}$', 'interpreter', 'latex');
ylabel('$\frac{FL^3}{EI}$', 'interpreter', 'latex');
legend('$w_{tip, FEM}$', '$w_{tip, ref}$', '$-u_{tip, FEM}$', '$-u_{tip, ref}$', 'interpreter', 'latex', 'location','best');
f1.Units = "centimeters";
f1.Position = [1 1 9.81, 8.01];
figureSaver(f1, "dispControl.pdf");