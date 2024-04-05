%Reads the data to plot the reconfiguration of a plate
%Compares with Gosselin's model (Gosselin et al., 2010)
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

linewidth = 2;
green = [75, 204, 40]/256;
blue = [10, 157, 255]/256;
red = [256,0,0]/256;
orange = [255, 201, 40]/256;

%% Load data
% Geometry
Lx = 50; % Length in direction 1 (before applying symmetry, if any)
Ly = 35; % ... in direction 2 (before applying symmetry, if any)
h = 1E-1; % Thickness

% Material properties
E = 68E3; % Young's modulus
nu = 0.; % Poisson's coefficient

%Load definition
B = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 1.225*(10^(-3))^3;

cd data;

%Meshes
% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));
%[coord, connec] = mesh_rectangle_Delaunay(Lx, Ly, 200, 10);

for i = 1:length(coord)
    if coord(i,1) == 50
        tip = i;
    end
end

list = dir;
CYs = [];
u = {};
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 8
        if item(1:8) == "solution"
            fprintf(item+"\n");
            res = open(item);
            u{end+1} = res.displacements;
            startind = find('Y' == item);
            CYs(end+1) = str2num(item(startind+2:end-4));
        end
    end
end

[sorted_CYs, indices] = sort(CYs);
sorted_u = u(indices);
CYs = [];
Rs = [];
map1 = linspace(blue(1), 220/256, 100);
map2 = linspace(blue(2), 220/256, 100);
map3 = linspace(blue(3), 220/256, 100);
map = [map1', map2', map3'];

for i = 1:length(sorted_CYs)
    x = sorted_u{:,i};
    u = matrix2vec(x);
    Rs(end+1) = get_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, 0);
    CYs(end+1) = sorted_CYs(i);
    if round(CYs(end),4) == 100.0
        fprintf("Reconfiguration = %f \t Vertical disp = %f Horizontal disp = %f", Rs(end), x(tip,3), x(tip,1));
    end
end

cd ..;

%% Reconfiguration figure

f1 = figure(1);
hold on;
plot(CYs, Rs, "lineStyle", '-', "color", blue, 'linewidth', linewidth)
model = table2array(readtable("CurveModel.csv"));
plot(model(:,1), model(:,2), '--k', 'linewidth', linewidth)
xlabel('$C_Y$')
ylabel('$\mathcal{R}$')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1E-1, 1E3])
ylim([0.088, 1]);
%axis("tight")
legend("FIRM", "Gosselin et al. (2010)", "location", "best");
f1.Units = "centimeters";
f1.Position = [1 1 11.5, 6.5];
figureSaver(f1, "PlateReconfigurationCurve.pdf");


%% 3D deformations (CY = 1, 10, 100)

ind1 = find(round(CYs,2) == 1);
u1 = sorted_u{ind1};
u1 = matrix2vec(u1);
f2 = figure(2);
[disp, ~] = get_displacements(u1);
def_nodes = get_deformed_nodes(coord, disp);
[nodes_angles, ~] = get_node_angles(def_nodes, connec);
display_deformation(connec, def_nodes, nodes_angles, true);
view(45,30)
colormap((map));

f2.Units = "centimeters";
f2.Position = [1 1 8.7168, 4.8114];
figureSaver(f2, "PlateReconfigurationCY1.png");

ind1 = find(round(CYs,2) == 10);
u1 = sorted_u{ind1};
u1 = matrix2vec(u1);
f3 = figure(3);
[disp, ~] = get_displacements(u1);
def_nodes = get_deformed_nodes(coord, disp);
[nodes_angles, ~] = get_node_angles(def_nodes, connec);
display_deformation(connec, def_nodes, nodes_angles, true);
view(45,30)
colormap((map));

f3.Units = "centimeters";
f3.Position = [1 1 8.7168, 4.8114];
figureSaver(f3, "PlateReconfigurationCY10.png");

ind1 = find(round(CYs,2) == 100);
u1 = sorted_u{ind1};
u1 = matrix2vec(u1);
f4 = figure(4);
[disp, ~] = get_displacements(u1);
def_nodes = get_deformed_nodes(coord, disp);
[nodes_angles, ~] = get_node_angles(def_nodes, connec);
display_deformation(connec, def_nodes, nodes_angles, true);
view(45,30)
colormap((map));

f4.Units = "centimeters";
f4.Position = [1 1 8.7168, 4.8114];
figureSaver(f4, "PlateReconfigurationCY100.png");



%% Deformation figure (CY = 1, 10, 100)
defmodel = table2array(readtable("DeformationModel.csv"));

x1=[]; y1=[];
x2=[]; y2=[];
x3=[]; y3=[];
for i = 1:length(CYs)
    x = sorted_u{:,i};
    u = matrix2vec(x);
    if round(CYs(i),2) == 1
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            x1(end+1) = def_nodes(j,1);
            y1(end+1) = def_nodes(j,3);
        end
    elseif round(CYs(i),2) == 10
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            x2(end+1) = def_nodes(j,1);
            y2(end+1) = def_nodes(j,3);
        end
    elseif round(CYs(i),2) == 100
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            x3(end+1) = def_nodes(j,1);
            y3(end+1) = def_nodes(j,3);
        end
    end
end
[x1s, indices] = sort(x1);
y1s = y1(indices);
[x2s, indices] = sort(x2);
y2s = y2(indices);
[x3s, indices] = sort(x3);
y3s = y3(indices);

f5 = figure(5);
hold on;
plot(x1s/Lx, y1s/Lx, "lineStyle", "-", "color", blue, 'linewidth', linewidth);
plot(defmodel(:,1), defmodel(:,2),'--k', 'linewidth', 3, 'linewidth', linewidth);
plot(x2s/Lx, y2s/Lx, "lineStyle", "-", "color", blue, 'linewidth', linewidth);
plot(defmodel(:,3), defmodel(:,4),'--k', 'linewidth', linewidth);
plot(x3s/Lx, y3s/Lx, "lineStyle", "-", "color", blue, 'linewidth', linewidth);
plot(defmodel(:,5), defmodel(:,6),'--k', 'linewidth', linewidth);
xlabel('x/L');
ylabel('z/L');
axis('tight', 'equal');
f5.Units = "centimeters";
f5.Position = [1 1 8, 8];
figureSaver(f5, "PlateReconfigurationDeformations.pdf");