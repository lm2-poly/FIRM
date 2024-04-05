%Reads the data to plot the reconfiguration of a slitted disk
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

map1 = linspace(green(1), 220/256, 100);
map2 = linspace(green(2), 220/256, 100);
map3 = linspace(green(3), 220/256, 100);
mapgreen = [map1', map2', map3'];

map1 = linspace(orange(1), 220/256, 100);
map2 = linspace(orange(2), 220/256, 100);
map3 = linspace(orange(3), 220/256, 100);
maporange = [map1', map2', map3'];

%% Load data
% Geometry
Re = 3.7;
Ri = 0.9;
N = 36;
tcut = 0.001;
h = 1E-2; % Thickness

% Material properties
B = 404E-3;
nu = 0.0;
E = B*12*(1-nu^2)/h^3; % Young's modulus

%Load definition
beta = Re/Ri;
D = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 1.225*(10^(-3))^3;

%Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

%% No shading
cd nonshaded

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
cd ..;

[sorted_CYs, indices] = sort(CYs);
sorted_u = u(indices);
CYs_noshade = [];
Rs_noshade = [];


for i = 1:length(sorted_CYs)
    x = sorted_u{:,i};
    u = matrix2vec(x);
    Rs_noshade(end+1) = get_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, 0);
    CYs_noshade(end+1) = sorted_CYs(i);
end

%% Deformation figure no shade (CY = 1, 10, 100)

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
            if abs(coord(j,2)) < 1E-3 && coord(j,1) > 0
                x1(end+1) = sqrt(def_nodes(j,1)^2 + def_nodes(j,2)^2) - Ri;
                y1(end+1) = def_nodes(j,3);
            end
        end
    elseif round(CYs(i),2) == 10
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            if abs(coord(j,2)) < 1E-3 && coord(j,1) > 0
                x2(end+1) = sqrt(def_nodes(j,1)^2 + def_nodes(j,2)^2) - Ri;
                y2(end+1) = def_nodes(j,3);
            end
        end
    elseif round(CYs(i),2) == 100
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            if abs(coord(j,2)) < 1E-3 && coord(j,1) > 0
                x3(end+1) = sqrt(def_nodes(j,1)^2 + def_nodes(j,2)^2) - Ri;
                y3(end+1) = def_nodes(j,3);
            end
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
plot(x1s/(Re-Ri), y1s/(Re-Ri), "lineStyle", "-", "color", green, 'linewidth', linewidth);
plot(x2s/(Re-Ri), y2s/(Re-Ri), "lineStyle", "-", "color", green, 'linewidth', linewidth);
plot(x3s/(Re-Ri), y3s/(Re-Ri), "lineStyle", "-", "color", green, 'linewidth', linewidth);
xlabel('$\frac{x}{R_e-R_i}$');
ylabel('$\frac{z}{R_e-R_i}$');
axis('tight', 'equal');
f5.Units = "centimeters";
f5.Position = [1 1 8, 8];

%% With shading

cd shaded

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
cd ..;

[sorted_CYs, indices] = sort(CYs);
sorted_u = u(indices);
CYs_shade = [];
Rs_shade = [];


for i = 1:length(sorted_CYs)
    x = sorted_u{:,i};
    u = matrix2vec(x);
    Rs_shade(end+1) = get_shaded_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, 0);
    CYs_shade(end+1) = sorted_CYs(i);
end

%% Deformation figure (CY = 1, 10, 100)
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
            if abs(coord(j,2)) < 1E-3 && coord(j,1) > 0
                x1(end+1) = sqrt(def_nodes(j,1)^2 + def_nodes(j,2)^2) - Ri;
                y1(end+1) = def_nodes(j,3);
            end
        end
    elseif round(CYs(i),2) == 10
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            if abs(coord(j,2)) < 1E-3 && coord(j,1) > 0
                x2(end+1) = sqrt(def_nodes(j,1)^2 + def_nodes(j,2)^2) - Ri;
                y2(end+1) = def_nodes(j,3);
            end
        end
    elseif round(CYs(i),2) == 100
        [disp, ~] = get_displacements(u);
        def_nodes = get_deformed_nodes(coord, disp);
        for j = 1:length(disp)
            if abs(coord(j,2)) < 1E-3 && coord(j,1) > 0
                x3(end+1) = sqrt(def_nodes(j,1)^2 + def_nodes(j,2)^2) - Ri;
                y3(end+1) = def_nodes(j,3);
            end
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
plot(x1s/(Re-Ri), y1s/(Re-Ri), "lineStyle", "-.", "color", green, 'linewidth', linewidth);
plot(x2s/(Re-Ri), y2s/(Re-Ri), "lineStyle", "-.", "color", green, 'linewidth', linewidth);
plot(x3s/(Re-Ri), y3s/(Re-Ri), "lineStyle", "-.", "color", green, 'linewidth', linewidth);
xlabel('$\frac{r-R_i}{R_e-R_i}$');
ylabel('$\frac{z}{R_e-R_i}$');
axis('tight', 'equal');
f5.Units = "centimeters";
f5.Position = [1 1 10 7];
figureSaver(f5, "DiskReconfigurationDeformations.pdf");

%% Plotting
model = table2array(readtable("modelDisk.csv"));

f1 = figure(2);
hold on;
grid on;
plot(CYs_noshade, Rs_noshade, '-', "color", green, 'linewidth', linewidth)
plot(CYs_shade, Rs_shade, '-.', "color", green, 'linewidth', linewidth)
plot(model(:,1), model(:,2), '--k')
xlabel('$C_Y$')
ylabel('$\mathcal{R}$')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1E-1, 1E3]);
legend("FIRM - no shade", "FIRM - with shade", "Gosselin et al. (2010)", 'location', 'best')
f1.Units = "centimeters";
f1.Position = [1 1 10 7];
figureSaver(f1, "SlittedDiskReconfigurationCurve.pdf");