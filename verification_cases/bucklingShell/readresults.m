%Reads the data to plot the buckling shell
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
R = 100E-3;
h = R/8.5; % Radius of the plate
theta = deg2rad(60);

Rh = R*sin(theta);
H = R - Rh/tan(theta);

% Material properties
E = 3102.75; % Young's modulus
nu = 0.3; % Poisson's coefficient

G = E/(2*(1+nu));

%Maximum load
p = 10;

% Mesh
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));
%mesh_size = 7.5E-3;
%[coord, connec] = mesh_disk(Rh, h, mesh_size, true);

for i = 1:length(coord)
    coord(i,3) = real(sqrt(R^2 - coord(i,1)^2 - coord(i,2)^2)) - R*cos(theta);
end

for i = 1:length(coord)
    if (coord(i,1)^2+coord(i,1)^2) < 0.01*Rh^2
        mid = i;
        break;
    end
end

%Comparison models
model = "GorissenEtAl.csv";
model = table2array(readtable(model));

cd arclength;
list = dir;
index = [];
ps = [];
u = {};
names = [];
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 5
        if item(end-3:end) == ".mat"
            fprintf(item+"\n");
            names(end+1) = str2num(item(1:end-4));
            res = open(item);
            ps(end+1) = p*res.lambda;
            u{end+1} = res.u;
        end
    end
end
cd ..;
[~, indices] = sort(names);
sorted_p = ps(indices);
sorted_u = u(indices);

dV = [0];

z = 2*R;

V0 = compute_volume(connec, coord, [0, 0, z]);
P = [0];
for i = 1:length(sorted_u)    
    u = sorted_u{i};
    [x, disp_norms] = get_displacements(u);
    def_nodes = get_deformed_nodes(coord, x);
    display_deformation(connec, def_nodes, disp_norms, true);
    pause(0.01);
    dV(end+1) = compute_volume(connec, def_nodes, [0, 0, z]) - V0;
    P(end+1) = sorted_p(i);
end

f1 = figure(2);
hold on;
xlabel("$\Delta V/R^3$");
ylabel("$P/\mu$");
scatter(model(:,1), model(:,2));
plot(dV/R^3, P/G, "-k");
text(0, 15E-3, "Note: Gorissen et al. use a different constitutive law for their materials, the agreement should therefore only be qualitative here");
legend("Gorissen et al., 2020", "FEM", "location", "best");
% f1.Units = "centimeters";
% f1.Position = [1 1 9, 5];
figureSaver(f1, "bucklingShell.pdf");