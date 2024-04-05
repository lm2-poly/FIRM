%Reads the data to plot the dead deformation of a plate
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
Lx = 10; %Length of the plate
Ly = 1; % Width of the plate
h = 10E-1; % Thickness

% Material properties
E = 1.2E6; % Young's modulus
nu = 0.; % Poisson's coefficient

%Defining the applied loads
I = Ly*h^3/12;
q0 = E*I/Lx^3;

% Mesh
%Reuse existing meshes
coord = table2array(readtable("data/nodes.csv"));
connec = table2array(readtable("data/connec.csv"));

%What node is the tip of the plate
for i = 1:length(coord)
    if coord(i,1) == Lx
        tip = i;
    end
end

cd data;
list = dir;
index = [];
Fs = [];
u = {};
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 8
        if item(1:8) == "solution"
            fprintf(item+"\n");
            index(end+1) = str2num(item(10:end-4));
            res = open(item);
            u{end+1} = res.displacements;
            Fsi = 0;
            for j = 1:length(coord(:,1))
                if coord(j,1) == 0
                    Fsi = Fsi - res.forces(j,3);
                end
            end
            %Applied pressure
            Fs(end+1) = Fsi/(Lx*Ly);
        end
    end
end
cd ..;
[~, indices] = sort(index);
sorted_Fs = Fs(indices);
sorted_u = u(indices);

w = [0];
u = [0];
F = [0];
for i = 1:length(sorted_Fs)    
    %Tip deflections
    F(end+1) = sorted_Fs(i);
    w(end+1) = sorted_u{i}(tip, 3);
    u(end+1) = -sorted_u{i}(tip, 1);
    fprintf("Load = %e\t Deflection = %e\n", F(end), w(end));
end

%% Deformation figure

%Comparison models
modelw = "modelw.csv";
modelw = table2array(readtable(modelw));

modelu = "modelu.csv";
modelu = table2array(readtable(modelu));

%Plotting the tip deflections (to compare with the article)
f1 = figure();
hold on;
plot(w/Lx, F*Lx^3/(E*I), 'blue');
scatter(modelw(:,1), modelw(:,2), 'blue');
plot(u/Lx,F*Lx^3/(E*I), 'red');
scatter(modelu(:,1), modelu(:,2), 'red');
xlabel('$\frac{\Delta}{L}$', 'interpreter', 'latex');
ylabel('$\frac{FL^3}{EI}$', 'interpreter', 'latex');
legend('$w_{tip, FEM}$', '$w_{tip, ref}$', '$-u_{tip, FEM}$', '$-u_{tip, ref}$', 'interpreter', 'latex', 'location','best');
xlim([0,1]);
ylim([0,40]);
f1.Units = "centimeters";
f1.Position = [1 1 9, 5];
figureSaver(f1, "deadDistributionCurve.pdf");