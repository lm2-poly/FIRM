%Reads the data to plot the extension of a kirigami sheet
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

% Geometry
cut_thickness = 1E-2;

t = 1;
Ls = 30;
dx = 5/2;
dy = 5;

Nx = 1;
Ny = 11;

L = 2*Ny*dy;
W = Nx*(Ls+2*dx);

% Material properties
E = 4E3; % Young's modulus
nu = 0.3; % Poisson's coefficient

cd IsobeEtAl;

coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

list = dir();
fs = [];
u = [];
index = [];

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
                if coord(j,2) == 0
                    Fsi = Fsi - res.forces(j,2);
                end
            end
            fs(end+1) = Fsi;
        end
    end
end
[~, indices] = sort(index);
sorted_fs = fs(indices);
sorted_u = u(indices);

cd ..;

fs = [];
ys = [];
for i = 1:length(sorted_fs)
    fs(end+1) = sorted_fs(i);
    x = sorted_u{i};
    ys(end+1) = max(abs(x(:,2)));
    def_nodes = get_deformed_nodes(coord, x);
    figure(1);
    display_deformation(connec, def_nodes, x(:,3), true);
    view(90,90)
    colormap([0,0,0]);
    pause(0.01);
end

model = table2array(readtable("IsobeData.xlsx"));

f1 = figure(2);
hold on;
grid off;
plot(ys, fs, '-k')
plot(model(:,1), model(:,2), '--k')
xlabel('$\Delta L$ [mm]')
ylabel('$F$ [N]')
legend("FEM", "Isobe's data", 'location', 'best')
axis("tight")
f1.Units = "centimeters";
f1.Position = [1 1 9.81, 8.01];
figureSaver(f1, "bendRead.pdf");