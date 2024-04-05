%Reads the data to plot the snapping of a hinged roof
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
L = 254;
R = 2540; % Radius of the plate
theta = 0.1;

% Material properties
E = 3102.75; % Young's modulus
nu = 0.3; % Poisson's coefficient

% Mesh
%Reuse existing meshes
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

%What node is the tip of the plate
for i = 1:length(coord)
    if coord(i,1) == 0 && coord(i,2) == L && ismember(i, connec(:))
        tip = i;
    end
end

loadnodes = [];

for i = 1:length(coord)
    if coord(i,1) == 0 && coord(i,2) == L
        for j = 1:length(connec)
            if ismember(i, connec(j,:))
                loadnodes(end+1) = i;
            end
        end
    end
end

%% h = 12.7

h = 12.7; % Thickness

cd h12.7;
list = dir;
index = [];
Fs = [];
u = {};
names = [];
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 5
        if item(end-3:end) == ".mat"
            fprintf(item+"\n");
            names(end+1) = str2num(item(1:end-4));
            res = open(item);
            Fs{end+1} = res.Fext*res.lambda;
            u{end+1} = res.u;
        end
    end
end
cd ..;
[~, indices] = sort(names);
sorted_F = Fs(indices);
sorted_u = u(indices);

w = [0];
P = [0];
for i = 1:length(sorted_u)    
    %Tip deflections
    F = sorted_F{i};
    Fx = get_displacements(F);
    loads = Fx(loadnodes, 3);
    P(end+1) = -4*sum(loads);
    u = sorted_u{i};
    x = get_displacements(u);
    w(end+1) = -x(tip, 3);
    fprintf("Load = %e\t Deflection = %e\n", P(end), w(end));
end

w1 = w; P1 = P;

%% h = 6.35
cd h6.35;
list = dir;
index = [];
Fs = [];
u = {};
names = [];
for i = 1:length(list)
    item = list(i).name;
    if length(item(:)) > 5
        if item(end-3:end) == ".mat"
            fprintf(item+"\n");
            names(end+1) = str2num(item(1:end-4));
            res = open(item);
            Fs{end+1} = res.Fext*res.lambda;
            u{end+1} = res.u;
        end
    end
end
cd ..;
[~, indices] = sort(names);
sorted_F = Fs(indices);
sorted_u = u(indices);

w = [0];
P = [0];
for i = 1:length(sorted_u)    
    %Tip deflections
    F = sorted_F{i};
    Fx = get_displacements(F);
    loads = Fx(loadnodes, 3);
    P(end+1) = -4*sum(loads);
    u = sorted_u{i};
    x = get_displacements(u);
    w(end+1) = -x(tip, 3);
    fprintf("Load = %e\t Deflection = %e\n", P(end), w(end));
end

w2 = w; P2 = P;

%% Deformation figure

allw = [w1, w2];
allP = [P1, P2];

%Comparison models
model1 = "h12.7.csv";
model1 = table2array(readtable(model1));
model2 = "h6.35.csv";
model2 = table2array(readtable(model2));

%Plotting the tip deflections (to compare with the article)
f1 = figure();
hold on;
xlabel("Displacements [mm]");
ylabel("Load [N]");
xlim([0, max(allw)]);
ylim([min(allP), max(allP)]);
plot(w1, P1, '-k');
scatter(model1(:,1), model1(:,2), 'ok');
plot(w2, P2, '--k');
scatter(model2(:,1), model2(:,2), 'sk');
legend('FEM: $h = 12.7$mm', 'Ref: $h = 12.7$mm', 'FEM: $h = 6.35$mm','Ref: $h = 6.35$mm', 'interpreter', 'latex', 'location','best');
f1.Units = "centimeters";
f1.Position = [1 1 9.81, 8.01];
figureSaver(f1, "HingedRoof.pdf");