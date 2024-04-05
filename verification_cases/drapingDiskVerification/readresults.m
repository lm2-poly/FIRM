%Reads the data to plot the draping disk
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

clear all; close all; clc;

Re = 70E-3;
Ri = 6E-3;
rho = 1.225;
beta = Re/Ri;
h = 75E-6; % Thickness
E = 4E9;
nu = 0.3;
B = E*h^3/(12*(1-nu^2));
CD = 1.0;

%% Mode C

cd C

coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

list = dir;
names = [];
CYs = [];
u = [];
Fext = [];
for i = 1:length(list)
    item =list(i).name;
    if length(list(i).name(:)) > 3
        if list(i).name(end-3:end) == ".mat"
            fprintf(item+"\n");
            res = open(item);
            names(end+1) = str2num(item(1:end-4));
            u(:,end+1) = res.u;
            Fext(:,end+1) = res.Fext;
            CYs(end+1) = res.lambda;
        end
    end
end

[~, indices] = sort(names);
sorted_u = u(:,indices);
sorted_Fext = Fext(:,indices);
sorted_CYs = CYs(indices);
Rs = [];
for i = 1:length(sorted_CYs)
    Cy = sorted_CYs(i);
    Rs(end+1) = get_reconfiguration_number(coord, connec, sorted_u(:,i), @(a) cos(a)^2, Ri);
end

fig = figure(11);
hold on;
grid on;
plot(sorted_CYs, Rs)
xlabel('$C_D C_Y$')
ylabel('$\mathcal{R}$')
ylim([0.1, 1.0])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis("tight")

cd ..

%% 3F

cd 3F

coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

list = dir;
names = [];
CYs = [];
u = [];
Fext = [];
for i = 1:length(list)
    item =list(i).name;
    if length(list(i).name(:)) > 3
        if list(i).name(end-3:end) == ".mat"
            fprintf(item+"\n");
            res = open(item);
            names(end+1) = str2num(item(1:end-4));
            u(:,end+1) = res.u;
            Fext(:,end+1) = res.Fext;
            CYs(end+1) = res.lambda;
        end
    end
end

[~, indices] = sort(names);
sorted_u = u(:,indices);
sorted_CYs = CYs(indices);
sorted_Fext = Fext(:,indices);
Rs = [];
for i = 1:length(sorted_CYs)
    Cy = sorted_CYs(i);
    Rs(end+1) = get_reconfiguration_number(coord, connec, sorted_u(:,i), @(a) cos(a)^2, Ri);
end

fig = figure(11);
hold on;
grid on;
plot(sorted_CYs, Rs)
xlabel('$C_D C_Y$')
ylabel('$\mathcal{R}$')
ylim([0.1, 1.0])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis("tight")

cd ..

%% 4F

cd 4F

coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

list = dir;
names = [];
CYs = [];
u = [];
Fext = [];
for i = 1:length(list)
    item =list(i).name;
    if length(list(i).name(:)) > 3
        if list(i).name(end-3:end) == ".mat"
            fprintf(item+"\n");
            res = open(item);
            names(end+1) = str2num(item(1:end-4));
            u(:,end+1) = res.u;
            Fext(:,end+1) = res.Fext;
            CYs(end+1) = res.lambda;
        end
    end
end

[~, indices] = sort(names);
sorted_u = u(:,indices);
sorted_CYs = CYs(indices);
sorted_Fext = Fext(:,indices);
Rs = [];
for i = 1:length(sorted_CYs)
    Cy = sorted_CYs(i);
    Rs(end+1) = get_reconfiguration_number(coord, connec, sorted_u(:,i), @(a) cos(a)^2, Ri);
end

fig = figure(11);
hold on;
grid on;
plot(sorted_CYs, Rs)
xlabel('$C_D C_Y$')
ylabel('$\mathcal{R}$')
ylim([0.1, 1.0])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis("tight")

cd ..

%% 5F

cd 5F

coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

list = dir;
names = [];
CYs = [];
u = [];
Fext = [];
for i = 1:length(list)
    item =list(i).name;
    if length(list(i).name(:)) > 3
        if list(i).name(end-3:end) == ".mat"
            fprintf(item+"\n");
            res = open(item);
            names(end+1) = str2num(item(1:end-4));
            u(:,end+1) = res.u;
            Fext(:,end+1) = res.Fext;
            CYs(end+1) = res.lambda;
        end
    end
end

[~, indices] = sort(names);
sorted_u = u(:,indices);
sorted_CYs = CYs(indices);
sorted_Fext = Fext(:,indices);
Rs = [];
for i = 1:length(sorted_CYs)
    Cy = sorted_CYs(i);
    Rs(end+1) = get_reconfiguration_number(coord, connec, sorted_u(:,i), @(a) cos(a)^2, Ri);
end

fig = figure(11);
hold on;
grid on;
plot(sorted_CYs, Rs)
xlabel('$C_D C_Y$')
ylabel('$\mathcal{R}$')
ylim([0.1, 1.0])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
axis("tight")

legend("C", "3F", "4F", "5F")

cd ..