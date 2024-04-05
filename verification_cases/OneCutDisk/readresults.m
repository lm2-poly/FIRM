%Reads the data to plot the reconfiguration of a single slit disk
%Compares with Schouveiler & Boudaoud's model and data (Schouveiler & Boudaoud., 2006)
%L. Schouveiler and A. Boudaoud, The rolling up of sheets in a steady flow, J. Fluid Mech. 10.1017/S0022112006000851
%(2006)
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

map1 = linspace(orange(1), 220/256, 100);
map2 = linspace(orange(2), 220/256, 100);
map3 = linspace(orange(3), 220/256, 100);
maporange = [map1', map2', map3'];

%% Parameters
% Geometry
Re = 100;
Ri = 1;
tcut = 0.1;
h = 0.1; % Thickness

% Material properties
B = 1.92;
nu = 0.3;
E = B*12*(1-nu^2)/h^3; % Young's modulus

%Load definition
D = E*h^3/(12*(1-nu^2));
CD = 2.0;
rho = 997*(10^(-3))^3;

%% Read data
coord = table2array(readtable("nodes.csv"));
connec = table2array(readtable("connec.csv"));

bestangle = pi;
for i = 1:length(coord)
    if coord(i,1)^2 + coord(i,2)^2 == Re^2
        if abs(atan2(coord(i,2), coord(i,1)) - pi) < bestangle
            masternode = i;
            bestangle = abs(atan2(coord(i,2), coord(i,1)) - pi);
        end
    end
end

%% No contact no shade

cd data
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
CYs_none = [];
Rs_none = [];
alphas_none = [];
map1 = linspace(orange(1), 220/256, 100);
map2 = linspace(orange(2), 220/256, 100);
map3 = linspace(orange(3), 220/256, 100);
map = [map1', map2', map3'];

for i = 1:length(sorted_CYs)
    x = sorted_u{:,i};
    def_nodes = coord+x(:,1:3);
    u = matrix2vec(x);
    Rs_none(end+1) = get_reconfiguration_number(coord, connec, u, @(a) cos(a)^2, Ri);
    CYs_none(end+1) = sorted_CYs(i);
    alphas_none(end+1) = abs(atan(def_nodes(masternode,3)/sqrt(def_nodes(masternode,1)^2 + def_nodes(masternode,2)^2)));
end

%% Rigid

res = open("rigid.mat");
rigid = res.ans;
x = rigid(:,1); y= rigid(:,2);
ys = @(q,x) q*x.^2;
q = nlinfit(x, y, ys, 2);
CD = q/(0.5*997*pi*(Re/1000)^2);

%% Plot angles
model = table2array(readtable("SchouveilerModel.csv"));
expdata = table2array(readtable("angles.csv"));

f1 = figure(3);
hold on;
grid on;

plot(CYs_none, alphas_none, '-', "color", orange, 'linewidth', linewidth)
plot(model(:,1)*CD, model(:,2), '--k', 'linewidth', linewidth)
scatter(expdata(:,1)*CD, expdata(:,2), "ok")
xlabel('$C_Y$')
ylabel('$\alpha$ [rad]')
set(gca, 'XScale', 'log')
xlim([1E-1, 100]);

f1.Units = "centimeters";
f1.Position = [1 1 9.5 7];
figureSaver(f1, "RollingSheetAngleCurve.pdf");

%% Plot reconfiguration
f1 = figure(2);
hold on;
grid on;

plot(CYs_none, Rs_none, '-', "color", orange, 'linewidth', linewidth)
plot(model(:,1), model(:,3), '--k', 'linewidth', linewidth)
res = open("flex1.mat");
res = res.ans;
scatter(res(:,1)*CD, res(:,2)/CD, "ok")
res = open("flex2.mat");
res = res.ans;
scatter(res(:,1)*CD, res(:,2)/CD, "ok")
xlabel('$C_Y$')
ylabel('$\mathcal{R}$')
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
xlim([1E-1, 100]);

f1.Units = "centimeters";
f1.Position = [1 1 9.5 7];
figureSaver(f1, "RollingSheetReconfigurationCurve.pdf");