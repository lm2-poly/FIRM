%Initial file to run to use the current implementation.
%Adds all the required paths to the current one and makes the default
%interpreter Latex. This code is based on Caselli and Bisegna (2013) 
% (to add to path the adequate elements) 
% as well as Max Bartholdt's answer (to change interpreter to latex)
%(https://www.mathworks.com/matlabcentral/answers/183311-setting-default-interpreter-to-latex)

%Author: Danick Lamoureux
%Creation date: 2023-05-20

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function importAll()
    clear all; close all; clc;
    
    %Adding the different directories to the current path
    path(".\", path);
    path('.\core_element\shell', path);
    path('.\corotational_filter\common', path);
    path('.\corotational_filter\shell', path);
    path('.\constitutive_law\2D', path);
    path('.\mesh', path);
    path('.\export2vtk', path);
    path('.\FEsolver', path);
    path('.\boundary_conditions', path);
    path('.\nonlinear', path);
    path('.\analysis', path);
    path('.\loads', path);
    
    % This script changes all interpreters from text to latex. (See Max
    % Bartholdt's answer)
    list_factory = fieldnames(get(groot,'factory'));
    index_interpreter = find(contains(list_factory,'Interpreter'));
    for i = 1:length(index_interpreter)
        default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
        set(groot, default_name,'latex');
    end
end