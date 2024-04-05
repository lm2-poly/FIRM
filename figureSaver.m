%General template to save figures
%Author: Danick Lamoureux

%Copyright (C) 2024 Danick Lamoureux

%This program is free software: you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation, either version 3 of the License, or any later version.

%This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

%You should have received a copy of the GNU General Public License along with 
% this program. If not, see https://www.gnu.org/licenses/.

function figureSaver(fig, filename)
    grid off;
    set(gca,'FontSize',12)
    set(fig, 'color', 'white');    
    exportgraphics(gcf, filename,...
        'ContentType','vector',...
        'BackgroundColor','none', ...
        'Resolution',600)
end