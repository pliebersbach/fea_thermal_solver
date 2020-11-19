clear
close all
clc

t = load('solution.mat');
x = load('points.mat');
mesh = load('mesh.mat');

%contourf(x(:,1),x(:,2),t)
mesh = mesh+1;
figure('units','normalized','position', [0 0.33 .3 .3]);
patch('Faces',mesh,'Vertices',x,'FaceVertexCData',t,'FaceColor','interp');
colormap('jet')
colorbar
title('Temperature Distribution (^{\circ}K)')
set(gca, 'fontsize', 16);
saveas(gcf, 'Temp.png')