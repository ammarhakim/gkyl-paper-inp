%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grant Johnson
% 8/12/2024
% Description: Computes and plots geodesics for a sphere a
% and a hyperboliod with refelcting BC's.
%%% Note: Part code was generated via CHATGPT, GJ (8/12/2024) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

% I.C.s for geodescis
ic_sphere_1 = [1, pi/2, -pi/2, 0.8, 1.0, 5];
ic_sphere_2 = [1, pi/2, -pi/2, 1.5, 1.0, 5.65];
ic_hyperbolic_1 = [-pi/2, -0.1, 1.5, 0.2, 6];
ic_hyperbolic_2 = [-pi/2, -0.5, 1.5, 0.9, 10];

% Figure setup
figure;
set(gcf, 'Units', 'inches', 'Position', [1, 1, 12, 3])

% Sphere ICs 1
subplot(1,4,1);
hold on;
title('Sphere - Initial Conditions 1', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Arial');
geodesic_sphere(ic_sphere_1(1),ic_sphere_1(2),ic_sphere_1(3),ic_sphere_1(4),ic_sphere_1(5),ic_sphere_1(6));
text(-1.5, -1.5, 1.5, '(a)', 'FontSize', 12);
hold off;

% Sphere ICs 2
subplot(1,4,2);
hold on;
title('Sphere - Initial Conditions 2', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Arial');
geodesic_sphere(ic_sphere_2(1),ic_sphere_2(2),ic_sphere_2(3),ic_sphere_2(4),ic_sphere_2(5),ic_sphere_2(6));
text(-1.5, -1.5, 1.5, '(b)', 'FontSize', 12); 
hold off;

% Hyperbolic ICs 1
subplot(1,4,3);
hold on;
title('Hyperbolic - Initial Conditions 1', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Arial');
geodesic_hyperbolic(ic_hyperbolic_1(1),ic_hyperbolic_1(2),ic_hyperbolic_1(3),ic_hyperbolic_1(4),ic_hyperbolic_1(5));
text(-1.5^1.85, -1.5^1.85, 1.5^2, '(c)', 'FontSize', 12);
hold off;

% Hyperbolic ICs 2
subplot(1,4,4);
hold on;
title('Hyperbolic - Initial Conditions 2', 'FontSize', 12, 'FontWeight', 'normal', 'FontName', 'Arial');
geodesic_hyperbolic(ic_hyperbolic_2(1),ic_hyperbolic_2(2),ic_hyperbolic_2(3),ic_hyperbolic_2(4),ic_hyperbolic_2(5));
text(-1.5^1.85, -1.5^1.85, 1.5^2, '(d)', 'FontSize', 12);
hold off;

% Save the figure as a png
exportgraphics(gcf, 'fig1_geodesic_plots.png', 'Resolution', 300);
