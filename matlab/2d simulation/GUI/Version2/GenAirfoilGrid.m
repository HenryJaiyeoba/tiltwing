function [xrow, yrow] = GenAirfoilGrid()
%to create 3D images of conceptual UAV

%SHAPE PARAMETERS:
%wing

%create airfoil
ang = pi*(-1:0.02:1);
rad = exp((1*ang).^2);

x0vec = rad .* cos(ang);
y0vec = rad .* sin(ang);
B0max = max(abs(x0vec));
xrow = x0vec / B0max + 0.5;
yrow = y0vec / B0max;








