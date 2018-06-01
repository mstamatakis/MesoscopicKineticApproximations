clear all
clc

[x,y] = meshgrid(-10:0.5:10,-10:0.5:10);

z = sin(sqrt(x.^2+y.^2))./sqrt(x.^2+y.^2);

surf(x,y,z)