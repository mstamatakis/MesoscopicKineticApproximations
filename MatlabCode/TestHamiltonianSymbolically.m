clear all 
clc

syms H0 H1 h g1 g2 g3 c1 c2 c3 c4 p1 p2 p3 p4 p5 p6
syms z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 z11 z12 z13 z14 z15 z16 z17 z18 z19

nn = 6;

% epsl = zeros(1,19);
% epsl(13) = 1;
% epsl(8) = 1;

epsl = ones(1,19);
zz = [z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 z11 z12 z13 z14 z15 z16 z17 z18 z19];

% HamiltonianK3NNC2(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,zz) 
HamiltonianK3NNC3(H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,p4,p5,p6,zz) 

Hamiltonian('K3NNC3',H0,H1,h,nn,g1,g2,g3,c1,c2,c3,c4,p1,p2,p3,p4,p5,p6,zz)