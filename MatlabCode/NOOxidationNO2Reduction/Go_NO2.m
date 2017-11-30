function out1 = Go_NO2(Temp) % Gibbs free energy in eV

if 298 <= Temp && Temp < 1200
    A = 16.10857;
    B = 75.89525;
    C = -54.38740;
    D = 14.30777;
    E = 0.239423;
    F = 26.17464;
    G = 240.5386;
    H = 33.09502;
elseif 1200 <= Temp && Temp <= 6000
    A = 56.82541;
    B = 0.738053;
    C =  -0.144721;
    D = 0.009777;
    E = -5.459911;
    F = 2.846456;
    G = 290.5056;
    H = 33.09502;
else
    error('Temperature out of range!')
end

% Gibbs free energy in J/mol
out1 = Ho(A,B,C,D,E,F,G,H,Temp) - Temp*So(A,B,C,D,E,F,G,H,Temp);
% out1 = DGo(A,B,C,D,E,F,G,H,Temp);

% 1 eV = 96.48530749925793 kJ/mol
out1 = out1/96.48530749925793;

end


