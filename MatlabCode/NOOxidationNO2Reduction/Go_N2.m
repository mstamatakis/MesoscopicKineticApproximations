function out1 = Go_N2(Temp) % Gibbs free energy in eV

if 100 <= Temp && Temp < 500
    A = 28.98641;
    B = 1.853978;
    C = -9.647459;
    D = 16.63537;
    E = 0.000117;
    F = -8.671914;
    G = 226.4168;
    H = 0.0;
elseif 500 <= Temp && Temp < 2000
    A = 19.50583;
    B = 19.88705;
    C = -8.598535;
    D = 1.369784;
    E = 0.527601;
    F = -4.935202;
    G = 212.3900;
    H = 0.0;
elseif 2000 <= Temp && Temp <= 6000
    A = 35.51872;
    B = 1.128728;
    C = -0.196103;
    D = 0.014662;
    E = -4.553760;
    F = -18.97091;
    G = 224.9810;
    H = 0.0;
else
    error('Temperature out of range!')
end

% Gibbs free energy in J/mol
out1 = Ho(A,B,C,D,E,F,G,H,Temp) - Temp*So(A,B,C,D,E,F,G,H,Temp);
% out1 = DGo(A,B,C,D,E,F,G,H,Temp);

% 1 eV = 96.48530749925793 kJ/mol
out1 = out1/96.48530749925793;

end
