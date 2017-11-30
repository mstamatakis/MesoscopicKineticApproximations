function out1 = Go_NO(Temp) % Gibbs free energy in eV

if 298 <= Temp && Temp < 1200
    A = 23.83491;
    B = 12.58878;
    C = -1.139011;
    D = -1.497459;
    E = 0.214194;
    F = 83.35783;
    G = 237.1219;
    H = 90.29114;
elseif 1200 <= Temp && Temp <= 6000
    A = 35.99169;
    B = 0.957170;
    C = -0.148032;
    D = 0.009974;
    E = -3.004088;
    F = 73.10787;
    G = 246.1619;
    H = 90.29114;
else
    error('Temperature out of range!')
end

% Gibbs free energy in J/mol
out1 = Ho(A,B,C,D,E,F,G,H,Temp) - Temp*So(A,B,C,D,E,F,G,H,Temp);
% out1 = DGo(A,B,C,D,E,F,G,H,Temp);

% 1 eV = 96.48530749925793 kJ/mol
out1 = out1/96.48530749925793;

end
