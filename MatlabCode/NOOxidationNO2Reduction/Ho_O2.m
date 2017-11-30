function out1 = Ho_O2(Temp) % Enthalpy in eV

if 100 <= Temp && Temp < 700
    A = 31.32234;
    B = -20.23531;
    C = 57.86644;
    D = -36.50624;
    E = -0.007374;
    F = -8.903471;
    G = 246.7945;
    H = 0.0;
elseif 700 <= Temp && Temp < 2000
    A = 30.03235;	
    B = 8.772972;	
    C = -3.988133;
    D = 0.788313;
    E = -0.741599;
    F = -11.32468;	
    G = 236.1663;	
    H = 0.0;
elseif 2000 <= Temp && Temp <= 6000
    A = 20.91111;
    B = 10.72071;
    C = -2.020498;
    D = 0.146449;
    E = 9.245722;
    F = 5.337651;
    G = 237.6185;
    H = 0.0;
else
    error('Temperature out of range!')
end

% Enthalpy in J/mol
out1 = Ho(A,B,C,D,E,F,G,H,Temp);

% 1 eV = 96.48530749925793 kJ/mol
out1 = out1/96.48530749925793;

end
