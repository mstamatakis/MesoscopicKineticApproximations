function out1 = Ho(A,B,C,D,E,F,G,H,Temp) % Enthalpy in kJ/mol

t = Temp/1000;

out1 = A*t + B*t^2/2 + C*t^3/3 + D*t^4/4 - E/t + F;%- H;

end
