function out1 = So(A,B,C,D,E,F,G,H,Temp) % Entropy in kJ/mol/K

t = Temp/1000;

out1 = A*log(t) + B*t + C*t^2/2 + D*t^3/3 - E/(2*t^2) + G; % in J

out1 = out1/1000; % convert to kJ

end
