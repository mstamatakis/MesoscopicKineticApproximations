function out1 = Go(A,B,C,D,E,F,G,H,Temp) % Gibbs energy in kJ/mol

t = Temp/1000;

out1 = A*(t-t*log(t)) - B*t^2/2 - C*t^3/6 - D*t^4/12 - E/(2*t) - G*t + F;

end
