clear all
clc
close all

global hplanck clight kboltz
addpath('./NOOxidationNO2Reduction/')

Temp = 480; % Kelvin
nFigs = 3; % 3-Figure CE - up to 1NN interactions
runK2NNC2 = true;

% One test call to EnergeticsAndKineticsTlnPr (also sets global parameters 
% and populates the gas species energy vector for the given temperature)
[EnrgGasSpecies,QvibO,ZPEO,...
    ApreO2ads,PERatioO2adsdes,EaO2ads,ProxFacO2ads,...
    ApreOdiff,PERatioOdiff,EaOdiff,ProxFacOdiff,ApreOdiff_sc,...
    ApreNOoxi,PERatioNOoxiNO2red,EaNOoxi,ProxFacNOoxi,ApreNOoxi_sc] = ...
    EnergeticsAndKineticsTlnPr(480,3,-4);

% Adsorption energy of oxygen: the value of -1.200 was hard-coded in the KMC energetics input creator
% DEadsO = 0.6*-1.200;
% DEadsO = 0.8*-1.200;
DEadsO = -1.200;

% Parameter H1: the adsorption energy of oxygen, with vibrational and 
% ZPE corrections, referenced with respect to NO2 and NO ideal gases at 1
% atm pressure
H1 = DEadsO+ZPEO-kboltz*Temp*log(QvibO); % -1.200; 

h = 0.300; % (also hard-coded in the KMC energetics input creator)
% h = 0.3*0.300;
Pres = 1; % bar
yO2 = 0.1;
nn = 6;
logyNO2yNO_equil = 1/2*log(Pres*yO2) + ...
    -(Go_NO2(Temp) - Go_NO(Temp) - 1/2*Go_O2(Temp))/(kboltz*Temp);
% logyNO2yNO_equil_compar = 1/2*log(Pres*yO2) + ...
%     -(EnrgGasSpecies(3) - EnrgGasSpecies(2) - 1/2*EnrgGasSpecies(1))/(kboltz*Temp);
% logyNO2yNOrange = -4; 
% logyNO2yNOrange = logyNO2yNO_equil;
logyNO2yNOrange = -4:0.5:8;
% logyNO2yNOrange = -2:0.02:6;
% logyNO2yNOrange = -4:1:8;

options = optimoptions('fsolve',...
    'TolFun',1e-14,...
    'TolX',1e-14,...
    'Display','iter-detailed');
    ...'Display','off');

%% Test the residuals function

i = 10;
logyNO2yNO = logyNO2yNOrange(i);

beta = 1/(kboltz*Temp);

muOstar = kboltz*Temp*logyNO2yNO + Go_NO2(Temp) - Go_NO(Temp);

% ApproxIdent = 'BP';
% CorrelLHS = {1};
% CorrelRHS = {2};
% nsites = 7;

ApproxIdent = 'BPEC';
CorrelLHS = {1,[1 2]};
CorrelRHS = {2,[2 3]};
nsites = 7;

H0 = 0;

gcguess = [0.25 0.10];


[PFcnBP,gcsolnBP,SThetaBP,SCorrelBP,SCorrelUnoccBP, ...
    SAvgekO2adsBP,SAvgekO2desBP] = ...
    SolveForThetaCorrelRates(beta,Temp,EaO2ads,DEadsO,ProxFacO2ads,0,H1,...
    muOstar,h,nn,ZPEO,QvibO,ApproxIdent,...
    CorrelLHS,CorrelRHS,...
    nsites,gcguess,options);


gc = gcguess;
clc
AllStt = AllStates(7);
AllStt = AllStt(:,7:-1:1);
AllEnergs = HamiltonianBPEC(H0,H1,h,nn,gc(1),gc(2),AllStt);
disp(' ');
for i = 1:2^7
    disp(sprintf('%1.0f  %20.16f',i,AllEnergs(i)));
end
disp(['gc parameter = ' sprintf('%0.16e ',gc)]);

outResid = Residual(ApproxIdent,CorrelLHS,CorrelRHS,nsites,beta,H0,H1,h,nn,muOstar,gc);

disp(['Residual = ' sprintf('%0.16e ',outResid)]);

return

ApproxIdent = 'MF';
CorrelLHS = {};
CorrelRHS = {};
nsites = 1;

H0 = 0;

[PFcnMF,gcsolnMF,SThetaMF,SCorrelMF,SCorrelUnoccMF, ...
    SAvgekO2adsMF,SAvgekO2desMF] = ...
    SolveForThetaCorrelRates(beta,Temp,EaO2ads,DEadsO,ProxFacO2ads,0,H1,...
    muOstar,h,nn,ZPEO,QvibO,ApproxIdent,...
    CorrelLHS,CorrelRHS,...
    nsites,gcguess,options);



ApproxIdent = 'K2NNC1';
nsites = 13; % number of sites in the approximation's cluster
CorrelLHS = {1,1,[1 2],[1 2]}; % left hand side in correlation equations
CorrelRHS = {2,8,[2 3],[2 8]}; % right hand side in correlation equations
gcguess = [0.067 0.22 -0.0084 -0.038]; % guess values for parameters (row vector)
H0 = 0; % Hamiltonian's constant (for renormalisations)

[PFcnK2NNC1,gcsolnK2NNC1,SThetaK2NNC1,SCorrelK2NNC1,SCorrelUnoccK2NNC1, ...
    SAvgekO2adsK2NNC1,SAvgekO2desK2NNC1] = ...
    SolveForThetaCorrelRates(beta,Temp,EaO2ads,DEadsO,ProxFacO2ads,H0,H1,...
    muOstar,h,nn,ZPEO,QvibO,ApproxIdent,...
    CorrelLHS,CorrelRHS,...
    nsites,gcguess,options);

disp(['Solution = ' sprintf('\n ') sprintf('%0.16e \n',gcsolnK2NNC1)]);

g1 = gcguess;

disp(' ');
disp(['g1 parameter = ' sprintf('\n ') sprintf('%0.16e \n',g1)]);

outResid = Residual(ApproxIdent,CorrelLHS,CorrelRHS,nsites,beta,H0,H1,h,nn,muOstar,gcguess);

disp(['Residual = ' sprintf('\n ') sprintf('%0.16e \n',outResid)]);


















