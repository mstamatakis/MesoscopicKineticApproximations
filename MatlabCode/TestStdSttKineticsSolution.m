clc

if 1
clear all

global hplanck clight kboltz

addpath('./NOOxidationNO2Reduction/')

Temp = 680; % Kelvin
nFigs = 3; % 3-Figure CE - up to 1NN interactions
compareMC = true;

% One test call to EnergeticsAndKineticsTlnPr (also sets global parameters 
% and populates the gas species energy vector for the given temperature)
[EnrgGasSpecies,QvibO,ZPEO,...
    ApreO2ads,PERatioO2adsdes,EaO2ads,ProxFacO2ads,...
    ApreOdiff,PERatioOdiff,EaOdiff,ProxFacOdiff,ApreOdiff_sc,...
    ApreNOoxi,PERatioNOoxiNO2red,EaNOoxi,ProxFacNOoxi,ApreNOoxi_sc] = ...
    EnergeticsAndKineticsTlnPr(Temp,3,-4);

beta = 1/(kboltz*Temp);

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
% logyNO2yNOrange = -4; %-1.65; 
% logyNO2yNOrange = logyNO2yNO_equil;
% logyNO2yNOrange = -4:0.1:8;
% logyNO2yNOrange = -2:0.02:6;
% logyNO2yNOrange = -4:0.5:8;
% logyNO2yNOrange = -4:1.0:8;
logyNO2yNOrange = -4;

%% Solver options

options = optimoptions('fsolve',...
    'TolFun',1e-14,...
    'TolX',1e-14,...
    'Display','iter-detailed');
    ...'Display','off');

%% Define all the approximations
% For each approximation the user has to supply a Hamiltonian function
% and a set of commands mainly defining the equations used to solve for the
% effective field parameters. The name of the approximation specifies the
% Hamiltonian function called; this has to be defined in the function
% contained in the file "Hamiltonian.m".

kAppr = 0;
ApproxName = {};

% The first approximation is the mean-field and is treated in a special way
% since correlations are neglected
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'MF'; % name of the approximation
nsites{kAppr} = 1; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {}; % left hand side in correlation equations
CorrelRHS{kAppr} = {}; % right hand side in correlation equations
gcguess{kAppr} = 0.3; % guess values for parameters (row vector)
H0(kAppr) = 0; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 0 0]; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'o'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the BP approximation without edge interactions
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BP'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2}; % right hand side in correlation equations
gcguess{kAppr} = 0.2332; % guess values for parameters (row vector)
H0(kAppr) = 0.3016; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [250 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'v'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '--'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the BPE approximation with edge interactions
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BPE'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2}; % right hand side in correlation equations
gcguess{kAppr} = 0.1803; % guess values for parameters (row vector)
H0(kAppr) = 0.3147; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [200 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = '>'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '--'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the BPEC approximation with edge interactions and an equation for
% the correlation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BPEC'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,[1 2]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,[2 3]}; % right hand side in correlation equations
gcguess{kAppr} = [0.1805   -0.0650]; % guess values for parameters (row vector)
H0(kAppr) = 0.3147; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [180 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = '^'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '--'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the K2NNC1 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC1'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8]}; % right hand side in correlation equations
gcguess{kAppr} = [0.0669    0.2188   -0.0084   -0.0380]; % guess values for parameters (row vector)
H0(kAppr) = 0.6568; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 's'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the K2NNC2 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC2'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2],[1 8],[1 8]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8],[2 4],[8 9]}; % right hand side in correlation equations
gcguess{kAppr} = [0.0866    0.3923    0.0723    0.0012   -0.0442   -0.1952]; % guess values for parameters (row vector)
H0(kAppr) = 0.7100; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 100 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'd'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the K2NNC3 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC3'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2],[1 8],[1 8],[2 5]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8],[2 4],[8 9],[4 8]}; % right hand side in correlation equations
gcguess{kAppr} = [0.0956    0.3970    0.0604   -0.0036   -0.0474   -0.1970   -0.0063]; % guess values for parameters (row vector)
H0(kAppr) = 0.6927; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 200 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'p'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size

% Next the K3NNC2 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K3NNC2'; % name of the approximation
nsites{kAppr} = 19; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14]}; % right hand side in correlation equations
gcguess{kAppr} = [0.1165    0.2745    0.3094    0.0524   -0.0310    0.1431   -0.0441    0.0694   -0.1377   -0.1548]; % guess values for parameters (row vector)
H0(kAppr) = 1.0713; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [200 200 0]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'h'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size

% Next the K3NNC3 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K3NNC3'; % name of the approximation
nsites{kAppr} = 19; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8],[1 14],[1 14],[ 1 14]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14],[2  5],[2  9],[14 15]}; % right hand side in correlation equations
gcguess{kAppr} = [0.1730    0.2695    0.1752   -0.0576   -0.0223   -0.0027    0.0029   -0.0080   -0.1360   -0.0806   -0.0579   -0.0771    0.1049]; % guess values for parameters (row vector)
H0(kAppr) = 1.0915; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [255 200 0]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'o'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size

%% Initialise arrays for data on coverage, correlations, and rates

% The following are intended to be indexed by 
% i -> logyNO2yNOrange and 
% j -> approximation type
Theta = [];
CorrelOccup = [];
CorrelUnocc = [];
AvgekO2ads = [];
AvgekO2des = [];

%% Solving
mOstarRange = [];

i = 1;
k = 1; % Approximation

logyNO2yNO = logyNO2yNOrange(i);

[EnrgGasSpecies,QvibO,ZPEO,...
    ApreO2ads,PERatioO2adsdes,EaO2ads,ProxFacO2ads,...
    ApreOdiff,PERatioOdiff,EaOdiff,ProxFacOdiff,ApreOdiff_sc,...
    ApreNOoxi,PERatioNOoxiNO2red,EaNOoxi,ProxFacNOoxi,ApreNOoxi_sc] = ...
    EnergeticsAndKineticsTlnPr(Temp,nFigs,logyNO2yNO);
% Undo any stiffness scalings:
% ApreNOoxi = ApreNOoxi/ApreNOoxi_sc;
% ApreNOoxi_sc = 1;
% ApreOdiff = ApreOdiff/ApreOdiff_sc;
% ApreOdiff_sc = 1;

% We care about the steady state behaviour only: calculate the coverage and
% correlations at that chemical potential of adsorbed oxygen, which is known
% assuming quasi-equilibrated NO oxidation / NO2 reduction

% beta = 1/(kboltz*Temp);

% Adsorbed atomic oxygen chemical potential referenced with respect to
% molecular NO2 and NO ideal gases at 1 atm
muOstar = kboltz*Temp*logyNO2yNO + Go_NO2(Temp) - Go_NO(Temp);
mOstarRange = [mOstarRange muOstar];

disp(['*** Solving ' ApproxName{k} ' for logyNO2yNO = ' num2str(logyNO2yNO)])

[PFcn,gcsoln,STheta,SCorrel,SCorrelUnocc,...
    SAvgekO2ads,SAvgekO2des] = ...
    SolveForThetaCorrelRates(beta,Temp,...
    EaO2ads,DEadsO,ProxFacO2ads,H0(k),H1,...
    muOstar,h,nn,ZPEO,QvibO,ApproxName{k},...
    CorrelLHS{k},CorrelRHS{k},...
    nsites{k},gcguess{k},options);

Theta(k,i) = STheta; % coverage
CorrelOccup(k,i) = SCorrel; % correlation of occupied sites
CorrelUnocc(k,i) = SCorrelUnocc; % correlation of unoccupied sites
AvgekO2ads(k,i) = SAvgekO2ads;
AvgekO2des(k,i) = SAvgekO2des;

gcsolnCell = num2cell(gcsoln);
% Solution has now converged, calculate correlations

end

%% Tests for residual when solving for the coverage...

% ResidualTheta(ApproxName{k},CorrelLHS{k},CorrelRHS{k},nsites{k},beta,...
%     H0(k),H1,h,nn,STheta,muOstar,gcsoln)

RfcnT = @(gc) ResidualTheta(ApproxName{k},CorrelLHS{k},CorrelRHS{k},nsites{k},beta,...
    H0(k),H1,h,nn,STheta,gc(1),gc(2:end));
gcsolnT = fsolve(RfcnT,[muOstar gcguess{k}],options);
gcsolnTCell = num2cell(gcsolnT);

disp('Difference between two solutions:');
[muOstar gcsoln].' - gcsolnT.'


% function [PFcn,gcsolnT,muOstar,CorrelOccup,CorrelUnocc, ...
%     AvgekO2ads,AvgekO2des] = ...
%     SolveForMuCorrelRates(beta,Temp,EaO2ads,DEadsO,ProxFacO2ads,H0,H1,...
%     Theta,h,nn,ZPEO,QvibO,ApproxIdent,...
%     CorrelLHS,CorrelRHS,nsites,muOstarguess,gcguess,options)


% [PFcn,gcsolnT,muOstar,CorrelOccup,CorrelUnocc, ...
%     AvgekO2ads,AvgekO2des] = ...
%     SolveForMuCorrelRates(beta,Temp,...
%     EaO2ads,DEadsO,ProxFacO2ads,H0(k),H1,...
%     STheta,h,nn,ZPEO,QvibO,ApproxName{k},...
%     CorrelLHS{k},CorrelRHS{k},...
%     nsites{k},muOstar,gcguess{k},options);

% STheta = 0.32;

[PFcn,gcsolnT,muOstar,CorrelOccup,CorrelUnocc, ...
    AvgekO2ads,AvgekO2des,AvgekNOoxi,AvgekNO2re] = ...
    SolveForMuCorrelRatesExtd(beta,Temp,...
    EaO2ads,DEadsO,ProxFacO2ads,...
    EaNOoxi,ProxFacNOoxi,...
    H0(k),H1,...
    STheta,h,nn,ZPEO,QvibO,ApproxName{k},...
    CorrelLHS{k},CorrelRHS{k},nsites{k},muOstar,gcguess{k},options);

yNO2 = 0.1*exp(-18);
yNO = yNO2*exp(-logyNO2yNO);

RNOoxi = ApreNOoxi*yNO*Pres*AvgekNOoxi.*STheta;
RNO2re = ApreNOoxi/PERatioNOoxiNO2red*yNO2*Pres*AvgekNO2re.*(1-STheta);
RO2ads = 2*nn*ApreO2ads*yO2*Pres*AvgekO2ads.*CorrelUnocc;
RO2des = 2*nn*ApreO2ads/PERatioO2adsdes*AvgekO2des.*CorrelOccup;
NetRate = RO2ads-RO2des;

% return

yNO2 = 0.1*exp(-18);
yNO = yNO2*exp(-logyNO2yNO);

ThetaRange = STheta + [-0.05:0.005:0.05];
oRes = [];
for ThetaRi = ThetaRange
    oRes = [oRes StdSttKineticsResidual(beta,Temp,Pres,yO2,...
        ApreO2ads,PERatioO2adsdes,EaO2ads,DEadsO,ProxFacO2ads,...
        ApreNOoxi,PERatioNOoxiNO2red,EaNOoxi,ProxFacNOoxi,...
        yNO2,yNO,H0(k),H1,...
        ThetaRi,h,nn,ZPEO,QvibO,ApproxName{k},...
        CorrelLHS{k},CorrelRHS{k},nsites{k},muOstar,gcguess{k},options)];
end
plot(ThetaRange,oRes,'o')

return
