if 1
clear all
clc

global hplanck clight kboltz

addpath('./NOOxidationNO2Reduction/')

Temp = 680; % Kelvin
nFigs = 3; % 3-Figure CE - up to 1NN interactions
 compareMC = true;
%compareMC = false;

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
logyNO2yNOrange = [-4 -3.5 -3 -2.5 -2 -1.5 -1 -0.5 0 0.2 0.5 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8];
%-4:0.5:8;
% logyNO2yNOrange = -4:0.5:-2;
% logyNO2yNOrange = -4:1.0:8;
% logyNO2yNOrange = -4;

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

%if 0
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
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size
%end
% Next the BP approximation without edge interactions
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BP'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2}; % right hand side in correlation equations
gcguess{kAppr} = 0.2332; % guess values for parameters (row vector)
H0(kAppr) = 0.3016; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [200 100 0]/255; % FOR PLOTTING: color
%[250 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = '^'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the BPE approximation with edge interactions
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BPE'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2}; % right hand side in correlation equations
gcguess{kAppr} = 0.1803; % guess values for parameters (row vector)
H0(kAppr) = 0.3147; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 0 255]/255; % FOR PLOTTING: color
%[200 200 0]/255; % FOR PLOTTING: color
%[200 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'd'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 7; % FOR PLOTTING: marker size
%if 0
% Next the BPEC approximation with edge interactions and an equation for
% the correlation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BPEC'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,[1 2]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,[2 3]}; % right hand side in correlation equations
gcguess{kAppr} = [0.1805   -0.0650]; % guess values for parameters (row vector)
H0(kAppr) = 0.3147; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [255 200 0]/255; % FOR PLOTTING: color
%[255 200 0]/255; % FOR PLOTTING: color
%[180 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'o'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

%end

if 0
    
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

% Next the K2NNC2 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC2T1'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2],[1 8],[1 8],[1 2 3]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8],[2 4],[8 9],[2 3 8]}; % right hand side in correlation equations
gcguess{kAppr} = [0.0866    0.3923    0.0723    0.0012   -0.0442   -0.1952   -0.1840]; % guess values for parameters (row vector)
H0(kAppr) = 0.7100; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 150 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = '<'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

end

if 0

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

% if 0
    
% % Next the K3NN approximation
% kAppr = kAppr + 1; % index of the approximation
% ApproxName{kAppr} = 'K3NN'; % name of the approximation
% nsites{kAppr} = 19; % number of sites in the approximation's cluster
% CorrelLHS{kAppr} = {1, 1,[1 2],[1  2],[1  2]}; % left hand side in correlation equations
% CorrelRHS{kAppr} = {8,14,[2 8],[2 14],[8 14]}; % right hand side in correlation equations
% gcguess{kAppr} = [0.1165    0.2745    0.0524   -0.0310    0.1431]; % guess values for parameters (row vector)
% H0(kAppr) = 1.0713; % Hamiltonian's constant (for renormalisations)
% ApproxColor{kAppr} = [200 100 0]/255; % FOR PLOTTING: color
% ApproxMarkr{kAppr} = '^'; % FOR PLOTTING: marker
% ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
% ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size

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
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% % Next the K3NNC3b approximation. This is similar to the above, but the
% % equation for the 3NN interactions uses as a reference (uncorrected) the
% % 2-5 correlation, instead of the 1-14.
% kAppr = kAppr + 1; % index of the approximation
% ApproxName{kAppr} = 'K3NNC3b'; % name of the approximation
% nsites{kAppr} = 19; % number of sites in the approximation's cluster
% CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8],[2  5],[2  5],[2  5]}; % left hand side in correlation equations
% CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14],[1 14],[2  9],[14 15]}; % right hand side in correlation equations
% gcguess{kAppr} = [0.1730    0.2695    0.1752   -0.0576   -0.0223   -0.0027    0.0029   -0.0080   -0.1360   -0.0806   -0.0579   -0.0771    0.1049]; % guess values for parameters (row vector)
% H0(kAppr) = 1.0915; % Hamiltonian's constant (for renormalisations)
% ApproxColor{kAppr} = [255 100 0]/255; % FOR PLOTTING: color
% ApproxMarkr{kAppr} = '+'; % FOR PLOTTING: marker
% ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
% ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size

% Next the K3NNC4 approximation.
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K3NNC4'; % name of the approximation
nsites{kAppr} = 19; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8],[2  5],[2  5],[2   5],[2 10],[2 10]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14],[1 14],[2  9],[14 15],[2 16],[9 14]}; % right hand side in correlation equations
gcguess{kAppr} = [0.2293    0.3057    0.1836   -0.0881   -0.0395   -0.0183   -0.0149   -0.0306   -0.1508   -0.0803   -0.0579   -0.1019    0.1096   -0.0063    0.0125]; % guess values for parameters (row vector)
H0(kAppr) = 1.0915; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [180 100 0]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'x'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size

% % Next the K3NNC5 approximation.
% kAppr = kAppr + 1; % index of the approximation
% ApproxName{kAppr} = 'K3NNC5'; % name of the approximation
% nsites{kAppr} = 19; % number of sites in the approximation's cluster
% CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8],[2  5],[2  5],[2   5],[2 10],[2 10],[ 8 11]}; % left hand side in correlation equations
% CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14],[1 14],[2  9],[14 15],[2 16],[9 14],[14 16]}; % right hand side in correlation equations
% gcguess{kAppr} = [0.1730    0.2695    0.1752   -0.0576   -0.0223   -0.0027    0.0029   -0.0080   -0.1360   -0.0806   -0.0579   -0.0771    0.1049 0.0 0.0 0.0]; % guess values for parameters (row vector)
% H0(kAppr) = 1.0915; % Hamiltonian's constant (for renormalisations)
% ApproxColor{kAppr} = [50 100 0]/255; % FOR PLOTTING: color
% ApproxMarkr{kAppr} = 'd'; % FOR PLOTTING: marker
% ApproxLinSpec{kAppr} = '--'; % FOR PLOTTING: line type
% ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
end
%% Initialise arrays for data on coverage, correlations, and rates

% The following are intended to be indexed by 
% i -> logyNO2yNOrange and 
% j -> approximation type
Theta = [];
CorrelOccup = [];
CorrelUnocc = [];
AvgekO2ads = [];
AvgekO2des = [];

%% Initialise arrays for extrapolation of solution guesses and H0 

% Degree of the interpolation for the continuation scheme
ndeg = 3;

for k = 1:kAppr
    for j = 1:ndeg+1
        H0v{k}(j) = NaN;
        gcguessv{k}(:,j) = NaN(length(gcguess{k}),1);
    end
    disp(['*** Guess vector and H0v for ' ApproxName{k}])
    disp(gcguessv{k});
    disp(H0v{k});
    disp(gcguess{k});
    disp(H0(k));
%     pause
end

% return

%% Entering main loop
mOstarRange = [];

for i = 1:length(logyNO2yNOrange);
        
    logyNO2yNO = logyNO2yNOrange(i);
    
    [EnrgGasSpecies,QvibO,ZPEO,...
        ApreO2ads,PERatioO2adsdes,EaO2ads,ProxFacO2ads,...
        ApreOdiff,PERatioOdiff,EaOdiff,ProxFacOdiff,ApreOdiff_sc,...
        ApreNOoxi,PERatioNOoxiNO2red,EaNOoxi,ProxFacNOoxi,ApreNOoxi_sc] = ...
        EnergeticsAndKineticsTlnPr(Temp,nFigs,logyNO2yNO);
    % Undo any stiffness scalings:
    ApreNOoxi = ApreNOoxi/ApreNOoxi_sc;
    ApreNOoxi_sc = 1;
    ApreOdiff = ApreOdiff/ApreOdiff_sc;
    ApreOdiff_sc = 1;
    
    % We care about the steady state behaviour only: calculate the coverage and
    % correlations at that chemical potential of adsorbed oxygen, which is known
    % assuming quasi-equilibrated NO oxidation / NO2 reduction
    
    beta = 1/(kboltz*Temp);

    % Adsorbed atomic oxygen chemical potential referenced with respect to
    % molecular NO2 and NO ideal gases at 1 atm
    muOstar = kboltz*Temp*logyNO2yNO + Go_NO2(Temp) - Go_NO(Temp);
    mOstarRange = [mOstarRange muOstar];

    for k = 1:kAppr

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

        [H0v,H0,gcguessv,gcguess] = ...
            InitialGuessesContinuation(ndeg,i,k,beta,logyNO2yNOrange,PFcn,gcsoln,...
        H0v,H0,gcguessv,gcguess);
    end
    
end

% Evaluation of rates

for i = 1:kAppr
    RO2ads(i,:) = 2*nn*ApreO2ads*yO2*Pres*AvgekO2ads(i,:).*CorrelUnocc(i,:);
    RO2des(i,:) = 2*nn*ApreO2ads/PERatioO2adsdes*AvgekO2des(i,:).*CorrelOccup(i,:);
    NetRate(i,:) = RO2ads(i,:)-RO2des(i,:);
    indic(i,:) = sign(NetRate(i,:));
    indicFwd(i,:) = indic(i,:) > 0;
    indicRev(i,:) = indic(i,:) < 0;
end

else

clear all
clc

load('SecondTestK3NNC3_Corrected1NNz1NNCorrel.mat')
% ApproxMSize(5) = 4; % FOR PLOTTING: marker size
% ApproxMSize(8) = 6; % FOR PLOTTING: marker size
% ApproxMSize(9) = 6; % FOR PLOTTING: marker size

% load('TempLast.mat')
% ApproxColor{kAppr} = [180 100 0]/255; % FOR PLOTTING: color
% ApproxMarkr{kAppr} = 'x'; % FOR PLOTTING: marker
% compareMC = true;

% load('FirstTestK3NNC4_T680K.mat')

end

%% PLOTTING COMMANDS FOLLOW

close all

FntSz = 20;
 FntNm = 'Times New Roman';
%FntNm = 'Calibri';
margx0 = 0.16;
margx1 = 0.01;
margy0 = 0.35;
margy1 = 0.1;

if compareMC
    if Temp == 480
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\KineticMonteCarlo\Postprocessing\logNO2NO_vs_TOF_nFig03_T480_NOMINAL.txt');
 %       M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\KineticMonteCarlo_z1AfterSystStiff\SimulationSet02\EnsemblesCreateTlogPrRevers\logNO2NO_vs_TOF_nFigs03_T480.txt');
        M=load('logNO2NO_vs_TOF_nFigs03_T480.txt')
    elseif Temp == 680
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\KineticMonteCarlo\Postprocessing\logNO2NO_vs_TOF_nFig03_T680_NOMINAL.txt');
%        M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\KineticMonteCarlo_z1AfterSystStiff\SimulationSet11\Postprocessing\logNO2NO_vs_TOF_nFigs03_T680.txt');
         M=load('logNO2NO_vs_TOF_nFigs03_T680.txt')
    end
    M = M(1:length(logyNO2yNOrange),:);
end

%% Reaction rates vs log(yNO2/yNO)

logconvc = 1/log(10);
% logconvc = 1;

posfig1 = [178   230  900   570];
figure('Color','w','Position',posfig1)
set(gcf,'InvertHardCopy','on','Color',[235 235 235]/255)

legendprop = {};
hplot = [];

hold on

NetRateGCMC = [];
if compareMC
    indicFwdKMC = M(:,3) > 0;
    indicRevKMC = M(:,3) < 0;
    hplot = [hplot plot(logconvc*M(indicFwdKMC,1),logconvc*M(indicFwdKMC,2), ...
        '-.sr','MarkerFaceColor','r','MarkerSize',8)]; % 12
    plot(logconvc*M(indicRevKMC,1),logconvc*M(indicRevKMC,2),'-.sr','MarkerSize',8) % 10
%    plot(logconvc*0.2,7.949052e-01,'>','Color',[0 0 0],'MarkerFaceColor','None','MarkerSize',6)
%    plot(logconvc*0.2,-2.573523e-01,'^','Color',[200 100 0]/255,'MarkerFaceColor','None','MarkerSize',6)
%    plot(logconvc*0.2,1.931566e-01,'d','Color',[0 0 255]/255,'MarkerFaceColor','None','MarkerSize',7)
%    plot(logconvc*0.2,1.786939e-01,'o','Color',[255 200 0]/255,'MarkerFaceColor','None','MarkerSize',6)
   % plot(4.342945e-02,logconvc*(-4.1295),'o','Color',[255 0 0]/255,'MarkerFaceColor','None','MarkerSize',8)
    NetRateGCMC = exp(M(:,2).');
    legendprop = {legendprop{:} 'KMC'};
end

% for k = 1:kAppr % k = [1 2 3 4 5 6 8 9 10]; % 1:kAppr
% for k = [1 2 3 4 5 6 8 9 10]; % 1:kAppr
for k = 1:kAppr
    hplot = [hplot plot(logconvc*logyNO2yNOrange(indicFwd(k,:)),logconvc*log(NetRate(k,indicFwd(k,:))),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor',ApproxColor{k},'MarkerSize',ApproxMSize(k))];
    plot(logconvc*logyNO2yNOrange(indicRev(k,:)),logconvc*log(-NetRate(k,indicRev(k,:))),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor','None','MarkerSize',ApproxMSize(k))
    legendprop = {legendprop{:} ApproxName{k}};
end
absAllNetRates = log(abs([NetRateGCMC; NetRate]));
% plot([1 1]*logyNO2yNO_equil,[min(min(absAllNetRates)) max(max(absAllNetRates))],'--',...
%     'Color',[255 210 0]/255)
plot([1 1]*logconvc*logyNO2yNO_equil,[-20 10],'--',...
    'Color','k')

hold off
box on

set(gca,...
    'FontName',FntNm,...
    'FontSize',FntSz,...
    ...'Ylim',[min(min(absAllNetRates)) max(max(absAllNetRates))],...
    'Ylim',logconvc*[-18 11],...
    'Xlim',logconvc*[-4.5 8.5],...
    ...'Ylim',[-15 10],...
    ...'YTick',[-18:4:6],...
    ...'Xlim',[0.95*min(ylab) 1.05*max(ylab)],...
    ...'YScale','log',...
    ...'YTick',logspace(-8,8,17),...
    'YMinorTick','on',...
    'XMinorTick','on',...
    'XScale','lin',...
    ...'XTick',[-4:1:8],...
    'Color','w',...
    'Layer','Bottom')

set(gca,'Position',[margx0 margy0 1-margx0-margx1 1-margy0-margy1]);

xlabel('log(P_{NO2}/P_{NO})',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'Position', [1.0, -9.8, 0],...  
    'FontSize',FntSz+6);

ylabel('log(TOF)',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',FntSz+6);
 if ~isempty(hplot)
     legend(hplot,legendprop,...
         'FontName',FntNm,...
         'FontSize',16,'Location','SouthEast')
     legend('boxoff')
 end
text(0.0,6.0,'T=680 K','FontSize',30)
%text(-1.8,5.5,'(a'')','FontSize',30) 
text(-1.8,5.9,'(d)','FontSize',30)
%text(4.342945e-02,5.372262e-01,'\Delta','FontSize',12)
%text(4.342945e-02,-5.170725e-01,'>','FontSize',12)
%text(4.342945e-02,-7.172230e-02,'>','FontSize',12)
%text(4.342945e-02,-8.573918e-02,'o','FontSize',20)
%text(4.342945e-02,-8.573918e-02,'o','FontSize',17,'Color',[255 200 0]/255)
%plot(4.342945e-02,-8.573918e-02,'o','FontSize',6)
eval(['print -dpng TOFvsRatioBPT680']);
eval(['print -depsc Fig5d']);


%% Coverages vs log(yNO2/yNO)

posfig2 = [178   230  900   570];
figure('Color','w','Position',posfig2)
set(gcf,'InvertHardCopy','on','Color',[235 235 235]/255)

legendprop = {};
hplot = [];

hold on

if compareMC
    indicFwdKMC = M(:,3) > 0;
    indicRevKMC = M(:,3) < 0;
    hplot = [hplot plot(logconvc*M(indicFwdKMC,1),M(indicFwdKMC,5),'-.sr',...
        'MarkerFaceColor','r','MarkerSize',8)]; % 12
    plot(logconvc*M(indicRevKMC,1),M(indicRevKMC,5),'-.sr','MarkerSize',8) % 10
    legendprop = {legendprop{:} 'KMC'};
end

% for k = 1:kAppr % k = [1 2 3 4 5 6 8 9 10]; % 1:kAppr
% for k = [1 2 3 4 5 6 8 9 10]; % 1:kAppr
for k = 1:kAppr
    hplot = [hplot plot(logconvc*logyNO2yNOrange(indicFwd(k,:)),Theta(k,indicFwd(k,:)),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor',ApproxColor{k},'MarkerSize',ApproxMSize(k))];
    plot(logconvc*logyNO2yNOrange(indicRev(k,:)),Theta(k,indicRev(k,:)),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor','None','MarkerSize',ApproxMSize(k))
    legendprop = {legendprop{:} ApproxName{k}};
end
plot([1 1]*logconvc*logyNO2yNO_equil,[-20 10],'--',...
    'Color','k')

hold off
box on

set(gca,...
    'FontName',FntNm,...
    'FontSize',FntSz,...
    'Ylim',[0 1],...
    'Xlim',logconvc*[-4.5 8.5],...
    'YMinorTick','on',...
    'XMinorTick','on',...
    'XScale','lin',...
    'XTick',[-4:1:8],...
    'YTick',[0:0.2:1],...
    'Color','w',...
    'Layer','Bottom')

% set(gca,'Position',[ ...
%     margx0*posfig1(3)/posfig2(3) ...
%     margy0*posfig1(4)/posfig2(4) ...
%     (1-margx0-margx1)*posfig1(3)/posfig2(3) ...
%     1-margy0-margy1]);
% 
set(gca,'Position',[margx0 margy0 1-margx0-margx1 1-margy0-margy1]);
xlabel('log(P_{NO2}/P_{NO})',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'Position', [1.0, -0.15, 0],...
    'FontSize',FntSz+6);

ylabel('O^{*} Coverage',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',FntSz+6);
if ~isempty(hplot)
    hleg2 = legend(hplot,legendprop,...
        ...'FontName',FntNm,...
        'FontName',FntNm,...
        'FontSize',16,'Location','NorthWest'); %'NorthWest')
    legend('boxoff')
end

%hlegpos2 = get(hleg2,'Position');
%set(gca,'Position',[ ...
%    margx0*posfig1(3)/posfig2(3) ...
%    margy0*posfig1(4)/posfig2(4) ...
%    (1-margx0-margx1)*posfig1(3)/posfig2(3) ...
%    (1-margy0-margy1)*posfig1(4)/posfig2(4) ...
%    ]);

text(-1.8,1.1,'(b)','FontSize',30)
eval(['print -dpng ThetaVsNO2NORatioBPT680']);
eval(['print -depsc ThetaVsNO2NORatioBPT680']);

return

%% Correlations vs log(yNO2/yNO)
%posfig3 = [178   230  900   570];
%figure('Color','w','Position',posfi3)
figure('Color','w','Position',[178   230  1000   570])
set(gcf,'InvertHardCopy','on','Color',[235 235 235]/255)
%figure('Color','w','Position',[178   230  1000   570])

legendprop = {};
hplot = [];

hold on

% if compareMC
%     indicFwdKMC = M(:,3) > 0;
%     indicRevKMC = M(:,3) < 0;
%     hplot = [hplot plot(M(indicFwdKMC,1),M(indicFwdKMC,2),'-.sr','MarkerFaceColor','r','MarkerSize',8)];
%     plot(M(indicRevKMC,1),M(indicRevKMC,2),'-.sr','MarkerSize',8)
%     NetRateGCMC = exp(M(:,2).');
%     legendprop = {legendprop{:} 'Kinetic MC'};
% end

for k = 1:kAppr
    hplot = [hplot plot(Theta(k,indicFwd(k,:)),CorrelOccup(k,indicFwd(k,:)),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor',ApproxColor{k},'MarkerSize',ApproxMSize(k))];
    plot(Theta(k,indicRev(k,:)),CorrelOccup(k,indicRev(k,:)),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor','None','MarkerSize',ApproxMSize(k));
    legendprop = {legendprop{:} ApproxName{k}};
end
% plot([1 1]*muOstar_qequil,[min(min(CorrelOccup)) max(max(CorrelOccup))],'--',...
%     'Color',[255 210 0]/255)

hold off
box on

set(gca,...
    'FontName',FntNm,...
    'FontSize',FntSz,...
    ...'Ylim',[min(min(absAllNetRates)) max(max(absAllNetRates))],...
    'Ylim',[-0.05 1.05],...
    'Xlim',[-0.05 1.05],...
    ...'Ylim',[-15 10],...
    ...'XTick',[-4:2:8],...
    ...'YTick',[-18:2:6],...
    ...'Xlim',[0.95*min(ylab) 1.05*max(ylab)],...
    ...'YScale','log',...
    ...'YTick',logspace(-8,8,17),...
    'YMinorTick','on',...
    'XMinorTick','on',...
    'XScale','lin',...
    'Color','none',...
    'Layer','Bottom')

set(gca,'Position',[margx0 margy0 1-margx0-margx1 1-margy0-margy1]);

xlabel('\Theta',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',18);

ylabel('Correlation',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',18);
if ~isempty(hplot)
    legend(hplot,legendprop,...
        'FontName',FntNm,...
        'FontSize',16,'Location','NorthWest')
    legend('boxoff')
end

eval(['print -dpng CorrelVsTheta']);

return



