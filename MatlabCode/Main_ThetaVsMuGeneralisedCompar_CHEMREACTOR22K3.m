if 1
clear all
clc

global hplanck clight kboltz

addpath('./NOOxidationNO2Reduction/')

Temp = 480; % 480; % Kelvin
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
muOstar_qequil = Go_NO2(Temp) - Go_NO(Temp);
muOstarLeft = kboltz*Temp*(-4) + Go_NO2(Temp) - Go_NO(Temp);
muOstarRange = -1.5:0.02:1.0;
%:0.02:1.0;
% muOstarRange = muOstarLeft:-0.02:-1.5;

%% Solver options

options = optimoptions('fsolve',...
    'TolFun',1e-8,...
    'TolX',1e-8,...
    'Display','iter-detailed','MaxIter',800,'MaxFunEvals',100000);
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
%kAppr = kAppr + 1; % index of the approximation
%ApproxName{kAppr} = 'MF'; % name of the approximation
%nsites{kAppr} = 1; % number of sites in the approximation's cluster
%CorrelLHS{kAppr} = {}; % left hand side in correlation equations
%CorrelRHS{kAppr} = {}; % right hand side in correlation equations
%gcguess{kAppr} = 0*0.3; % guess values for parameters (row vector)
%H0(kAppr) = 0*0; % Hamiltonian's constant (for renormalisations)
%ApproxColor{kAppr} = [0 0 0]; % FOR PLOTTING: color
%ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
%ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
%ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
%ApproxMWith(kAppr)= 10;
if 0
% Next the BP approximation without edge interactions
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BP'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2}; % right hand side in correlation equations
gcguess{kAppr} = 0*0.2332; % guess values for parameters (row vector)
H0(kAppr) = 0*0.3016; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [250 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '--'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
ApproxMWith(kAppr)= 2;

% Next the BPE approximation with edge interactions
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BPE'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2}; % right hand side in correlation equations
gcguess{kAppr} = 0*0.1803; % guess values for parameters (row vector)
H0(kAppr) = 0*0.3147; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [200 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ':'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
ApproxMWith(kAppr)= 2;

% Next the BPEC approximation with edge interactions and an equation for
% the correlation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'BPEC'; % name of the approximation
nsites{kAppr} = 7; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,[1 2]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,[2 3]}; % right hand side in correlation equations
gcguess{kAppr} = 0*[0.1805   -0.0650]; % guess values for parameters (row vector)
H0(kAppr) = 0*0.3147; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [180 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-.'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
ApproxMWith(kAppr)= 2;

% Next the K2NNC1 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC1'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8]}; % right hand side in correlation equations
gcguess{kAppr} = 0*[0.0669    0.2188   -0.0084   -0.0380]; % guess values for parameters (row vector)
H0(kAppr) = 0*0.6568; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '--'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size

% Next the K2NNC2 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC2'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2],[1 8],[1 8]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8],[2 4],[8 9]}; % right hand side in correlation equations
gcguess{kAppr} = 0*[0.0866    0.3923    0.0723    0.0012   -0.0442   -0.1952]; % guess values for parameters (row vector)
H0(kAppr) = 0*0.7100; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 100 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ':'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size
%if 0
% Next the K2NNC3 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K2NNC3'; % name of the approximation
nsites{kAppr} = 13; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1,[1 2],[1 2],[1 8],[1 8],[2 5]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,[2 3],[2 8],[2 4],[8 9],[4 8]}; % right hand side in correlation equations
gcguess{kAppr} = 0*[0.0956    0.3970    0.0604   -0.0036   -0.0474   -0.1970   -0.0063]; % guess values for parameters (row vector)
H0(kAppr) = 0*0.6927; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 200 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-.'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
end
%if 0
% The first approximation is the mean-field and is treated in a special way
% since correlations are neglected
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'MF'; % name of the approximation
nsites{kAppr} = 1; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {}; % left hand side in correlation equations
CorrelRHS{kAppr} = {}; % right hand side in correlation equations
gcguess{kAppr} = 1*0.3; % guess values for parameters (row vector)
H0(kAppr) = 1*0; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [0 0 0]; % FOR PLOTTING: color
ApproxMarkr{kAppr} = ''; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size
ApproxMWith(kAppr)= 10;
%end
%if 0
% % Next the K3NNC1 approximation
 kAppr = kAppr + 1; % index of the approximation
 ApproxName{kAppr} = 'K3NNC1'; % name of the approximation
 nsites{kAppr} = 19; % number of sites in the approximation's cluster
 CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2]}; % left hand side in correlation equations
 CorrelRHS{kAppr} = {2,8,14,[2,3],[2,8],[2 14],[8 14]}; % right hand side in correlation equations
 gcguess{kAppr} = 0*[0.1165  0.2745   0.3094   0.0524   -0.0310   0.1431  -0.0441]; % guess values for parameters (row vector)
 H0(kAppr) = 1*1.0713; % Hamiltonian's constant (for renormalisations)
 ApproxColor{kAppr} = [200 100 0]/255; % FOR PLOTTING: color
 ApproxMarkr{kAppr} = '^'; % FOR PLOTTING: marker
 ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
 ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size
 ApproxMWith(kAppr)= 2;
%end
%if 0
% Next the K3NNC2 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K3NNC2'; % name of the approximation
nsites{kAppr} = 19; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14]}; % right hand side in correlation equations
%gcguess{kAppr} = 1*[0.1165    0.2745    0.3094    0.0524   -0.0310    0.1431   -0.0441    0.0694   -0.1377   -0.1548]; % guess values for parameters (row vector)
%gcguess{kAppr} =  1*[0.1730    0.2695    0.1752   -0.0576   -0.0223   -0.0027    0.0029   -0.0080   -0.1360   -0.0806]
gcguess{kAppr} = 1*[0.0866    0.3923    0.0723    0.0012   -0.0442   -0.1952  0.1165  0.2745   0.3094   0.0524];
H0(kAppr) = 1*1.0713; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [255 200 0]/255; % FOR PLOTTING: color
%[0 0 255]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'o'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = ''; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 6; % FOR PLOTTING: marker size
ApproxMWith(kAppr)= 2;
if 0
% Next the K3NNC3 approximation
kAppr = kAppr + 1; % index of the approximation
ApproxName{kAppr} = 'K3NNC3'; % name of the approximation
nsites{kAppr} = 19; % number of sites in the approximation's cluster
CorrelLHS{kAppr} = {1,1, 1,[1 2],[1 2],[1  2],[1  2],[1 8],[1 8],[1  8],[1 14],[1 14],[ 1 14]}; % left hand side in correlation equations
CorrelRHS{kAppr} = {2,8,14,[2 3],[2 8],[2 14],[8 14],[2 4],[8 9],[3 14],[2  5],[2  9],[14 15]}; % right hand side in correlation equations
gcguess{kAppr} = 0*[0.1730    0.2695    0.1752   -0.0576   -0.0223   -0.0027    0.0029   -0.0080   -0.1360   -0.0806   -0.0579   -0.0771    0.1049]; % guess values for parameters (row vector)
H0(kAppr) = 0*1.0915; % Hamiltonian's constant (for renormalisations)
ApproxColor{kAppr} = [255 200 0]/255; % FOR PLOTTING: color
ApproxMarkr{kAppr} = 'o'; % FOR PLOTTING: marker
ApproxLinSpec{kAppr} = '-'; % FOR PLOTTING: line type
ApproxMSize(kAppr) = 8; % FOR PLOTTING: marker size
ApproxMWith(kAppr)= 2;
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
%     disp(['*** Guess vector and H0v for ' ApproxName{k}])
%     disp(gcguessv{k});
%     disp(H0v{k});
%     disp(gcguess{k});
%     disp(H0(k));
%     pause
end

for k = 1:kAppr
    Theta(k,1:length(muOstarRange)) = NaN; % coverage
    CorrelOccup(k,1:length(muOstarRange)) = NaN; % correlation of occupied sites
    CorrelUnocc(k,1:length(muOstarRange)) = NaN; % correlation of unoccupied sites
    AvgekO2ads(k,1:length(muOstarRange)) = NaN;
    AvgekO2des(k,1:length(muOstarRange)) = NaN;
    gcsolnv{k} = NaN(length(gcguess{k}),length(muOstarRange));
end

% return

%% Entering main loop

for i = 1:length(muOstarRange);
    
    muOstar = muOstarRange(i);
    logyNO2yNO = (muOstar - Go_NO2(Temp) + Go_NO(Temp))/(kboltz*Temp);
    
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

    for k = 1:kAppr

        disp(['*** Solving ' ApproxName{k} ' for logyNO2yNO = ' num2str(logyNO2yNO) ' for muOstar = ' num2str(muOstar)])
            
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

        gcsolnv{k}(:,i) = gcsoln;
        
        [H0v,H0,gcguessv,gcguess] = ...
            InitialGuessesContinuation(ndeg,i,k,beta,muOstarRange,PFcn,gcsoln,...
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

% clear all
% clc
% load('ThetaVsMu_AllApproxMF_K3NNC3.mat')
% ApproxMSize(9) = 5; % FOR PLOTTING: marker size
% compareMC = true;

clear all
clc
load('ThetaVsMu_AllApproxMF_K3NNC3_T480.mat')
% load('ThetaVsMu_AllApproxMF_K3NNC3_T800.mat')
compareMC = true;

end

% return

%% PLOTTING COMMANDS FOLLOW

% close all

FntSz = 20;
FntNm = 'Times New Roman';
% margx0 = 0.08;
% margx1 = 0.01;
% margy0 = 0.15;
% margy1 = 0.05;
margx0 = 0.16;
margx1 = 0.01;
margy0 = 0.25;
margy1 = 0.05;

if compareMC
    if Temp == 480
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_384x288_FullRange\results_output.txt');
%  M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_096x072_FullRange_T480_GrndIniStat\results_output.txt');
   % M = load('/Users/mpineda/Desktop/metropolis/metropolis/MetropolisLatticeGasFortranRestartable/OutputLegion/Run_096x072_FullRange_T480_GrndIniStat_CONT1/results_output.txt');
     M=load('results_output.txt')
   %%M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_096x072_FullRange_T480_GrndIniStat_CONT1\results_output.txt');        
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_576x432_FullRange_T480_GrndIniStat\results_output.txt');
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_576x432_FullRange_T480_IniStatTiled4x4_48hr\results_output.txt');
    elseif Temp == 680
        M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\KineticMonteCarlo\Postprocessing\logNO2NO_vs_TOF_nFig03_T680_NOMINAL.txt');
    elseif Temp == 800
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_384x288_FullRange_T800\results_output.txt');
%         M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_576x432_FullRange_T800\results_output.txt');
        M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_576x432_FullRange_T800_IniStatTiled4x4_48hr\results_output.txt');
    elseif Temp == 1000
        M = load('E:\Academic\Research\CEng_Lectureship_UCL\Computations\ReducedKineticModels\MetropolisLatticeGasFortranRestartable\MetropolisLatticeGasFortranRestartable\OutputLegion\Run_384x288_FullRange_T1000\results_output.txt');
    end
    M(:,1) = M(:,1) - (-1.20 - H1);
end

%% Coverage vs muO*

posfig1 = [178   330  720   570];

figure('Color','w','Position',posfig1);
set(gcf,'InvertHardCopy','on','Color',[235 235 235]/255)

legendprop = {};
hplot = [];

Di = 3;

hold on

if compareMC
    hplot = [hplot plot(M(1:2*Di:end,1),M(1:2*Di:end,2),'-.sr','MarkerFaceColor','r','MarkerSize',5)];
    legendprop = {legendprop{:} 'MMC'};
end

for k = 1:kAppr
    hplot = [hplot plot(muOstarRange(1:Di:end),Theta(k,1:Di:end),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor',ApproxColor{k},'MarkerSize',ApproxMSize(k)),...
        'Linewidth',2.0];
    legendprop = {legendprop{:} ApproxName{k}};
end
%plot([1 1]*muOstar_qequil,[min(min(Theta)) max(max(Theta))],'--',...
%    'Color',[255 210 0]/255)

hold off
box on

set(gca,...
    'FontName',FntNm,...
    'FontSize',FntSz,...
    ...'Ylim',[min(min(absAllNetRates)) max(max(absAllNetRates))],...
    'Ylim',[-0.05 1.05],...
    'Xlim',[-1.53 1.03],...
    ...'Ylim',[-15 10],...
    ...'XTick',[-4:2:8],...
    ...'YTick',[-18:2:6],...
    ...'Xlim',[0.95*min(ylab) 1.05*max(ylab)],...
    ...'YScale','log',...
    ...'YTick',logspace(-8,8,17),...
    'YMinorTick','on',...
    'XMinorTick','on',...
    'XScale','lin',...
    'Color','w',...
    'Layer','Bottom')

set(gca,'Position',[margx0 margy0 1-margx0-margx1 1-margy0-margy1]);

xlabel('\mu_{O*}',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',FntSz+6);

ylabel('O^{*} Coverage',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',FntSz+6);
 if ~isempty(hplot)
     legend(hplot,legendprop,...
         'FontName',FntNm,...
         'FontSize',16,'Location','NorthWest')
     legend('boxoff')
 end
%text(-1.715,1.09,'(e)','FontSize',30)
text(-1.78,1.09,'(c)','FontSize',30)
eval(['print -dpng ThetaVsMuK3']);
eval(['print -depsc Fig4c']);
%eval(['print -dmeta ThetaVsMu']);


%% Correlation vs muO*

%posfig2 = [178   130  1000   570];
posfig2 = [178   330  720   570];

figure('Color','w','Position',posfig2)
set(gcf,'InvertHardCopy','on','Color',[235 235 235]/255)

legendprop = {};
hplot = [];

hold on

Di = 3;

if compareMC
    hplot = [hplot plot(M(1:Di:end,2),M(1:Di:end,4),'-.sr','MarkerFaceColor','r','MarkerSize',5)];
    legendprop = {legendprop{:} 'MMC'};
end

for k = 1:kAppr
    hplot = [hplot plot(Theta(k,1:Di:end),CorrelOccup(k,1:Di:end),...
        [ApproxLinSpec{k} ApproxMarkr{k}],'Color',ApproxColor{k},...
        'MarkerFaceColor',ApproxColor{k},'MarkerSize',ApproxMSize(k)),...
        'Linewidth',10];
    legendprop = {legendprop{:} ApproxName{k}};
end
% plot([1 1]*muOstar_qequil,[min(min(CorrelOccup)) max(max(CorrelOccup))],'--',...
%     'Color',[255 210 0]/255)

hold off
box on

set(gca,...
    'FontName',FntNm,...
    'FontSize',FntSz,...
    'Ylim',[-0.05 1.05],...
    'Xlim',[-0.02 1.02],...
    ...'Ylim',[min(min(absAllNetRates)) max(max(absAllNetRates))],...
    ...'Ylim',[-0.05e-4 4.05e-4],...
    ...'Ylim',[-0.05e-4 5.05e-3],...
    ...'Xlim',[-0.05 1.0/3+0.05],...
    ...'Ylim',[-15 10],...
    ...'XTick',[-4:2:8],...
    ...'YTick',[-18:2:6],...
    ...'Xlim',[0.95*min(ylab) 1.05*max(ylab)],...
    ...'YScale','log',...
    ...'YTick',logspace(-8,8,17),...
    'YMinorTick','on',...
    'XMinorTick','on',...
    'XScale','lin',...
    'Color','w',...
    'Layer','Bottom')

set(gca,'Position',[margx0 margy0 1-margx0-margx1 1-margy0-margy1]);

xlabel('O^{*} Coverage',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',FntSz+6);

ylabel('Correlation',...
    'FontName',FntNm,...
    'FontWeight','Bold',...
    'FontSize',FntSz+6);
if ~isempty(hplot)
    hleg2 = legend(hplot,legendprop,...
        'FontName',FntNm,...
        'FontSize',16,'Location','NorthWest');
    legend('boxoff')
end
hlegpos2 = get(hleg2,'Position');
set(gca,'Position',[ ...
    margx0*posfig1(3)/posfig2(3) ...
    margy0*posfig1(4)/posfig2(4) ...
    (1-margx0-margx1)*posfig1(3)/posfig2(3) ...
    (1-margy0-margy1)*posfig1(4)/posfig2(4) ...
    ]);
%text(-0.1,1.09,'(f)','FontSize',30)
text(-0.125,1.09,'(f)','FontSize',30)
eval(['print -dpng CorrelVsThetaK3']);
eval(['print -depsc Fig4f']);
%eval(['print -dmeta CorrelVsTheta']);

return



