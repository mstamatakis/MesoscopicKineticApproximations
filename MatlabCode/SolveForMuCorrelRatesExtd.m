function [PFcn,gcsolnT,muOstarT,CorrelOccup,CorrelUnocc, ...
    AvgekO2ads,AvgekO2des,AvgekNOoxi,AvgekNO2re] = ...
    SolveForMuCorrelRatesExtd(beta,Temp,...
    EnrgGasSpecies,EaO2ads,DEadsO,ProxFacO2ads,EaNOoxi,ProxFacNOoxi,H0,H1,...
    Theta,h,nn,ZPEO,QvibO,ApproxIdent,...
    CorrelLHS,CorrelRHS,nsites,muOstarguess,gcguess,options)

global kboltz

% First solve for the parameters of the approximation and the chemical
% potential
RfcnT = @(gc) ResidualTheta(ApproxIdent,CorrelLHS,CorrelRHS,nsites,beta,...
    H0,H1,h,nn,Theta,gc(1),gc(2:end));
% muOstarguess
% gcguess
gcsolnT = fsolve(RfcnT,[muOstarguess gcguess],options);
muOstarT = gcsolnT(1);
gcsolnT = gcsolnT(2:end);
gcsolnTCell = num2cell(gcsolnT);

if isempty(CorrelLHS) % MF approximation treated specially (trivial evaluations given Theta)
    PFcn = exp(-beta*H0) + exp(-beta*(H0 + DEadsO + nn*h*Theta-muOstarT));
    DAMF = 2*(DEadsO + nn*h*Theta); % free energy difference between final and initial state of O2 adsorption
    CorrelOccup = Theta^2; % correlation of occupied sites
    CorrelUnocc = 1-2*Theta+CorrelOccup; % correlation of unoccupied sites
    AvgekO2ads = exp(-beta*(EaO2ads+ProxFacO2ads*(DAMF-(2*DEadsO+h))));
    AvgekO2des = exp(-beta*(EaO2ads-(2*DEadsO+h)-(1-ProxFacO2ads)*(DAMF-2*DEadsO)));
    
    DAoxi0  = EnrgGasSpecies(3) - (EnrgGasSpecies(2) + DEadsO);
    DAoxiMF = EnrgGasSpecies(3) - (EnrgGasSpecies(2) + DEadsO + nn*h*Theta);
    Eaoxi = max([0, DAoxiMF,EaNOoxi-ProxFacNOoxi*nn*h*Theta]);
    Eared = max([0,-DAoxiMF,EaNOoxi-DAoxi0+(1-ProxFacNOoxi)*nn*h*Theta]);
    AvgekNOoxi = exp(-beta*(Eaoxi));
    AvgekNO2re = exp(-beta*(Eared));
%     disp(['DErxnMF = ' num2str(DAoxiMF) ' ' num2str(Eaoxi-Eared)]);
    return
end

% Initialisation of variables
PFcn = 0;
PFcnEps0Occup = 0;
PFcnEps0Unocc = 0;
PFcnEps0Eps1Occup = 0;
PFcnEps0Eps1Unocc = 0;
AvgekO2ads = 0;
AvgekO2des = 0;

% Create the array containing all occupancy states
AllStatesArray = AllStates(nsites);

% Calculate all quantities of interest using vectorised expressions
EAllMicroStates = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnTCell{:},...
    AllStatesArray);

PFcnTerm = exp(-beta*(EAllMicroStates - muOstarT*sum(AllStatesArray,2)));

PFcn = sum(PFcnTerm);
PFcnEps0Occup = sum(AllStatesArray(:,1).*PFcnTerm);
PFcnEps0Unocc = sum((1-AllStatesArray(:,1)).*PFcnTerm);
PFcnEps0Eps1Occup = sum(all(AllStatesArray(:,1:2) == 1,2).*PFcnTerm);
PFcnEps0Eps1Unocc = sum(all(AllStatesArray(:,1:2) == 0,2).*PFcnTerm);

% Desorption can happen from the following states
indxDesorp = find(AllStatesArray(:,1) == 1 & AllStatesArray(:,2) == 1).';
AllStatesDesorpIni = AllStatesArray(indxDesorp,:);
AllStatesDesorpFin = [repmat([0 0],length(indxDesorp),1) AllStatesDesorpIni(:,3:end)]; % Unoccupied state (final)

EAllStatesDesorpIni = EAllMicroStates(indxDesorp);
EmicrostateBothUnocc = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnTCell{:},AllStatesDesorpFin);

PFcnTermDesorp = PFcnTerm(indxDesorp);

DErxn = EAllStatesDesorpIni - EmicrostateBothUnocc -2*(ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the adsorption event
Eafwd = EaO2ads+ProxFacO2ads*(DErxn-2*DEadsO-h);
Earev = Eafwd - DErxn;
AvgekO2des = sum(PFcnTermDesorp.*exp(-beta*Earev));

% Adsorption can happen in the following states
indxAdsorp = find(AllStatesArray(:,1) == 0 & AllStatesArray(:,2) == 0).';
AllStatesAdsorpIni = AllStatesArray(indxAdsorp,:);
AllStatesAdsorpFin = [repmat([1 1],length(indxAdsorp),1) AllStatesAdsorpIni(:,3:end)]; % Occupied state (final)

EAllStatesAdsorpIni = EAllMicroStates(indxAdsorp);
EmicrostateBothOccup = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnTCell{:},AllStatesAdsorpFin);

PFcnTermAdsorp = PFcnTerm(indxAdsorp);

DErxn = EmicrostateBothOccup - EAllStatesAdsorpIni -2*(ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the adsorption event
Eafwd = EaO2ads+ProxFacO2ads*(DErxn-2*DEadsO-h);
Earev = Eafwd - DErxn;
AvgekO2ads = sum(PFcnTermAdsorp.*exp(-beta*Eafwd));

% NO oxidation can happen from the following states
indxNOoxi = find(AllStatesArray(:,1) == 1).';
AllStatesNOoxiIni = AllStatesArray(indxNOoxi,:);
AllStatesNOoxiFin = [zeros(length(indxNOoxi),1) AllStatesNOoxiIni(:,2:end)]; % Unoccupied state (final)

EAllStatesNOoxiIni = EAllMicroStates(indxNOoxi);
EmicrostateCentrUnocc = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnTCell{:},AllStatesNOoxiFin);

PFcnTermNOoxi = PFcnTerm(indxNOoxi);

DENOoxi0  = EnrgGasSpecies(3) - (EnrgGasSpecies(2) + DEadsO);
DErxn = EnrgGasSpecies(3) - EnrgGasSpecies(2) + EmicrostateCentrUnocc - EAllStatesNOoxiIni + (ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the NO oxidation event
Eafwd = max([zeros(size(DErxn)),DErxn,EaNOoxi+ProxFacNOoxi*(DErxn-DENOoxi0)],[],2);
Earev = Eafwd - DErxn;
AvgekNOoxi = sum(PFcnTermNOoxi.*exp(-beta*Eafwd));

% NO2 reduction can happen from the following states
indxNO2re = find(AllStatesArray(:,1) == 0).';
AllStatesNO2reIni = AllStatesArray(indxNO2re,:);
AllStatesNO2reFin = [ones(length(indxNO2re),1) AllStatesNO2reIni(:,2:end)]; % Occupied state (final)

EAllStatesNO2reIni = EAllMicroStates(indxNO2re);
EmicrostateCentrOccup = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnTCell{:},AllStatesNO2reFin);

PFcnTermNO2re = PFcnTerm(indxNO2re);

DENOoxi0  = EnrgGasSpecies(3) - (EnrgGasSpecies(2) + DEadsO);
DErxn = EnrgGasSpecies(3) - EnrgGasSpecies(2) + EAllStatesNO2reIni - EmicrostateCentrOccup + (ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the NO oxidation event
Eafwd = max([zeros(size(DErxn)),DErxn,EaNOoxi+ProxFacNOoxi*(DErxn-DENOoxi0)],[],2);
Earev = Eafwd - DErxn;
AvgekNO2re = sum(PFcnTermNO2re.*exp(-beta*Earev));

% disp(['DErxn = ' num2str(DErxn(37)) ' ' num2str(Eafwd(37)-Earev(37))]);

% Normalisations

% Theta = PFcnEps0Occup/PFcn;
CorrelOccup = PFcnEps0Eps1Occup/PFcn;
CorrelUnocc = PFcnEps0Eps1Unocc/PFcn;
AvgekO2ads = AvgekO2ads/PFcnEps0Eps1Unocc;
AvgekO2des = AvgekO2des/PFcnEps0Eps1Occup;
AvgekNOoxi = AvgekNOoxi/PFcnEps0Occup;
AvgekNO2re = AvgekNO2re/PFcnEps0Unocc;

end