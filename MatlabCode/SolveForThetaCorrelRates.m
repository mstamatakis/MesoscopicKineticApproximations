function [PFcn,gcsoln,Theta,CorrelOccup,CorrelUnocc, ...
    AvgekO2ads,AvgekO2des] = ...
    SolveForThetaCorrelRates(beta,Temp,EaO2ads,DEadsO,ProxFacO2ads,H0,H1,...
    muOstar,h,nn,ZPEO,QvibO,ApproxIdent,...
    CorrelLHS,CorrelRHS,nsites,gcguess,options)

global kboltz

% First solve for the parameters of the approximation
Rfcn = @(gc) Residual(ApproxIdent,CorrelLHS,CorrelRHS,nsites,beta,H0,H1,h,nn,muOstar,gc);
gcsoln = fsolve(Rfcn,gcguess,options);
gcsolnCell = num2cell(gcsoln);

if isempty(CorrelLHS) % MF approximation treated specially
    Theta = gcsoln;
    PFcn = exp(-beta*H0) + exp(-beta*(H0 + DEadsO + nn*h*Theta-muOstar));
    DAMF = 2*(DEadsO + nn*h*Theta); % free energy difference between final and initial state of O2 adsorption
    CorrelOccup = Theta^2; % correlation of occupied sites
    CorrelUnocc = 1-2*Theta+CorrelOccup; % correlation of unoccupied sites
    AvgekO2ads = exp(-beta*(EaO2ads+ProxFacO2ads*(DAMF-(2*DEadsO+h))));
    AvgekO2des = exp(-beta*(EaO2ads-(2*DEadsO+h)-(1-ProxFacO2ads)*(DAMF-2*DEadsO)));
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
EAllMicroStates = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnCell{:},...
    AllStatesArray);

PFcnTerm = exp(-beta*(EAllMicroStates - muOstar*sum(AllStatesArray,2)));

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
EmicrostateBothUnocc = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnCell{:},AllStatesDesorpFin);

PFcnTermDesorp = PFcnTerm(indxDesorp);

DErxn = EAllStatesDesorpIni - EmicrostateBothUnocc -2*(ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the desorption event
Eafwd = EaO2ads+ProxFacO2ads*(DErxn-2*DEadsO-h);
Earev = Eafwd - DErxn;
AvgekO2des = sum(PFcnTermDesorp.*exp(-beta*Earev));

% Adsorption can happen in the following states
indxAdsorp = find(AllStatesArray(:,1) == 0 & AllStatesArray(:,2) == 0).';
AllStatesAdsorpIni = AllStatesArray(indxAdsorp,:);
AllStatesAdsorpFin = [repmat([1 1],length(indxAdsorp),1) AllStatesAdsorpIni(:,3:end)]; % Occupied state (final)

EAllStatesAdsorpIni = EAllMicroStates(indxAdsorp);
EmicrostateBothOccup = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnCell{:},AllStatesAdsorpFin);

PFcnTermAdsorp = PFcnTerm(indxAdsorp);

DErxn = EmicrostateBothOccup - EAllStatesAdsorpIni -2*(ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the adsorption event
Eafwd = EaO2ads+ProxFacO2ads*(DErxn-2*DEadsO-h);
Earev = Eafwd - DErxn;
AvgekO2ads = sum(PFcnTermAdsorp.*exp(-beta*Eafwd));

% for i = indxAdsorp
%     stateV = AllStatesArray(i,:);
%     stateVBothOcc = [1 1 stateV(3:end)]; % Occupied state (final)
% 
%     Emicrostate = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnCell{:},stateV);
%     EmicrostateBothOcc = Hamiltonian(ApproxIdent,H0,H1,h,nn,gcsolnCell{:},stateVBothOcc);
%     
%     PFcnTerm = exp(-beta*(Emicrostate-muOstar*sum(stateV)));
% 
%     DErxn = EmicrostateBothOcc - Emicrostate -2*(ZPEO-kboltz*Temp*log(QvibO)); % DErxn for the adsorption event
%     Eafwd = EaO2ads+ProxFacO2ads*(DErxn-2*DEadsO-h);
%     Earev = Eafwd - DErxn;
%     AvgekO2ads = AvgekO2ads + PFcnTerm*exp(-beta*Eafwd);
% end

Theta = PFcnEps0Occup/PFcn;
CorrelOccup = PFcnEps0Eps1Occup/PFcn;
CorrelUnocc = PFcnEps0Eps1Unocc/PFcn;
AvgekO2ads = AvgekO2ads/PFcnEps0Eps1Unocc;
AvgekO2des = AvgekO2des/PFcnEps0Eps1Occup;

end