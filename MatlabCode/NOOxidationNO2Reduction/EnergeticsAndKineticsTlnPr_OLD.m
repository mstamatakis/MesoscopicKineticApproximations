function [EnrgGasSpecies,QvibO,ZPEO,...
    ApreO2ads,PERatioO2adsdes,EaO2ads,ProxFacO2ads,...
    ApreOdiff,PERatioOdiff,EaOdiff,ProxFacOdiff,ApreOdiff_sc,...
    ApreNOoxi,PERatioNOoxiNO2red,EaNOoxi,ProxFacNOoxi,ApreNOoxi_sc] = ...
    EnergeticsAndKineticsTlnPr(Temp,nFigs,logyNO2yNO)

global hplanck clight kboltz

hplanck = 4.135667516e-15; % eV*s
clight = 2.99792458e+10; % cm/s
% kboltz = 8.6173324e-5; % eV/K
kboltz = 1.380650321752056d-23*6.24150974d+18; % eV/K

% hplanck = 9.5306980570977e-14; % kcal/mol*s
% clight = 2.99792458e+10; % cm/s
% kboltz = 1.98587514214586e-3; % kcal/mol/K

ConvFPE = 6.24150964712042e-7; % Angstrom^3/eV = 6.24150964712042e-7*bar^-1
% ConvFPE = 1.43836378618515E-5; % Angstrom^3/(kcal/mol) = 1.43836378618515E-5*bar^-1

ElementNames = {'H','C','N','O'};
% 1 amu = 1.03642696948404e-28 eV*s^2/Angstrom^2
% 1 amu = 2.38845904951737e-27 kcal/mol*s^2/Angstrom^2
ElementAMass = [1.00794, 12.0107, 14.00672, 15.9994]*1.03642696948404e-28; % eV*s^2/Angstrom^2

NamesGasSpecies = {
    'O2'        % 1
    'NO'        % 2
    'NO2'       % 2
    };

GasSpeciesElemMatrix = [
    0 0 0 2;
    0 0 1 1;
    0 0 1 2;
    ];

EnrgGasSpecies = [
    Go_O2(Temp) - Go_O2(Temp)
    Go_NO(Temp) - 1/2*Go_O2(Temp) - 1/2*Go_N2(Temp)
    Go_NO2(Temp) - Go_O2(Temp) - 1/2*Go_N2(Temp)
    ]; % We use the O2 and N2 Gibbs free energies as reference

% EnrgGasSpecies = [
%     Go_O2(Temp)
%     Go_NO(Temp)
%     Go_NO2(Temp)
%     ];

% EnrgGasSpecies = [
%     Go_O2(Temp) - 2*Go_NO2(Temp) + 2*Go_NO(Temp)
%     Go_NO(Temp) - Go_NO(Temp)
%     Go_NO2(Temp) - Go_NO2(Temp)
%     ]; % We use the NO and NO2 Gibbs free energies as reference

% EnrgGasSpecies = [
%     Go_O2(Temp) - Go_O2(Temp)
%     Go_NO(Temp) - Go_NO(Temp)
%     Go_NO2(Temp) - Go_NO(Temp) - 1/2*Go_O2(Temp)
%     ]; % We use the O2 and NO Gibbs free energies as reference

GasSpeciesMMass = GasSpeciesElemMatrix*ElementAMass.';

VibFreqsO = [429 380 377];
QvibO = VibrationalPartitionFunction(VibFreqsO,Temp);
ZPEO = ZeroPointEnergy(VibFreqsO);

% EO2gasRef = - Go_O2(Temp) - 8.903471*0.01036427230133138; % eV/O2
EO2gasRef = - Go_O2(Temp); % eV/O2

% O2 adsorption step
Asite = 1/1.45; % Angstrom^2
kapa = 1/6; % orientation sticking coefficient
ApreO2ads = kapa*ConvFPE*Asite/sqrt(2*pi*GasSpeciesMMass(1)*kboltz*Temp);
% ApreO2des = kboltz*Temp/hplanck;
ApreO2des = ApreO2ads/QvibO^2*exp((EO2gasRef+2*ZPEO)/(kboltz*Temp));
PERatioO2adsdes = ApreO2ads/ApreO2des;
switch nFigs
    case 0
        EaO2ads =  0.020; % for conceptual model
    case 3
        EaO2ads =  0.020; % for 3FigCE
    case 4
        EaO2ads = -0.650; % for 4FigCE
    case 5
        EaO2ads = -0.780; % for 5FigCE
    case 7
        EaO2ads = -0.731; % for 7FigCE
    case 8
        EaO2ads = -0.472; % for 8FigCE
    case 12
        EaO2ads = -0.476; % for 12FigCE
    otherwise
        error('Invalid number of figures specified. Must be one of: 0, 3, 4, 5, 8 or 12.')
end
ProxFacO2ads = 1.000;

% O diffusion step
ApreOdiff = 0.1e-7*kboltz*Temp/hplanck;
PERatioOdiff = 1.0;
EaOdiff = 0.100;
ProxFacOdiff = 0.500;

% NO oxidation / NO2 reduction step
% ApreNOoxi = 2.328363491607911e+002*kboltz*Temp/hplanck/QvibO*exp((1/2*EO2gasRef+ZPEO)/(kboltz*Temp));
ApreNOoxi = 2e-2*kboltz*Temp/hplanck/QvibO*exp((1/2*EO2gasRef+ZPEO)/(kboltz*Temp));
PERatioNOoxiNO2red = 1.0/QvibO*exp((1/2*EO2gasRef+ZPEO)/(kboltz*Temp));
EaNOoxi = 0.000;
ProxFacNOoxi = 0.000;


Press = 1;
GasMolFracO2 = 0.1;
logyNO2yNO_equil = 1/2*log(Press*GasMolFracO2)-(2*Go_NO2(Temp) - 2*Go_NO(Temp) - Go_O2(Temp))/(2*kboltz*Temp);

if (Temp == 480)
    
    if logyNO2yNO < logyNO2yNO_equil % on the left of the equilibrium point
        
        ApreNOoxi_1 = 5.049191309179336e+016;
        ApreNOoxi_2 = 7.933538334782410e+008;
        logyNO2yNO_1 = -4.0;
        logyNO2yNO_2 = 8.0;
        
        ApreNOoxi_sc = (1/ApreNOoxi)*ApreNOoxi_1*(ApreNOoxi_2/ApreNOoxi_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));
    
        ApreOdiff_1 = 1.443672714377862e+004;
        ApreOdiff_2 = 5.968846757660124e-006;
        logyNO2yNO_1 = -4.0;
        logyNO2yNO_2 = 8.0;
        
        ApreOdiff_sc = (1/ApreOdiff)*ApreOdiff_1*(ApreOdiff_2/ApreOdiff_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));

    else % on the right of the equilibrium point
        
        ApreNOoxi_1 = 9.455328029780760e+011;
        ApreNOoxi_2 = 5.077819808750501e+012;
        logyNO2yNO_1 = 4.659193890237543;
        logyNO2yNO_2 = 8.285602152932793;
        
        ApreNOoxi_sc = (1/ApreNOoxi)*ApreNOoxi_1*(ApreNOoxi_2/ApreNOoxi_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));    

        ApreOdiff_1 = 1.595074953674987e-002;
        ApreOdiff_2 = 1.900194873843102e-002;
        logyNO2yNO_1 = 4.659193890237543;
        logyNO2yNO_2 = 8.285602152932793;
        
        ApreOdiff_sc = (1/ApreOdiff)*ApreOdiff_1*(ApreOdiff_2/ApreOdiff_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));
        
    end

elseif (Temp == 580)
    
    if logyNO2yNO < logyNO2yNO_equil % on the left of the equilibrium point
        
        ApreNOoxi_1 = 2.042273658443739e+018;
        ApreNOoxi_2 = 3.749186548374233e+009;
        logyNO2yNO_1 = -4.0;
        logyNO2yNO_2 = 8.0;
        
        ApreNOoxi_sc = (1/ApreNOoxi)*ApreNOoxi_1*(ApreNOoxi_2/ApreNOoxi_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));    

        ApreOdiff_1 = 5.167603364673340e+005;
        ApreOdiff_2 = 6.549875364009890e-005;
        logyNO2yNO_1 = -4.0;
        logyNO2yNO_2 = 8.0;
        
        ApreOdiff_sc = (1/ApreOdiff)*ApreOdiff_1*(ApreOdiff_2/ApreOdiff_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));

    else % on the right of the equilibrium point
        
        ApreNOoxi_1 = 1.108632992950299e+015;
        ApreNOoxi_2 =  5.742136515546386e+015;
        logyNO2yNO_1 = 2.0;
        logyNO2yNO_2 = 8.0;
        
        ApreNOoxi_sc = (1/ApreNOoxi)*ApreNOoxi_1*(ApreNOoxi_2/ApreNOoxi_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));    

        ApreOdiff_1 = 1.493997881096791e+002;
        ApreOdiff_2 = 7.738124195965912e+002;
        logyNO2yNO_1 = 2.0;
        logyNO2yNO_2 = 8.0;
        
        ApreOdiff_sc = (1/ApreOdiff)*ApreOdiff_1*(ApreOdiff_2/ApreOdiff_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));
        
    end
    
elseif (Temp >= 680)
    
    if logyNO2yNO < logyNO2yNO_equil % on the left of the equilibrium point
        
        ApreNOoxi_1 = 3.369750827903273e+019;
        ApreNOoxi_2 = 1.072301532145342e+011;
        logyNO2yNO_1 = -4.0;
        logyNO2yNO_2 = 8.0;
        
        ApreNOoxi_sc = (1/ApreNOoxi)*ApreNOoxi_1*(ApreNOoxi_2/ApreNOoxi_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));
    
        ApreOdiff_1 = 1.384282104932539e+007;
        ApreOdiff_2 = 1.736861633116972e-003;
        logyNO2yNO_1 = -4.0;
        logyNO2yNO_2 = 8.0;
        
        ApreOdiff_sc = (1/ApreOdiff)*ApreOdiff_1*(ApreOdiff_2/ApreOdiff_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));

    else % on the right of the equilibrium point
        
        ApreNOoxi_1 = 5.518055883663435e+016;
        ApreNOoxi_2 = 1.118055883663435e+018;
        logyNO2yNO_1 = 0.5;
        logyNO2yNO_2 = 8.0;
        
        ApreNOoxi_sc = (1/ApreNOoxi)*ApreNOoxi_1*(ApreNOoxi_2/ApreNOoxi_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));
    
        ApreOdiff_1 = 1.136184959367195e+003;
        ApreOdiff_2 = 1.136184959367195e+004;
        logyNO2yNO_1 = 0.5;
        logyNO2yNO_2 = 8.0;
        
        ApreOdiff_sc = (1/ApreOdiff)*ApreOdiff_1*(ApreOdiff_2/ApreOdiff_1)^ ...
            ((logyNO2yNO-logyNO2yNO_1)/(logyNO2yNO_2-logyNO2yNO_1));
        
    end
    
end

ApreNOoxi = ApreNOoxi*ApreNOoxi_sc;
ApreOdiff = ApreOdiff*ApreOdiff_sc;

return

end