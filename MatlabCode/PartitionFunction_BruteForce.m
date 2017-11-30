function [PFcn,PFcnEps1_Z,PFcnEps0_Z,PFcnPartial] = ...
    PartitionFunction_BruteForce(varargin)

% This function assumes arguments are given as:
% [Approximation Identifier],[Correlations Requested],...
%         nsites,beta,H0,H1,h,nn,mu,[Hamiltonian Parameters]

% The Correlations Requested are given by a single cell array with vectors
% specifying the sites appearing in the pertinent correlation.

% The Hamiltonian Parameters pertain to the effective interactions and are 
% passed on to the Hamiltonian function in the same order as given in the 
% calling function.

% In all implementations of these scripts the central site's index is 1

ApproxIdent = varargin{1}; 
CorrelCell = varargin{2};
nsites = varargin{3};
beta = varargin{4};
H0 = varargin{5};
H1 = varargin{6};
h = varargin{7};
nn = varargin{8};
mu = varargin{9};

% Initialisation of variables
PFcn = 0;
PFcnEps1_Z = 0;
PFcnEps0_Z = 0;
PFcnPartial = zeros(1,length(CorrelCell));

% Create the array containing all occupancy states
AllStatesArray = AllStates(nsites);

% Calculate all quantities of interest using vectorised expressions
EAllMicroStates = Hamiltonian(ApproxIdent,H0,H1,h,nn,varargin{10:end},...
    AllStatesArray);

PFcnTerm = exp(-beta*(EAllMicroStates - mu*sum(AllStatesArray,2)));

PFcn = sum(PFcnTerm);
PFcnEps1_Z = sum(AllStatesArray(:,1).*PFcnTerm);
PFcnEps0_Z = sum((1-AllStatesArray(:,1)).*PFcnTerm);
for j = 1:length(CorrelCell)
    PFcnPartial(j) = sum(all(AllStatesArray(:,CorrelCell{j}) == 1,2).*PFcnTerm);
end

return











