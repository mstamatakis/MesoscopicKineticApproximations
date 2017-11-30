function outResid = ResidualTheta(ApproxIdent,CorrelLHS,CorrelRHS,nsites,beta,...
    H0,H1,h,nn,Theta,muGuess,HamiltPars)

%% Generalised residual constructor. 


nEq = length(CorrelLHS);

if strcmp(ApproxIdent,'MF')
    if nEq == 0 % this is valid only for the MF approximation
        z = HamiltPars; % This Hamiltonian has only one parameter: the coverage
        outResid(1) = z - 1/(1+exp(beta*(H1-muGuess+h*nn*z)));
        outResid(2) = log(Theta) - log(z);
        return
    else
        error('For the mean-field (MF) approximation no equations can be specified.')
    end
end

if nEq == 0
    error('No equations specified. This is valid only for the mean-field (MF) approximation.')
end

if length(CorrelLHS) ~= length(CorrelRHS)
    error('CorrelLHS and CorrelRHS need to have the same number of elements.')
end

if length(CorrelLHS) ~= length(CorrelRHS)
    error('The number of equations does not seem to equal the number of unknowns.')
end

HamiltParsCell = num2cell(HamiltPars);

[PFcn,PFcnEps1_Z,PFcnEps0_Z,PFcnPartial] = ...
    PartitionFunction_BruteForce(ApproxIdent,...
    {CorrelLHS{:},CorrelRHS{:}},...
    nsites,beta,H0,H1,h,nn,muGuess,HamiltParsCell{:});

outResid = NaN(1,nEq);
for i = 1:nEq
%     % The equations to be solved look like:
%     % 1 - CorrelationLHS/CorrelationRHS = 0
%     outResid(i) = 1-PFcnPartial(i)/PFcnPartial(nEq+i);

    % The equations to be solved look like:
    % log(CorrelationLHS) - log(CorrelationRHS) = 0
    % This has improved convergence
    outResid(i) = log(PFcnPartial(i)) - log(PFcnPartial(nEq+i));
end

outResid(nEq+1) = log(Theta) - log(PFcnEps1_Z/PFcn);

end
