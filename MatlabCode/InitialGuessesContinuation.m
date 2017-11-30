function [H0v,H0,gcguessv,gcguess] = ...
    InitialGuessesContinuation(ndeg,i,k,beta,xVarRange,PFcn,gcsoln,...
        H0v,H0,gcguessv,gcguess)

% Create good initial guesses via the continuation/interpolation scheme
if i <= ndeg % First few iterations where we need to populate the guess vectors
    
    % The guess vector will be set to the latest solution (0th
    % order continuation scheme) and the H0 to minus the free energy
    
    H0v{k}(i) = 1/beta*(log(PFcn)) + H0(k);
    gcguessv{k}(:,i) = gcsoln.';
    
    H0(k) = H0v{k}(i);
    gcguess{k} = gcguessv{k}(:,i);
    
%     disp(gcguessv{k});
%     disp(H0v{k});
%     disp(gcguess{k});
%     disp(H0(k));
%     pause
    
else
    if i == length(xVarRange), return, end
    % Apply the Interpolation/Continuation scheme:
    
    % ** Store the new solution values
    H0v{k}(ndeg+1) = 1/beta*(log(PFcn)) + H0(k);
    gcguessv{k}(:,ndeg+1) = gcsoln.';
    
    % ** Create the new guess values via extrapolation
    pCoef = polyfit(xVarRange(i-ndeg:i),H0v{k},ndeg);
    H0(k) = polyval(pCoef,xVarRange(i+1));
    for j = 1:length(gcguess{k})
        pCoef = polyfit(xVarRange(i-ndeg:i),gcguessv{k}(j,:),ndeg);
        gcguess{k}(j) = polyval(pCoef,xVarRange(i+1));
    end
    
%     disp(gcguessv{k});
%     disp(H0v{k});
%     disp(gcguess{k});
%     disp(H0(k));
%     pause
    
    % ** Shift back
    for j = 1:ndeg
        H0v{k}(j) = H0v{k}(j+1);
        gcguessv{k}(:,j) = gcguessv{k}(:,j+1);
    end
end

end