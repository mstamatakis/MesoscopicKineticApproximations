function out1 = VibrationalPartitionFunction(vibfreqs,Temp)

global hplanck clight kboltz

indxs = vibfreqs > 0;

Theta = (hplanck*clight/kboltz)*vibfreqs;

TcharT = Theta(indxs)./Temp;

% out1 = prod(exp(-TcharT./2)./(1 - exp(-TcharT))); % This includes the ZPE
out1 = prod(1./(1 - exp(-TcharT))); % This does not include the ZPE

end