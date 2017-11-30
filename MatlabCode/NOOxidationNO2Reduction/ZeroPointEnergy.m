function out1 = ZeroPointEnergy(vibfreqs)

global hplanck clight

indxs = vibfreqs > 0;

out1 = 1/2*hplanck*clight*sum(vibfreqs(indxs));

return

end