function out1 = HamiltonianBP(H0,H1,h,nn,g1,epsl)

% Implementation of the Hamiltonian of Bether Peierls without edge interactions
% epsl is the state vector. The energy of the microstate is returned (out1)
% In this implementation the central site's index is 1

if nn ~= 6, error('This function is only valid for nn = 6!'), end
if size(epsl,2) ~= 7, error('The length of the state vector should be 7!'), end

out1 = H0 + H1*epsl(:,1) + ... % central site 1-body term
    (H1+g1)*sum(epsl(:,2:7),2) + ... % 1NN sites, 1-body + effective field
    h*epsl(:,1).*sum(epsl(:,2:7),2); % central-1NN 2-body terms

end