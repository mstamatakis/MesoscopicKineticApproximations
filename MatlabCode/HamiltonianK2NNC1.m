function out1 = HamiltonianK2NNC1(H0,H1,h,nn,g1,g2,c1,c2,epsl)

% Implementation of the Hamiltonian of Kikuchi with up to 2NN interactions 
% epsl is the state vector. The energy of the microstate is returned (out1)
% In this implementation the central site's index is 1

if nn ~= 6, error('This function is only valid for nn = 6!'), end
if size(epsl,2) ~= 13, error('The length of the state vector should be 13!'), end

out1 = H0 + H1*epsl(:,1) + ... % central site 1-body term
    (H1+g1)*sum(epsl(:,2:7),2) + ... % 1NN sites, 1-body + effective field
    (H1+g2)*sum(epsl(:,8:13),2) + ... % 2NN sites, 1-body + effective field
    h*epsl(:,1).*sum(epsl(:,2:7),2); % central-1NN 2-body terms

EdgeInter1 = 0;
for i = 2:6
    EdgeInter1 = EdgeInter1 + (h+c1)*epsl(:,i).*epsl(:,i+1);
end
EdgeInter1 = EdgeInter1 + (h+c1)*epsl(:,7).*epsl(:,2); % 1NN-1NN 2-body terms

EdgeInter12 = 0;
for i = 8:12
    ineigh1 = i-6;
    ineigh2 = ineigh1 + 1;
    EdgeInter12 = EdgeInter12 + ...
        (h+c2)*(epsl(:,i).*epsl(:,ineigh1) + epsl(:,i).*epsl(:,ineigh2));
end
EdgeInter12 = EdgeInter12 + ...
    (h+c2)*(epsl(:,13).*epsl(:,7) + epsl(:,13).*epsl(:,2)); % 1NN-2NN 2-body terms

out1 = out1 + EdgeInter1 + EdgeInter12;

end