function out1 = HamiltonianK3NN2corr(H0,H1,h,nn,g1,g2,epsl)

% Implementation of the Hamiltonian of Kikuchi with up to 2NN interactions 
% epsl is the state vector. The energy of the microstate is returned (out1)
% In this implementation the central site's index is 1

if nn ~= 6, error('This function is only valid for nn = 6!'), end
if size(epsl,2) ~= 19, error('The length of the state vector should be 19!'), end

out1 = H0 + H1*epsl(:,1) + ... % central site 1-body term
    (H1)*sum(epsl(:,2:7),2) + ... % 1NN sites, 1-body + effective field
    (H1+g1)*sum(epsl(:,8:13),2) + ... % 2NN sites, 1-body + effective field
    (H1+g2)*sum(epsl(:,14:19),2) + ... % 3NN sites, 1-body + effective field
    h*epsl(:,1).*sum(epsl(:,2:7),2); % central-1NN 2-body terms

EdgeInter1 = 0;
for i = 2:6
    EdgeInter1 = EdgeInter1 + h*epsl(:,i).*epsl(:,i+1);
end
EdgeInter1 = EdgeInter1 + h*epsl(:,7).*epsl(:,2); % 1NN-1NN 2-body terms

EdgeInter12 = 0;
for i = 8:12
    ineigh1 = i-6;
    ineigh2 = ineigh1 + 1;
    EdgeInter12 = EdgeInter12 + ...
        (h)*epsl(:,i).*(epsl(:,ineigh1) + epsl(:,ineigh2));
end
EdgeInter12 = EdgeInter12 + ...
    (h)*epsl(:,13).*(epsl(:,7) + epsl(:,2)); % 1NN-2NN 2-body terms

EdgeInter13 = 0;
for i = 14:19
    ineigh1 = i-12;
    EdgeInter13 = EdgeInter13 + (h)*epsl(:,i).*epsl(:,ineigh1);
end % 1NN-3NN 2-body terms

EdgeInter23 = 0;
for i = 15:19
    ineigh1 = i-7;
    ineigh2 = ineigh1+1;
    EdgeInter23 = EdgeInter23 + ...
        (h)*epsl(:,i).*(epsl(:,ineigh1) + epsl(:,ineigh2));
end 
EdgeInter23 = EdgeInter23 + ...
    (h)*epsl(:,14).*(epsl(:,13) + epsl(:,8)); % 2NN-3NN 2-body terms

% Sum up all terms:
out1 = out1 + EdgeInter1 + EdgeInter12 + EdgeInter13 + EdgeInter23;

end