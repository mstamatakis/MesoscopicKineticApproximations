function out1 = HamiltonianK2NNC2T1(H0,H1,h,nn,g1,g2,c1,c2,p1,p2,p3,epsl)

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

EdgeInterP11 = p1* ( ...
    epsl(:,2).*epsl(:,4) + ...
    epsl(:,2).*epsl(:,6) + ...
    epsl(:,3).*epsl(:,5) + ...
    epsl(:,3).*epsl(:,7) + ...
    epsl(:,4).*epsl(:,6) + ...
    epsl(:,5).*epsl(:,7) ...
    ); % 1NN-1NN 2nd nearest neighbor 2-body effective interactions

EdgeInterP22 = 0;
for i = 8:12
    EdgeInterP22 = EdgeInterP22 + p2*epsl(:,i).*epsl(:,i+1);
end
EdgeInterP22 = EdgeInterP22 +  p2*epsl(:,13).*epsl(:,8); % 2NN-2NN 2nd nearest neighbor 2-body effective interactions

TripleInter1 = p3*( ...
    epsl(:,2).*epsl(:,3).*epsl(:,8)+ ...
    epsl(:,3).*epsl(:,4).*epsl(:,9)+ ...
    epsl(:,4).*epsl(:,5).*epsl(:,10)+ ...
    epsl(:,5).*epsl(:,6).*epsl(:,11)+ ...
    epsl(:,6).*epsl(:,7).*epsl(:,12)+ ...
    epsl(:,7).*epsl(:,2).*epsl(:,13));

% Sum up all terms:
out1 = out1 + EdgeInter1 + EdgeInter12 + EdgeInterP11 + EdgeInterP22 + TripleInter1;

end