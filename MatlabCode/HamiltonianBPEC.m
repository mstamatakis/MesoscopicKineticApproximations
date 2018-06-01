function out1 = HamiltonianBPEC(H0,H1,h,nn,g1,c1,epsl)

% Implementation of the Hamiltonian of Bethe-Peiers with edge interactions 
% epsl is the state vector. The energy of the microstate is returned (out1)

if nn ~= 6, error('This function is only valid for nn = 6!'), end
if size(epsl,2) ~= 7, error('The length of the state vector should be 7!'), end

out1 = H0 + H1*epsl(:,1) + (H1+g1)*sum(epsl(:,2:nn+1),2) + h*epsl(:,1).*sum(epsl(:,2:nn+1),2);

EdgeInter = (h+c1)*epsl(:,2).*epsl(:,nn+1);
for i = 2:nn
    EdgeInter = EdgeInter + (h+c1)*epsl(:,i).*epsl(:,i+1);
end

out1 = out1 + EdgeInter;

end